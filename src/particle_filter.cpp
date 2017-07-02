/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;
 // Friend function for easy printing of particle
std::ostream& operator<< (std::ostream& os, const Particle& P) {
  os << "Particle #" << P.id << " x: " << P.x << " y: " << P.y << " theta: " << P.theta << " weight: " << P.weight << std::endl;
  return os;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  num_particles = 10;
  weights.reserve(num_particles);
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i) {
    //Add particle(+noise) to vector of particles
    particles.push_back(Particle(i, dist_x(gen), dist_y(gen), dist_theta(gen), 1.0 ));
  }

  //DEBUG print of all init particles
  std::for_each(particles.begin(), particles.end(), [](auto& p) {std::cout << p;} );
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);

  for (auto& p : particles) {
    //Motion Prediction
    //yaw rate ==0
    if (abs(yaw_rate) < 0.001) {
      p.x += velocity*cos(p.theta)*delta_t;
      p.y += velocity*sin(p.theta)*delta_t;
    }
    //yaw_rate > 0
    else {
      p.x += (velocity / yaw_rate) * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
      p.y += (velocity / yaw_rate) * (-cos(p.theta + yaw_rate*delta_t) + cos(p.theta));
    }
    p.theta += yaw_rate*delta_t;

    //Add gaussian noise
    p.x += dist_x(gen);
    p.y += dist_y(gen);
    p.theta += dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(double sensor_range, std::vector<LandmarkObs>& observations, Map& map) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//For every particle
	//Transform observation from vehicle to map
	//Find best match for the given observation

	for (int i = 0; i < num_particles; ++i) {
		double particle_theta = particles[i].theta;
		double particle_x = particles[i].x;
		double particle_y = particles[i].y;
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		//Transform observation from vehicle to map
		for (int obs = 0; obs < observations.size(); ++obs) {
			double obs_x = observations[obs].x;
			double obs_y = observations[obs].y;
			double obs_trans_x = obs_x*cos(particle_theta) + obs_y*sin(particle_theta) + particle_x;
			double obs_trans_y = obs_x*sin(particle_theta) - obs_y*cos(particle_theta) + particle_y;

			//For every observation, find the best landmark in the map
			//if the landmark is within sensor_range from current position, only check then
			
			double best_distance = std::numeric_limits<double>::max();
			int best_landmark_id = -1;
			bool found_landmark = false; 
			double distance = 0; //Distance of observation from landmark
			for (int l = 0; l < map.landmark_list.size(); ++l) {
				if (abs(map.landmark_list[l].x_f - particles[i].x) < sensor_range && abs(map.landmark_list[l].y_f - particles[i].y) < sensor_range) {
					distance = abs(map.landmark_list[l].x_f - obs_trans_x) + abs(map.landmark_list[l].y_f - obs_trans_y);
					if (distance < best_distance) {
						best_distance = distance;
						best_landmark_id = l;
						found_landmark = true;
					}
				}
			}
			if (found_landmark) {
				particles[i].sense_x.push_back(obs_trans_x);
				particles[i].sense_y.push_back(obs_trans_y);
				particles[i].associations.push_back(best_landmark_id);
			}	
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html


	dataAssociation(sensor_range, observations, map);

	//Common denominator used by all measurements;
	//We can use a simplified version of the calc because std_xy = 0, we only have std_x and std_y
	double denom = pow(2 * M_PI, 2.0) * std_landmark[0] * std_landmark[1];


	double particle_weight;
	for (int i = 0; i < num_particles; ++i) {
		particle_weight = 1.0;
		for (int j = 0; j < particles[i].associations.size(); ++j) {
			double d_x = particles[i].sense_x[j] - map.landmark_list[particles[i].associations[j]].x_f;
			double d_y = particles[i].sense_y[j] - map.landmark_list[particles[i].associations[j]].y_f;


			//Calculate weight of particle
			particle_weight *= d_x*d_x*std_landmark[1] * std_landmark[1] + d_y*d_y*std_landmark[0] * std_landmark[0];
			particle_weight /= denom;
		}
		particles[i].weight = particle_weight;
		weights[i] = particle_weight; // This weight vector is used for resampling
	}
		
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::vector<Particle> resampledParticles(num_particles);
	//We need to find a resampled particle for every particle
	double w_max = 2 * (*std::max_element(weights.begin(), weights.end()));
	std::discrete_distribution<double> distrib(weights.begin(), weights.end());
	for (int i = 0; i < num_particles; ++i) {
		resampledParticles[i] = particles[distrib(gen)];
	}
	
	//move the resampled data back to original
	particles = std::move(resampledParticles);
}

Particle ParticleFilter::SetAssociations(Particle& particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
