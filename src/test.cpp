#include <iostream>
#include <ctime>
#include <iomanip>
#include <random>
#include "particle_filter.h"

using namespace std;

void test1() {
  Particle P = Particle(0, 1, 2, 3, 4);
  cout << P << endl;
  Particle Q = P;
  cout << Q << endl;
}

void test2() {
  default_random_engine gen;
  normal_distribution<double> dist(0, 1);
  cout << dist(gen) << endl;
}

int runTests() {
  test1();
  test2();
  return 1;
}
