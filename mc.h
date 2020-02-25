#ifndef mc_header
#define mc_header

#include <cmath>
#include <random>
#include <vector>

#include "constants.h"
#include "math_extra.h"

using namespace std;

bool rosenbluthBending(vector<vector<double>>&);

vector<vector<double>> generateBendingSegments(vector<double>, int);
vector<vector<double>> generateBendingSegments(int);

double bendingEnergy (vector<double>, vector<double>);
double electricEnergy (vector<double>);
double bendingProbability (vector<double>, vector<double>);
double electricProbability (vector<double>);

#endif
