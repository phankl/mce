#ifndef analysis_header
#define analysis_header

#include <cmath>
#include <iostream>
#include <vector>

#include "constants.h"
#include "math_extra.h"

using namespace std;

vector<double> cosAverage(vector<vector<double>>);
vector<double> cosCorrelation(vector<vector<double>>);
vector<double> meanSquaredAngle(vector<vector<double>>);

double cosAverageTheory(double);
double cosCorrelationTheory(double);
double meanSquaredAngleTheory(double);

double energy(vector<vector<double>>);

#endif
