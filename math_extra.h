#ifndef math_extra_header
#define math_extra_header

#include <cmath>
#include <iostream>
#include <vector>

#include "constants.h"

using namespace std;

int doubleFactorial(int);
double risingFactorial(double, int);
double hypergeometric1f1(double, double, double);
double sinc(double);
double coth(double);
double dawson(double);

double randomSineGaussian(double);

vector<double> operator + (vector<double>, vector<double>);
vector<double> operator - (vector<double>, vector<double>);
vector<double> operator * (double, vector<double>);
double operator * (vector<double>, vector<double>);

vector<double> cross (vector<double>, vector<double>);
vector<double> rotate (vector<double>, vector<double>, double);

double length(vector<double>);
void normalise(vector<double>&);

#endif
