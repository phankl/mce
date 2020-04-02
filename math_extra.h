#ifndef math_extra_header
#define math_extra_header

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "constants.h"

using namespace std;

int doubleFactorial(int);
double risingFactorial(double, int);
double hypergeometric1f1(double, double, double);
double hypergeometric0f1(double, double);
double sinc(double);
double coth(double);
double dawson(double);

double randomSineGaussian(double);
double randomGaussianCosine(double);
double randomSineGaussianCosine(double);
void initRandomGaussianCosine(double);

vector<double> operator + (vector<double>, vector<double>);
vector<double> operator - (vector<double>, vector<double>);
vector<double> operator * (double, vector<double>);
double operator * (vector<double>, vector<double>);

vector<double> cross (vector<double>, vector<double>);
vector<double> rotate (vector<double>, vector<double>, double);

double length(vector<double>);
void normalise(vector<double>&);

double newton(function<double(double)>,double);

#endif
