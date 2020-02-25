#include "math_extra.h"

int doubleFactorial(int n) {
  int result = 1;
  for (; n > 0; n -= 2) result *= n;
  return result;
}

double risingFactorial (double x, int n) {
  double result = 1.0;
  for (int i = 0; i < n; i++) result *= x + i;
  return result;
}

double hypergeometric1f1 (double a, double b, double z) {
  double result = 0.0;
  double zPower = 1.0;
  for (int i = 0; i < 50; i++) {
    result += risingFactorial(a,i) / (risingFactorial(b,i) * risingFactorial(1,i)) * zPower;
    zPower *= z;
  }

  return result;
}

double sinc (double x) {
  if (x == 0.0) return 1.0;
  else return sin(x) / x;
}

double coth (double x) {
  double c = cosh(x);
  double s = sinh(x);
  if (isinf(c) && isinf(s)) return 1.0;
  return c / s;
}

double dawson (double x) {
  double xSq = x * x;

  double result = 0.0;
  double xPower = x;
  int sign = 1;
  int twoPower = 1;
  for (int i = 0; i < 20; i++) {
    result += sign * twoPower * xPower / doubleFactorial(2*i+1);
    xPower *= xSq;
    sign *= -1;
    twoPower *= 2;
  }

  return result;
}

double randomSineGaussian (double a) {
  uniform_real_distribution<double> uniform(0.0,1.0);
  double ng = -0.5 * (exp(-a*M_PI*M_PI) - 1.0) / a;
  double u1,u2;
  double ratio,x;
  
  // rejection sampling of sin(x)*exp(-a*x^2)

  do {
    u1 = uniform(rng);
    u2 = uniform(rng);

    // inverse transform sampling of x*exp(-a*x^2)

    x = sqrt(-log(1.0 - 2.0*a*ng*u1)/a);
    ratio = sinc(x);
  }
  while (u2 > ratio);
  if (uniform(rng) < 0.5) x *= -1;
  
  return x;
}

vector<double> operator + (vector<double> a, vector<double> b) {
  int size = a.size();
  vector<double> c(size);
  for (int i = 0; i < size; i++) c[i] = a[i] + b[i];
  return c;
}

vector<double> operator - (vector<double> a, vector<double> b) {
  int size = a.size();
  vector<double> c(size);
  for (int i = 0; i < size; i++) c[i] = a[i] - b[i];
  return c;
}

vector<double> operator * (double a, vector<double> b) {
  int size = b.size();
  vector<double> c(size);
  for (int i = 0; i < size; i++) c[i] = a * b[i];
  return c;
}

double operator * (vector<double> a, vector<double> b) {
  int size = a.size();
  double result = 0;
  for (int i = 0; i < size; i++) result += a[i] * b[i];
  return result;
}

vector<double> cross (vector<double> a, vector<double> b) {
  vector<double> c(3);
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return c;
}

vector<double> rotate (vector<double> r, vector<double> n, double theta) {
  double c = cos(theta);
  double s = sin(theta);
  return (1.0-c)*(n*r)*n + c*r - s*cross(n,r);
}

double length(vector<double> a) {
  return sqrt(a*a);
}

void normalise(vector<double> &a) {
  double len = length(a);
  int size = a.size();
  a = (1.0/len) * a;
}
