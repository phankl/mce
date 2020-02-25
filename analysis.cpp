#include "analysis.h"

vector<double> meanSquaredAngle (vector<vector<double>> configuration) {
  vector<double> msa(segmentNumber);
  if (field == 0.0) 
    for (int i = 0; i < segmentNumber; i++) 
      msa[i] = pow(configuration[i][0] - configuration[0][0],2);
  else
    for (int i = 0; i < segmentNumber; i++) 
      msa[i] = pow(configuration[i][0],2);

  return msa;
}

vector<double> cosAverage (vector<vector<double>> configuration) {
  vector<double> c(segmentNumber);
  for (int i = 0; i < segmentNumber; i++) c[i] = cos(configuration[i][0]);
  return c;
}

vector<double> cosCorrelation (vector<vector<double>> configuration) {
  vector<double> correlation(segmentNumber);
  if (dimension == 2)
    for (int i = 0; i < segmentNumber; i++)
      correlation[i] = cos(configuration[i][0] - configuration[0][0]);
  else if (dimension == 3) {
    double theta1 = configuration[0][0];
    double phi1 = configuration[0][1];
    double c1 = cos(theta1);
    double s1 = sin(theta1);
    for (int i = 0; i < segmentNumber; i++) {
      double theta2 = configuration[i][0];
      double phi2 = configuration[i][1];
      double c2 = cos(theta2);
      double s2 = sin(theta2);
      correlation[i] = c1*c2 + cos(phi1-phi2)*s1*s2;
    }
  }

  return correlation;
}

double meanSquaredAngleTheory (double s) {
  double b = susceptibility * pow(field,order);
  double alpha = sqrt(b*order/stiffness);

  if (field == 0.0) return (dimension-1)*s / (stiffness*beta);
  return (dimension-1)*cosh((chainLength-s)*alpha)*cosh(s*alpha) / (sinh(chainLength*alpha)*beta*sqrt(stiffness*b*order));
}

double cosAverageTheory (double s) {
  if (field == 0.0) return 0.0;
  double b = susceptibility * pow(field,order);
  double alpha = sqrt(b * order / stiffness);
  double z = -cosh((chainLength-s)*alpha) * cosh(s*alpha) / (2.0*sinh(chainLength*alpha)*beta*stiffness*alpha);
  return hypergeometric1f1(0.5*(dimension-1),0.5,z);
}

double cosCorrelationTheory (double s) {
  if (field == 0.0) {
    return exp(-0.5*s*(dimension-1)/(stiffness*beta));
  }
  else {
    if (dimension == 2) {
      double b = susceptibility * pow(field,order);
      return exp(-tanh(0.5*s*sqrt(b*order/stiffness)) / (beta*sqrt(stiffness*b*order)));
    }
  }
}

double energy (vector<vector<double>> configuration) {
  double result = electricEnergy(configuration[0]);
  for (int i = 1; i < segmentNumber; i++) {
    result += bendingEnergy(configuration[i-1],configuration[i]);
    result += electricEnergy(configuration[i]);
  }
  return result;
}
