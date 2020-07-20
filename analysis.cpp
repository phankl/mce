#include "analysis.h"

vector<double> meanSquaredAngle (vector<vector<double>> configuration) {
  vector<double> msa(segmentNumber);
  if (energyMode == "square") {
    for (int i = 0; i < segmentNumber; i++) {
      msa[i] = 0.0;
      for (int j = 0; j < dimension-1; j++)
        msa[i] += pow(configuration[i][j],2);
    }
  }
  else if (energyMode == "cosine") {
    if (dimension == 2) {
    for (int i = 0; i < segmentNumber; i++) 
      msa[i] = pow(configuration[i][0],2);
    }
  }

  return msa;
}

vector<double> cosAverage (vector<vector<double>> configuration) {
  vector<double> c(segmentNumber);
  if (energyMode == "square") {
    for (int i = 0; i < segmentNumber; i++) {
      double thetaSquared = 0.0;
      for (int j = 0; j < dimension-1; j++) 
        thetaSquared += pow(configuration[i][j],2);
      c[i] = hypergeometric0f1(0.5*(dimension-1),-0.25*(dimension-1)*thetaSquared);
    }
  }
  else if (energyMode == "cosine") {
    for (int i = 0; i < segmentNumber; i++) 
      c[i] = fabs(cos(configuration[i][0]));
  }
  return c;
}

vector<double> cosCorrelation (vector<vector<double>> configuration) {
  vector<double> correlation(segmentNumber);
  if (energyMode == "square") {
    for (int i = 0; i < segmentNumber; i++) {
      double thetaSquared = 0.0;
      for (int j = 0; j < dimension-1; j++) 
        thetaSquared += pow(configuration[i][j]-configuration[0][j],2);
      correlation[i] = hypergeometric0f1(0.5*(dimension-1),-0.25*(dimension-1)*thetaSquared);
    } 
  }
  else if (energyMode == "cosine") {
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
  }

  return correlation;
}

double extension (vector<vector<double>> configuration) {
  double z = 0.0;
  for (int i = 0; i < segmentNumber; i++)
    z += cos(configuration[i][0]);

  return fabs(z) / segmentNumber;
}

double meanSquaredAngleTheory (double s) {
  double alpha = sqrt(b*order/stiffness);

  if (field == 0.0) return (dimension-1)*s / (stiffness*beta);
  return (dimension-1)*cosh((chainLength-s)*alpha)*cosh(s*alpha) / (sinh(chainLength*alpha)*beta*sqrt(stiffness*b*order));
}

double cosAverageTheory (double s) {
  if (field == 0.0) return 0.0;
  double alpha = sqrt(b * order / stiffness);
  double z = -cosh((chainLength-s)*alpha) * cosh(s*alpha) / (2.0*sinh(chainLength*alpha)*beta*stiffness*alpha);
  return exp((dimension-1)*z);
}

double cosSquaredAverageTheory (double s) {
  if (field == 0.0) return 0.0;
  double alpha = sqrt(b * order / stiffness);
  double z = -2.0 * cosh((chainLength-s)*alpha) * cosh(s*alpha) / (sinh(chainLength*alpha)*beta*stiffness*alpha);
  if (dimension == 2) return 0.5 * (1.0 + exp(z));
  else return hypergeometric1f1(0.5*(dimension-2.0),dimension-2.0,(dimension-1)*z);
}

double cosCorrelationTheory (double s) {
  if (field == 0.0) {
    return exp(-0.5*s*(dimension-1)/(stiffness*beta));
  }
  else {
    double arg = exp(2.0*chainLength*alpha) + 3.0*exp(2.0*s*alpha) - exp(3.0*s*alpha) -  3.0*exp((2.0*chainLength + s)*alpha);
    arg *= exp(-2.0*s*alpha) * (exp(s*alpha) - 1.0)/ (4.0*beta*sqrt(stiffness*b*order)*(exp(2.0*chainLength*alpha) - 1.0));
    return exp((dimension-1.0)*arg);
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
