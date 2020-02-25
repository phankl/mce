#include "mc.h"

bool rosenbluthBending (vector<vector<double>> &oldConfiguration) {
  vector<vector<double>> newConfiguration(segmentNumber,vector<double>(dimension-1));
  vector<double> newWeights(trialNumber,0.0);
  double acceptance = 1.0;

  for (int i = 0; i < segmentNumber; i++) {
    // new configuration

    vector<vector<double>> newSegments;
    if (i == 0) newSegments = generateBendingSegments(trialNumber);
    else newSegments = generateBendingSegments(newConfiguration[i-1],trialNumber);

    double newWeightSum = 0.0;
    for (int j = 0; j < trialNumber; j++) {
      newWeights[j] = electricProbability(newSegments[j]);
      newWeightSum += newWeights[j];
    }

    // pick new segment

    discrete_distribution<int> discrete(newWeights.begin(),newWeights.end());
    newConfiguration[i] = newSegments[discrete(rng)];

    // old configuration

    vector<vector<double>> oldSegments;
    if (i == 0) oldSegments = generateBendingSegments(trialNumber-1);
    else oldSegments = generateBendingSegments(oldConfiguration[i-1],trialNumber-1);

    double oldWeightSum = electricProbability(oldConfiguration[i]);
    for (int j = 0; j < trialNumber-1; j++) 
      oldWeightSum += electricProbability(oldSegments[j]);

    acceptance *= newWeightSum / oldWeightSum;
  }

  // accept or discard new configuration

  if (acceptance > 1.0) acceptance = 1.0;
  discrete_distribution<int> accept({1.0-acceptance,acceptance});
  if (accept(rng)) {
    oldConfiguration = newConfiguration;
    return true;
  }
  else return false;
}

vector<vector<double>> generateBendingSegments (int k) {
  uniform_real_distribution<double> pi(0.0,M_PI);
  uniform_real_distribution<double> twoPiCentred(-M_PI,M_PI);

  vector<vector<double>> segments(k,vector<double>(dimension-1));
  
  // assume angle relative to electric field does not exceed pi
  
  if (dimension == 2)
    for (int i = 0; i < k; i++)
      segments[i][0] = twoPiCentred(rng);
  else if (dimension == 3)
    for (int i = 0; i < k; i++) {
      segments[i][0] = twoPiCentred(rng);
      segments[i][1] = pi(rng);
    }

  return segments;
}

vector<vector<double>> generateBendingSegments (vector<double> previousSegment, int k) {
  double a = 0.5 * beta * stiffness / (segmentLength * (dimension-1));
  double sigma = sqrt(0.5/a);

  normal_distribution<double> gaussian(0.0,sigma);

  vector<vector<double>> segments(k,vector<double>(dimension-1));

  if (dimension == 2)
    for (int i = 0; i < k; i++) {
      double relativeAngle = gaussian(rng);
      segments[i][0] = previousSegment[0] + relativeAngle;
    }
  else if (dimension == 3) {
  }

  return segments;
}

double bendingEnergy (vector<double> segment1, vector<double> segment2) {
  double relativeAngle;
  
  if (dimension == 2) {
    double theta1 = segment1[0];
    double theta2 = segment2[0];
    relativeAngle = theta2 - theta1;
  }
  else if (dimension == 3) {
    double theta1 = segment1[0];
    double phi1 = segment1[1];
    double theta2 = segment2[0];
    double phi2 = segment2[1];
    double c1 = cos(theta1);
    double s1 = sin(theta1);
    double c2 = cos(theta2);
    double s2 = sin(theta2);

    // Assume relative angle does not exceed pi
    relativeAngle = acos(c1*c2 + cos(phi1-phi2)*s1*s2);
  }

  return 0.5 * stiffness / segmentLength * pow(relativeAngle,2);
}

double electricEnergy (vector<double> segment) {
  if (energyMode == "cosine") 
    return -susceptibility * segmentLength * pow(cos(segment[0])*field,order);
  else if (energyMode == "square")
    return -b * segmentLength * (1.0 - 0.5*order*pow(segment[0],2));
}

double bendingProbability (vector<double> segment1, vector<double> segment2) {
  double energy = bendingEnergy(segment1,segment2);
  return exp(-beta*energy);
}

double electricProbability (vector<double> segment) {
  double energy = electricEnergy(segment);
  return exp(-beta*energy);
}
