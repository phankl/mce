#include "mc.h"

bool rosenbluthBending (vector<vector<double>> &oldConfiguration) {
  vector<vector<double>> newConfiguration(segmentNumber,vector<double>(dimension-1));
  vector<double> newWeights(trialNumber,0.0);

  newConfiguration[0] = generateBendingSegments();
  double acceptance = electricProbability(newConfiguration[0]) / electricProbability(oldConfiguration[0]);
  for (int i = 1; i < segmentNumber; i++) {
    // new configuration

    vector<vector<double>> newSegments;
    newSegments = generateBendingSegments(newConfiguration[i-1],trialNumber);

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
    oldSegments = generateBendingSegments(oldConfiguration[i-1],trialNumber-1);

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

vector<double> generateBendingSegments () {
  uniform_real_distribution<double> piCentred(-0.5*M_PI,0.5*M_PI);
  uniform_real_distribution<double> twoPiCentred(-M_PI,M_PI);
  uniform_real_distribution<double> uniform(0.0,1.0);
  uniform_real_distribution<double> twoPi(0.0,2.0*M_PI);

  if (field == 0.0) {
    vector<double> segment(dimension-1,0.0);
    return segment;
  }
  else if (energyMode == "square") {
    double sigma = 1.0 / sqrt(beta*order*susceptibility*field*segmentLength);
    normal_distribution<double> gaussian(0.0,sigma);
    vector<double> segment(dimension-1);
    for (int i = 0; i < dimension-1; i++) segment[i] = gaussian(rng);
    return segment;
  }
  else if (energyMode == "cosine") {
    //if (dimension == 2) return {randomGaussianCosine(segmentLength*beta*b)};
    if (dimension == 2) return {twoPiCentred(rng)};
    //else if (dimension == 3) return {randomSineGaussianCosine(segmentLength*beta*b),piCentred(rng)};
    else if (dimension == 3) return {acos(2.0*uniform(rng)-1),twoPi(rng)};
  }
}

vector<vector<double>> generateBendingSegments (vector<double> previousSegment, int k) {
  double a = 0.5 * beta * stiffness / segmentLength;
  double sigma = sqrt(0.5/a);

  normal_distribution<double> gaussian(0.0,sigma);
  uniform_real_distribution<double> piCentred(-0.5*M_PI,0.5*M_PI);

  vector<vector<double>> segments(k,vector<double>(dimension-1));

  if (energyMode == "square") {
    for (int i = 0; i < k; i++)
      for (int j = 0; j < dimension-1; j++)
        segments[i][j] = previousSegment[j] + gaussian(rng);
  }
  else if (energyMode == "cosine") {
    if (dimension == 2) 
      for (int i = 0; i < k; i++) segments[i][0] = previousSegment[0] + gaussian(rng);
    else if (dimension == 3) {
      double previousTheta = previousSegment[0];
      double previousPhi = previousSegment[1];
      vector<double> previousVector({
        sin(previousTheta)*cos(previousPhi),
        sin(previousTheta)*sin(previousPhi),
        cos(previousTheta)
      });
      vector<double> zAxis({0.0,0.0,1.0});
      vector<double> rotationAxis = cross(previousVector,zAxis);
      normalise(rotationAxis);
      double rotationAngle = acos(previousVector[2]);

      int previousWindings;
      if (previousTheta < 0) previousWindings = ceil(previousTheta/M_PI);
      else if (previousTheta > 0) previousWindings = floor(previousTheta/M_PI);
      else previousWindings = 0;

      for (int i = 0; i < k; i++) {
        double relativeTheta = randomSineGaussian(a);
        double relativePhi = piCentred(rng);
        vector<double> newVector({
          sin(relativeTheta)*cos(relativePhi),
          sin(relativeTheta)*sin(relativePhi),
          cos(relativeTheta)
        });
        if (previousTheta != 0.0) newVector = rotate(newVector,rotationAxis,rotationAngle);
       
        double newTheta = acos(newVector[2]);
        if (newVector[0] < 0) newTheta *= -1.0;
        //newTheta += previousWindings*M_PI;
        double newPhi = atan(newVector[1]/newVector[0]);
        
        segments[i][0] = newTheta;
        segments[i][1] = newPhi;
      }
    }
  }

  return segments;
}

double bendingEnergy (vector<double> segment1, vector<double> segment2) { 
  if (energyMode == "square") {
    double thetaSquared = 0.0;
    for (int i = 0; i < dimension-1; i++) 
      thetaSquared += pow(segment1[i]-segment2[i],2);
    return 0.5 * stiffness / segmentLength * thetaSquared;
  }
  else if (energyMode == "cosine")
    if (dimension == 2) {
      double relativeAngle = segment1[0] - segment2[0];
      return 0.5 * stiffness / segmentLength * pow(relativeAngle,2);
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
      double relativeAngle = acos(c1*c2 + cos(phi1-phi2)*s1*s2);
      return 0.5 * stiffness / segmentLength * pow(relativeAngle,2);
    }
}

double electricEnergy (vector<double> segment) {
  if (energyMode == "square") {
    double thetaSquared = 0.0;
    for (int i = 0; i < dimension-1; i++)
      thetaSquared += pow(segment[i],2);
    return -b * segmentLength * (1.0 - 0.5*order*thetaSquared);
  }
  else if (energyMode == "cosine") 
    return -susceptibility * segmentLength * pow(cos(segment[0])*field,order);
}

double bendingProbability (vector<double> segment1, vector<double> segment2) {
  double energy = bendingEnergy(segment1,segment2);
  return exp(-beta*energy);
}

double electricProbability (vector<double> segment) {
  double energy = electricEnergy(segment);
  return exp(-beta*energy);
}
