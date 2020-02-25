#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "analysis.h"
#include "constants.h"
#include "io.h"
#include "math_extra.h"
#include "mc.h"

using namespace std;

int main (int argc, char* argv[]) {

  ofstream xyzFile(xyzFileName);
  ofstream cosFile(cosFileName);
  ofstream correlationFile(correlationFileName);
  ofstream msaFile(msaFileName);
  ofstream cosTheoryFile(cosTheoryFileName);
  ofstream correlationTheoryFile(correlationTheoryFileName);
  ofstream msaTheoryFile(msaTheoryFileName);

  vector<vector<double>> configuration(segmentNumber,vector<double>(dimension-1,0.0));
  if (dumpXYZ) writeXYZ(configuration,0,xyzFile);

  int acceptedSteps = 0;
  double energyAverage = 0.0;
  vector<double> cosAvg(segmentNumber);
  vector<double> correlation(segmentNumber);
  vector<double> msa(segmentNumber);

  for (long int i = 0; i < steps; i++) {
    if (rosenbluthMode == "bending" && rosenbluthBending(configuration)) acceptedSteps++;
    if (dumpXYZ) writeXYZ(configuration,i+1,xyzFile);
    if (analyse) {
      cosAvg = cosAvg + cosAverage(configuration);
      correlation = correlation + cosCorrelation(configuration);
      msa = msa + meanSquaredAngle(configuration);
      energyAverage += energy(configuration);
    }
  }

  cout << "Accepted steps: " << acceptedSteps << "/" << steps << endl;
  if (analyse) {
    energyAverage /= steps;
    cosAvg = (1.0/steps) * cosAvg;
    correlation = (1.0/steps) * correlation;
    msa = (1.0/steps) * msa;

    cout << "Energy average: " << energyAverage << endl;
    
    writeFile(cosAvg,cosFile);
    writeFile(correlation,correlationFile);
    writeFile(msa,msaFile);
    writeTheoryFile(cosAverageTheory,cosTheoryFile);
    writeTheoryFile(cosCorrelationTheory,correlationTheoryFile);
    writeTheoryFile(meanSquaredAngleTheory,msaTheoryFile);
  }

  xyzFile.close();
  cosFile.close();
  correlationFile.close();
  msaFile.close();
  cosTheoryFile.close();
  correlationTheoryFile.close();
  msaTheoryFile.close();

  return 0;
}
