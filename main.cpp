#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>

#include "analysis.h"
#include "constants.h"
#include "io.h"
#include "math_extra.h"
#include "mc.h"

using namespace std;

int main (int argc, char* argv[]) {

  // MPI setup
  int rank;
  int nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  initRandomGaussianCosine(beta*segmentLength*b);

  /*
  int localBinNumber = 1000;
  int samples = 100000000;
  double a = beta*segmentLength*b;
  double binWidth = 2.0*M_PI / localBinNumber;
  double cosine = 0.0;
  vector<int> bins(localBinNumber);
  for (int i = 0; i < samples; i++) {
    double x = randomSineGaussianCosine(a);
    int bin = floor((x+M_PI)/binWidth);
    cosine += cos(x);
    bins[bin]++;
  }

  ofstream random("random.dat");
  for (int i = 0; i < localBinNumber; i++)
    random << (double(i) + 0.5)*binWidth - M_PI << " " << double(bins[i])/samples/binWidth << endl;
  random.close();
  cout << cosine / samples << endl;
  */
  /*
  long int testNumber = 1000000;
  double testTheta = -0.1*M_PI;
  double testPhi = -1.3*M_PI;
  vector<double> testSegment({testTheta,testPhi});
  vector<vector<double>> generatedSegments = generateBendingSegments(testSegment,testNumber);
  vector<double> mean(3,0.0);
  for (long int i = 0; i < testNumber; i++) {
    double theta = generatedSegments[i][0];
    double phi = generatedSegments[i][1];
    vector<double> position({
      sin(theta)*cos(phi),
      sin(theta)*sin(phi),
      cos(theta)
    });
    mean = mean + position;
  }
  mean = (1.0/testNumber) * mean;
  cout << "Previous segment: " << sin(testTheta)*cos(testPhi) << " " << sin(testTheta)*sin(testPhi) << " " << cos(testTheta) << endl;
  cout << "Mean: " << mean[0] << " " << mean[1] << " " << mean[2] << endl;
  */
  ofstream xyzFile;

  vector<vector<double>> configuration(segmentNumber,vector<double>(dimension-1,0.0));
  if (rank == 0 && dumpXYZ) {
    xyzFile.open(xyzFileName);
    writeXYZ(configuration,0,xyzFile);
  }
  
  long int acceptedSteps = 0;
  double energyAverage = 0.0;
  vector<double> cosAvg(segmentNumber);
  vector<double> correlation(segmentNumber);
  vector<double> msa(segmentNumber);

  // split loop over MPI procs

  long int end = steps / nproc;
  if (rank < steps % nproc) end++;
  for (long int i = 0; i < end; i++) {
    if (rosenbluthMode == "bending" && rosenbluthBending(configuration)) acceptedSteps++;
    if (dumpXYZ && rank == 0) writeXYZ(configuration,i+1,xyzFile);
    if (analyse) {
      for (int j = 0; j < segmentNumber; j++) {
        //cout << configuration[j][0] << endl;
      }
      cosAvg = cosAvg + cosAverage(configuration);
      correlation = correlation + cosCorrelation(configuration);
      msa = msa + meanSquaredAngle(configuration);
      energyAverage += energy(configuration);
    }
    if (rank == 0) {
      cout << fixed << setprecision(2);
      cout << "Progress: " << double(i)/end*100.0 << "%\r" << flush;
    }
  }
  if (rank == 0) cout << "Progress: " << 100.0  << "%" << endl;

  // combine results

  long int masterAcceptedSteps;
  double masterEnergyAverage;
  vector<double> masterCosAvg(segmentNumber);
  vector<double> masterCorrelation(segmentNumber);
  vector<double> masterMsa(segmentNumber);
  MPI_Reduce(&acceptedSteps,&masterAcceptedSteps,1,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&energyAverage,&masterEnergyAverage,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(cosAvg.data(),masterCosAvg.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(correlation.data(),masterCorrelation.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(msa.data(),masterMsa.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (analyse && rank == 0) {
    masterEnergyAverage /= steps;
    masterCosAvg = (1.0/steps) * masterCosAvg;
    masterCorrelation = (1.0/steps) * masterCorrelation;
    masterMsa = (1.0/steps) * masterMsa;

    cout << "Accepted steps: " << masterAcceptedSteps << "/" << steps << endl;
    cout << "Energy average: " << masterEnergyAverage << endl;  
    
    ofstream cosFile(cosFileName);
    ofstream correlationFile(correlationFileName);
    ofstream msaFile(msaFileName);
    ofstream cosTheoryFile(cosTheoryFileName);
    ofstream correlationTheoryFile(correlationTheoryFileName);
    ofstream msaTheoryFile(msaTheoryFileName);
    
    writeFile(masterCosAvg,cosFile);
    writeFile(masterCorrelation,correlationFile);
    writeFile(masterMsa,msaFile);
    writeTheoryFile(cosAverageTheory,cosTheoryFile);
    writeTheoryFile(cosCorrelationTheory,correlationTheoryFile);
    writeTheoryFile(meanSquaredAngleTheory,msaTheoryFile);
    
    cosFile.close();
    correlationFile.close();
    msaFile.close();
    cosTheoryFile.close();
    correlationTheoryFile.close();
    msaTheoryFile.close();
    xyzFile.close();
  }

  MPI_Finalize();

  return 0;
}
