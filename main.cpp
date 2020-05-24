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

  ofstream xyzFile;

  vector<vector<double>> configuration(segmentNumber,vector<double>(dimension-1,0.0));
  if (rank == 0 && dumpXYZ) {
    xyzFile.open(xyzFileName);
    if (endCriterion == "steps") writeXYZ(configuration,0,xyzFile);
  }
  
  long int acceptedSteps = 0;
  long int totalSteps = 0;
  double energyAverage = 0.0;
  vector<double> cosAvg(segmentNumber);
  vector<double> cosSquaredAvg(segmentNumber);
  vector<double> correlation(segmentNumber);
  vector<double> msa(segmentNumber);
  
  double energyTemp = energy(configuration);
  vector<double> cosAvgTemp = cosAverage(configuration);
  vector<double> correlationTemp = cosCorrelation(configuration);
  vector<double> msaTemp = meanSquaredAngle(configuration);
  vector<double> cosSquaredAvgTemp(segmentNumber);
  for (int i = 0; i < segmentNumber; i++)
    cosSquaredAvgTemp[i] = cosAvgTemp[i] * cosAvgTemp[i];

  // split loop over MPI procs

  long int end = steps / nproc;
  if (rank < steps % nproc) end++;
  if (endCriterion == "steps")
    for (long int i = 0; i < end; i++) {
      bool accepted = false;
      if (rosenbluthMode == "bending" && rosenbluthBending(configuration)){
        acceptedSteps++;
        accepted = true;
      }
      if (dumpXYZ && rank == 0) writeXYZ(configuration,i+1,xyzFile);
      if (analyse) {
        if (accepted) {
          cosAvgTemp = cosAverage(configuration);
          correlationTemp = cosCorrelation(configuration);
          msaTemp = meanSquaredAngle(configuration);
          energyTemp = energy(configuration);
          for (int j = 0; j < segmentNumber; j++)
            cosSquaredAvgTemp[j] = cosAvgTemp[j] * cosAvgTemp[j];
        }
        cosAvg = cosAvg + cosAvgTemp;
        cosSquaredAvg = cosSquaredAvg + cosSquaredAvgTemp;
        correlation = correlation + correlationTemp;
        msa = msa + msaTemp;
        energyAverage += energyTemp;
      }
      if (rank == 0) {
        cout << fixed << setprecision(2);
        cout << "Progress: " << double(i)/end*100.0 << "%\r" << flush;
      }
    }
  else
    for (; acceptedSteps < end; totalSteps++) {
      bool accepted = false;
      if (rosenbluthMode == "bending" && rosenbluthBending(configuration)){
        acceptedSteps++;
        accepted = true;
      }
      if (dumpXYZ && rank == 0 && accepted) writeXYZ(configuration,acceptedSteps,xyzFile);
      if (analyse) {
        if (accepted) {
          cosAvgTemp = cosAverage(configuration);
          correlationTemp = cosCorrelation(configuration);
          msaTemp = meanSquaredAngle(configuration);
          energyTemp = energy(configuration);
          for (int j = 0; j < segmentNumber; j++)
            cosSquaredAvgTemp[j] = cosAvgTemp[j] * cosAvgTemp[j];
        }
        cosAvg = cosAvg + cosAvgTemp;
        cosSquaredAvg = cosSquaredAvg + cosSquaredAvgTemp;
        correlation = correlation + correlationTemp;
        msa = msa + msaTemp;
        energyAverage += energyTemp;
      }
      if (rank == 0) {
        cout << fixed << setprecision(2);
        cout << "Progress: " << double(acceptedSteps)/end*100.0 << "%\r" << flush;
      }
    }
    
  if (rank == 0) cout << "Progress: " << 100.0  << "%" << endl;

  // combine results

  long int masterAcceptedSteps;
  long int masterTotalSteps;
  double masterEnergyAverage;
  vector<double> masterCosAvg(segmentNumber);
  vector<double> masterCosSquaredAvg(segmentNumber);
  vector<double> masterCorrelation(segmentNumber);
  vector<double> masterMsa(segmentNumber);
  MPI_Reduce(&acceptedSteps,&masterAcceptedSteps,1,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&totalSteps,&masterTotalSteps,1,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&energyAverage,&masterEnergyAverage,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(cosAvg.data(),masterCosAvg.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(cosSquaredAvg.data(),masterCosSquaredAvg.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(correlation.data(),masterCorrelation.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(msa.data(),masterMsa.data(),segmentNumber,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (endCriterion == "steps") masterTotalSteps = steps;
  
  if (analyse && rank == 0) {
    masterEnergyAverage /= masterTotalSteps;
    masterCosAvg = (1.0/masterTotalSteps) * masterCosAvg;
    masterCosSquaredAvg = (1.0/masterTotalSteps) * masterCosSquaredAvg;
    masterCorrelation = (1.0/masterTotalSteps) * masterCorrelation;
    masterMsa = (1.0/masterTotalSteps) * masterMsa;

    double orderParameter = 0.0;
    for (int i = 0; i < segmentNumber; i++)
      orderParameter += masterCosSquaredAvg[i];
    orderParameter /= segmentNumber;
    if (dimension == 2) orderParameter = 2.0*orderParameter - 1.0;
    else if (dimension == 3) orderParameter = 0.5*(3.0*orderParameter - 1.0);

    cout << fixed << setprecision(4);
    cout << "Accepted steps: " << masterAcceptedSteps << "/" << masterTotalSteps 
         << " = " << double(masterAcceptedSteps)/double(masterTotalSteps) << endl;
    cout << "Energy average: " << masterEnergyAverage << endl;
    cout << "Order parameter: " << orderParameter << endl;
    
    ofstream cosFile(cosFileName);
    ofstream cosSquaredFile(cosSquaredFileName);
    ofstream correlationFile(correlationFileName);
    ofstream msaFile(msaFileName);
    ofstream cosTheoryFile(cosTheoryFileName);
    ofstream cosSquaredTheoryFile(cosSquaredTheoryFileName);
    ofstream correlationTheoryFile(correlationTheoryFileName);
    ofstream msaTheoryFile(msaTheoryFileName);
    
    writeFile(masterCosAvg,cosFile);
    writeFile(masterCosSquaredAvg,cosSquaredFile);
    writeFile(masterCorrelation,correlationFile);
    writeFile(masterMsa,msaFile);
    writeTheoryFile(cosAverageTheory,cosTheoryFile);
    writeTheoryFile(cosSquaredAverageTheory,cosSquaredTheoryFile);
    writeTheoryFile(cosCorrelationTheory,correlationTheoryFile);
    writeTheoryFile(meanSquaredAngleTheory,msaTheoryFile);
    
    cosFile.close();
    cosSquaredFile.close();
    correlationFile.close();
    msaFile.close();
    cosTheoryFile.close();
    cosSquaredTheoryFile.close();
    correlationTheoryFile.close();
    msaTheoryFile.close();
    xyzFile.close();
  }

  MPI_Finalize();

  return 0;
}
