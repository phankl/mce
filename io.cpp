#include "io.h"

void writeXYZ (vector<vector<double>> configuration, int step, ofstream &file) {
  file << segmentNumber + 1 << endl;
  file << "Atoms. Timestep: " << step << endl;

  vector<double> position(dimension);
  vector<double> segment(dimension);

  // write initial position at origin

  file << "1 ";
  for (int i = 0; i < dimension; i++) file << "0.0 ";
  if (dimension == 2) file << "0.0";
  file << endl;

  // compute all subsequent positions and write to file

  for (int i = 0; i < segmentNumber; i++) {
    if (dimension == 2) {
      double theta = configuration[i][0];
      segment[0] = sin(theta);
      segment[1] = cos(theta);
    }
    else if (dimension == 2) {
      double theta = configuration[i][0];
      double phi = configuration[i][1];
      segment[0] = sin(theta) * cos(phi);
      segment[1] = sin(theta) * sin(phi);
      segment[2] = cos(theta);
    }
    position = position + segment;  
    
    file << "1 ";
    for (int j = 0; j < dimension; j++) file << position[j] << " ";
    if (dimension == 2) file << "0.0";
    file << endl;
  }
}

void writeFile (vector<double> data, ofstream &file) {
  for (int i = 0; i < segmentNumber; i++) {
    double s = i * segmentLength;
    file << s << " " << data[i] << endl;
  }
}

void writeTheoryFile (function<double(double)> theory, ofstream &file) {
  for (int i = 0; i < segmentNumber; i++) {
    double s = (double(i) + 0.5) * segmentLength;
    file << s << " " << theory(s) << endl;
  }
}
