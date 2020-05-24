#ifndef header_constants
#define header_constants

#include <fstream>
#include <limits>
#include <random>
#include <string>

using namespace std;

extern const int dimension;
extern const long int steps;

extern const int segmentNumber;
extern const double chainLength;
extern const double segmentLength;

extern const int order;
extern const double stiffness;
extern const double susceptibility;
extern const double field;
extern const double beta;

extern const double b;
extern const double alpha;

extern const int binNumber;
extern const double newtonDelta;
extern const double newtonEpsilon;

extern vector<double> gaussianCosineValues;
extern double gaussianCosineBound;

extern const int trialNumber;

extern mt19937_64 rng;

extern const string rosenbluthMode;
extern const string energyMode;
extern const string endCriterion;

extern const bool dumpXYZ;
extern const bool analyse;

extern const string xyzFileName;
extern const string cosFileName;
extern const string cosSquaredFileName;
extern const string msaFileName;
extern const string cosTheoryFileName;
extern const string cosSquaredTheoryFileName;
extern const string msaTheoryFileName;
extern const string correlationFileName;
extern const string correlationTheoryFileName;

#endif
