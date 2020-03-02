#include "constants.h"

const int dimension = 3;
const long int steps = 1000000;

const int segmentNumber = 100;
const double chainLength = 1.0;
const double segmentLength = chainLength / segmentNumber;

const int order = 2;
const double stiffness = 1.0;
const double susceptibility = 1.0;
const double field = 1.0;
const double beta = 1.0;

const double b = susceptibility * pow(field, order);
const double alpha = sqrt(b * order / stiffness);

const int binNumber = 100;
const double newtonDelta = 1.0e-8;
const double newtonEpsilon = 1.0e-8;

vector<double> gaussianCosineValues(binNumber);
double gaussianCosineBound = 0.0;

const int trialNumber = 1;

random_device noise;
uniform_int_distribution<int> seed(0,numeric_limits<int>::max());
mt19937_64 rng(seed(noise));

const string rosenbluthMode = "bending";
const string energyMode = "cosine";

const bool dumpXYZ = false;
const bool analyse = true;

const string xyzFileName = "mc.xyz";
const string cosFileName = "cos.dat";
const string msaFileName = "msa.dat";
const string cosTheoryFileName = "cosTheory.dat";
const string msaTheoryFileName = "msaTheory.dat";
const string correlationFileName = "correlation.dat";
const string correlationTheoryFileName = "correlationTheory.dat";
