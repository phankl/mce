#include "constants.h"

const int dimension = 2;
const long int steps = 100000;

const int segmentNumber = 100;
const double chainLength = 1.0;
const double segmentLength = chainLength / segmentNumber;

const int order = 1;
const double stiffness = 1.0;
const double susceptibility = 1.0;
const double field = 1.0;
const double beta = 1.0;

const double b = susceptibility * pow(field, order);
const double alpha = sqrt(b * order / stiffness);

const int trialNumber = 1;

mt19937_64 rng(125555);

const string rosenbluthMode = "bending";
const string energyMode = "square";

const bool dumpXYZ = true;
const bool analyse = true;

const string xyzFileName = "mc.xyz";
const string cosFileName = "cos.dat";
const string msaFileName = "msa.dat";
const string cosTheoryFileName = "cosTheory.dat";
const string msaTheoryFileName = "msaTheory.dat";
const string correlationFileName = "correlation.dat";
const string correlationTheoryFileName = "correlationTheory.dat";
