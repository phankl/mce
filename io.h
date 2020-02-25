#ifndef io_header
#define io_header

#include <fstream>
#include <functional>
#include <vector>

#include "constants.h"
#include "math_extra.h"

using namespace std;

void writeXYZ(vector<vector<double>>, int, ofstream&);
void writeFile(vector<double>, ofstream&);
void writeTheoryFile(function<double(double)>, ofstream&);

#endif
