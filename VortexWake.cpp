#include "FileIO.h"
#include "Util.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
using namespace std;

int main() {
  string filename("Flow0010.00000.plt");
  std::vector<std::string> variables = {"x", "y", "u", "v", "W_z"};
  std::vector<int> N = {985, 4497, 1};
  std::vector<std::vector<double>> data;
  int isdouble;
  std::map<int, int> vm;
  InputTec360_binary(filename, N, data, isdouble);

  OutputTec360_binary("field.plt", variables, N, data, isdouble);

  return 0;
}