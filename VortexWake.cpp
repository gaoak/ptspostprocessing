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
  std::vector<std::string> variables = {"x", "y", "u", "v", "p", "W_z"};
  std::vector<DataPack>  zones(3);
  zones[0].N = {985, 4497, 1};
  zones[1].N = {101, 1, 1};
  zones[2].N = {101, 1, 1};
  std::vector<std::vector<double>> dataBody;
  int isdouble = 0;
  std::map<int, int> vm;

  InputTec360_FSILBM2D(filename, zones);

  for(size_t z = 0; z < zones.size(); ++z) {
    OutputTec360_binary("field" + to_string(z) + ".plt", variables, zones[z].N, zones[z].data, isdouble);
  }

  return 0;
}