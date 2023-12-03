#include "FileIO.h"
#include "Util.h"
#include "IncFlow.h"
#include "LineData.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {
  string filenameformat;
  if(argc>1) {
    filenameformat = argv[1];
  }
  IncFlow fieldData;
  fieldData.InputData("Points_R00000_T000.000000.plt");
  std::vector<double> ux = fieldData.GetPhys(fieldData.GetPhysID("u_x"));
  std::vector<double> uy = fieldData.GetPhys(fieldData.GetPhysID("u_y"));
  std::vector<double> uz = fieldData.GetPhys(fieldData.GetPhysID("u_z"));
  std::vector<double> vx = fieldData.GetPhys(fieldData.GetPhysID("v_x"));
  std::vector<double> vy = fieldData.GetPhys(fieldData.GetPhysID("v_y"));
  std::vector<double> vz = fieldData.GetPhys(fieldData.GetPhysID("v_z"));
  std::vector<double> wx = fieldData.GetPhys(fieldData.GetPhysID("w_x"));
  std::vector<double> wy = fieldData.GetPhys(fieldData.GetPhysID("w_y"));
  std::vector<double> wz = fieldData.GetPhys(fieldData.GetPhysID("w_z"));
  std::vector<double> Lu = fieldData.GetPhys(fieldData.GetPhysID("Lu"));
  std::vector<double> Lu = fieldData.GetPhys(fieldData.GetPhysID("Lv"));
  std::vector<double> Lu = fieldData.GetPhys(fieldData.GetPhysID("Lw"));
  std::vector<double> vorz(ux.size(), 0);
  for(size_t i=0; i<ux.size(); ++i) {
    vorz[i] = ;
  }
  return 0;
}