#include "FileIO.h"
#include "Util.h"
#include "LBMData.h"
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

void processField(IncFlow &baseflow) {
  vector<int> N = baseflow.GetN();
  int Np = N[0] * N[1] * N[2];
  vector<double> ux = baseflow.GetPhys(baseflow.GetPhysID("u"));
  vector<double> u(Np, 0.);
  for(int i=0; i<Np; ++i) u[i] = ux[i] - 1.;
  vector<double> v = baseflow.GetPhys(baseflow.GetPhysID("v"));
  vector<double> p = baseflow.GetPhys(baseflow.GetPhysID("p"));
  vector<double> Txx = baseflow.GetPhys(baseflow.GetPhysID("Txx"));
  vector<double> Txy = baseflow.GetPhys(baseflow.GetPhysID("Txy"));
  vector<double> Tyy = baseflow.GetPhys(baseflow.GetPhysID("Tyy"));
  vector<double> kux(Np, 0.), pu(Np, 0.), tauum(Np, 0.);
  for(int i=0; i<Np; ++i) {
    kux[i] = 0.5 * (u[i] * u[i] + v[i] * v[i]) * ux[i];
    pu[i] = p[i] * u[i];
    tauum[i] = - (Txx[i]*u[i] + Txy[i]*v[i]);
  }
  // integrate
  vector<double> intkux(N[0], 0.);
  vector<double> intpu(N[0], 0.);
  vector<double> inttauum(N[0], 0.);
  double dy = baseflow.GetDetx()[1];
  for(int i=0; i<N[0]; ++i) {
    for(int j=0; j<N[1]; ++j) {
      int ind = Index(N,{i,j,0});
      intkux[i] += kux[ind];
      intpu[i] += pu[ind];
      inttauum[i] += tauum[ind];
    }
  }
  for(int i=0; i<N[0]; ++i) {
    intkux[i] *= dy;
    intpu[i] *= dy;
    inttauum[i] *= dy;
  }
  ofstream ofile("surface.dat");
  ofile << "variables = x, kflux, press, visc\n";
  for(int i=0; i<N[0]; ++i) {
    ofile << baseflow.GetCoordValue(0, i) << " " << intkux[i] << " " << intpu[i] << " " << inttauum[i] << "\n";
  }
  ofile.close();
}

int main(int argc, char* argv[]) {
  string flowfilename;
  if(argc<2) {
    cout << "1 input files requred: flow field" << std::endl;
    return argc - 5;
  }
  flowfilename = argv[1];
  IncFlow baseflow;
  baseflow.InputData(flowfilename);
  printf("processing file %s\n", argv[1]);
  processField(baseflow);
  printf("Finished.\n");
  return 0;
}