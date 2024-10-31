#include "FileIO.h"
#include "Util.h"
#include "IncFlow.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  string inputformat("interp%d.plt");
  int startn = 0, endn=0;
  if (argc>2) {
    startn = std::stoi(argv[1]);
    endn   = std::stoi(argv[2]);
  }
  IncFlow flow;
  int ncount = 0;
  std::vector<vector<double>> PhysMean(4);
  for(int nf = startn; nf <= endn; ++nf) {
    char buffer[100];
    sprintf(buffer, inputformat.c_str(), nf);
    string filename(buffer);
    flow.InputData(filename);
    std::vector<int> N = flow.GetN();
    std::vector<double> u = flow.GetPhys(flow.GetPhysID("u"));
    std::vector<double> v = flow.GetPhys(flow.GetPhysID("v"));
    std::vector<double> w = flow.GetPhys(flow.GetPhysID("w"));
    std::vector<double> p = flow.GetPhys(flow.GetPhysID("p"));
    std::vector<vector<double>> phys;
    phys.push_back(u);
    phys.push_back(v);
    phys.push_back(w);
    phys.push_back(p);
    ShiftIndex<double>(N, phys, 1);
    if((int)PhysMean[0].size()<N[2]) {
      for(size_t n=0; n<PhysMean.size(); ++n) {
        PhysMean[n].resize(N[2], 0.);
      }
    }
    int nseg = N[0] * N[1];
    for(int i=0; i<N[2]; ++i) {
        int offset = nseg * i;
        for(int j=0; j<nseg; ++j) {
        for(size_t n=0; n<PhysMean.size(); ++n) {
            PhysMean[n][i] += phys[n][offset+j];
        }
        }
    }
    ncount += nseg;
  }
  for(size_t n=0; n<PhysMean.size(); ++n) {
    for(size_t i=0; i<PhysMean[n].size(); ++i) {
      PhysMean[n][i] /= ncount;
    }
  }
  double dy = flow.GetDetx()[1];
  double mu = 1. / 2857.;
  double tau = mu * (PhysMean[0][1] - PhysMean[0][0] ) / dy;
  double utau = sqrt(tau);
  double delta = mu / utau;
  ofstream meanfile("mean.dat");
  meanfile << "variables = yplus, uplus, y, u, v, w, p\n";
  for(size_t i=0; i<PhysMean[0].size(); ++i) {
    double y = i * dy;
    meanfile << y / delta << " " << PhysMean[0][i] / utau << " " << y << " " << PhysMean[0][i] << " " << PhysMean[1][i] << " " << PhysMean[2][i] << " " << PhysMean[3][i] << endl;
  }
  meanfile.close();
  return 0;
}