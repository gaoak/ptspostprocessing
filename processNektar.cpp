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

std::string LoadPhiFields(std::string phifilename, IncFlow &baseflow) {
    IncFlow phifile;
    phifile.InputData(phifilename);
    for(int i=0; i<phifile.GetNumPhys();++i) {
      baseflow.AddPhysics(phifile.GetPhysVarName(i), phifile.GetPhys(i));
    }
    if(baseflow.GetPhysID("W_z")==-1) {
      cout << "compute vorticity." << endl;
      baseflow.CalculateVorticity();
      baseflow.OutputData("combinedfile.plt");
    }

    std::vector<double> W_zfield    = baseflow.GetPhys(baseflow.GetPhysID("W_z"));
    std::vector<double> ufield    = baseflow.GetPhys(baseflow.GetPhysID("u"));
    std::vector<double> vfield    = baseflow.GetPhys(baseflow.GetPhysID("v"));
    std::vector<double> Phi1xfield = baseflow.GetPhys(baseflow.GetPhysID("phi1_x"));
    std::vector<double> Phi1yfield = baseflow.GetPhys(baseflow.GetPhysID("phi1_y"));

    std::vector<double> integrand(W_zfield.size());

    for(size_t i=0; i<W_zfield.size(); ++i) {
        integrand[i] = W_zfield[i]*(-vfield[i]*Phi1xfield[i] + ufield[i]*Phi1yfield[i]);
    }
    double volumeforce = baseflow.Integrate(integrand);
    char buffer[1000];
    sprintf(buffer, "%25.14f ", volumeforce);
    return std::string(buffer);
}

int main(int argc, char* argv[]) {
  string flowfilename, phifilename, bndfilename0, bndfilename1;
  if(argc<5) {
    cout << "4 input files requred: phi field[uniform grid, integration field], LBM field, boundary line with x y a and n, boundary line with phi" << std::endl;
    return argc - 5;
  }
  flowfilename = argv[1];
  phifilename = argv[2];
  bndfilename0 = argv[3]; // x, y, nx, ny, ax, ay
  bndfilename1 = argv[4]; // x, y, z, phi0
  // load uniform grid of phi file
  IncFlow baseflow;
  baseflow.InputData(flowfilename);
  // work flow 1, load LBM file and compute Q, Dxx, Dxy, Dyy, Laplace u, Laplace v
  std::string volumeres = LoadPhiFields(phifilename, baseflow);
  
  printf("Volume force; friction force; acceleration force; viscous pressure force\n");
  printf("RESULT %s\n", volumeres.c_str());
  return 0;
}