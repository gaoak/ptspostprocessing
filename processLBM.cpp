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

std::string LoadLBMFields(std::string lbmfilename, IncFlow &baseflow) {
    std::vector<std::string> variables = {"x", "y", "p", "u", "v", "W_z"};
    std::vector<std::vector<int>> Ns = {{985, 1, 4497}, {101, 1, 1}, {101, 1, 1}};
    LBMData lbmfile(lbmfilename, Ns);
    std::vector<std::vector<double>> u1;
    lbmfile.Interpolation(baseflow, u1);
    baseflow.AddPhysics(variables[3], u1[3]);
    baseflow.AddPhysics(variables[4], u1[4]);
    baseflow.CalculateVorticity();
    baseflow.CalculateForcePartition2D();
    //baseflow.OutputData("combinedfile.plt");

    std::vector<double> Qfield    = baseflow.GetPhys(baseflow.GetPhysID("Q"));
    std::vector<double> Phi0field = baseflow.GetPhys(baseflow.GetPhysID("phi0"));
    for(size_t i=0; i<Qfield.size(); ++i) {
        Qfield[i] *= Phi0field[i];
    }
    double volumeforce = 2.0 * baseflow.Integrate(Qfield);
    char buffer[1000];
    sprintf(buffer, "%25.14f ", volumeforce);
    return std::string(buffer);
}

std::string LoadBoundaryFields(std::string bndfilename0, std::string bndfilename1, IncFlow &baseflow, double nu) {
    LineData boundfield, tmpfield;
    std::vector<std::string> bndvars{"x", "y", "nx", "ny", "ax", "ay"};
    boundfield.InputData(bndfilename0, bndvars, true);
    //add phi fields
    std::vector<std::string>  tmpvars{"x", "y", "z", "phi0", "phi1", "phi2"};
    tmpfield.InputData(bndfilename1, tmpvars, true);
    std::vector<double> phi0 = tmpfield.GetPhys(tmpfield.GetPhysID("phi0"));
    //add flow info
    std::map<int,double> addfield;
    addfield[baseflow.GetPhysID("Dxx")] = 0.;
    addfield[baseflow.GetPhysID("Dxy")] = 0.;
    addfield[baseflow.GetPhysID("Dyy")] = 0.;
    addfield[baseflow.GetPhysID("uLaplace")] = 0.;
    addfield[baseflow.GetPhysID("vLaplace")] = 0.;
    boundfield.InterpolateFrom(baseflow, addfield);
    //boundfield.OutputData("test_" + bndfilename0);
    // do integration
    std::vector<double> nx = boundfield.GetPhys(boundfield.GetPhysID("nx"));
    std::vector<double> ny = boundfield.GetPhys(boundfield.GetPhysID("ny"));
    std::vector<double> ax = boundfield.GetPhys(boundfield.GetPhysID("ax"));
    std::vector<double> ay = boundfield.GetPhys(boundfield.GetPhysID("ay"));
    std::vector<double> Dxx = boundfield.GetPhys(boundfield.GetPhysID("Dxx"));
    std::vector<double> Dxy = boundfield.GetPhys(boundfield.GetPhysID("Dxy"));
    std::vector<double> uLaplace = boundfield.GetPhys(boundfield.GetPhysID("uLaplace"));
    std::vector<double> vLaplace = boundfield.GetPhys(boundfield.GetPhysID("vLaplace"));
    std::vector<double> acceleration(phi0.size(), 0);
    std::vector<double> friction(phi0.size(), 0);
    std::vector<double> diffusion(phi0.size(), 0);
    for(size_t i=0; i<phi0.size(); ++i) {
        acceleration[i] = -(ax[i] * nx[i] + ay[i] * ny[i]) * phi0[i];
        friction[i] = 2. * nu * (Dxx[i] * nx[i] + Dxy[i] * ny[i]);
        diffusion[i] = nu * (uLaplace[i] * nx[i] + vLaplace[i] * ny[i]) * phi0[i];
    }
    double acceforce = boundfield.Integrate(acceleration);
    double fricforce = boundfield.Integrate(friction);
    double diffforce = boundfield.Integrate(diffusion);
    char buffer[1000];
    sprintf(buffer, "%25.14f %25.14f %25.14f", fricforce, acceforce, diffforce);
    return std::string(buffer);
}

int main(int argc, char* argv[]) {
  string phifilename, lbmfilename, bndfilename0, bndfilename1;
  if(argc<5) {
    cout << "4 input files requred: phi field[uniform grid, integration field], LBM field, boundary line with x y a and n, boundary line with phi" << std::endl;
    return argc - 5;
  }
  phifilename = argv[1];
  lbmfilename = argv[2];
  bndfilename0 = argv[3]; // x, y, nx, ny, ax, ay
  bndfilename1 = argv[4]; // x, y, z, phi0
  // load uniform grid of phi file
  IncFlow baseflow;
  baseflow.InputData(phifilename);
  // work flow 1, load LBM file and compute Q, Dxx, Dxy, Dyy, Laplace u, Laplace v
  std::string volumeres = LoadLBMFields(lbmfilename, baseflow);
  // work flow 2, load boundary file x, y, nx, ny, ax, ay
  double nu = 1/200.;
  std::string bndres = LoadBoundaryFields(bndfilename0, bndfilename1, baseflow, nu);
  
  printf("Volume force; friction force; acceleration force; viscous pressure force\n");
  printf("RESULT %s %s\n", volumeres.c_str(), bndres.c_str());
  return 0;
}