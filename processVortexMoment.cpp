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

std::string VortexMoment(IncFlow &baseflow) {
    std::vector<double> W_zfield    = baseflow.GetPhys(baseflow.GetPhysID("W_z"));
    std::vector<double> integrandx(W_zfield.size());
    std::vector<double> integrandy(W_zfield.size());
    std::vector<double> integrandx2y2(W_zfield.size());

    for(size_t i=0; i<W_zfield.size(); ++i) {
      double x = baseflow.GetCoordValue(0, i) ;
      double y = baseflow.GetCoordValue(1, i) ;
      if(y<0) {
        W_zfield[i] = 0.;
      }
      integrandx[i] = W_zfield[i]*x;
      integrandy[i] = W_zfield[i]*y;
    }
    double gamma = baseflow.Integrate(W_zfield);
    double centerx = baseflow.Integrate(integrandx) / gamma;
    double centery = baseflow.Integrate(integrandy) / gamma;
    for(size_t i=0; i<W_zfield.size(); ++i) {
      double x = baseflow.GetCoordValue(0, i)  - centerx;
      double y = baseflow.GetCoordValue(1, i) - centery ;
      integrandx2y2[i] = W_zfield[i]*(x*x+y*y);
    }
    double sigma = sqrt(baseflow.Integrate(integrandx2y2) / gamma);
    char buffer[1000];
    sprintf(buffer, "%25.14f %25.14f %25.14f %25.14f ", gamma, centerx, centery, sigma);
    return std::string(buffer);
}

int main(int argc, char* argv[]) {
  string flowfilename;
  if(argc<2) {
    cout << "1 input files requred: flow field with vorticity" << std::endl;
    return argc - 2;
  }
  flowfilename = argv[1];
  // load uniform grid of phi file
  IncFlow baseflow;
  baseflow.InputData(flowfilename);
  string res = VortexMoment(baseflow);
  printf("Volume force; friction force; acceleration force; viscous pressure force\n");
  printf("%s %s\n", flowfilename.c_str(), res.c_str());
  return 0;
}