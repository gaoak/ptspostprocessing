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

void SelectPoint(IncFlow & fieldData, int &maxid) {
  std::vector<double> uy = fieldData.GetPhys(fieldData.GetPhysID("u_y"));
  std::vector<double> vx = fieldData.GetPhys(fieldData.GetPhysID("v_x"));
  std::vector<int> Np = fieldData.GetN();
  maxid = 0;
  double maxv = 0, vorz;
  for(int p1=0; p1<Np[1]; ++p1) {
    for(int p0=0; p0<Np[0]; ++p0) {
      std::vector<int> tmpi = {p0, p1, 1};
      int ind = Index(Np, tmpi);
      vorz = fabs(vx[ind] - uy[ind]);
      if(vorz>maxv) {
        maxv = vorz;
        maxid = ind;
      }
    }
  }
}

void GetVorticityTerms(IncFlow & fieldData, int maxid, std::vector<double> & vorticity, std::vector<double> &mBomega, std::vector<double> & curla) {
  double ux = fieldData.GetPhys(fieldData.GetPhysID("u_x"))[maxid];
  double uy = fieldData.GetPhys(fieldData.GetPhysID("u_y"))[maxid];
  double uz = fieldData.GetPhys(fieldData.GetPhysID("u_z"))[maxid];
  double vx = fieldData.GetPhys(fieldData.GetPhysID("v_x"))[maxid];
  double vy = fieldData.GetPhys(fieldData.GetPhysID("v_y"))[maxid];
  double vz = fieldData.GetPhys(fieldData.GetPhysID("v_z"))[maxid];
  double wx = fieldData.GetPhys(fieldData.GetPhysID("w_x"))[maxid];
  double wy = fieldData.GetPhys(fieldData.GetPhysID("w_y"))[maxid];
  double wz = fieldData.GetPhys(fieldData.GetPhysID("w_z"))[maxid];
  double Lu = fieldData.GetPhys(fieldData.GetPhysID("Lu"))[maxid];
  double Lu = fieldData.GetPhys(fieldData.GetPhysID("Lv"))[maxid];
  double Lu = fieldData.GetPhys(fieldData.GetPhysID("Lw"))[maxid];
  vorticity.resize(3);
  mBomega.resize(3);
  curla.resize(3);
  double nu = 1/800.;
  curla[0] = nu*Lu;
  curla[1] = nu*Lv;
  curla[2] = nu*Lw;
  vorticity[0] = wy - vz;
  vorticity[1] = uz - wx;
  vorticity[2] = vx - uy;
  double theta = ux + vy + wz;
  mBomega[0] = (ux - theta) * vorticity[0] +
               (uy        ) * vorticity[1] +
               (uz        ) * vorticity[2];
  mBomega[1] = (vx        ) * vorticity[0] +
               (vy - theta) * vorticity[1] +
               (vz        ) * vorticity[2];
  mBomega[2] = (wx        ) * vorticity[0] +
               (wy        ) * vorticity[1] +
               (wz - theta) * vorticity[2];
}

void GetCoordsTerms(IncFlow & fieldData, int maxid, std::vector<double> &Xcoords, std::vector<double> &Ycoords, std::vector<double> &Zcoords) {
  Xcoords.resize(6, 0.);
  Ycoords.resize(6, 0.);
  Zcoords.resize(6, 0.);
  double x = fieldData.GetPhys(fieldData.GetPhysID("x"));
  double y = fieldData.GetPhys(fieldData.GetPhysID("y"));
  double z = fieldData.GetPhys(fieldData.GetPhysID("z"));
  std::vector<int> maxindex;
  invIndex(Np, maxid, maxindex);
  std::vector<double> Xcoords(6, 0.), Ycoords(6, 0.), Zcoords(6, 0.);
  for(int d=0; d<3; ++d) {
    std::vector<int> index = maxindex;
    index[d] -= 1;
    int im = Index(Np, index);
    Xcoords[d*2] = x[im]; Ycoords[d*2] = y[im]; Zcoords[d*2] = z[im];
    index = maxindex;
    index[d] += 1;
    int ip = Index(Np, index);
    Xcoords[2*d+1] = x[ip]; Ycoords[2*d+1] = y[ip]; Zcoords[2*d+1] = z[ip];
  }
}

int main(int argc, char* argv[]) {
  string filenameformat;
  if(argc>1) {
    filenameformat = argv[1];
  }
  IncFlow fieldData;
  fieldData.InputData("Points_R00000_T000.000000.plt");
  int maxid = 0;
  SelectPoint(fieldData, maxid);

  return 0;
}