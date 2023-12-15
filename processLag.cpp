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

void GetVorticityTerms(IncFlow & fieldData, int maxid, std::vector<double> & vorticity, std::vector<double> &Bsurf, std::vector<double> &mBomega, std::vector<double> & curla) {
  double nu = 1/800.;
  vorticity.resize(3);
  mBomega.resize(3);
  curla.resize(3);
  Bsurf.resize(9);
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
  double Lv = fieldData.GetPhys(fieldData.GetPhysID("Lv"))[maxid];
  double Lw = fieldData.GetPhys(fieldData.GetPhysID("Lw"))[maxid];
  curla[0] = nu*Lu;
  curla[1] = nu*Lv;
  curla[2] = nu*Lw;
  vorticity[0] = wy - vz;
  vorticity[1] = uz - wx;
  vorticity[2] = vx - uy;
  double theta = ux + vy + wz;
  Bsurf[0] = theta - ux;
  Bsurf[1] =       - uy;
  Bsurf[2] =       - uz;
  Bsurf[3] =       - vx;
  Bsurf[4] = theta - vy;
  Bsurf[5] =       - vz;
  Bsurf[6] =       - wx;
  Bsurf[7] =       - wy;
  Bsurf[8] = theta - wz;
  for(int d=0; d<3; ++d) {
    int i= 3*d;
    mBomega[d] = - Bsurf[i] * vorticity[0] - Bsurf[i+1] * vorticity[1] - Bsurf[i+2] * vorticity[2];
  }
}

void GetCoords(IncFlow & fieldData, int maxid, std::vector<double> &coords, std::vector<double> &Xcoords, std::vector<double> &Ycoords, std::vector<double> &Zcoords) {
  coords.resize(3, 0.);
  Xcoords.resize(6, 0.);
  Ycoords.resize(6, 0.);
  Zcoords.resize(6, 0.);
  std::vector<int> Np = fieldData.GetN();
  std::vector<int> maxindex;
  invIndex(Np, maxid, maxindex);
  for(int d=0; d<3; ++d) {
    std::vector<int> index = maxindex;
    int i = Index(Np, index);
    coords[d] = fieldData.GetCoordValue(d, i);
    index[d] -= 1;
    int im = Index(Np, index);
    Xcoords[d*2] = fieldData.GetCoordValue(0, im);
    Ycoords[d*2] = fieldData.GetCoordValue(1, im);
    Zcoords[d*2] = fieldData.GetCoordValue(2, im);
    index = maxindex;
    index[d] += 1;
    int ip = Index(Np, index);
    Xcoords[2*d+1] = fieldData.GetCoordValue(0, ip);
    Ycoords[2*d+1] = fieldData.GetCoordValue(1, ip);
    Zcoords[2*d+1] = fieldData.GetCoordValue(2, ip);
  }
}

void GetGradCoords(std::vector<double> &dx, std::vector<double> &Xcoords, std::vector<double> &Ycoords, std::vector<double> &Zcoords, std::vector<double> &F) {
  F.resize(9, 0.);
  for(int d=0; d<3; ++d) {
    F[0+d] = (Xcoords[2*d+1] - Xcoords[2*d]) / (2. * dx[d]);
    F[3+d] = (Ycoords[2*d+1] - Ycoords[2*d]) / (2. * dx[d]);
    F[6+d] = (Zcoords[2*d+1] - Zcoords[2*d]) / (2. * dx[d]);
  }
}

int main(int argc, char* argv[]) {
  IncFlow fieldData;
  string filenameformat("Points");
  if(argc>1) {
    filenameformat = argv[1];
  }
  filenameformat += "_R%05d_T%010.6lf.plt";
  char buff[100];
  double time = 0., dt = 0.001953125;
  std::vector<double> dx0;
    int maxid = 0;
  for(int n=0; n<=2048; ++n) {
    time = dt *n;
    sprintf(buff, filenameformat.c_str(), 0, time);
    string filename(buff);
    fieldData.InputData(filename);
    if(n==0) {
      dx0 = fieldData.GetDetx();
      SelectPoint(fieldData, maxid);
    }
    cout << time << " ";
    std::vector<double>  vorticity;
    std::vector<double>  Bsurf;
    std::vector<double>  mBomega;
    std::vector<double>  curla;
    GetVorticityTerms(fieldData, maxid, vorticity, Bsurf, mBomega, curla);
    std::vector<double> coords;
    std::vector<double> Xcoords;
    std::vector<double> Ycoords;
    std::vector<double> Zcoords;
    GetCoords(fieldData, maxid, coords, Xcoords, Ycoords, Zcoords);
    std::vector<double> F;
    GetGradCoords(dx0, Xcoords, Ycoords, Zcoords, F);
    for(auto v:coords) {
      cout << v << " ";
    }
    cout << -9999 << " ";
    for(auto v:vorticity) {
      cout << v << " ";
    }
    cout << -9999 << " ";
    for(auto v:mBomega) {
      cout << v << " ";
    }
    cout << -9999 << " ";
    for(auto v:curla) {
      cout << v << " ";
    }
    cout << -9999 << " ";
    for(auto v:Bsurf) {
      cout << v << " ";
    }
    cout << -9999 << " ";
    for(auto v:F) {
      cout << v << " ";
    }
    cout << endl;
  }
  return 0;
}