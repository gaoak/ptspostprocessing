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

void CalculatePhi(int m, map<string, vector<double>> &phi)
{
    if (phi.find("x") == phi.end() || phi.find("y") == phi.end())
    {
        cout << "x and y are missing when calculating phi." << endl;
        exit(-1);
    }
    vector<double> x, y;
    x            = phi["x"];
    y            = phi["y"];
    int np       = x.size();
    phi["px"]    = vector<double>(np, 0.);
    phi["px_x"]  = vector<double>(np, 0.);
    phi["px_y"]  = vector<double>(np, 0.);
    phi["px_xx"] = vector<double>(np, 0.);
    phi["px_yy"] = vector<double>(np, 0.);
    phi["px_xy"] = vector<double>(np, 0.);
    phi["px_yx"] = phi["px_xy"];
    phi["py"]    = vector<double>(np, 0.);
    phi["py_x"]  = vector<double>(np, 0.);
    phi["py_y"]  = vector<double>(np, 0.);
    phi["py_xx"] = vector<double>(np, 0.);
    phi["py_yy"] = vector<double>(np, 0.);
    phi["py_xy"] = vector<double>(np, 0.);
    phi["py_yx"] = phi["py_xy"];
    double a0 = 0.5;
    double c0 = std::pow(a0, m + 1) / m;
    double c1 = c0 * m;
    double c2 = c0 * m * (m + 1);
    for (int i = 0; i < np; ++i)
    {
        double r     = std::hypot(x[i], y[i]);
        double t     = std::atan2(x[i], y[i]);
        double ir    = 1. / r;
        double irm   = std::pow(ir, m);
        double irm1  = irm * ir;
        double irm2  = irm1 * ir;
        double k0    = c0 * irm;
        phi["px"][i]    = -k0 * cos(m * t);
        phi["py"][i]    = -k0 * sin(m * t);
        double k1    = c1 * irm1;
        phi["px_x"][i]  = k1 * cos((m + 1) * t);
        phi["px_y"][i]  = k1 * sin((m + 1) * t);
        phi["py_x"][i]  = k1 * sin((m + 1) * t);
        phi["py_y"][i]  = -k1 * cos((m + 1) * t);
        double k2    = c2 * irm2;
        phi["px_xx"][i] = -k2 * cos((m + 2) * t);
        phi["px_xy"][i] = -k2 * sin((m + 2) * t);
        phi["px_yy"][i] = k2 * cos((m + 2) * t);
        phi["py_xx"][i] = -k2 * sin((m + 2) * t);
        phi["py_xy"][i] = k2 * cos((m + 2) * t);
        phi["py_yy"][i] = k2 * sin((m + 2) * t);
    }
}

void CalculatePhi0(map<string, vector<double>> &phi)
{
    if (phi.find("x") == phi.end() || phi.find("y") == phi.end())
    {
        cout << "x and y are missing when calculating phi." << endl;
        exit(-1);
    }
    vector<double> x, y;
    x            = phi["x"];
    y            = phi["y"];
    int np       = x.size();
    phi["px"]    = vector<double>(np, 0.);
    phi["px_x"]  = vector<double>(np, 0.);
    phi["px_y"]  = vector<double>(np, 0.);
    phi["px_xx"] = vector<double>(np, 0.);
    phi["px_yy"] = vector<double>(np, 0.);
    phi["px_xy"] = vector<double>(np, 0.);
    phi["px_yx"] = phi["px_xy"];
    phi["py"]    = vector<double>(np, 0.);
    phi["py_x"]  = vector<double>(np, 0.);
    phi["py_y"]  = vector<double>(np, 0.);
    phi["py_xx"] = vector<double>(np, 0.);
    phi["py_yy"] = vector<double>(np, 0.);
    phi["py_xy"] = vector<double>(np, 0.);
    phi["py_yx"] = phi["py_xy"];
    double a0 = 0.5;
    for (int i = 0; i < np; ++i)
    {
        double r     = std::hypot(x[i], y[i]);
        double t     = std::atan2(x[i], y[i]);
        double ir    = 1. / r;
        double ir2   = ir * ir;
        double k0    = a0;
        phi["px"][i]    = k0 * log(r);
        double k1    = a0 * ir;
        phi["px_x"][i]  = k1 * cos(t);
        phi["px_y"][i]  = k1 * sin(t);
        double k2    = a0 * ir2;
        phi["px_xx"][i] = -k2 * cos(2. * t);
        phi["px_xy"][i] = -k2 * sin(2. * t);
        phi["px_yy"][i] = k2 * cos(2. * t);
    }
}

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
  ofile << std::scientific << std::setprecision(20);
  ofile << "variables = x, kflux, press, visc\n";
  for(int i=0; i<N[0]; ++i) {
    ofile << baseflow.GetCoordValue(0, i) << " " << intkux[i] << " " << intpu[i] << " " << inttauum[i] << "\n";
  }
  ofile.close();
}

void processField2(IncFlow &baseflow, int m) {
  double mu = 1. / 100.;
  map<string, vector<double>> phi;
  if(m==0) {
    CalculatePhi0(phi);
  } else {
    CalculatePhi(m, phi);
  }
  vector<int> N = baseflow.GetN();
  int Np = N[0] * N[1] * N[2];
  vector<double> u = baseflow.GetPhys(baseflow.GetPhysID("u"));
  vector<double> v = baseflow.GetPhys(baseflow.GetPhysID("v"));
  vector<double> p = baseflow.GetPhys(baseflow.GetPhysID("p"));
  vector<double> ux = baseflow.GetPhys(baseflow.GetPhysID("ux"));
  vector<double> uy = baseflow.GetPhys(baseflow.GetPhysID("uy"));
  vector<double> vx = baseflow.GetPhys(baseflow.GetPhysID("vx"));
  vector<double> vy = baseflow.GetPhys(baseflow.GetPhysID("vy"));
  vector<double> om(Np, 0.), lx(Np, 0.), ly(Np, 0.), Txx(Np, 0), Txy(Np, 0.), Tyy(Np, 0.);
  vector<double> taux(Np, 0.), tauy(Np, 0.);
  for(int i=0; i<Np; ++i) {
    om[i] = vx[i] -  uy[i];
    Txx[i] = mu * ux[i] * 2.;
    Txy[i] = mu * (uy[i] + vx[i]);
    Tyy[i] = mu * vy[i] * 2.;
    taux[i] = Txx[i];
    tauy[i] = Txy[i];
    lx[i] = -om[i] * v[i];
    ly[i] = om[i] * u[i];
  }
  vector<double> fwtx(Np, 0.), fwty(Np, 0.);
  vector<double> fwsx(Np, 0.), fwsy(Np, 0.);
  vector<double> lwx(Np, 0.), lwy(Np, 0.);
  for(int i=0; i<Np; ++i) {
    lwx[i] = lx[i] * phi["px_x"][i] + ly[i] * phi["px_y"][i];
    lwy[i] = lx[i] * phi["py_x"][i] + ly[i] * phi["py_y"][i];
    fwtx[i] = u[i] * phi["px_x"][i] + v[i] * phi["px_y"][i];
    fwty[i] = u[i] * phi["py_x"][i] + v[i] * phi["py_y"][i];
    fwsx[i] = (ux[i]*u[i]+vx[i]*v[i]) * phi["px_x"][i] + (uy[i]*u[i]+vy[i]*v[i]) * phi["px_y"][i] + Txx[i]*phi["px_xx"][i] + 2.*Txy[i]*phi["px_xy"][i] + Tyy[i] *phi["px_yy"][i];
    fwsy[i] = (ux[i]*u[i]+vx[i]*v[i]) * phi["py_x"][i] + (uy[i]*u[i]+vy[i]*v[i]) * phi["py_y"][i] + Txx[i]*phi["py_xx"][i] + 2.*Txy[i]*phi["py_xy"][i] + Tyy[i] *phi["py_yy"][i];
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
  ofile << std::scientific << std::setprecision(20);
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
  processField2(baseflow);
  printf("Finished.\n");
  return 0;
}