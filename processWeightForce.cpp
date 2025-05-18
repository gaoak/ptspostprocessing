#include "FileIO.h"
#include "IncFlow.h"
#include "LBMData.h"
#include "LineData.h"
#include "Util.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

using namespace std;

double sum(vector<int> &N, vector<double> &data, int is, int neli, int js, int nelj) {
  vector<double> sumy(nelj+1, 0.);
  for(int j=0; j<=nelj; ++j) {
    int offset = Index(N, {is, j+js, 0});
    sumy[j] = 0.5 * (data[offset] + data[offset + neli]);
    for(int i=1; i<neli; ++i) {
      sumy[j] += data[offset + i];
    }
  }
  double sum = 0.5 * (sumy[0] + sumy[nelj]);
  for(int j=1; j<nelj; ++j) {
    sum += sumy[j];
  }
  return sum;
}

// Nt, x number elements X y number elements
// 256 x 256
void CalculateWeight(IncFlow &wfield, vector<int> Nt, vector<vector<double>> &data, vector<vector<double>> &bnd) {
  // Number of grid points, 4097 x 4097
  vector<double> dx = wfield.GetDetx();
  vector<int> N = wfield.GetN();
  vector<int> dN = N;
  int Npts = 1;
  double area = 1.;
  for(size_t i=0; i<N.size(); ++i) {
    Npts *= Nt[i];
    dN[i] = (N[i] - 1) / Nt[i];
    area *= dx[i];
  }
  // volume weight
  data.resize(4);
  for(auto &it : data) {
    it.resize(Npts, 0.);
  }
  vector<double> u = wfield.GetPhys(wfield.GetPhysID("u"));
  vector<double> v = wfield.GetPhys(wfield.GetPhysID("v"));
  vector<double> muLu = wfield.GetPhys(wfield.GetPhysID("muLu"));
  vector<double> muLv = wfield.GetPhys(wfield.GetPhysID("muLv"));
  for(int j=0; j<Nt[1]; ++j) {
    int js = j * dN[1];
    for(int i=0; i<Nt[0]; ++i) {
      int is = i * dN[0];
      int id = i + j * Nt[0];
      data[0][id] = -sum(N, u, is, dN[0], js, dN[1]) * area;
      data[1][id] = -sum(N, v, is, dN[0], js, dN[1]) * area;
      data[2][id] = sum(N, muLu, is, dN[0], js, dN[1]) * area;
      data[3][id] = sum(N, muLv, is, dN[0], js, dN[1]) * area;
    }
  }
  // bnd weight
  double mu = 1./100.;
  vector<double> ux = wfield.GetPhys(wfield.GetPhysID("u_x"));
  vector<double> uy = wfield.GetPhys(wfield.GetPhysID("u_y"));
  vector<double> vx = wfield.GetPhys(wfield.GetPhysID("v_x"));
  vector<double> vy = wfield.GetPhys(wfield.GetPhysID("v_y"));
  vector<double> Txx(ux.size(), 0.), Txy(ux.size(), 0.), Tyy(ux.size(), 0.);
  for (size_t i = 0; i < ux.size(); ++i) {
    Txx[i] = mu * ux[i] * 2.;
    Txy[i] = mu * (uy[i] + vx[i]);
    Tyy[i] = mu * vy[i] * 2.;
  }
  int Nbndp = 2 * (Nt[0] + Nt[1]);
  bnd.resize(2);
  for(auto &it : bnd) {
    it.resize(Nbndp, 0.);
  }
  int of = 0;
  for(int i=0; i<Nt[0]; ++i) {
    int idu = i * dN[0];
    bnd[0][i+of] = 0.5 * (Txy[idu] + Txy[idu + dN[0]]);
    bnd[1][i+of] = 0.5 * (Tyy[idu] + Tyy[idu + dN[0]]);
    for(int k=1; k<dN[0]; ++k) {
      bnd[0][i+of] += Txy[idu + k];
      bnd[1][i+of] += Tyy[idu + k];
    }
    bnd[0][i+of] *= dx[0];
    bnd[1][i+of] *= dx[0];
  }
  of = Nt[0];
  for(int i=0; i<Nt[0]; ++i) {
    int idu = i * dN[0] + (N[1] - 1) * N[0];
    bnd[0][i+of] = 0.5 * (Txy[idu] + Txy[idu + dN[0]]);
    bnd[1][i+of] = 0.5 * (Tyy[idu] + Tyy[idu + dN[0]]);
    for(int k=1; k<dN[0]; ++k) {
      bnd[0][i+of] += Txy[idu + k];
      bnd[1][i+of] += Tyy[idu + k];
    }
    bnd[0][i+of] *= -dx[0];
    bnd[1][i+of] *= -dx[0];
  }
  of = Nt[0] * 2;
  for(int j=0; j<Nt[1]; ++j) {
    int idu = j * dN[1] * N[0];
    bnd[0][j+of] = 0.5 * (Txx[idu] + Txx[idu + N[0] * dN[1]]);
    bnd[1][j+of] = 0.5 * (Txy[idu] + Txy[idu + N[0] * dN[1]]);
    for(int k=1; k<dN[1]; ++k) {
      bnd[0][j+of] += Txx[idu + N[0] * k];
      bnd[1][j+of] += Txy[idu + N[0] * k];
    }
    bnd[0][j+of] *= dx[1];
    bnd[1][j+of] *= dx[1];
  }
  of = Nt[0] * 2 + Nt[1];
  for(int j=0; j<Nt[1]; ++j) {
    int idu = N[0] - 1 + j * dN[1] * N[0];
    bnd[0][j+of] = 0.5 * (Txx[idu] + Txx[idu + N[0] * dN[1]]);
    bnd[1][j+of] = 0.5 * (Txy[idu] + Txx[idu + N[0] * dN[1]]);
    for(int k=1; k<dN[1]; ++k) {
      bnd[0][j+of] += Txx[idu + N[0] * k];
      bnd[1][j+of] += Txy[idu + N[0] * k];
    }
    bnd[0][j+of] *= -dx[1];
    bnd[1][j+of] *= -dx[1];
  }
}

void GetFlow(IncFlow &ffield, vector<int> Nt, vector<vector<double>> &data, vector<vector<double>> &bnd) {
  vector<int> N = ffield.GetN();
  vector<int> dN = N;
  int Npts = 1;
  for(size_t i=0; i<N.size(); ++i) {
    Npts *= Nt[i];
    dN[i] = (N[i] - 1) / Nt[i];
  }
  // volume data
  data.resize(6);
  for(auto &it : data) {
    it.resize(Npts, 0.);
  }
  vector<double> x = ffield.GetCoord(ffield.GetCoordID("x"));
  vector<double> y = ffield.GetCoord(ffield.GetCoordID("y"));
  vector<double> ax = ffield.GetPhys(ffield.GetPhysID("ax"));
  vector<double> ay = ffield.GetPhys(ffield.GetPhysID("ay"));
  vector<double> u = ffield.GetPhys(ffield.GetPhysID("u"));
  vector<double> v = ffield.GetPhys(ffield.GetPhysID("v"));
  for(int j=0; j<Nt[1]; ++j) {
    int js = j * dN[1] + dN[1] / 2;
    for(int i=0; i<Nt[0]; ++i) {
      int idu = i * dN[0] + dN[0] / 2 + js * N[0];
      int id = i + j * Nt[0];
      data[0][id] = x[idu];
      data[1][id] = y[idu];
      data[2][id] = ax[idu];
      data[3][id] = ay[idu];
      data[4][id] = u[idu];
      data[5][id] = v[idu];
    }
  }
  // boundary data
  int Nbndp = 2 * (Nt[0] + Nt[1]);
  bnd.resize(4);
  for(auto &it : bnd) {
    it.resize(Nbndp, 0.);
  }
  int of = 0;
  for(int i=0; i<Nt[0]; ++i) {
    int idu = i * dN[0] + dN[0] / 2;
    bnd[0][i+of] = x[idu];
    bnd[1][i+of] = y[idu];
    bnd[2][i+of] = u[idu];
    bnd[3][i+of] = v[idu];
  }
  of = Nt[0];
  for(int i=0; i<Nt[0]; ++i) {
    int idu = i * dN[0] + dN[0] / 2 + (N[1] - 1) * N[0];
    bnd[0][i+of] = x[idu];
    bnd[1][i+of] = y[idu];
    bnd[2][i+of] = u[idu];
    bnd[3][i+of] = v[idu];
  }
  of = Nt[0] * 2;
  for(int j=0; j<Nt[1]; ++j) {
    int idu = (j * dN[1] + dN[1] / 2) * N[0];
    bnd[0][j+of] = x[idu];
    bnd[1][j+of] = y[idu];
    bnd[2][j+of] = u[idu];
    bnd[3][j+of] = v[idu];
  }
  of = Nt[0] * 2 + Nt[1];
  for(int j=0; j<Nt[1]; ++j) {
    int idu = N[0] - 1 + (j * dN[1] + dN[1] / 2) * N[0];
    bnd[0][j+of] = x[idu];
    bnd[1][j+of] = y[idu];
    bnd[2][j+of] = u[idu];
    bnd[3][j+of] = v[idu];
  }
}

void processField(string wfieldname, string flowfieldname, vector<int> Nt) {
  IncFlow wfield;
  IncFlow flowfield;
  wfield.InputData(wfieldname);
  flowfield.InputData(flowfieldname);
  vector<vector<double>> wdata, wbnd, fdata, fbnd;
  CalculateWeight(wfield, Nt, wdata, wbnd);
  GetFlow(flowfield, Nt, fdata, fbnd);
  double sum = 0., sum1 = 0., sum0 = 0.;
  for(size_t i=0; i<wdata[0].size(); ++i) {
    for(int k=0; k<4; ++k) {
      sum0 += wdata[k][i] * fdata[2 + k][i];
    }
  }
  for(size_t i=0; i<wbnd[0].size(); ++i) {
    for(int k=0; k<2; ++k) {
      sum1 += wbnd[k][i] * fbnd[2 + k][i];
    }
  }
  sum = sum0 + sum1;
  cout << std::scientific << std::setprecision(20);
  cout << "SUMFORCE " << sum << " " << sum0 << " " << sum1 << endl;
  vector<string> volstr = {"x", "y", "ax", "ay", "u", "v"};
  OutputTec360_binary("volume.plt", volstr, Nt, fdata, 0);
  vector<string> bndstr = {"x", "y", "u", "v"};
  OutputTec360_ascii("bnd.dat", bndstr, fbnd);
}

int main(int argc, char *argv[]) {
  if (argc < 5) {
    cout << "2 input files requred: weight, flow" << std::endl;
    return -argc;
  }
  string weightfilename = argv[1];
  string flowfilename = argv[2];
  int M = stoi(argv[3]);
  int N = stoi(argv[4]);
  printf("processing file %s\n", argv[1]);
  vector<int> Nt = {M, N, 1};
  processField(weightfilename, flowfilename, Nt);
  printf("Finished.\n");
  return 0;
}
