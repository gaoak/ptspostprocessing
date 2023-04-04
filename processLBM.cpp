#include "FileIO.h"
#include "Util.h"
#include "LBMData.h"
#include "IncFlow.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <iostream>
using namespace std;
double Integration(const std::vector<double> &dx, const std::vector<double> &dy, const std::vector<double> &data);
int ExtractDxDy(const std::vector<DataPack> & zones, std::vector<double> &dx, std::vector<double> &dy);

double ProcessWakeData(const std::vector<DataPack> & zones, std::string message) {
    int Np = zones[0].data[0].size();
    int Idx = 0, Idy = 1, Idu = 3, Idv = 4, Idvor = 5;
    double xmin = -5, xmax = 25., ymin = -5., ymax = 5.;
    double x0 = zones[1].data[0][0];
    double y0 = zones[1].data[1][0];
    std::vector<double> dx, dy;
    ExtractDxDy(zones, dx, dy);
    xmin += x0;
    xmax += x0;
    ymin += y0;
    ymax += y0;
    xmin = max(xmin, zones[0].data[0][0]);
    ymin = max(ymin, zones[0].data[1][0]);
    xmax = min(xmax, zones[0].data[0][Np-1]);
    ymax = min(ymax, zones[0].data[1][Np-1]);
    double Area = (ymax - ymin) * (xmax - xmin);
    std::vector<double> enstrophy(Np, 0.);
    std::vector<double> kinetic(Np, 0.);
    for(int i=0; i<Np; ++i) {
        double x = zones[0].data[Idx][i];
        double y = zones[0].data[Idy][i];
        double u = zones[0].data[Idu][i];
        double v = zones[0].data[Idv][i];
        double W_z = zones[0].data[Idvor][i];
        if(xmin < x && x < xmax && ymin < y && y < ymax) {
            enstrophy[i] = 0.5 * W_z * W_z;
            kinetic[i] = 0.5 * (u*u + v*v);
        }
    }
    double totalenstrophy = Integration(dx, dy, enstrophy) / Area;
    double totalkinetic = Integration(dx, dy, kinetic) / Area;
    std::cout << scientific << totalenstrophy << " " << totalkinetic << " " << message << std::endl;
    return 0;
}

int ExtractDxDy(const std::vector<DataPack> & zones, std::vector<double> &dx, std::vector<double> &dy) {
    std::vector<int> N = zones[0].N;
    int Idx = 0, Idy = 1;
    dx.resize(N[0], 0.);
    dy.resize(N[1], 0.);
    for(int i=0; i<N[0]-1; ++i) {
        dx[i] = zones[0].data[Idx][i+1] - zones[0].data[Idx][i];
    }
    dx[N[0]-1] = 0.;
    for(int j=0; j<N[1]-1; ++j) {
        dy[j] = zones[0].data[Idy][(j+1)*N[0]] - zones[0].data[Idy][j*N[0]];
    }
    dy[N[1]-1] = 0.;
    return 0;
}

double Integration(const std::vector<double> &dx, const std::vector<double> &dy, const std::vector<double> &data) {
    double sum = 0.;
    int offset = 0;
    for(size_t j=0; j<dy.size();++j) {
        for(size_t i=0; i<dx.size(); ++i, ++offset) {
            if(i && j) {
                sum += (dx[i-1] + dx[i]) * (dy[j-1] + dy[j]) * data[offset];
            }
        }
    }
    return sum * 0.25;
}

std::string processMessage(std::string message) {
    if(message.size()==0) {
        return message;
    }
    string cholder = "KD.";
    vector<size_t> holder(cholder.size(), 0);
    size_t underline = message.find('_');
    if(underline == string::npos) {
        return "";
    }
    for(size_t i=0; i<cholder.size(); ++i) {
        holder[i] = message.find(cholder[i], underline);
        if(holder[i] == string::npos) {
            return "";
        }
    }
    string res = message.substr(holder[0]+1, holder[1] - holder[0]-1) + " ";
    res += message.substr(holder[1]+2, holder[2] - holder[1]-2) + " ";
    if(message[0]=='v' || message[0]=='V') {
        res += "override";
    }
    return res;
}

int main(int argc, char* argv[]) {
  string phifilename, lbmfilename, bndfilename0, bndfilename1;
  if(argc<5) {
    cout << "4 input files requred: phi field[uniform grid, integration field], LBM field, boundary line with x y a and n, boundary line with phi" << std::endl;
    return argc - 5;
  }
  phifilename = argv[1];
  lbmfilename = argv[2];
  bndfilename0 = argv[3];
  bndfilename1 = argv[4];
  // load uniform grid of phi file
  IncFlow baseflow;
  baseflow.InputData(phifilename);
  // load LBM file
  std::vector<std::string> variables = {"x", "y", "p", "u", "v", "W_z"};
  std::vector<std::vector<int>> Ns = {{985, 1, 4497}, {101, 1, 1}, {101, 1, 1}};
  LBMData lbmfile(lbmfilename, Ns);
  std::vector<std::vector<double>> u1;
  lbmfile.Interpolation(baseflow, u1);
  baseflow.AddPhysics(variables[3], u1[3]);
  baseflow.AddPhysics(variables[4], u1[4]);
  baseflow.CalculateVorticity();
  //baseflow.OutputData("combinedfile.plt");
  // load boundary file x, y, nx, ny, ax, ay
  std::vector<std::vector<double> > bnddata0;
  InputPoints_ascii(bndfilename0, bnddata0);
  ;
  // get volume integral
  std::vector<double> Qfield    = baseflow.GetPhys(baseflow.GetPhysID("Q"));
  std::vector<double> Phi0field = baseflow.GetPhys(baseflow.GetPhysID("Phi0"));
  for(size_t i=0; i<Qfield.size(); ++i) {
    Qfield[i] *= Phi0field[i];
  }
  double volumeforce = 2.0 * baseflow.Integrate(Qfield);
  return 0;
}

int LoadPhiField(int argc, char* argv[]) {
  string filename("0.plt");
  string message;
  if(argc>1) {
    filename = argv[1];
  }
  if(argc>2) {
    message = argv[2];
  }
  StructuredData file;
  file.InputData(filename);
  file.OutputData("test_"+filename);
  return 0;
}