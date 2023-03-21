#include "FileIO.h"
#include "Util.h"
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

double ProcessWakeData(const std::vector<DataPack> & zones) {
    int Np = zones[0].data[0].size();
    int Idx = 0, Idy = 1, Idu = 2, Idv = 3, Idvor = 5;
    double xmin = -5, xmax = 25., ymin = -5., ymax = 5.;
    double x0 = zones[1].data[0][0];
    double y0 = zones[1].data[1][0];
    double Area = (ymax - ymin) * (xmax - xmin);
    xmin += x0;
    xmax += x0;
    ymin += y0;
    ymax += y0;
    std::vector<double> dx, dy;
    ExtractDxDy(zones, dx, dy);
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
    std::cout << "Value: " << totalenstrophy << ", " << totalkinetic << std::endl;
    return 0;
}

int ExtractDxDy(const std::vector<DataPack> & zones, std::vector<double> &dx, std::vector<double> &dy) {
    std::vector<int> N = zones[0].N;
    int Idx = 0, Idy = 1;
    dx.resize(N[0]);
    dy.resize(N[1]);
    for(int i=0; i<N[0]; ++i) {
        dx[i] = zones[0].data[Idx][i+1] - zones[0].data[Idx][i];
    }
    for(int j=0; j<N[1]; ++j) {
        dy[j] = zones[0].data[Idy][(j+1)*N[0]] - zones[0].data[Idy][j*N[0]];
    }
    return 0;
}

double Integration(const std::vector<double> &dx, const std::vector<double> &dy, const std::vector<double> &data) {
    double sum = 0.;
    int offset = 0;
    for(size_t i=0; i<dx.size(); ++i) {
        for(size_t j=0; j<dy.size();++j, ++offset) {
            if(i && j) {
                sum += (dx[i-1] + dx[i]) * (dy[j-1] + dy[j]) * data[offset];
            }
        }
    }
    return sum * 0.25;
}

int main(int argc, char* argv[]) {
  string filename("0.plt");
  if(argc>1) {
    filename = argv[1];
  }
  std::vector<std::string> variables = {"x", "y", "u", "v", "p", "W_z"};
  std::vector<DataPack>  zones(3);
  zones[0].N = {985, 1, 4497};
  zones[1].N = {101, 1, 1};
  zones[2].N = {101, 1, 1};
  std::vector<std::vector<double>> dataBody;
  int isdouble = 0;
  std::map<int, int> vm;

  InputTec360_FSILBM2D(filename, zones);
  ShiftIndex<double>(zones[0].N, zones[0].data, 1);

  for(size_t z = 0; z < zones.size(); ++z) {
    OutputTec360_binary("field" + to_string(z) + ".plt", variables, zones[z].N, zones[z].data, isdouble);
  }
  ProcessWakeData(zones);

  return 0;
}