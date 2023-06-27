#include "FileIO.h"
#include "Util.h"
#include "LBMData.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <iostream>
using namespace std;
double average(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const double averangex[], const double averangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data);
double Integration(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const std::vector<double> &dx, const std::vector<double> &dy, const double integrangex[], const double integrangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data);
int ExtractDxDy(const std::vector<DataPack> & zones, std::vector<double> &dx, std::vector<double> &dy);
int Extractxyuvw(const std::vector<DataPack> & zones, std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &v,std::vector<std::vector<double> > &W_z);

double ProcessWakeData(const std::vector<DataPack> & zones, std::string message) {
    int Np = zones[0].data[0].size();
    std::vector<int> N = zones[0].N;
    int Nx = N[0];
    int Ny = N[1];
    // calculate zone [-20,20]*[0,20], modify as needed
    double xmin = -20., xmax = 20., ymin = 0., ymax = 20.;
    double Nu = 0.00625;
    double x0 = 0;
    double y0 = 0;
    xmin += x0;
    xmax += x0;
    ymin += y0;
    ymax += y0;
    xmin = max(xmin, zones[0].data[0][0]);
    ymin = max(ymin, zones[0].data[1][0]);
    xmax = min(xmax, zones[0].data[0][Np-1]);
    ymax = min(ymax, zones[0].data[1][Np-1]);
    // average range [-5,5]*[1,1], integrate range [-5,5]*[1,20],modify as needed
    double averangex[2]= {-5.0,5.0}, averangey[2]= {1.0,1.0};
    double integrangex[2]= {-5.0,5.0}, integrangey[2]= {1.0,20.0};
    double Area = (ymax - 1.0) * (5.0 - (-5.0));
    std::vector<std::vector<double> > x(Ny,vector<double>(Nx,0)); 
    std::vector<std::vector<double> > y(Ny,vector<double>(Nx,0)); 
    std::vector<std::vector<double> > u(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > v(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > W_z(Ny,vector<double>(Nx,0));
    std::vector<double> dx, dy;
    ExtractDxDy(zones, dx, dy);
    std::vector<std::vector<double> > enstrophy(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > kinetic(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > velocity_u(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > velocity_v(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > vorticity(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > boundaryvortexflux(Ny,vector<double>(Nx,0));
    double totalenstrophy = 0., totalkinetic = 0., avervelocity_u = 0., avervelocity_v = 0., avervorticity = 0., totalboundaryvortexflux = 0.;
    fstream outfile;
    outfile.open("out.txt", ios::app);
    if (!outfile.is_open())
    {
        cout << "打开文件失败" << endl;
        return -1;
    }
    Extractxyuvw(zones, x, y, u, v, W_z);
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(xmin < x[j][i] && x[j][i] < xmax && ymin < y[j][i] && y[j][i] < ymax) {
                enstrophy[j][i] = 0.5 * W_z[j][i] * W_z[j][i];
                kinetic[j][i] = 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                velocity_u[j][i]  = u[j][i];
                velocity_v[j][i]  = v[j][i];
                vorticity[j][i] = W_z[j][i];
                boundaryvortexflux[j][i] = Nu * (W_z[j+1][i]-W_z[j-1][i])/(dy[j-1] + dy[j]);
            }
        }
    }
    totalenstrophy = Integration(x, y, dx, dy, integrangex, integrangey, Nx, Ny, enstrophy) / Area;
    totalkinetic = Integration(x, y, dx, dy, integrangex, integrangey, Nx, Ny, kinetic) / Area;
    avervelocity_u = average(x, y, averangex, averangey, Nx, Ny, velocity_u);
    avervelocity_v = average(x, y, averangex, averangey, Nx, Ny, velocity_v);
    avervorticity = average(x, y, averangex, averangey, Nx, Ny, vorticity);
    totalboundaryvortexflux = average(x, y, averangex, averangey, Nx, Ny, boundaryvortexflux);
    std::cout << scientific << "totalenstrophy" << " " << "totalkinetic" << " " << "avervelocity_u" << " " << "avervelocity_v" << " " << "avervorticity" << " " << "totalboundaryvortexflux" << message << std::endl;
    std::cout << scientific << totalenstrophy << " " << totalkinetic << " " << avervelocity_u << " " << avervelocity_v << " " << avervorticity << " " << totalboundaryvortexflux <<  message << std::endl;
    outfile << scientific << totalenstrophy << " " << totalkinetic << " " << avervelocity_u << " " << avervelocity_v << " " << avervorticity << " " << totalboundaryvortexflux << endl;
    outfile.close();
    return 0;
}

int Extractxyuvw(const std::vector<DataPack> & zones, std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &v,std::vector<std::vector<double> > &W_z) {
    int Np = zones[0].data[0].size();
    std::vector<int> N = zones[0].N;
    int Nx = N[0];
    int i,j;
    int Idx = 0, Idy = 1, Idu = 3, Idv = 4, Idvor = 5;
    for(int n=0; n<Np; ++n) {
        i = n % Nx;
        j = n / Nx;
        x[j][i] = zones[0].data[Idx][n];
        y[j][i] = zones[0].data[Idy][n];
        u[j][i] = zones[0].data[Idu][n];
        v[j][i] = zones[0].data[Idv][n];
        W_z[j][i] = zones[0].data[Idvor][n];
    }
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

double average(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const double averangex[], const double averangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data) {
    double sum = 0.;
    int offset = 0;
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(averangey[0]-0.01 < y[j][i] && y[j][i] < averangey[1]+0.01 && averangex[0]-0.01 < x[j][i] && x[j][i] < averangex[1]+0.01){
                sum += data[j][i];
                offset += 1;
            }
        }
    }
    return sum/offset;
}

double Integration(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const std::vector<double> &dx, const std::vector<double> &dy, const double integrangex[], const double integrangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data) {
    double sum = 0.;
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(integrangey[0]-0.01 < y[j][i] && y[j][i] < integrangey[1]+0.01 && integrangex[0]-0.01 < x[j][i] && x[j][i] < integrangex[1]+0.01){
                if( 0 < j && j < Ny-1 && 0 < i && i < Nx-1) {
                    sum += (dx[i-1] + dx[i]) * (dy[j-1] + dy[j]) * data[j][i];
                }
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
  string filename("0.plt");
  string message;
  if(argc>1) {
    filename = argv[1];
  }
  if(argc>2) {
    message = argv[2];
  }
  std::vector<std::string> variables = {"x", "y", "p", "u", "v", "W_z"};
  std::vector<DataPack>  zones(1);
  zones[0].N = {997, 1, 1997};
//   zones[1].N = {51, 1, 1};
  std::vector<std::vector<double> > dataBody;
  std::map<int, int> vm;

  InputTec360_FSILBM2D(filename, zones);
  ShiftIndex<double>(zones[0].N, zones[0].data, 1);

  message = processMessage(message);
  ProcessWakeData(zones, message);

  return 0;
}