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
//attention!!! only window have <io.h>, the existence of <io.h> in unix/Linux needs to be verified
#include <io.h>
using namespace std;
double average(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const double averangex[], const double averangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data);
double DisplaThick_LineIntegration(const std::vector<std::vector<double> > &y, const std::vector<double> &dy, const double &Ny, const double &u_inf, const std::vector<double> &data, const double &beam_y_end);
double MomentumThick_LineIntegration(const std::vector<std::vector<double> > &y, const std::vector<double> &dy, const double &Ny, const double &u_inf, const std::vector<double> &data, const double &beam_y_end);
double EnergyThick_LineIntegration(const std::vector<std::vector<double> > &y, const std::vector<double> &dy, const double &Ny, const double &u_inf, const std::vector<double> &data, const double &beam_y_end);
double LineIntegration(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const std::vector<double> &dx, const std::vector<double> &dy, const double integrangex[], const double integrangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data);
double SurfaceIntegration(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const std::vector<double> &dx, const std::vector<double> &dy, const double integrangex[], const double integrangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data);
int ExtractDxDy(const std::vector<DataPack> & zones, std::vector<double> &dx, std::vector<double> &dy);
int Extractxypuvw(const std::vector<DataPack> & zones, std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, std::vector<std::vector<double> > &p, std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &v,std::vector<std::vector<double> > &W_z);

double ProcessWakeData(const std::vector<DataPack> & zones, std::string message, const double &u_inf, const double &beam_y_end, const double &volumforce) {
    int Np = zones[0].data[0].size();
    std::vector<int> N = zones[0].N;
    int Nx = N[0];
    int Ny = N[1];

    //  modify as needed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // calculate zone [-20,20]*[0,20], modify as needed
    double xmin = -20., xmax = 20., ymin = 0., ymax = 20.;
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
    double Nu = 0.00625;
    double Mu = 0.00625;
    double Rho = 1.0;
    // Basic information integration range
    double averangex[2]= {-5.0,5.0}, averangey[2]= {1.0,1.0};
    double integrangex[2]= {-5.0,5.0}, integrangey[2]= {1.0,20.0};
    double sampflowrangex[2]= {0.0,0.0},sampflowrangey[2]= {0.0,20.0};
    // Power integration range
    double disspation_integrangex[2]= {-5.0,5.0}, disspation_integrangey[2]= {0.0,1.0};
    double presspower_averangex[2]= {5.0,5.0}, presspower_averangey[2]= {0.0,1.0};
    double kineticpower_averangex[2]= {-5.0,5.0}, kineticpower_averangey[2]= {1.0,1.0};
    double Uref = 0.0625;
    double Pref = 0.0002441406;
    double Fref = 0.0019531250;
    double Eref = 0.0039062500;
    double Area = (ymax - 1.0) * (5.0 - (-5.0));
    //  modify as needed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::vector<std::vector<double> > x(Ny,vector<double>(Nx,0)); 
    std::vector<std::vector<double> > y(Ny,vector<double>(Nx,0)); 
    std::vector<std::vector<double> > p(Ny,vector<double>(Nx,0)); 
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
    std::vector<std::vector<double> > velocity_ugradient(Ny,vector<double>(Nx,0));
    double totalenstrophy = 0., totalkinetic = 0., avervelocity_u = 0., avervelocity_v = 0., avervorticity = 0., averboundaryvortexflux = 0., avervelocity_ugradient = 0.;
    Extractxypuvw(zones, x, y, p, u, v, W_z);

    fstream BasicInfor;
    BasicInfor.open("BasicInfor.txt", ios::app);
    if (!BasicInfor.is_open())
    {
        cout << "打开文件失败" << endl;
        return -1;
    }
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(xmin < x[j][i] && x[j][i] < xmax && ymin < y[j][i] && y[j][i] < ymax) {
                enstrophy[j][i] = 0.5 * W_z[j][i] * W_z[j][i];
                kinetic[j][i] = 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                velocity_u[j][i]  = u[j][i];
                velocity_v[j][i]  = v[j][i];
                vorticity[j][i] = W_z[j][i];
                boundaryvortexflux[j][i] = Nu * (W_z[j+1][i]-W_z[j-1][i])/(dy[j-1] + dy[j]);
                velocity_ugradient[j][i] = (u[j+1][i]-u[j-1][i])/(dy[j-1] + dy[j]);
            }
        }
    }
    totalenstrophy = SurfaceIntegration(x, y, dx, dy, integrangex, integrangey, Nx, Ny, enstrophy) / Area;
    totalkinetic = SurfaceIntegration(x, y, dx, dy, integrangex, integrangey, Nx, Ny, kinetic) / Area;
    avervelocity_u = average(x, y, averangex, averangey, Nx, Ny, velocity_u);
    avervelocity_v = average(x, y, averangex, averangey, Nx, Ny, velocity_v);
    avervorticity = average(x, y, averangex, averangey, Nx, Ny, vorticity);
    averboundaryvortexflux = average(x, y, averangex, averangey, Nx, Ny, boundaryvortexflux);
    avervelocity_ugradient = average(x, y, averangex, averangey, Nx, Ny, velocity_ugradient);
    std::cout << scientific << "totalenstrophy" << " " << "totalkinetic" << " " << "avervelocity_u" << " " << "avervelocity_v" << " " << "avervorticity" << " " << "averboundaryvortexflux" << " " << "avervelocity_ugradient" << message << std::endl;
    std::cout << scientific << totalenstrophy << " " << totalkinetic << " " << avervelocity_u << " " << avervelocity_v << " " << avervorticity << " " << averboundaryvortexflux << " " << avervelocity_ugradient <<  message << std::endl;
    BasicInfor << scientific << totalenstrophy << " " << totalkinetic << " " << avervelocity_u << " " << avervelocity_v << " " << avervorticity << " " << averboundaryvortexflux << " " << avervelocity_ugradient << endl;
    BasicInfor.close();

    // Output flow field point data
    fstream SampFlowPoint;
    SampFlowPoint.open("SampFlowPoint.txt", ios::app);
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(sampflowrangey[0]-0.005 < y[j][i] && y[j][i] < sampflowrangey[1]+0.005 && sampflowrangex[0]-0.005 < x[j][i] && x[j][i] < sampflowrangex[1]+0.005){
                SampFlowPoint << x[j][i] << " " << y[j][i] << " " << p[j][i] << " " << u[j][i] << " " << v[j][i] << " " << W_z[j][i] << " " << enstrophy[j][i] << " " << kinetic[j][i] << " " << boundaryvortexflux[j][i] << endl;
            }
        }
    }
    SampFlowPoint.close();

    // Calculation of boundary layer displacement thickness, momentum thickness and energy thickness
    std::vector<double> aver_u(Ny,0);
    double displathick = 0.,momenthick = 0., energythick = 0.;
    fstream BounLayThick;
    BounLayThick.open("BounLayThick.txt", ios::app);
    for(int j=0; j<Ny; ++j) {
        double tmp[2]= {y[j][0],y[j][0]};
        aver_u[j]=average(x, y, averangex, tmp, Nx, Ny, velocity_u);
    }
    displathick = DisplaThick_LineIntegration(y, dy, Ny, u_inf, aver_u, beam_y_end);
    momenthick = MomentumThick_LineIntegration(y, dy, Ny, u_inf, aver_u, beam_y_end);
    energythick = EnergyThick_LineIntegration(y, dy, Ny, u_inf, aver_u, beam_y_end);
    std::cout << scientific << "displathick"  << " " << "momenthick" << " " << "energythick" <<  message << std::endl;
    std::cout << scientific << displathick  << " " << momenthick << " " << energythick <<  message << std::endl;
    BounLayThick << scientific << displathick << " " << momenthick << " " << energythick << endl;
    BounLayThick.close();

    // Calculate power
    std::vector<std::vector<double> > disspation(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > presspower(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > kineticpower(Ny,vector<double>(Nx,0));
    std::vector<std::vector<double> > volumforcepower(Ny,vector<double>(Nx,0));
    double totaldisspation = 0., totalpresspower = 0., totalkineticpower = 0., totalvolumforcepower = 0.;
    fstream Power;
    Power.open("Power.txt", ios::app);
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(xmin < x[j][i] && x[j][i] < xmax && ymin < y[j][i] && y[j][i] < ymax) {
                disspation[j][i] = Mu * ((u[j+1][i]-u[j-1][i])/(dy[j-1] + dy[j])) * ((u[j+1][i]-u[j-1][i])/(dy[j-1] + dy[j])) * Uref*Uref/Pref;
                presspower[j][i] = u[j][i] * p[j][i] / Rho * 0.5*Uref*Uref*Uref/Pref;
                kineticpower[j][i] = Nu * ((0.5 * u[j+1][i]*u[j+1][i])-(0.5 * u[j-1][i]*u[j-1][i]))/(dy[j-1] + dy[j]) * Eref*Eref/Pref;
                volumforcepower[j][i] = volumforce * u[j][i] * Uref*Fref/Pref;
            }
        }
    }
    totaldisspation = SurfaceIntegration(x, y, dx, dy, disspation_integrangex, disspation_integrangey, Nx, Ny, disspation);
    totalpresspower = LineIntegration(x, y, dx, dy, presspower_averangex, presspower_averangey, Nx, Ny, presspower);
    totalkineticpower = LineIntegration(x, y, dx, dy, kineticpower_averangex, kineticpower_averangey, Nx, Ny, kineticpower);
    presspower_averangex[0] = -presspower_averangex[0];
    presspower_averangex[1] = -presspower_averangex[1];
    totalpresspower = totalpresspower + LineIntegration(x, y, dx, dy, presspower_averangex, presspower_averangey, Nx, Ny, presspower);
    kineticpower_averangey[0] = kineticpower_averangey[0]-1;
    kineticpower_averangey[1] = kineticpower_averangey[1]-1;
    totalkineticpower = totalkineticpower + LineIntegration(x, y, dx, dy, kineticpower_averangex, kineticpower_averangey, Nx, Ny, kineticpower);
    totalvolumforcepower = SurfaceIntegration(x, y, dx, dy, disspation_integrangex, disspation_integrangey, Nx, Ny, volumforcepower);
    std::cout << scientific << "totaldisspation"  << " " << "totalpresspower" << " " << "totalkineticpower" << " " << "totalvolumforcepower" <<  message << std::endl;
    std::cout << scientific << totaldisspation  << " " << totalpresspower << " " << totalkineticpower << " " << totalvolumforcepower <<  message << std::endl;
    Power << scientific << totaldisspation  << " " << totalpresspower << " " << totalkineticpower << " " << totalvolumforcepower << endl;
    Power.close();
    
    return 0;
}

int Extractxypuvw(const std::vector<DataPack> & zones, std::vector<std::vector<double> > &x, std::vector<std::vector<double> > &y, std::vector<std::vector<double> > &p, std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &v,std::vector<std::vector<double> > &W_z) {
    int Np = zones[0].data[0].size();
    std::vector<int> N = zones[0].N;
    int Nx = N[0];
    int i,j;
    int Idx = 0, Idy = 1, Idp = 2, Idu = 3, Idv = 4, Idvor = 5;
    for(int n=0; n<Np; ++n) {
        i = n % Nx;
        j = n / Nx;
        x[j][i] = zones[0].data[Idx][n];
        y[j][i] = zones[0].data[Idy][n];
        p[j][i] = zones[0].data[Idp][n];
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
            if(averangey[0]-0.005 < y[j][i] && y[j][i] < averangey[1]+0.005 && averangex[0]-0.005 < x[j][i] && x[j][i] < averangex[1]+0.005){
                sum += data[j][i];
                offset += 1;
            }
        }
    }
    return sum/offset;
}

double DisplaThick_LineIntegration(const std::vector<std::vector<double> > &y, const std::vector<double> &dy, const double &Ny, const double &u_inf, const std::vector<double> &data, const double &beam_y_end) {
    double sum = 0.;
    for(int j=0; j<Ny; ++j) {
        if( 0 < j && j < Ny-1 ) {
            if(beam_y_end < y[j][0]) {
            sum += (dy[j-1] + dy[j]) * (1-data[j]/u_inf);
            }
        }
    }
    return sum * 0.5;
}

double MomentumThick_LineIntegration(const std::vector<std::vector<double> > &y, const std::vector<double> &dy, const double &Ny, const double &u_inf, const std::vector<double> &data, const double &beam_y_end) {
    double sum = 0.;
    for(int j=0; j<Ny; ++j) {
        if( 0 < j && j < Ny-1 ) {
            if(beam_y_end < y[j][0]) {
                sum += (dy[j-1] + dy[j]) * (data[j]/u_inf) * (1-data[j]/u_inf);
            }
        }
    }
    return sum * 0.5;
}

double EnergyThick_LineIntegration(const std::vector<std::vector<double> > &y, const std::vector<double> &dy, const double &Ny, const double &u_inf, const std::vector<double> &data, const double &beam_y_end) {
    double sum = 0.;
    for(int j=0; j<Ny; ++j) {
        if( 0 < j && j < Ny-1 ) {
            if(beam_y_end < y[j][0]) {
                sum += (dy[j-1] + dy[j]) * (data[j]/u_inf) * (1-(data[j]*data[j])/(u_inf*u_inf));
            }
        }
    }
    return sum * 0.5;
}

double LineIntegration(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const std::vector<double> &dx, const std::vector<double> &dy, const double integrangex[], const double integrangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data) {
    double sum = 0.;
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i) {
            if(integrangey[0]-0.005 < y[j][i] && y[j][i] < integrangey[1]+0.005 && integrangex[0]-0.005 < x[j][i] && x[j][i] < integrangex[1]+0.005) {
                if( 0 < j && j < Ny-1 && 0 < i && i < Nx-1) {
                    if(integrangex[0]==integrangex[1]) {
                        sum += (dy[j-1] + dy[j]) * data[j][i];
                    }else if(integrangey[0]==integrangey[1]){
                        sum += (dx[j-1] + dx[j]) * data[j][i];
                    }
                }
            }
        }
    }
    return sum * 0.5;
}

double SurfaceIntegration(const std::vector<std::vector<double> > &x, const std::vector<std::vector<double> > &y, const std::vector<double> &dx, const std::vector<double> &dy, const double integrangex[], const double integrangey[], const double &Nx, const double &Ny, const std::vector<std::vector<double> > &data) {
    double sum = 0.;
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i){
            if(integrangey[0]-0.005 < y[j][i] && y[j][i] < integrangey[1]+0.005 && integrangex[0]-0.005 < x[j][i] && x[j][i] < integrangex[1]+0.005){
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

void getFileNames(string path, vector<string>& files) {
	intptr_t hFile = 0;
	//file information
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*.plt").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if(strcmp(fileinfo.name , "Flow0000.00000.plt" ) !=0 )//search each file except "Flow0000.00000.plt"
			{
				files.push_back(fileinfo.name);
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
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
//   zones[0].N = {997, 1, 1997};
//   zones[1].N = {51, 1, 1};
  std::vector<std::vector<double> > dataBody;
  std::map<int, int> vm;

  vector<string> fileNames;
  string path("ptspostprocessing\\build");//choose directory
  getFileNames(path, fileNames);
  
// Set the far-field velocity and calculation range of boundary layer loss thickness at different times
// If you don't take into account the loss of the boundary layer thickness you can ignore this
  double t = -0.78125;
//   Y-axis coordinates at the end of beam. Be used to determine the calculation range of boundary layer loss thickness
  double beam_y_end[66] ={0.44381, 0.43872, 0.43552, 0.43372, 0.43354, 0.43491, 0.43829, 0.44374, 0.45132, 0.4609, 0.47298, 0.48792, 0.50713, 0.53125, 0.56191, 0.60025, 0.64965, 0.7134, 0.79821, 0.90358, 0.99839, 0.98758, 0.86612, 0.73716, 0.63271, 0.56496, 0.52729, 0.50351, 0.48613, 0.47182, 0.46031, 0.45124, 0.44422, 0.43889, 0.43531, 0.43334, 0.43337, 0.43517, 0.43875, 0.44405, 0.45123, 0.4605, 0.4727, 0.48803, 0.50752, 0.53163, 0.56196, 0.6, 0.6494, 0.71332, 0.79823, 0.90349, 0.99825, 0.98761, 0.86604, 0.73737, 0.63322, 0.56535, 0.52733, 0.50329, 0.48593, 0.47185, 0.4606, 0.45157, 0.44431, 0.4387};
  int BeamInstantY = 0.;

  for (const auto &ph : fileNames) {
    filename = ph;
    std::cout << filename << endl;
    zones[0].N = {997, 1, 1997};
    InputTec360_FSILBM2D(filename, zones);
    ShiftIndex<double>(zones[0].N, zones[0].data, 1);

    message = processMessage(message);

    // If you don't take into account the loss of the boundary layer thickness you can ignore this
    double u_inf = cos(2 * 3.141592653589793 * 0.02 * t);
    double volumforce = -0.2513274122871834624 * sin(2 * 3.141592653589793 * 0.02 * t);
    // -0.2513274122871834624=-4.908738521234052e-4/0.0019531250, 0.0019531250 is Fref, 0.02 = (1/0.00125)/Tref
    if(u_inf == 0) {
        u_inf = 0.000001;
    }

    ProcessWakeData(zones, message, u_inf, beam_y_end[BeamInstantY], volumforce);

    // If you don't take into account the loss of the boundary layer thickness you can ignore this
    BeamInstantY += 1;
    t += 0.78125;
    
  }
  return 0;
}