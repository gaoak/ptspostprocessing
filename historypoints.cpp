#include<fstream>
#include<map>
#include<cmath>
#include<algorithm>
#include<set>
#include<fstream>
#include "FileIO.h"
#include "Util.h"
using namespace std;

int InputHistoryPointsLine(const std::string filename, double stime, double etime,
                std::vector<int> &N, std::vector<std::vector<double> > &data) {
    N.clear();
    std::ifstream file(filename.c_str());
    if(!file.is_open()) {
        printf("error: unable to open file %s\n", filename.c_str());
        exit(-1);
    }
    char buffer[1000];
    // read header
    file.getline(buffer, sizeof(buffer));
    // read coordinate
    std::vector<double> z;
    std::vector<double> value;
    file.getline(buffer, sizeof(buffer));
    while(buffer[0]=='#') {
        parserDouble(buffer, value);
        z.push_back(value[3]);
        file.getline(buffer, sizeof(buffer));
    }
    N.push_back(z.size());
    // decide data size
    parserDouble(buffer, value);
    data.resize(1+value.size());
    // read physical values
    int index = 0;
    while(!file.eof()) {
        parserDouble(buffer, value);
        if (value[0]<stime) continue;
        if (value[0]>etime) break;
        data[0].push_back(z[index % z.size()]);
        for(size_t i=0; i<value.size(); ++i) {
            data[i+1].push_back(value[i]);
        }
        ++index;
        file.getline(buffer, sizeof(buffer));
    }
    N.push_back(data[0].size() / N[0]);
    file.close();
    return z.size();
}

int main(int argc, char *argv[]) {
    // filename, start time, end time
    if (argc<4) {
        printf("argc is %d [<4]\n", argc);
    }
    string pointfilename(argv[1]);
    double stime = 0.05;
    sscanf(argv[2], "%lf", &stime);
    double etime = 0.05;
    sscanf(argv[3], "%lf", &etime);

    std::vector<std::vector<double> > data; // size [nvar + 1][all data], data[0] is z location
    std::vector<int> N; // [number of points][number of time steps]
    InputHistoryPointsLine(pointfilename, stime, etime, N, data);
    std::vector<std::string> variables = {"z", "t", "u", "v", "w", "p"};
    OutputTec360_binary("history.plt", variables, N, data, 0);

    //do spectral
    std::vector<int> Nsp = {N[0]/2+1, N[1]};
    std::vector<std::string> vsp = {"beta", "t", "w"};
    std::vector<std::vector<double> > datasp(3); //beta, time, spectral
    for(size_t i=0; i<datasp.size(); ++i)
    {
        datasp[i].resize(Nsp[0]*Nsp[1]);
    }
    int offset0 = 0, offset1 = 0;
    double Zlen = data[0][N[0]-1] - data[0][0];
    for (int i=0; i<N[1]; ++i)
    {
        getRealPowerSpectral(&(data[4][offset0]), &(datasp[2][offset1]),
            &(datasp[0][offset1]),Zlen, N[0]-1);
        for(int j=0; j<Nsp[0]; ++j)
        {
            datasp[1][j+offset1] = data[1][j+offset0];
        }
        offset0 += N[0];
        offset1 += Nsp[0];
    }
    OutputTec360_binary("historyspectral.plt", vsp, Nsp, datasp, 0);

    //do summation
    ofstream energy("energy.dat"); // spanwise L2 normal at each time
    energy << "variables = ";
    for(size_t i=0; i<variables.size(); ++i)
    {
        energy << variables[i] << " ";
    }
    energy << "\n";
    int offset = 0;
    for (int i=0; i<N[1]; ++i)
    {
        energy << data[0][offset] << " " << data[1][offset] << " ";
        for(size_t v=2; v<variables.size(); ++v)
        {
            double sum = 0.;
            for(int j=0; j<N[0]-1; ++j)
            {
                sum += data[v][j+offset] * data[v][j+offset];
            }
            sum = sqrt(sum / (N[0]-1));
            energy << sum << " ";
        }
        energy << "\n";
        offset += N[0];
    }
    energy.close();
    return 0;
}