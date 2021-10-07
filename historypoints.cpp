#include<fstream>
#include<map>
#include<cmath>
#include<algorithm>
#include<set>
#include "FileIO.h"
#include "Util.h"
using namespace std;

int InputHistoryPointsLine(const std::string filename, double stime, double etime,
                std::vector<int> &N, std::vector<std::vector<double> > &data) {
    N.clear();
    std::ifstream file(filename.c_str());
    if(!file.is_open()) {
        printf("error: unable to open file %s\n", filename.c_str());
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

    std::vector<std::vector<double> > data;
    std::vector<int> N;
    InputHistoryPointsLine(pointfilename, stime, etime, N, data);
    std::vector<std::string> variables = {"z", "t", "u", "v", "w", "p"};
    OutputTec360_binary("history.plt", variables, N, data, 0);
    return 0;
}