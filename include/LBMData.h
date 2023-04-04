#ifndef HEADERFILELBMDATA_H
#define HEADERFILELBMDATA_H
#include<string>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<map>
#include<tuple>
#include "Util.h"

struct DataPack {
    std::vector<int> N;
    std::vector<std::vector<double> > data;
};

/***
 * structured data on the uniform grid
***/
class LBMData {
public:
    LBMData(const std::string filename, const std::vector<std::vector<int>> &N);
    int Interpolation(std::vector<std::vector<double>> &x1, std::vector<std::vector<double>> &u1);
protected:
    int Interpolation(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &u,
        std::vector<std::vector<double>> &x1, std::vector<std::vector<double>> &u1);
    int Extractxy();
    std::vector<DataPack> m_zones;
    std::vector<std::vector<double>> m_x;
};
#endif