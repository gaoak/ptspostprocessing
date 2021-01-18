#ifndef STRUCTUREDDATA_H
#define STRUCTUREDDATA_H
#include<string>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<map>
#include "Tecplotwraper.h"
#include "Util.h"

class StructuredData {
public:
    StructuredData(const std::vector<int> &N, const std::vector<double> &range);
    int OutputCSV(std::string filename);
    int OutputTec360(std::string filename);
    int LoadCSV(std::string filename);
    int Smoothing(double sigma, const std::vector<double> &idata, std::vector<double> &odata);
    int Diff(int dir, const std::vector<double> &idata, std::vector<double> &odata, int order = 2);
    int m_Np;
    std::vector<int> m_N;
    std::vector<double> m_range;
    std::vector<double> m_dx;
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::string m_varList;
private:
    int GenPoints();
    int ParserCSVHeader(const char * header);
};
#endif