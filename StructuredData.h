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
    int Smoothing(double sigma, std::vector<std::vector<double> > &odata);
    int Diff(std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &du, int dir, int order);
    int GetTotPoints();
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::string m_varList;
private:
    int GenPoints();
    int ParserCSVHeader(const char * header);
    std::vector<int> m_N;        //dimension 3
    std::vector<double> m_range; //dimension 6
    std::vector<double> m_dx;    //dimension 3
    int m_Np;
};
#endif