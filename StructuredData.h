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
    int GetTotPoints();
    int AddPhysics(std::string var, void * func);
    int Smoothing(double sigma, std::vector<int> &field, bool inplace = true);
    int MaskBoundary(double sigma, std::vector<int> &field, std::map<int, double> def);
    int Diff(std::vector<int > &field, int dir, int order = 2);
    double GetPhysNorm(int f, int p);
    double GetPhysValue(int f, int i);
    double GetCoordValue(int f, int i);
    int GetNumPhys();
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::string m_varList;
private:
    int GenPoints();
    int ParserCSVHeader(const char * header);
    int Smoothing(double sigma, std::vector<std::vector<double> > &odata);
    int Diff(std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &du, int dir, int order);
    std::vector<int> m_N;        //dimension 3
    std::vector<double> m_range; //dimension 6
    std::vector<double> m_dx;    //dimension 3
    int m_Np;
};
#endif