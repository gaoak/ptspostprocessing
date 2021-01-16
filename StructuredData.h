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
    StructuredData(int i, int j, int k, double x0, double x1,
                   double y0, double y1, double z0, double z1);
    int OutputCSV(std::string filename);
    int OutputTec360(std::string filename);
    int LoadCSV(std::string filename);
    int Smoothing(double l, const std::vector<double> &idata, std::vector<double> &odata);
    int Diff(int dir, const std::vector<double> &idata, std::vector<double> &odata, int order = 2);
    int m_i;
    int m_j;
    int m_k;
    int m_Np;
    double m_dx;
    double m_dy;
    double m_dz;
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::string m_varList;
private:
    int GenPoints(double x0, double x1, double y0, double y1, double z0, double z1);
    int ArrayIndex(int i, int j, int k);
    int ParserCSVHeader(const char * header);
};
#endif