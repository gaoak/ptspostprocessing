#ifndef STRUCTUREDDATA_H
#define STRUCTUREDDATA_H
#include "TECIO.h"
#include "MASTER.h"
#include<string>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
int OutputTec360(std::string filename, std::string variables,
                 int i, int j, int k, std::vector<void*> data,
                 int isdouble = 1,
                 int debug = 1,
                 int filetype = 0,
                 int fileformat = 0);

void parserDouble(const char * cstr, std::vector<double> & value);

class StructuredData {
public:
    StructuredData(int i, int j = 1, int k = 1);
    int OutputCSV(std::string filename);
    int OutputTec360(std::string filename);
    int LoadCSV(std::string filename);
    int GenPoints(double x0, double x1, double y0, double y1, double z0, double z1);
    int m_i;
    int m_j;
    int m_k;
    int m_Np;
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::string m_varList;
private:
    int ArrayIndex(int i, int j, int k);
    int ParserCSVHeader(const char * header);
};
#endif