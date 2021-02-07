#ifndef STRUCTUREDDATA_H
#define STRUCTUREDDATA_H
#include<string>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<map>
#include<tuple>
#include "Util.h"

class CoordSystem {
public:
    CoordSystem();
    CoordSystem(const std::vector<double> &range);
    CoordSystem(const std::vector<double> &range, const std::vector<std::vector<double> > & axis);
    std::vector<double> m_o;
    std::vector<std::vector<double> > m_e;
    void ToPhysCoord(std::vector<double> &x) const;
    void ToCompCoord(std::vector<double> &x) const;
};

class StructuredData {
public:
    StructuredData(const std::vector<int> &N, const std::vector<double> &range);
    StructuredData(const std::vector<int> &N, const std::vector<double> &range,
                   const std::vector<std::vector<double> > & axis);
    StructuredData();
    int OutputData(std::string filename, const bool info = true);
    int InputData(std::string filename, const bool info = true);
    int GetTotPoints();
    int AddPhysics(std::string var, void * func);
    int AddPhysics(std::string var, const std::vector<double> &data);
    int RemovePhysics(int i);
    int Smoothing(double sigma, std::vector<int> &field, bool inplace = true);
    int Smoothing(double sigma, std::vector<std::vector<double> > &odata);
    int MaskBoundary(double sigma, std::vector<int> &field, std::map<int, double> def);
    int Diff(std::vector<int > &field, int dir, int order = 2);
    double GetPhysNorm(int f, int p);
    int ExtractPlane(const std::vector<double> &data, std::pair<int, int> plane,
                     const std::vector<int> &range, std::vector<int> & N,
                     std::vector<double> &odata);
    double GetPhysValue(int f, int i);
    double GetCoordValue(int f, int i);
    int GetNumPhys();
    int CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                        const std::vector<int> &skip, std::map<int, double> &field,
                        const StructuredData & big);
    int CopyToSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                        const std::vector<int> &skip, std::map<int, double> &field,
                        StructuredData & small);
    int InterpolateFrom(const StructuredData & origin, std::map<int,double> field);
    int InterpolatePoint (const std::vector<double> & x, std::map<int, double> field,
                        std::map<int, double> &value) const;
    void clear();
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::vector<std::string> m_vars;
protected:
    int GetInterpDimension() const;
    int ParserCSVHeader(const char * header);
    int GenPoints(const std::vector<double> &range);
    int ResetAxis();
    int ReSetNp();
    int Diff(std::vector<std::vector<double> > &u,
             std::vector<std::vector<double> > &du, int dir, int order);
    std::vector<int> m_N;        //dimension 3
    int m_Np;
    CoordSystem m_axis;
    std::vector<double> m_dx;    //dimension 3
};
#endif