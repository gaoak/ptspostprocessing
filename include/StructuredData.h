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
    int ShuffleIndex(std::map<int, int> ReIndex, std::vector<int> dir,
            std::map<int, int> pm);
    int ResetAxis();
    std::vector<double> GetRange();
    inline std::vector<int> GetN() {return m_N;}
    inline int GetTotPoints() {return m_Np;}
    inline std::string GetPhysVarName(int i) {return m_vars[(int)m_x.size()+i];}
    inline double GetPhysValue(int f, int i) {return m_phys[f][i];}
    inline double GetCoordValue(int f, int i) {return m_x[f][i];}
    inline void SetPhysValue(double v, int f, int i) {m_phys[f][i] = v;}
    inline void SetCoordValue(double v, int f, int i) {m_x[f][i] = v;};
    inline int GetNumPhys() {return m_phys.size();}
protected:
    int GetInterpDimension() const;
    int ParserCSVHeader(const char * header);
    int GenPoints(const std::vector<double> &range);
    int ReSetNp();
    int Diff(std::vector<std::vector<double> > &u,
             std::vector<std::vector<double> > &du, int dir, int order);
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::vector<int> m_N;        //dimension 3
    std::vector<std::string> m_vars;
    int m_Np;
    CoordSystem m_axis;
    std::vector<double> m_dx;    //dimension 3
};
#endif