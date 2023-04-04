#ifndef HEADERFILESTRUCTUREDDATA_H
#define HEADERFILESTRUCTUREDDATA_H
#include<string>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<map>
#include<tuple>
#include "Util.h"

/***
 * coordinate systems, origin, basis vectors
 * range (x0, lengthx, y0, lengthy, z0, lengthz)
***/
class CoordSystem {
public:
    CoordSystem();
    CoordSystem(const std::vector<double> &range);
    CoordSystem(const std::vector<double> &range, const std::vector<std::vector<double> > & axis);
    CoordSystem(const CoordSystem & coord);
    CoordSystem& operator=(const CoordSystem & coord);
    std::vector<double> m_o;
    std::vector<std::vector<double> > m_e;
    void ToPhysCoord(std::vector<double> &x) const;
    void ToCompCoord(std::vector<double> &x) const;
};

/***
 * structured data on the uniform grid
***/
class StructuredData {
public:
    StructuredData(const std::vector<int> &N, const std::vector<double> &range);
    StructuredData(const std::vector<int> &N, const std::vector<double> &range,
                   const std::vector<std::vector<double> > & axis);
    StructuredData();
    StructuredData(const StructuredData & data);
    StructuredData& operator=(const StructuredData & data);
    int OutputData(std::string filename, const bool info = true);
    int InputData(std::string filename, const bool info = true);
    int AddPhysics(std::string var, void func());
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
    double Integrate(const std::vector<double> &data);
    void clear();
    int ShuffleIndex(std::map<int, int> ReIndex, std::vector<int> dir,
            std::map<int, int> pm);
    int ResetAxis();
    int HasField(std::string var);
    std::vector<double> GetRange();
    inline std::vector<int> GetN() {return m_N;}
    inline int GetTotPoints() {return m_Np;}
    inline std::string GetPhysVarName(int i) {return m_vars[(int)m_x.size()+i];}
    inline std::vector<double> & GetPhys(int i) {return m_phys[i]; }
    inline double GetPhysValue(size_t f, int i)  {return f<m_phys.size() ? m_phys[f][i] : 0.;}
    inline double GetCoordValue(size_t f, int i) {return f<m_x.size()    ? m_x[f][i]    : 0.;}
    inline void SetPhysValue(double v, size_t f, int i)  {if(f<m_phys.size()) m_phys[f][i] = v;}
    inline void SetCoordValue(double v, size_t f, int i) {if(f<m_x.size()   ) m_x[f][i]    = v;};
    inline int GetNumPhys() {return m_phys.size();}
    inline int GetNumCoords() {return m_x.size();}
    inline int GetVelocityDimension() {return m_velocityDim; }
protected:
    void UpdateVelocityDimension();
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
    int m_velocityDim;
};
#endif