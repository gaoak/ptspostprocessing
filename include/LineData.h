#ifndef HEADERFILELINEDATA_H
#define HEADERFILELINEDATA_H
#include<string>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<map>
#include<tuple>
#include "Util.h"
#include "StructuredData.h"

class LineData {
public:
    LineData(){};
    int OutputData(std::string filename, const bool info = true);
    int InputData(std::string filename, std::vector<std::string> &vars, const bool info = true);
    int AddPhysics(std::string var, const std::vector<double> &data);
    int RemovePhysics(int i);
    int InterpolateFrom(const StructuredData & origin, std::map<int,double> field);
    double Integrate(const std::vector<double> &data);
    int GetPhysID(std::string v); //get the physics id from the given variable name, -1 is returned if not found
    inline int GetTotPoints() {return m_Np;}
    inline std::vector<double> & GetPhys(int i) {return m_phys[i]; }
protected:
    std::vector<std::vector<double> > m_x; // i first, then j, k last
    std::vector<std::vector<double> > m_phys; //physics fields
    std::vector<std::string> m_vars;
    int m_Np;
};
#endif