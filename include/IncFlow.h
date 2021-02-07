#ifndef INCFLOW_H
#define INCFLOW_H
#include "StructuredData.h"
#include "Body.h"
#include<set>
class IncFlow : public StructuredData {
public:
    IncFlow();
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            std::string bodyname, std::vector<double> param);
    int OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot, const std::vector<double> &omega);
    int TransformCoord(const std::vector<double> &x0);
    int ExtractCore(const double sigma, std::vector<std::vector<double> > & cores, std::set<int> &searched,
                    std::vector<double> &inputcenter, const std::vector<int> &vf, const int field,
                    const int direction, const bool stoponwall, const double threshold);
    int ExtractCore(const double sigma, std::vector<std::vector<double> > & cores, std::set<int> &searched,
        std::vector<double> &inputcenter, const std::vector<int> &vf,
        const int field = -4, const bool stoponwall = true, const double threshold = 0.);
    int ExtractVortexParam2Dplane(const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
        std::vector<double> &planevorticity, std::vector<double> &radius, double &circulation);
    int CalculateVorticity(int order = 2);
    std::pair<int, std::vector<int> > GetProceedDirection(const std::vector<double> &vor, double sign);
    int CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                        const std::vector<int> &skip, std::map<int, double> &field,
                        const IncFlow & big);
    int GetSubdomainRange(const std::vector<int> &center, double radius, std::vector<int> &range);
    int SearchOneCoreXYZplane(
        std::vector<int> &intcenter, std::vector<double> &physcenter, std::vector<double> &info,
        const std::vector<int> &v, const double range, const bool ismax);
protected:
    Body m_body;
};

enum vortexcore{
    x = 0,
    y = 1,
    z = 2,
    r1 = 3,
    r2 = 4,
    g = 5,
    size = 6
};
#endif