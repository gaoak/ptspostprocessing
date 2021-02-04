#ifndef INCFLOW_H
#define INCFLOW_H
#include "StructuredData.h"
#include "Body.h"
class IncFlow : public StructuredData {
public:
    IncFlow();
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            std::string bodyname, std::vector<double> param);
    int OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot, const std::vector<double> &omega);
    int TransformCoord(const std::vector<double> &x0);
    int ExtractCore(const double sigma, std::vector<std::vector<double> > & cores,
        std::vector<std::vector<double> > & radius, std::vector<double> &circulation,
        std::vector<double> &inputcenter, const std::vector<int> &vf,
        const int field, const int direction, const bool stoponwall, const double threshold );
    int ExtractCore(const double sigma, std::vector<std::vector<double> > & cores,
        std::vector<std::vector<double> > & radius, std::vector<double> &circulation,
        std::vector<double> &inputcenter, const std::vector<int> &vf,
        const int field = -4, const bool stoponwall = true, const double threshold = 0.);
    int ExtractVortexParam2Dplane(const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
        std::vector<double> &planevorticity, std::vector<double> &radius, double &circulation);
    int ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
        std::vector<double> &data, std::vector<double> &core, bool ismax = true);
    int CalculateVorticity(int order = 2);
    std::pair<int, std::vector<int> > GetProceedDirection(const std::vector<double> &vor, double sign);
    int PurgeDifferentSign(const std::vector<int> &N, const std::vector<double> &v, std::vector<double> &data, double sign);
    int CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                        const std::vector<int> &skip, const IncFlow & big);
    int GetSubdomainRange(const std::vector<int> &center, double radius, std::vector<int> &range);
protected:
    Body m_body;
};
#endif