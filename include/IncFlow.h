#ifndef INCFLOW_H
#define INCFLOW_H
#include "StructuredData.h"
#include "Body.h"
#include<set>

enum VortexMethod {
    XYZPlane,
    PerpendicularPlane,
    VortexMethodSize
};

enum VortexExtractionStopReason {
    StopRepeatPoint,
    StopThreshold,
    StopReachWall,
    StopOutDomain,
    StopError,
    StopMaxTry,
    StopSuccess,
    StopRefineError,
    VortexExtractionStopReasonSize
};

enum CoreInfo {
    Corex,
    Corey,
    Corez,
    CoreR1,
    CoreR2,
    CoreGamma,
    CoreVectx,
    CoreVecty,
    CoreVectz,
    CoreCriterion,
    CoreInfoSize
};

class IncFlow : public StructuredData {
public:
    IncFlow();
    IncFlow(const std::vector<int> &N, const std::vector<double> &range);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            std::string bodyname, std::vector<double> param);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            const std::vector<std::vector<double> > & axis,
            std::string bodyname, std::vector<double> param);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            const std::vector<std::vector<double> > & axis);
    IncFlow(const IncFlow & flow);
    int OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot,
            const std::vector<double> &omega);
    int TransformCoord(const std::vector<double> &x0);
    int ExtractCoreByPoint(std::vector<std::vector<double> > & cores, std::set<int> &searched,
            std::vector<double> &inputcenter, const std::vector<int> vf,
            const int field = -4, const bool stoponwall = true, const double threshold = 0.,
            const double walltol = 1.e-6);
    VortexExtractionStopReason RefineCore(
            std::vector<std::vector<double> > & cores, const std::vector<int> vf, const int field);
    int CalculateVorticity(int order = 2);
    int InterpolateFrom(const IncFlow & origin, std::map<int,double> field);
    int CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
            const std::vector<int> &skip, std::map<int, double> &field, const IncFlow & big);
protected:
    int ExtractVortexParam2Dplane(
            const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
            std::vector<double> &planevorticity, std::vector<double> &radius, double &circulation);
    VortexExtractionStopReason ExtractCoreByPointDirection(
            std::vector<std::vector<double> > & cores, std::set<int> &searched,
            std::vector<double> &inputcenter, const std::vector<int> vf, const int field,
            const int direction, const bool stoponwall, const double threshold, const double walltol);
    int SearchOneCoreXYZplane(
            std::vector<int> &intcenter, std::vector<double> &physcenter, std::vector<double> &info,
            const std::vector<int> &v, const double range, const bool ismax, int dir = -1);
    int SearchOneCorePerpendicular(
            std::vector<int> &intcenter, std::vector<double> &physcenter, std::vector<double> &info,
            const std::vector<int> &v, const double range, const bool ismax);
    int GetSubdomainRange(const std::vector<int> &center, double radius, std::vector<int> &range);
    std::pair<int, std::vector<int> > GetProceedDirection(const std::vector<double> &vor, double sign);
    Body m_body;
};
#endif