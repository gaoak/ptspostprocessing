#ifndef HEADERFILEINCFLOW_H
#define HEADERFILEINCFLOW_H
#include "StructuredData.h"
#include "Body.h"
#include<set>
#define SLICERESOLUTION 52

enum VortexMethod {
    XYZPlane,
    VorticityLine,
    VortexMethodSize
};

enum VortexExtractionStopReason {
    StopRepeatPoint,
    StopThreshold,
    StopReachWall,
    StopOutDomain,
    StopParamsError,
    StopPerpendicularError,
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
    IncFlow(std::string bodyname, std::vector<double> param);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            std::string bodyname, std::vector<double> param);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            const std::vector<std::vector<double> > & axis,
            std::string bodyname, std::vector<double> param);
    IncFlow(const std::vector<int> &N, const std::vector<double> &range,
            const std::vector<std::vector<double> > & axis);
    IncFlow(const IncFlow & flow);
    IncFlow& operator=(const IncFlow & flow);
    int OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot,
            const std::vector<double> &omega);
    int TransformCoord(const std::vector<double> &x0);
    int ExtractCoreByPoint(std::vector<std::vector<double> > & cores, std::set<int> &searched,
            std::vector<double> &inputcenter, const std::vector<int> vf,
            const int field = -4, const bool stoponwall = true, const double threshold = 0.,
            const double walltol = 1.e-6, VortexMethod vm = XYZPlane);
    VortexExtractionStopReason RefineCore(
            std::vector<std::vector<double> > & cores, const std::vector<int> vf, const int field);
    int CalculateVorticity(int order = 2);
    int CalculateVorticity2D(int order = 2);
    int CalculateVorticity3D(int order = 2);
    int CalculateForcePartition2D(int order = 2);
    int InterpolateFrom(const IncFlow & origin, std::map<int,double> field);
    int CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
            const std::vector<int> &skip, std::map<int, double> &field, const IncFlow & big);
    int Extract2DVortex(std::vector<std::vector<int>> &intcenters,
            std::vector<std::vector<double>> &physcenters,
            std::vector<std::vector<double>> &info, const std::vector<int> &v,
            const std::pair<int, int> &plane, const double threshold);
    bool IsInBody(std::vector<double> p, double tol);
    int MaskExpData();
protected:
    int ExtractVortexParam2Dplane(
            const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
            std::vector<double> &planevorticity, std::vector<double> &radius, double &circulation);
    VortexExtractionStopReason ExtractCoreByPointDirection(
            std::vector<std::vector<double> > & cores, std::set<int> &searched,
            std::vector<double> &inputcenter, const std::vector<int> vf, const int field, const int direction,
            const bool stoponwall, const double threshold, const double walltol, VortexMethod vm);
    VortexExtractionStopReason ExtractCoreByPointDirectionXYZ(
            std::vector<std::vector<double> > & cores, std::set<int> &searched,
            std::vector<double> &inputcenter, const std::vector<int> vf, const int field,
            const int direction, const bool stoponwall, const double threshold, const double walltol,
            const bool onestep = false);
    VortexExtractionStopReason ExtractCoreByPointDirectionVorticityLine(
            std::vector<std::vector<double> > & cores, std::set<int> &searched,
            std::vector<double> &inputcenter, const std::vector<int> vf, const int field, const int direction,
            const bool stoponwall, const double threshold, const double walltol);
    int SearchOneCoreXYZplane(
            std::vector<int> &intcenter, std::vector<double> &physcenter, std::vector<double> &info,
            const std::vector<int> &v, const double range, const bool ismax, int dir = -1);
    int SearchAllCoreXYZplane(
            std::vector<std::vector<int>> &intcenters,
            std::vector<std::vector<double>> &physcenters, std::vector<std::vector<double>> &info, const std::vector<int> &v,
            const bool ismax, const std::pair<int, int> plane, const double threshold);
    int SearchOneCorePerpendicular(
            std::vector<double> &physcenter, std::vector<double> &info,
            const std::vector<int> &v, const double range, const bool ismax);
    int GetSubdomainRange(const std::vector<int> &center, double radius, std::vector<int> &range);
    std::pair<int, std::vector<int> > GetProceedDirectionInt(const std::vector<double> &vor, double sign);
    std::pair<int, std::vector<double> > GetProceedDirection(const std::vector<double> &vor, double sign);
    Body m_body;
};
#endif