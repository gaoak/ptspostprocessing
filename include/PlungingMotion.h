#ifndef HEADERFILEPLUNGINGDATA_H
#define HEADERFILEPLUNGINGDATA_H
#include<string>
#include "IncFlow.h"
class PlungingMotion {
public:
    PlungingMotion(std::string dataconfigue);
    int ProcessCFDWingData(int dir = 1);
    int ProcessEXPWingData(int dir = 1);
    int ProcessCFDAirfoilData(int dir = 1);
    int Dumppoints();
    std::string GetInFileName(int n);
    std::string GetOutFileName(int n);
    std::string GetVortexCoreFileName(int n);
protected:
    int ProcessFiniteWingData(IncFlow &flow, int n);
    int ProcessAirfoilData(IncFlow &flow, int n);
    int ProcessSmoothing(IncFlow &flow, double sigma);
    int ProcessVorticity(IncFlow &flow);
    int ProcessVortexCore(IncFlow &flow, int n, double sigma, std::vector<std::vector<double> > &cores,
            bool outfield = false);
    int GenerateFileSeries();
    int TransformBathCoord(IncFlow &flow, int n);
    int TransformBathVarsName(IncFlow &flow);
    int DoSpanwiseAverage(IncFlow &flow, int n, bool outfield);
    int MaskExpData(IncFlow &flow);
    double GetFilePhase(int n);
    double PlungingVelocity(double phase, double phi);
    double PlungingLocation(double phase, double phi);
    int Resampling(IncFlow &flow);
    std::vector<int> m_N;
    std::vector<double> m_range;
    double m_k;
    double m_A;
    double m_phi;
    std::vector<int> m_file;
    std::vector<int> m_fileseries;
    std::vector<double> m_phase;
    std::string m_inputformat;
    std::string m_outputformat;
    std::string m_airfoil;
    std::string m_vortexcorefileformat;
    double m_AoA;
    std::vector<double> m_span;
    std::vector<double> m_sigma;
    std::vector<int> m_vortexcoreVar;
    std::vector<double> m_initcenter;
    int m_stoponwall;
    double m_threshold;
    int m_calculateVorticityQ;
    int m_translation;
    VortexMethod m_vortexmethod;
    int m_processVortexCoreCount;
};
#endif