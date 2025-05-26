#ifndef HEADERFILEFlappingDATA_H
#define HEADERFILEFlappingDATA_H
#include<string>
#include "IncFlow.h"
class FlappingMotion {
public:
    FlappingMotion(std::string dataconfigue);
    int ProcessCFDWingData(int dir = 1);
    int ProcessEXPWingData(int dir = 1);
    int ProcessCFDAirfoilData(int dir = 1);
    int TransformFld(int dir = 1);
    void TransformBodyVel(const std::vector<double> &bodyVel, const double omega, const double theta);
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
    double GetFilePhase(int n);
    double FlappingAngularVelocity(double phase, double phi);
    double FlappingAngularAmp(double phase, double phi);
    int Resampling(IncFlow &flow);
    std::vector<int> m_N;
    std::vector<double> m_range;
    double m_f;
    double m_Amp;
    double m_phi;

    std::vector<int> m_file;
    std::vector<int> m_fileseries;
    std::vector<double> m_phase;
    std::vector<double> m_pivot;
    std::vector<std::vector<double>> m_bodyVel;
    std::vector<std::vector<double>> m_bodyVel2;
    std::vector<std::vector<double>> m_bodyLocation;
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
    int m_num;
    VortexMethod m_vortexmethod;
    int m_processVortexCoreCount;
};
#endif