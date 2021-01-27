#ifndef PLUNGINGDATA_H
#define PLUNGINGDATA_H
#include<string>
#include "IncFlow.h"
class PlungingMotion {
public:
    PlungingMotion(std::string dataconfigue);
    int ProcessFlowData(int dir = 1);
    int Dumppoints();
    double GetFilePhase(int n);
    int OutputVortexCore(std::string filename, IncFlow &flow);
    double PlungingVelocity(double phase, double phi);
    double PlungingLocation(double phase, double phi);
    std::string GetInFileName(int n);
    std::string GetOutFileName(int n);
protected:
    std::vector<int> m_N;
    std::vector<double> m_range;
    double m_k;
    double m_A;
    double m_phi;
    std::vector<int> m_file;
    std::vector<double> m_phase;
    std::string m_inputformat;
    std::string m_outputformat;
    std::string m_airfoil;
    double m_AoA;
    double m_sigma;
    int m_vortexcoreVar;
    std::vector<double> m_initcenter;
    int m_stoponwall;
    int m_initDirection;
    double m_threshold;
};
#endif