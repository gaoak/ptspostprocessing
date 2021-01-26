#ifndef PLUNGINGDATA_H
#define PLUNGINGDATA_H
#include<string>
class PlungingMotion {
public:
    PlungingMotion(std::string dataconfigue);
    int ProcessFlowData();
    int Dumppoints();
    double GetFilePhase(int n);
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
};
#endif