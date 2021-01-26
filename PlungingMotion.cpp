#include<fstream>
#include<map>
#include<cmath>
#include "Util.h"
#include "PlungingMotion.h"
#include "IncFlow.h"

PlungingMotion::PlungingMotion(std::string dataconfigue) {
    std::ifstream conf(dataconfigue.c_str());
    if(!conf.is_open()) {
        printf("error: unable to open configue file %s\n", dataconfigue.c_str());
    }
    char buffer[1000];
    std::map<std::string, std::string> param;
    conf.getline(buffer, sizeof(buffer));
    while(!conf.eof()) {
        std::vector<std::string> p;
        parserString(buffer, p);
        if(p.size()>1) {
            param[p[0]] = p[1];
        }
        conf.getline(buffer, sizeof(buffer));
    }
    if(param.count("k")) {
        m_k = StringToDouble(param["k"]);
    }
    if(param.count("A")) {
        m_A = StringToDouble(param["A"]);
    }
    if(param.count("phi")) {
        m_phi = StringToDouble(param["phi"]);
    }
    if(param.count("phi")) {
        m_phi = StringToDouble(param["phi"]);
    }
    if(param.count("input")) {
        m_inputformat = param["input"];
    }
    if(param.count("output")) {
        m_outputformat = param["output"];
    }
    if(param.count("filesnumber")) {
        parserUInt(param["filesnumber"].c_str(), m_file);
    }
    if(param.count("phase")) {
        parserDouble(param["phase"].c_str(), m_phase);
    }
    if(param.count("body")) {
        m_airfoil = param["body"];
    }
    if(param.count("AoA")) {
        m_AoA = StringToDouble(param["AoA"].c_str());
        m_AoA = m_AoA / 180. * M_PI;
    }
    if(param.count("vortexcorevar")) {
        m_vortexcoreVar = StringToDouble(param["vortexcorevar"].c_str());
    }
    if(param.count("N")) {
        parserUInt(param["N"].c_str(), m_N);
    }
    if(param.count("range")) {
        parserDouble(param["range"].c_str(), m_range);
    }
    if(param.count("sigma")) {
        m_sigma = StringToDouble(param["sigma"].c_str());
    }
    if(param.count("initcenter")) {
        parserDouble(param["initcenter"].c_str(), m_initcenter);
    }
}

double PlungingMotion::GetFilePhase(int n) {
    return m_phase[0] + m_phase[1] * (n - m_file[0]) / m_file[1];
}

std::string PlungingMotion::GetInFileName(int n) {
    char buffer[100];
    sprintf(buffer, m_inputformat.c_str(), n);
    std::string res(buffer);
    return res;
}

std::string PlungingMotion::GetOutFileName(int n) {
    char buffer[100];
    sprintf(buffer, m_outputformat.c_str(), n);
    std::string res(buffer);
    return res;
}

double PlungingMotion::PlungingVelocity(double phase, double phi) {
    return -m_k*m_A*sin(phase*2.*M_PI + phi);
}

double PlungingMotion::PlungingLocation(double phase, double phi) {
    return 0.5*m_A*cos(phase*2.*M_PI + phi);
    //0.5 A cos(2 k t + phi)
}

int PlungingMotion::Dumppoints() {
    int count  = 0;
    for(int n=m_file[0]; n<m_file[2]; n+=m_file[1]) {
        IncFlow flow(m_N, m_range, m_airfoil, {m_AoA});
        flow.OutputCSV(GetOutFileName(n)+".csv");
        ++count;
    }
    return count;
}

int PlungingMotion::ProcessFlowData() {
    int count  = 0;
    std::vector<double> center = m_initcenter;
    for(int n=m_file[0]; n<m_file[2]; n+=m_file[1]) {
        IncFlow flow(m_N, m_range, m_airfoil, {m_AoA});
        flow.LoadCSV(GetInFileName(n)+".csv");
        double v0 = PlungingVelocity(GetFilePhase(n), m_phi);
        flow.OverWriteBodyPoint({0., v0, 0.}, {0., 0., 0.}, {0., 0., 0.});
        //vorticity Q
        //u, v, w, p, W_x, W_y, W_z, Q
        std::vector<int> field = {0,1,2};
        flow.Smoothing(m_sigma, field);
        flow.CalculateVorticity();
        //vortex
        std::vector<std::vector<double> > cores;
        std::vector<std::vector<double> > radius;
        std::vector<double> circulation;
        flow.ExtractCore(m_sigma, cores, radius, circulation, center, -2, m_vortexcoreVar);
        std::string filename = "core" + GetOutFileName(n) + TECPLOTEXT;
        std::ofstream ofile(filename.c_str());
        ofile << "variables = x,y,z,radius1,radius2,Gamma" << std::endl;
        for(int i=0; i<cores.size(); ++i) {
            ofile << cores[i][0] << " "
                << cores[i][1] << " "
                << cores[i][2] << " "
                << radius[i][0] << " "
                << radius[i][1] << " "
                << circulation[i] << "\n";
        }
        ofile.close();
        field = {7};
        flow.Smoothing(m_sigma, field);
        flow.OutputTec360(GetOutFileName(n) + TECPLOTEXT);
    }
    return count;
}