#include<fstream>
#include<map>
#include<cmath>
#include<algorithm>
#include "Util.h"
#include "PlungingMotion.h"

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
    } else {
        m_k = 0.;
    }
    if(param.count("A")) {
        m_A = StringToDouble(param["A"]);
    } else {
        m_A = 0.;
    }
    if(param.count("phi")) {
        m_phi = StringToDouble(param["phi"]);
    } else {
        m_phi = 0.;
    }
    if(param.count("input")) {
        m_inputformat = param["input"];
    } else {
        m_inputformat = "unsetinputfile%d";
    }
    if(param.count("output")) {
        m_outputformat = param["output"];
    } else {
        m_outputformat = "unsetoutputfile%d";
    }
    if(param.count("filesnumber")) {
        parserUInt(param["filesnumber"].c_str(), m_file);
    }
    if(param.count("phase")) {
        parserDouble(param["phase"].c_str(), m_phase);
    }
    if(param.count("body")) {
        m_airfoil = param["body"];
    } else {
        m_airfoil = "0000";
    }
    if(param.count("AoA")) {
        m_AoA = StringToDouble(param["AoA"].c_str());
        m_AoA = m_AoA / 180. * M_PI;
    } else {
        m_AoA = 0.;
    }
    if(param.count("span")) {
        parserDouble(param["span"].c_str(), m_span);
    } else {
        m_span.push_back(std::numeric_limits<double>::min());
        m_span.push_back(std::numeric_limits<double>::max());
    }
    if(param.count("threshold")) {
        m_threshold = StringToDouble(param["threshold"].c_str());
    } else {
        m_threshold = 0.;
    }
    if(param.count("vortexcorevar")) {
        parserInt(param["vortexcorevar"].c_str(), m_vortexcoreVar);
    } else {
        m_vortexcoreVar = std::vector<int>(4,0);
    }
    if(param.count("stoponwall")) {
        m_stoponwall = myRound<double>(StringToDouble(param["stoponwall"].c_str()));
    } else {
        m_stoponwall = 0;
    }
    if(param.count("translation")) {
        m_translation = myRound<double>(StringToDouble(param["translation"].c_str()));
    } else {
        m_translation = 0;
    }
    if(param.count("calculateVorticityQ")) {
        m_calculateVorticityQ = myRound<double>(StringToDouble(param["calculateVorticityQ"].c_str()));
    } else {
        m_calculateVorticityQ = 0;
    }
    if(param.count("intialdirection")) {
        m_initDirection = myRound<double>(StringToDouble(param["intialdirection"].c_str()));
    }
    if(param.count("N")) {
        parserUInt(param["N"].c_str(), m_N);
    }
    if(param.count("range")) {
        parserDouble(param["range"].c_str(), m_range);
    }
    if(param.count("sigma")) {
        m_sigma = StringToDouble(param["sigma"].c_str());
    } else {
        m_sigma = -1.;
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

int PlungingMotion::OutputVortexCore(std::string filename, IncFlow &flow) {
    if(m_vortexcoreVar.size()<4 || m_vortexcoreVar[3]==0) {
        return -1;
    }
    for(int i=0; i<3; ++i) {
        if(m_vortexcoreVar[i]<=0 || m_vortexcoreVar[i]>flow.GetNumPhys()) {
            return -1;
        }
    }
    std::vector<std::vector<double> > cores;
    std::vector<std::vector<double> > radius;
    std::vector<double> circulation;
    flow.ExtractCore(m_sigma, cores, radius, circulation, m_initcenter, m_vortexcoreVar,
        m_vortexcoreVar[3], m_initDirection, m_stoponwall, m_threshold);
    if(cores.size()==0) {
        printf("no vortex core found with threshold %f\n", m_threshold);
        return 0;
    }
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
    return cores.size();
}

int PlungingMotion::Dumppoints() {
    int count  = 0;
    for(int n=m_file[0]; n<m_file[2]; n+=m_file[1]) {
        IncFlow flow(m_N, m_range, m_airfoil, {m_AoA});
        flow.OutputData(GetOutFileName(n));
        ++count;
    }
    return count;
}

int PlungingMotion::ProcessFlowData(int dir) {
    std::vector<int> filen;
    for(int n=m_file[0]; n<m_file[2]; n+=m_file[1]) {
        filen.push_back(n);
    }
    if(dir<0) {
        std::reverse(filen.begin(), filen.end());
    }
    int count  = 0;
    for(int k=0; k<filen.size(); ++k) {
        int n = filen[k];
        IncFlow flow(m_N, m_range, m_airfoil, {m_AoA, m_span[0], m_span[1]});
        flow.InputData(GetInFileName(n));
        if(m_airfoil.compare("0000")!=0) {
            double v0 = PlungingVelocity(GetFilePhase(n), m_phi);
            flow.OverWriteBodyPoint({0., v0, 0.}, {0., 0., 0.}, {0., 0., 0.});
        }
        if(m_translation) {
            double h0 = PlungingLocation(GetFilePhase(n), m_phi);
            flow.TransformCoord({0., h0, 0.});
        }
        if(m_sigma>0.) {
            std::vector<int> field;
            for(int i=0; i<flow.GetNumPhys(); ++i) {
                field.push_back(i);
            }
            flow.Smoothing(m_sigma, field);
        }
        //vorticity Q
        //u, v, w, p, W_x, W_y, W_z, Q
        if(m_calculateVorticityQ) {
            flow.CalculateVorticity();
        }
        //vortex
        std::string filename("vortexcore");
        filename += std::to_string(n) + ".dat";
        OutputVortexCore(filename, flow);
        flow.OutputData(GetOutFileName(n));
    }
    return count;
}