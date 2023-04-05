#include "LineData.h"
#include "FileIO.h"
#include <set>

int LineData::OutputData(std::string filename) {
    std::vector<std::vector<double> > data = m_x;
    for(size_t i=0; i<m_phys.size(); ++i) {
        data.push_back(m_phys[i]);
    }
    OutputTec360_ascii(filename, m_vars, data);
    return 0;
}
int LineData::InputData(std::string filename, std::vector<std::string> &vars, const bool closed) {
    std::vector<std::vector<double> > data;
    InputPoints_ascii(filename, data);
    m_x.clear();
    m_phys.clear();
    std::set<std::string> coordsname;
    coordsname.insert("x");
    coordsname.insert("y");
    coordsname.insert("z");
    for(size_t i=0; i<vars.size(); ++i) {
        if(coordsname.find(vars[i])!=coordsname.end()) {
            m_x.push_back(data[i]);
            m_vars.push_back(vars[i]);
        }
    }
    for(size_t i=0; i<vars.size(); ++i) {
        if(coordsname.find(vars[i])==coordsname.end()) {
            m_phys.push_back(data[i]);
            m_vars.push_back(vars[i]);
        }
    }
    if(m_x.size()) {
        m_Np = m_x[0].size();
    } else {
        m_Np = 0;
    }
    updateIntegrationWeight(closed);
    return (int) m_phys.size();
}
int LineData::AddPhysics(std::string var, const std::vector<double> &data) {
    m_vars.push_back(var);
    m_phys.push_back(data);
    if((int)data.size() != m_Np) {
        m_phys[(int)m_phys.size()-1].resize(m_Np);
    }
    return (int)data.size();
}
int LineData::RemovePhysics(int i) {
    m_phys.erase(m_phys.begin() + i);
    m_vars.erase(m_vars.begin() + i);
    return -1;
}

int LineData::InterpolateFrom(const StructuredData & origin, std::map<int,double> field) {
    std::vector<double> x(origin.GetNumCoords(), 0.);
    std::vector<std::vector<double>> data(field.size());
    for(int i=0; i<m_Np; ++i) {
        for(size_t d=0; d<m_x.size(); ++d) {
            x[d] = m_x[d][i];
        }
        std::map<int, double> value;
        origin.InterpolatePoint (x, field, value);
        int d = 0;
        for(auto it : value) {
            data[d++].push_back(it.second);
        }
    }
    int d = 0;
    for(auto it : field) {
        AddPhysics(origin.GetPhysVarName(it.first), data[d++]);
    }
    return data.size();
}

int LineData::updateIntegrationWeight(const bool closed) {
    m_weight.resize(m_Np);
    std::vector<double> seg(m_Np, 0.);
    for(int i=0; i<m_Np; ++i) {
        double sum = 0.;
        int i1 = i + 1;
        if(i1 == m_Np) i1 = 0;
        for(size_t d=0; d<m_x.size(); ++d) {
            sum += (m_x[d][i1]-m_x[d][i])*(m_x[d][i1]-m_x[d][i]);
        }
        seg[i] = sqrt(sum);
    }
    m_weight[0] = 0.5 * seg[0];
    m_weight[m_Np-1] = 0.5 * seg[m_Np-2];
    for(int i=1; i<m_Np-1; ++i) {
        m_weight[i] = 0.5 * (seg[i-1] + seg[i]);
    }
    if(closed) {
        m_weight[0] += 0.5 * seg[m_Np - 1];
        m_weight[m_Np-1] += 0.5 * seg[m_Np - 1];
    }
    return 0;
}

double LineData::Integrate(const std::vector<double> &data) {
    double sum = 0.;
    for(int i=0; i<m_Np; ++i) {
        sum += data[i] * m_weight[i];
    }
    return sum;
}
int LineData::GetPhysID(std::string v) {
        for(size_t i=0; i<m_vars.size(); ++i) {
        if(v==m_vars[i]) {
            return (int)i - (int)m_x.size();
        }
    }
    printf("Field %s not found.\n", v.c_str());
    return -1;
}