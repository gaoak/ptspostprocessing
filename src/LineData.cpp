#include "LineData.h"

int LineData::OutputData(std::string filename, const bool info) {
    ;
}
int LineData::InputData(std::string filename, std::vector<std::string> &vars, const bool info) {
    
    InputPoints_ascii(filename, );
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
    ;
}
double LineData::Integrate(const std::vector<double> &data) {
    ;
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