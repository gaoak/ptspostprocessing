#include<string>
#include<vector>
#include<cmath>
#include "StructuredData.h"
#include "Dataprocessing.h"

StructuredData::StructuredData(const std::vector<int> &N, const std::vector<double> &range) {
    m_N = N;
    for(int i=N.size(); i<3; ++i) {
        m_N.push_back(1);
    }
    m_range = range;
    for(int i=range.size(); i<3; ++i) {
        m_range.push_back(0.);
        m_range.push_back(0.);
    }
    m_Np = 1.;
    if(N.size()==0) {
        printf("error: 0 points\n");
        m_Np = 0;
    } else {
        for(int i=0; i<m_N.size(); ++i) {
            m_Np *= m_N[i];
        }
    }
    for(int i=0; i<3; ++i) {
        m_x.push_back(std::vector<double>(m_Np, 0.));
    }
    m_varList = "x, y, z";
    GenPoints();
}

int StructuredData::GetTotPoints() {
    return m_Np;
}

int StructuredData::OutputCSV(std::string filename) {
    std::ofstream file(filename.c_str());
    file << "# " << m_varList << "\n";
    file << std::scientific << std::setprecision(16);
    for(int index=0; index<m_Np; ++index) {
        file << m_x[0][index] << "," << m_x[1][index] << "," << m_x[2][index];
        for(int v=0; v<m_phys.size(); ++v) {
            file << "," << m_phys[v][index];
        }
        file << "\n";
    }
    return m_Np + 1;
}

int StructuredData::ParserCSVHeader(const char * header) {
    int i = 0;
    for(; header[i]=='#' || header[i]==' '; ++i);
    m_varList = header + i;
    int vcount = -2;
    for(;header[i];++i) {
        vcount += header[i]==',';
    }
    return vcount;
}

int StructuredData::LoadCSV(std::string filename) {
    std::ifstream file(filename.c_str());
    char buffer[1000];
    std::vector<double> value;
    file.getline(buffer, sizeof(buffer));
    int vcount = ParserCSVHeader(buffer);
    m_phys.clear();
    for(int i=0; i<vcount; ++i) {
        m_phys.push_back(std::vector<double>(m_Np, 0.));
    }

    file.getline(buffer, sizeof(buffer));
    int index = 0;
    while(!file.eof()) {
        parserDouble(buffer, value);
        for(int i=0; i<vcount; ++i) {
            m_phys[i][index] = value[3+i];
        }
        ++index;
        file.getline(buffer, sizeof(buffer));
    }
    return index;
}

int StructuredData::GenPoints() {
    for(int k=0; k<m_N[2]; ++k) {
        for(int j=0; j<m_N[1]; ++j) {
            for(int i=0; i<m_N[0]; ++i) {
                std::vector<int> ind(3);
                ind[0] = i; ind[1] = j; ind[2] = k;
                int index = Index(m_N, ind);
                if(m_N[0]>1) {
                    m_dx[0] = (m_range[1] - m_range[0])/(m_N[0]-1);
                    m_x[0][index] = m_range[0] + i*m_dx[0];
                } else {
                    m_dx[0] = std::nan("1");
                    m_x[0][index] = m_range[0];
                }
                if(m_N[1]>1) {
                    m_dx[1] = (m_range[3] - m_range[2])/(m_N[1]-1);
                    m_x[1][index] = m_range[2] + j*m_dx[1];
                } else {
                    m_dx[1] = std::nan("1");
                    m_x[1][index] = m_range[2];
                }
                if(m_N[2]>1) {
                    m_dx[2] = (m_range[5] - m_range[4])/(m_N[2]-1);
                    m_x[2][index] = m_range[4] + k*m_dx[2];
                } else {
                    m_dx[2] = std::nan("1");
                    m_x[2][index] = m_range[4];
                }
            }
        }
    }
    return m_Np;
}

int StructuredData::OutputTec360(std::string filename) {
    int isdouble = 1;
    std::vector<void*> data;
    for(int i=0; i<m_x.size(); ++i) {
        data.push_back((void*) (m_x[i].data()) );
    }
    for(int i=0; i<m_phys.size(); ++i) {
        data.push_back((void*) (m_phys[i].data()) );
    }
    ::OutputTec360(filename, m_varList, m_N[0], m_N[1], m_N[2], data, isdouble);
    return m_x.size() + m_phys.size();
}

int StructuredData::Smoothing(double sigma, std::vector<std::vector<double> > &odata) {
    KernelSmooth kernel;
    std::vector<double> sigma_dx(3);
    for(int i=0; i<3; ++i) {
        sigma_dx[i] = sigma / m_dx[i];
    }
    return kernel.DoSmooth(sigma_dx, m_N, odata);
}

int StructuredData::Diff(std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &du, int dir, int order) {
    if(m_N[dir]==0) {
        printf("error: cannot calculate finite difference in 1 data, dir %d\n", dir);
        return 0;
    }
    Derivative der;
    if(dir) {
        ShiftIndex(m_N, u, -dir);
    }
    der.Diff(m_N, u, du, m_dx[dir], order);
    if(dir) {
        ShiftIndex(m_N, u, dir);
        ShiftIndex(m_N, du, dir);
    }
}