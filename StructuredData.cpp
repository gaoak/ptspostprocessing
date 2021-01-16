#include<string>
#include<vector>
#include<cmath>
#include "StructuredData.h"

StructuredData::StructuredData(int i, int j, int k, double x0, double x1,
                               double y0, double y1, double z0, double z1) {
    m_i = i;
    m_j = j;
    m_k = k;
    m_Np = m_i * m_j * m_k;
    if(m_Np==0) {
        printf("error: 0 points\n");
    }
    for(int i=0; i<3; ++i) {
        m_x.push_back(std::vector<double>(m_Np, 0.));
    }
    m_varList = "x, y, z";
    GenPoints(x0, x1, y0, y1, z0, z1);
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

int StructuredData::GenPoints(double x0, double x1, double y0, double y1, double z0, double z1) {
    for(int k=0; k<m_k; ++k) {
        for(int j=0; j<m_j; ++j) {
            for(int i=0; i<m_i; ++i) {
                int index = ArrayIndex(i, j, k);
                if(m_i>1) {
                    m_dx = (x1 - x0)/(m_i-1);
                    m_x[0][index] = x0 + i*m_dx;
                } else {
                    m_dx = std::nan("1");
                    m_x[0][index] = x0;
                }
                if(m_j>1) {
                    m_dy = (y1 - y0)/(m_j-1);
                    m_x[1][index] = y0 + j*m_dy;
                } else {
                    m_dy = std::nan("1");
                    m_x[1][index] = y0;
                }
                if(m_k>1) {
                    m_dz = (z1 - z0)/(m_k-1);
                    m_x[2][index] = z0 + k*m_dz;
                } else {
                    m_dz = std::nan("1");
                    m_x[2][index] = z0;
                }
            }
        }
    }
    return m_Np;
}

int StructuredData::ArrayIndex(int i, int j, int k) {
    if(i<0 || i>=m_i || j<0 || j>=m_j || k<0 || k>=m_k) {
        printf("error: ilegal index (%d,%d,%d) for array size (%d,%d,%d)\n", i, j, k, m_i, m_j, m_k);
    }
    return i + j*m_i + k*m_i*m_j;
}



int StructuredData::OutputTec360(std::string filename)
{
    int isdouble = 1;
    std::vector<void*> data;
    for(int i=0; i<m_x.size(); ++i) {
        data.push_back((void*) (m_x[i].data()) );
    }
    for(int i=0; i<m_phys.size(); ++i) {
        data.push_back((void*) (m_phys[i].data()) );
    }
    ::OutputTec360(filename, m_varList, m_i, m_j, m_k, data, isdouble);
    return m_x.size() + m_phys.size();
}