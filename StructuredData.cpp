#include<string>
#include<vector>
#include<cmath>
#include "StructuredData.h"
#include "Dataprocessing.h"
#include "FileIO.h"

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
    m_dx = std::vector<double>(3, 1.);
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
    m_vars.push_back("x");
    m_vars.push_back("y");
    m_vars.push_back("z");
    GenPoints();
}

int StructuredData::GetTotPoints() {
    return m_Np;
}

int StructuredData::GetNumPhys() {
    return m_phys.size();
}

int StructuredData::AddPhysics(std::string var, void * func) {
    double(*physfun)(std::vector<double>) = (double(*)(std::vector<double>))func;
    m_vars.push_back(var);
    m_phys.push_back(std::vector<double>(m_Np));
    for(int i=0; i<m_Np; ++i) {
        std::vector<double> p;
        for(int k=0; k<m_x.size(); ++k) {
            p.push_back(m_x[k][i]);
        }
        for(int k=0; k<m_phys.size(); ++k) {
            p.push_back(m_phys[k][i]);
        }
        m_phys[m_phys.size()-1][i] = physfun(p);
    }
    return m_Np;
}

int StructuredData::AddPhysics(std::string var, const std::vector<double> &data) {
    m_vars.push_back(var);
    m_phys.push_back(data);
    if(data.size() != m_Np) {
        m_phys[m_phys.size()-1].resize(m_Np);
    }
    return data.size();
}

double StructuredData::GetPhysNorm(int f, int p) {
    double sum = 0.;
    for(int i=0; i<m_Np; ++i) {
        if(p>0) {
            sum += std::pow(std::fabs(m_phys[f][i]), p);
        } else {
            sum = std::max(sum, std::fabs(m_phys[f][i]));
        }
    }
    return sum;
}

double StructuredData::GetPhysValue(int f, int i) {
    return m_phys[f][i];
}

double StructuredData::GetCoordValue(int f, int i) {
    return m_x[f][i];
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

int StructuredData::OutputData(std::string filename) {
    int isdouble = 1;
    std::vector<std::vector<double> > data;
    for(int i=0; i<m_x.size(); ++i) {
        data.push_back(m_x[i]);
    }
    for(int i=0; i<m_phys.size(); ++i) {
        data.push_back(m_phys[i]);
    }
    std::string ext = filename.substr(filename.size()-4, 4);
    if(0 == ext.compare(".dat")) {
        OutputTec360_ascii(filename, m_vars, m_N, data, isdouble);
    } else if(0 == ext.compare(".plt")) {
        OutputTec360_binary(filename, m_vars, m_N, data, isdouble);
    } else if(0 == ext.compare(".csv")) {
        OutputCSV(filename, m_vars, m_N, data);
    } else {
        printf("error: unsupported tecplot file type %s\n", filename.c_str());
        return -1;
    }
    printf("output file %s\n", filename.c_str());
    return m_x.size() + m_phys.size();
}

int StructuredData::InputData(std::string filename) {
    std::vector<int> N;
    std::vector<std::vector<double> > data;
    int isdouble;
    std::string ext = filename.substr(filename.size()-4, 4);
    if(0 == ext.compare(".plt")) {
        InputTec360_binary(filename, m_vars, N, data, isdouble);
    } else if(0 == ext.compare(".csv")) {
        InputCSV(filename, m_vars, m_N, data, isdouble);
    } else {
        printf("error: unsupported tecplot file type %s\n", filename.c_str());
        return -1;
    }
    m_x.clear();
    m_phys.clear();
    for(int i=0; i<m_vars.size(); ++i) {
        if(i<3) {
            m_x.push_back(data[i]);
        } else {
            m_phys.push_back(data[i]);
        }
    }
    printf("Read file %s\n", filename.c_str());
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

int StructuredData::ExtractPlane(const std::vector<double> &data, std::pair<int, int> plane, std::vector<int> & N, std::vector<double> &odata) {
    odata.clear();
    N.resize(2);
    int dir = plane.first;
    if(dir%3==0) {
        N[0] = m_N[1];
        N[1] = m_N[2];
        for(int k=0; k<m_N[2]; ++k) {
            for(int j=0; j<m_N[1]; ++j) {
                std::vector<int> ind = {plane.second, j, k};
                odata.push_back(data[Index(m_N, ind)]);
            }
        }
    }
    if(dir%3==1) {
        N[0] = m_N[2];
        N[1] = m_N[0];
        for(int i=0; i<m_N[0]; ++i) {
            for(int k=0; k<m_N[2]; ++k) {
                std::vector<int> ind = {i, plane.second, k};
                odata.push_back(data[Index(m_N, ind)]);
            }
        }
    }
    if(dir%3==2) {
        N[0] = m_N[0];
        N[1] = m_N[1];
        for(int j=0; j<m_N[1]; ++j) {
            for(int i=0; i<m_N[0]; ++i) {
                std::vector<int> ind = {i, j, plane.second};
                odata.push_back(data[Index(m_N, ind)]);
            }
        }
    }
}

int StructuredData::Diff(std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &du, int dir, int order) {
    if(m_N[dir]==0) {
        printf("error: cannot calculate finite difference in 1 data, dir %d\n", dir);
        return 0;
    }
    Derivative der;
    std::vector<int> N = m_N;
    if(dir) {
        ShiftIndex(N, u, -dir);
    }
    der.Diff(N, u, du, m_dx[dir], order);
    if(dir) {
        std::vector<int> tN = N;
        ShiftIndex<double>(N, u, dir);
        ShiftIndex<double>(tN, du, dir);
    }
}

int StructuredData::Smoothing(double sigma, std::vector<int> &field, bool inplace) {
    std::vector<std::vector<double> > data;
    for(int i=0; i<field.size(); ++i) {
        data.push_back(m_phys[field[i]]);
    }
    Smoothing(sigma, data);
    if(!inplace) {
        for(int i=0; i<field.size(); ++i) {
            std::string var("S(");
            var += m_vars[i+3] + ")";
            m_vars.push_back(var);
            m_phys.push_back(data[i]);
        }
    } else {
        for(int i=0; i<field.size(); ++i) {
            m_phys[field[i]] = data[i];
        }
    }
}
int StructuredData::Diff(std::vector<int > &field, int dir, int order) {
    std::vector<std::vector<double> > u;
    std::vector<std::vector<double> > du;
    for(int i=0; i<field.size(); ++i) {
        u.push_back(m_phys[field[i]]);
        du.push_back(m_phys[field[i]]);
    }
    Diff(u, du, dir, order);
    for(int i=0; i<field.size(); ++i) {
        std::string var = m_vars[i+3] + "_";
        var += 'x' + dir;
        m_vars.push_back(var);
        m_phys.push_back(du[i]);
    }
}

int StructuredData::MaskBoundary(double sigma, std::vector<int> &field, std::map<int, double> def) {
    std::vector<int> numbers(3, -1);
    std::vector<int> index(3);
    for(int i=0; i<3; ++i) {
        if(m_N[i]>1) {
            numbers[i] = myRound<double>(sigma/m_dx[i]);
        }
    }
    for(int i=0; i<m_Np; ++i) {
        invIndex(m_N, i, index);
        for(int k=0; k<3; ++k) {
            if(m_N[k]==1) continue;
            if(index[k]<numbers[k] || m_N[k]-1-index[k]<numbers[k]) {
                for(auto f=field.begin(); f!=field.end(); ++f) {
                    m_phys[*f][i] = def[*f];
                }
                break;
            }
        }
    }
}