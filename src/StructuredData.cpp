#include<string>
#include<vector>
#include<cmath>
#include "StructuredData.h"
#include "Dataprocessing.h"
#include "FileIO.h"

StructuredData::StructuredData(const std::vector<int> &N, const std::vector<double> &range) {
    m_N = N;
    m_range = range;
    ReSetNp();
    ReSetDx();
    GenPoints();
}

StructuredData::StructuredData(){
    m_N.clear();
    m_range.clear();
    ReSetNp();
    ReSetDx();
    GenPoints();
}

int StructuredData::GenPoints() {
    m_x.resize(m_N.size());
    for(int i=0; i<m_N.size(); ++i) {
        m_x[i].resize(m_Np);
    }
    if(m_Np==0) {
        return 0;
    }
    std::vector<int> N = m_N;
    N.push_back(1);
    N.push_back(1);
    N.push_back(1);
    for(int k=0; k<N[2]; ++k) {
        int tmp2 = k * N[0] * N[1];
        for(int j=0; j<N[1]; ++j) {
            int tmp1 = j * N[0];
            for(int i=0; i<N[0]; ++i) {
                int index = i + tmp1 + tmp2;
                std::vector<int> tmpi = {i,j,k};
                for(int d=0; d<m_N.size(); ++d) {
                    if(m_N[d]>1) {
                        m_x[d][index] = m_range[2*d] + tmpi[d]*m_dx[d];
                    } else {
                        m_x[d][index] = m_range[2*d];
                    }
                }
            }
        }
    }
    if(m_N.size()>0) {
        m_vars.push_back("x");
    }
    if(m_N.size()>1) {
        m_vars.push_back("y");
    }
    if(m_N.size()>2) {
        m_vars.push_back("z");
    }
    return m_Np;
}

int StructuredData::GetTotPoints() {
    return m_Np;
}

int StructuredData::GetNumPhys() {
    return m_phys.size();
}

int StructuredData::ReSetNp() {
    if(m_N.size()==0) {
        m_Np = 0;
    } else {
        m_Np = 1.;
        for(int i=0; i<m_N.size(); ++i) {
            m_Np *= m_N[i];
        }
    }
    return m_Np;
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
    return sum/m_Np;
}

double StructuredData::GetPhysValue(int f, int i) {
    return m_phys[f][i];
}

double StructuredData::GetCoordValue(int f, int i) {
    return m_x[f][i];
}

int StructuredData::ReSetDx() {
    m_dx.resize(m_N.size());
    for(int i=0; i<m_N.size(); ++i) {
        if(m_N[i]>1) {
            m_dx[i] = (m_range[2*i+1] - m_range[2*i])/(m_N[i]-1);
        } else {
            m_dx[i] = std::nan("1");
        }
    }
    return m_N.size();
}

int StructuredData::OutputData(std::string filename, const bool info) {
    std::clock_t c_start = std::clock();
    int isdouble = TEC360USEDOUBLE;
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
    if(info) {
        std::clock_t c_end = std::clock();
        double time_elapsed_ms = (c_end-c_start) * 1. / CLOCKS_PER_SEC;
        printf("output file %s, cpu time %fs\n", filename.c_str(), time_elapsed_ms);
    }
    return m_x.size() + m_phys.size();
}

int StructuredData::InputData(std::string filename, const bool info) {
    std::clock_t c_start = std::clock();
    std::vector<std::vector<double> > data;
    int isdouble;
    std::string ext = filename.substr(filename.size()-4, 4);
    if(0 == ext.compare(".plt")) {
        InputTec360_binary(filename, m_vars, m_N, m_range, data, isdouble);
    } else if(0 == ext.compare(".csv")) {
        InputCSV(filename, m_vars, m_N, m_range, data, isdouble);
    } else {
        printf("error: unsupported tecplot file type %s\n", filename.c_str());
        return -1;
    }
    ReSetNp();
    ReSetDx();
    //copy data
    m_x.clear();
    m_phys.clear();
    for(int i=0; i<m_vars.size(); ++i) {
        if(i<3) {
            m_x.push_back(data[i]);
        } else {
            m_phys.push_back(data[i]);
        }
    }
    if(info) {
        std::clock_t c_end = std::clock();
        double time_elapsed_ms = (c_end-c_start) * 1. / CLOCKS_PER_SEC;
        printf("Read file %s, cpu time %fs\n", filename.c_str(), time_elapsed_ms);
    }
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

int StructuredData::ExtractPlane(const std::vector<double> &data, std::pair<int, int> plane,
                                 std::vector<int> & N, std::vector<double> &odata) {
    N.resize(2);
    int dir = plane.first;
    int N01 = m_N[0] * m_N[1];
    if(dir%3==0) {
        N[0] = m_N[1];
        N[1] = m_N[2];
        odata.resize(N[0]*N[1]);
        int count = 0;
        for(int k=plane.second; k<m_Np; k+=N01) {
            int jmax = N01 + k;
            for(int j=k; j<jmax; j+=m_N[0]) {
                odata[count++] = data[j];
            }
        }
    }
    if(dir%3==1) {
        N[0] = m_N[2];
        N[1] = m_N[0];
        odata.resize(N[0]*N[1]);
        int count = 0;
        int tmp1 = plane.second * m_N[0];
        for(int i=0; i<m_N[0]; ++i) {
            for(int k= i + tmp1; k<m_Np; k+=N01) {
                odata[count++] = data[k];
            }
        }
    }
    if(dir%3==2) {
        N[0] = m_N[0];
        N[1] = m_N[1];
        odata.resize(N[0]*N[1]);
        int count = 0;
        int tmp2 = plane.second * N01;
        int imax = tmp2 + N01;
        for(int i=tmp2; i<imax; ++i) {
            odata[count++] = data[i];
        }
    }
    return N[0] * N[1];
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
    int res = der.Diff(N, u, du, m_dx[dir], order);
    if(dir) {
        std::vector<int> tN = N;
        ShiftIndex<double>(N, u, dir);
        ShiftIndex<double>(tN, du, dir);
    }
    return res;
}

int StructuredData::CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &rawNe,
                                    const std::vector<int> &skip, const StructuredData & big) {
    if(Ns.size()<big.m_N.size()) {
        printf("error: subdomain has a different dimension\n");
        return -1;
    }
    if(big.m_Np==0) {
        printf("error: origial domain has no element\n");
        clear();
        return -1;
    }
    std::vector<int> Ne = rawNe;
    m_N.resize(big.m_N.size());
    m_range.resize(big.m_range.size());
    for(int i=0; i<big.m_N.size(); ++i) {
        if(Ne[i] - Ns[i]<=0) {
            printf("error: subdomain have zero elements in %d direction\n", i);
        }
        m_N[i] = (Ne[i] - Ns[i])/skip[i];
        if((Ne[i] - Ns[i])%skip[i]) {
            m_N[i] += 1;
        }
        Ne[i] = Ns[i] + skip[i] * (m_N[i] - 1);
        m_range[2*i] = big.m_range[2*i] + big.m_dx[i] * Ns[i];
        m_range[2*i + 1] = big.m_range[2*i] + big.m_dx[i] * Ne[i];
    }
    ReSetNp();
    ReSetDx();
    m_vars = big.m_vars;
    //copy data
    m_x.resize(big.m_x.size());
    m_phys.resize(big.m_phys.size());
    for(int i=0; i<m_x.size(); ++i) {
        m_x[i].resize(m_Np);
    }
    for(int i=0; i<m_phys.size(); ++i) {
        m_phys[i].resize(m_Np);
    }
    //copy data
    if(Ns.size()==1) {
        int count = 0;
        for(int i=Ns[0]; i<=Ne[0]; i+=skip[0]) {
            for(int p=0; p<big.m_x.size(); ++p) {
                m_x[p][count] = big.m_x[p][i];
            }
            for(int p=0; p<big.m_phys.size(); ++p) {
                m_phys[p][count] = big.m_phys[p][i];
            }
            ++count;
        }
    } else if(Ns.size()==2) {
        int count = 0;
        int jstart = Ns[1] * big.m_N[0];
        int jmax = Ne[1] * big.m_N[0];
        int jskip = skip[1] * big.m_N[0];
        for(int j=jstart; j<=jmax; j+=jskip) {
            int imax = j + Ne[0];
            for(int i=Ns[0]+j; i<=imax; i+=skip[0]) {
                for(int p=0; p<big.m_x.size(); ++p) {
                    m_x[p][count] = big.m_x[p][i];
                }
                for(int p=0; p<big.m_phys.size(); ++p) {
                    m_phys[p][count] = big.m_phys[p][i];
                }
                ++count;
            }
        }
    } else if(Ns.size()==3) {
        int count = 0;
        int N01 = big.m_N[0] * big.m_N[1];
        int kstart = Ns[2] * N01;
        int kmax = Ne[2] * N01;
        int kskip = skip[2] * N01;
        for(int k=kstart; k<=kmax; k+=kskip) {
            int jstart = Ns[1] * big.m_N[0];
            int jmax = Ne[1] * big.m_N[0];
            int jskip = skip[1] * big.m_N[0];
            for(int j=jstart; j<=jmax; j+=jskip) {
                int imax = k + j + Ne[0];
                for(int i=Ns[0]+j + k; i<=imax; i+=skip[0]) {
                    for(int p=0; p<big.m_x.size(); ++p) {
                        m_x[p][count] = big.m_x[p][i];
                    }
                    for(int p=0; p<big.m_phys.size(); ++p) {
                        m_phys[p][count] = big.m_phys[p][i];
                    }
                    ++count;
                }
            }
        }
    } else {
        printf("error unsupported dimension %d\n", (int)Ns.size());
        return 0;
    }
    return m_Np;
}

void StructuredData::clear() {
    m_x.clear();
    m_phys.clear();
    m_vars.clear();
    m_N.clear();
    m_range.clear();
    m_dx.clear();
    m_Np = 0;
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
    return field.size() * m_Np;
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
    return field.size() * m_Np;
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
    return field.size() * m_Np;
}