#include<string>
#include<vector>
#include<cmath>
#include<set>
#include "StructuredData.h"
#include "Dataprocessing.h"
#include "FileIO.h"

CoordSystem::CoordSystem() {
    m_o = {0.,0.,0.};
    m_e.push_back({1., 0., 0.});
    m_e.push_back({0., 1., 0.});
    m_e.push_back({0., 0., 1.});
}

CoordSystem::CoordSystem(const std::vector<double> &o)
    : CoordSystem() {
    for(int i=0; i<(int)o.size() && i<3; ++i) {
        m_o[i] = o[i];
    }
}

CoordSystem::CoordSystem(const std::vector<double> &o,
    const std::vector<std::vector<double> > & axis)
    : CoordSystem() {
    for(int i=0; i<(int)o.size() && i<3; ++i) {
        m_o[i] = o[i];
    }
    for(int i=0; i<(int)axis.size() && i<3; ++i) {
        for(int j=0; j<(int)axis[i].size() && j<3; ++j) {
            m_e[i][j] = axis[i][j];
        }
        NormalizeVect(m_e[i]);
    }
}

CoordSystem& CoordSystem::operator=(const CoordSystem & coord) {
    m_o = coord.m_o;
    m_e = coord.m_e;
    return *this;
}

CoordSystem::CoordSystem(const CoordSystem & coord) {
    *this = coord;
}

void CoordSystem::ToPhysCoord(std::vector<double> &x) const {
    std::vector<double> tmp = x;
    for(int i=tmp.size(); i<3; ++i) {
        tmp.push_back(0.);
    }
    for(int i=0; i<(int)x.size(); ++i) {
        x[i] = m_o[i] + tmp[0] * m_e[0][i] + tmp[1] * m_e[1][i] + tmp[2] * m_e[2][i];
    }
}

void CoordSystem::ToCompCoord(std::vector<double> &x) const {
    std::vector<double> tmp(3, 0.);
    int N = std::min(3, (int)x.size());
    for(int i=0; i<N; ++i) {
        tmp[i] = x[i] - m_o[i];
    }
    for(int i=0; i<(int)x.size() && i<3; ++i) {
        x[i] = DotVect(m_e[i], tmp);
    }
}

StructuredData::StructuredData(const std::vector<int> &N, const std::vector<double> &range)
    : StructuredData(N, range, std::vector<std::vector<double> >()) {
}

StructuredData::StructuredData(const std::vector<int> &N, const std::vector<double> &range,
                   const std::vector<std::vector<double> > & axis) {
    std::vector<double> tmp;
    for(int i=0; i<(int)range.size(); i+=2) {
        tmp.push_back(range[i]);
    }
    m_axis = CoordSystem(tmp, axis);
    m_N = N;
    ReSetNp();
    GenPoints(range);
    m_velocityDim = -1;
}

StructuredData::StructuredData()
    : m_axis() {
    m_velocityDim = -1;
}

StructuredData& StructuredData::operator=(const StructuredData & data) {
    m_x = data.m_x;
    m_phys = data.m_phys;
    m_N = data.m_N;
    m_vars = data.m_vars;
    m_Np = data.m_Np;
    m_axis = data.m_axis;
    m_dx = data.m_dx;
    return *this;
}

StructuredData::StructuredData(const StructuredData & data) {
    *this = data;
}

int StructuredData::GenPoints(const std::vector<double> &range) {
    m_x.resize(m_N.size());
    for(int i=0; i<(int)m_N.size(); ++i) {
        m_x[i].resize(m_Np);
    }
    if(m_Np==0) {
        return 0;
    }
    m_dx.resize(m_N.size());
    for(int i=0; i<(int)m_N.size(); ++i) {
        if(m_N[i]>1) {
            m_dx[i] = range[2*i+1] /(m_N[i]-1);
        } else if(m_N[i]==1) {
            m_dx[i] = 1.;
        } else {
            m_dx[i] = std::nan("1");
        }
    }
    std::vector<int> N = m_N;
    std::vector<double> x(3, 0.);
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
                for(int d=0; d<(int)m_N.size(); ++d) {
                    if(m_N[d]>1) {
                        x[d] = tmpi[d]*m_dx[d];
                    } else {
                        x[d] = 0.;
                    }
                }
                m_axis.ToPhysCoord(x);
                for(int d=0; d<(int)m_N.size(); ++d) {
                    m_x[d][index] = x[d];
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

int StructuredData::RemovePhysics(int i) {
    m_phys.erase(m_phys.begin() + i);
    m_vars.erase(m_vars.begin() + i);
    return -1;
}

int StructuredData::ReSetNp() {
    if(m_N.size()==0) {
        m_Np = 0;
    } else {
        m_Np = 1.;
        for(int i=0; i<(int)m_N.size(); ++i) {
            m_Np *= m_N[i];
        }
    }
    return m_Np;
}

int StructuredData::AddPhysics(std::string var, void func()) {
    double(*physfun)(std::vector<double>) = (double(*)(std::vector<double>))func;
    m_vars.push_back(var);
    int nphys = m_phys.size();
    m_phys.resize(nphys + 1); m_phys[nphys].resize(m_Np);
    for(int i=0; i<m_Np; ++i) {
        std::vector<double> p;
        for(int k=0; k<(int)m_x.size(); ++k) {
            p.push_back(m_x[k][i]);
        }
        for(int k=0; k<(int)m_phys.size(); ++k) {
            p.push_back(m_phys[k][i]);
        }
        m_phys[(int)m_phys.size()-1][i] = physfun(p);
    }
    return m_Np;
}

int StructuredData::AddPhysics(std::string var, const std::vector<double> &data) {
    m_vars.push_back(var);
    m_phys.push_back(data);
    if((int)data.size() != m_Np) {
        m_phys[(int)m_phys.size()-1].resize(m_Np);
    }
    return (int)data.size();
}

int StructuredData::HasField(std::string var) {
    for(size_t i=0; i<m_vars.size(); ++i) {
        if (var==m_vars[i]) {
            return 1;
        }
    }
    return 0;
}

double StructuredData::GetPhysNorm(int f, int p) {
    double sum = NormVect(m_phys[f], p);
    return sum/m_Np;
}

int StructuredData::GetPhysID(std::string v) {
    for(size_t i=0; i<m_vars.size(); ++i) {
        if(v==m_vars[i]) {
            return (int)i - (int)m_x.size();
        }
    }
    printf("Field %s not found.\n", v.c_str());
    return -1;
}

int StructuredData::OutputData(std::string filename, const bool info) {
    std::clock_t c_start = std::clock();
    int isdouble = TEC360USEDOUBLE;
    std::vector<std::vector<double> > data;
    for(size_t i=0; i<m_x.size(); ++i) {
        data.push_back(m_x[i]);
    }
    for(size_t i=0; i<m_phys.size(); ++i) {
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

int StructuredData::ResetAxis() {
    std::vector<double> e0, e1, e2;
    std::vector<double> x0 = {m_x[0][0], m_x[1][0], GetCoordValue(2, 0)};
    int i0 = 1;
    e0 = {m_x[0][i0]-x0[0], m_x[1][i0]-x0[1], GetCoordValue(2, i0)-x0[2]};
    if(m_N.size()>1 && m_N[1]>1) {
        int i1 = Index(m_N, {0, 1, 0});
        e1 = {m_x[0][i1]-x0[0], m_x[1][i1]-x0[1], GetCoordValue(2, i1)-x0[2]};
    } else {
        int imax = FindAbsMax(3, e0.data());
        e1 = {1.,1.,1.};
        e1[imax] = 0.;
        e1 = CrossVect(e0, e1);
    }
    if(m_N.size()>2 && m_N[2]>1) {
        int i2 = Index(m_N, {0, 0, 1});
        e2 = {m_x[0][i2]-x0[0], m_x[1][i2]-x0[1], GetCoordValue(2, i2)-x0[2]};
    } else {
        e2 = CrossVect(e0, e1);
    }
    std::vector<std::vector<double> > e;
    e.push_back(e0);
    e.push_back(e1);
    e.push_back(e2);
    m_axis = CoordSystem(x0, e);
    m_dx.resize(m_N.size());
    for(int i=0; i<(int)m_N.size(); ++i) {
        if(m_N[i]>1) {
            m_dx[i] = std::sqrt(DotVect(e[i], e[i]));
        } else if(m_N[i]==1) {
            m_dx[i] = 1.;
        } else {
            m_dx[i] = std::nan("1");
        }
    }
    return 0;
}

int StructuredData::ShuffleIndex(std::map<int, int> ReIndex, std::vector<int> dir,
        std::map<int, int> pm) {
    //only shuffle index, not change variable order
    //place ReIndex[i] to i
    int dim = m_x.size();
    std::map<int, int> RI;
    std::set<int> sec, fir;
    bool needshuff = false;
    int count = 0;
    for(auto it=ReIndex.begin(); it!=ReIndex.end(); ++it) {
        if(it->second<0 || it->second>=dim || sec.find(it->second)!=sec.end()) {
            continue;
        }
        if(count!=it->second) {
            needshuff = true;
        }
        fir.insert(count);
        sec.insert(it->second);
        RI[count++] = it->second;
    }
    if(!needshuff || fir!=sec || (int)fir.size()!=dim) {
        return 0;
    }
    std::vector<std::vector<double> > x = m_x;
    std::vector<std::vector<double> > phys = m_phys;
    std::vector<int> N = m_N;
    std::vector<std::string> vars = m_vars;
    N.push_back(1); N.push_back(1); N.push_back(1);
    std::vector<int> targ(3, 0);
    for(int d=0; d<dim; ++d) {
        m_N[RI[d]] = N[d];
    }
    for(int k=0; k<N[2]; ++k) {
        int tmp2 = k * N[0] * N[1];
        for(int j=0; j<N[1]; ++j) {
            int tmp1 = j * N[0] + tmp2;
            for(int i=0; i<N[0]; ++i) {
                int index = i + tmp1;
                std::vector<int> orig = {i, j, k};
                for(int d=0; d<dim; ++d) {
                    if(dir[d]<0) {
                        targ[RI[d]] = N[d]-1 - orig[d];
                    } else {
                        targ[RI[d]] = orig[d];
                    }
                }
                int newind = Index(m_N, targ);
                for(size_t d=0; d<m_x.size(); ++d) {
                    m_x[RI[d]][newind] = x[d][index];
                }
                for(size_t d=0; d<m_phys.size(); ++d) {
                    m_phys[pm[d]][newind] = phys[d][index];
                }
            }
        }
    }
    ResetAxis();
    return (m_x.size() + m_phys.size() ) * m_Np;
}

void StructuredData::UpdateVelocityDimension() {
    m_velocityDim = 0;
    std::set<std::string> velocity;
    velocity.insert("u");
    velocity.insert("v");
    velocity.insert("w");
    for(size_t i=0; i<m_phys.size(); ++i) {
        if(velocity.find(m_vars[i+m_x.size()])!=velocity.end()) {
            ++m_velocityDim;
        }
    }
}

int StructuredData::InputData(std::string filename, const bool info) {
    std::clock_t c_start = std::clock();
    std::vector<std::vector<double> > data;
    int isdouble;
    std::string ext = filename.substr(filename.size()-4, 4);
    int ncoor = -1;
    std::map<int, int> vm;
    std::vector<std::string> vars;
    if(0 == ext.compare(".plt")) {
        ncoor = InputTec360_binary(filename, vars, m_N, data, isdouble, vm);
    } else if(0 == ext.compare(".csv")) {
        ncoor = InputCSV(filename, vars, m_N, data, isdouble, vm);
    } else {
        printf("error: unsupported tecplot file type %s\n", filename.c_str());
        return -1;
    }
    ReSetNp();
    //copy data
    m_x.clear();
    m_phys.resize((int)vars.size()-ncoor);
    m_vars.resize(vars.size());
    for(int i=0; i<(int)vars.size(); ++i) {
        m_vars[vm[i]] = vars[i];
        if(vm[i]<ncoor) {
            m_x.push_back(data[i]);
        } else {
            m_phys[vm[i]-ncoor] = data[i];
        }
    }
    //reset axis
    ResetAxis();
    std::vector<int> dir(m_N.size(), 1);
    std::map<int, int> pm;
    for(int i=0; i<(int)m_phys.size(); ++i) {
        pm[i] = i;
    }
    ShuffleIndex(vm, dir, pm);
    if(info) {
        std::clock_t c_end = std::clock();
        double time_elapsed_ms = (c_end-c_start) * 1. / CLOCKS_PER_SEC;
        printf("Read file %s, cpu time %fs\n", filename.c_str(), time_elapsed_ms);
    }
    UpdateVelocityDimension();
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
                                 const std::vector<int> &Range, std::vector<int> & N,
                                 std::vector<double> &odata) {
    //range = {imin, imax, jmin, jmax, kmin, kmax}
    std::vector<int> range = Range;
    for(int i=0; i<3; ++i) {
        if(plane.first==i) {
            continue;
        }
        if(range[2*i] < 0) {
            range[2*i] = 0;
        }
        if(range[2*i+1]>=m_N[i]) {
            range[2*i+1] = m_N[i] - 1;
        }
        if(range[2*i]>range[2*i+1]) {
            odata.clear();
            return 0;
        }
    }
    N.resize(2);
    int dir = plane.first;
    int N01 = m_N[0] * m_N[1];
    if(dir%3==0) {
        N[0] = range[3] - range[2] + 1;
        N[1] = range[5] - range[4] + 1;
        odata.resize(N[0]*N[1]);
        int count = 0;
        int kstart = range[4] * N01 + plane.second;
        int kend = range[5] * N01 + plane.second;
        int tmp2 = m_N[0] * range[2];
        int tmp3 = m_N[0] * range[3];
        for(int k=kstart; k<=kend; k+=N01) {
            int jmax = tmp3 + k;
            for(int j= tmp2 + k; j<=jmax; j+=m_N[0]) {
                odata[count++] = data[j];
            }
        }
    }
    if(dir%3==1) {
        N[0] = range[5] - range[4] + 1;
        N[1] = range[1] - range[0] + 1;
        odata.resize(N[0]*N[1]);
        int count = 0;
        int tmp4 = plane.second * m_N[0] + N01 * range[4];
        int tmp5 = plane.second * m_N[0] + N01 * range[5];
        for(int i=range[0]; i<=range[1]; ++i) {
            int kmax = i + tmp5;
            for(int k= i + tmp4; k<=kmax; k+=N01) {
                odata[count++] = data[k];
            }
        }
    }
    if(dir%3==2) {
        N[0] = range[1] - range[0] + 1;
        N[1] = range[3] - range[2] + 1;
        odata.resize(N[0]*N[1]);
        int count = 0;
        int tmp2 = plane.second * N01;
        int jstart = range[2] * m_N[0] + tmp2;
        int jend = range[3] * m_N[0] + tmp2;
        for(int j=jstart; j<=jend; j+=m_N[0]) {
            int imax = range[1] + j;
            for(int i=range[0] + j; i<=imax; ++i) {
                odata[count++] = data[i];
            }
        }
    }
    return N[0] * N[1];
}

int StructuredData::Diff(std::vector<std::vector<double> > &u, std::vector<std::vector<double> > &du, int dir, int order) {
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

int StructuredData::GetInterpDimension() const {
    int dim = m_N.size();
    if(dim == 3 && m_N[2]==1) {
        dim = 2;
    }
    if(dim==2 && m_N[1]==1) {
        dim = 1;
    }
    return dim;
 }

//Evaluate values of field[physics index, default value] at one point x
int StructuredData::InterpolatePoint(const std::vector<double> & x, std::map<int, double> field,
                                    std::map<int, double> &value) const {
    value = field;
    if(x.size() < m_N.size()) {
        printf("error incorrect input for interpolation\n");
        return 0;
    }
    int dim = GetInterpDimension();
    for(int i=0; i<dim; ++i) {
        if(m_N[i] <= 1) {
            printf("error cannot do interpolation on boundary\n");
        }
    }
    std::vector<int> index(m_N.size(), 0);
    std::vector<double> xi = x;
    m_axis.ToCompCoord(xi);
    double geomtol = 1E-12;
    for(int i=0; i<dim; ++i) {
        xi[i] /= m_dx[i];
        if(xi[i]>=m_N[i]-1.+geomtol || xi[i] < -geomtol) {
            return 0;
        }
        index[i] = xi[i];
        xi[i] -= index[i];
        if(index[i]==m_N[i]-1) {
            index[i] -= 1;
            xi[i] += 1.;
        }
    }
    std::vector<int> ind;
    if(dim>=1) {
        ind.push_back(Index(m_N, {index[0]+0,index[1]+0,index[2]+0}));
        ind.push_back(Index(m_N, {index[0]+1,index[1]+0,index[2]+0}));
    }
    if(dim>=2) {
        ind.push_back(Index(m_N, {index[0]+0,index[1]+1,index[2]+0}));
        ind.push_back(Index(m_N, {index[0]+1,index[1]+1,index[2]+0}));
    }
    if(dim>=3) {
        ind.push_back(Index(m_N, {index[0]+0,index[1]+0,index[2]+1}));
        ind.push_back(Index(m_N, {index[0]+1,index[1]+0,index[2]+1}));
        ind.push_back(Index(m_N, {index[0]+0,index[1]+1,index[2]+1}));
        ind.push_back(Index(m_N, {index[0]+1,index[1]+1,index[2]+1}));
    }
    std::vector<double> w;
    Interpolation::CalcWeight(xi, dim, w);
    std::vector<double> stencil(w.size());
    for(auto it=field.begin(); it!=field.end(); ++it) {
        int f = it->first;
        for(int i=0; i<(int)w.size(); ++i) {
            stencil[i] = m_phys[f][ind[i]];
        }
        value[f] = DotVect(stencil, w);
    }
    return value.size();
}

int StructuredData::CopyToSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                        const std::vector<int> &skip, std::map<int, double> &field,
                        StructuredData & small) {
    return small.CopyAsSubDomain(Ns, Ne, skip, field, *this);
}

std::vector<double> StructuredData::GetRange() {
    std::vector<double> res(6);
    for(size_t i=0; i<m_N.size(); ++i) {
        res[2*i] = m_axis.m_o[i];
        res[2*i+1] = m_dx[i] * (m_N[i] - 1);
    }
    return res;
}

int StructuredData::InterpolateFrom(const StructuredData & origin, std::map<int,double> rawfield) {
    if(m_N.size() != origin.m_N.size()) {
        printf("error dimension mismatch in interpolation\n");
        return 0;
    }
    int ncoor = origin.m_x.size();
    m_vars.resize(ncoor);
    std::map<int,double> field;
    for(auto it=rawfield.begin(); it!=rawfield.end(); ++it) {
        if(it->first<(int)origin.m_phys.size()) {
            field[it->first] = it->second;
            m_vars.push_back(origin.m_vars[ncoor + it->first]);
        }
    }
    m_phys.resize(field.size());
    for(size_t i=0; i<field.size(); ++i) {
        m_phys[i].resize(m_Np);
    }
    std::map<int, double> value;
    int count = 0;
    for(int i=0; i<m_Np; ++i) {
        std::vector<double> x;
        for(int d=0; d<(int)m_N.size(); ++d) {
            x.push_back(m_x[d][i]);
        }
        count += origin.InterpolatePoint(x, field, value);
        int v = 0;
        for(auto it=value.begin(); it!=value.end(); ++it) {
            m_phys[v++][i] = it->second;
        }
    }
    return count;
}

int StructuredData::CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &rawNe,
                                    const std::vector<int> &skip, std::map<int, double> &rawfield,
                                    const StructuredData & big) {
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
    for(int i=0; i<(int)big.m_N.size(); ++i) {
        if(Ne[i] - Ns[i]<=0) {
            printf("error: subdomain have zero elements in %d direction\n", i);
        }
        m_N[i] = (Ne[i] - Ns[i])/skip[i];
        if((Ne[i] - Ns[i])%skip[i]) {
            m_N[i] += 1;
        }
        Ne[i] = Ns[i] + skip[i] * (m_N[i] - 1);
        m_axis.m_o[i] = big.m_axis.m_o[i] + big.m_dx[i] * Ns[i];
    }
    m_axis.m_e = big.m_axis.m_e;
    ReSetNp();
    m_vars.clear();
    int ncoor = big.m_x.size();
    for(int i=0; i<ncoor; ++i) {
        m_vars.push_back(big.m_vars[i]);
    }
    std::map<int, double> field;
    for(auto it=rawfield.begin(); it!=rawfield.end(); ++it) {
        if(it->first < (int)big.m_phys.size()) {
            field[it->first] = it->second;
            m_vars.push_back(big.m_vars[ncoor + it->first]);
        }
    }
    //copy data
    std::map<int, int> fm;
    int countf = 0;
    for(auto it=field.begin(); it!=field.end(); ++it) {
        fm[countf++] = it->first;
    }
    m_x.resize(big.m_x.size());
    m_phys.resize(fm.size());
    for(size_t i=0; i<m_x.size(); ++i) {
        m_x[i].resize(m_Np);
    }
    for(size_t i=0; i<m_phys.size(); ++i) {
        m_phys[i].resize(m_Np);
    }
    //copy data
    if(Ns.size()==1) {
        int count = 0;
        for(int i=Ns[0]; i<=Ne[0]; i+=skip[0]) {
            for(int p=0; p<(int)big.m_x.size(); ++p) {
                m_x[p][count] = big.m_x[p][i];
            }
            for(int p=0; p<(int)m_phys.size(); ++p) {
                m_phys[p][count] = big.m_phys[fm[p]][i];
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
                for(int p=0; p<(int)big.m_x.size(); ++p) {
                    m_x[p][count] = big.m_x[p][i];
                }
                for(int p=0; p<(int)m_phys.size(); ++p) {
                    m_phys[p][count] = big.m_phys[fm[p]][i];
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
                    for(int p=0; p<(int)big.m_x.size(); ++p) {
                        m_x[p][count] = big.m_x[p][i];
                    }
                    for(int p=0; p<(int)m_phys.size(); ++p) {
                        m_phys[p][count] = big.m_phys[fm[p]][i];
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
    m_dx.clear();
    m_Np = 0;
    m_velocityDim = -1;
}

int StructuredData::Smoothing(double sigma, std::vector<int> &field, bool inplace) {
    std::vector<std::vector<double> > data;
    for(int i=0; i<(int)field.size(); ++i) {
        data.push_back(m_phys[field[i]]);
    }
    Smoothing(sigma, data);
    if(!inplace) {
        for(int i=0; i<(int)field.size(); ++i) {
            std::string var("S(");
            var += m_vars[i+3] + ")";
            m_vars.push_back(var);
            m_phys.push_back(data[i]);
        }
    } else {
        for(int i=0; i<(int)field.size(); ++i) {
            m_phys[field[i]] = data[i];
        }
    }
    return field.size() * m_Np;
}

int StructuredData::Diff(std::vector<int > &field, int dir, int order) {
    std::vector<std::vector<double> > u;
    std::vector<std::vector<double> > du;
    for(int i=0; i<(int)field.size(); ++i) {
        u.push_back(m_phys[field[i]]);
        du.push_back(m_phys[field[i]]);
    }
    Diff(u, du, dir, order);
    for(int i=0; i<(int)field.size(); ++i) {
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

/***
 * volume integral of data over the whole domain, using rectangle rule
***/
double StructuredData::Integrate(const std::vector<double> &data) {
    double sum = 0.;
    for(size_t j=0; j<data.size();++j) {
        sum += data[j];
    }
    for(size_t j=0; j<m_dx.size(); ++j) {
        if(m_N[j]==1) continue;
        sum *= m_dx[j];
    }
    return sum;
}