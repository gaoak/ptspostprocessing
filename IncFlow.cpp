#include<tuple>
#include<set>
#include "IncFlow.h"
#include "Util.h"

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range) 
        :StructuredData(N, range) {
}

int IncFlow::CalculateVorticity(int order) {
    std::vector<std::vector<double> > u, ux(3), uy(3), uz(3);
    u.push_back(m_phys[0]);
    u.push_back(m_phys[1]);
    u.push_back(m_phys[2]);
    for(int i=0; i<3; ++i) {
        ux[i].resize(m_Np);
        uy[i].resize(m_Np);
        uz[i].resize(m_Np);
    }
    Diff(u, ux, 0, order);
    Diff(u, uy, 1, order);
    Diff(u, uz, 2, order);
    int id;
    //W_x
    m_varList += ",W_x,W_y,W_z";
    id = m_phys.size();
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    for(int i=0; i<m_Np; ++i) {
        m_phys[id  ][i] = uy[2][i] - uz[1][i];
        m_phys[id+1][i] = uz[0][i] - ux[2][i];
        m_phys[id+2][i] = ux[1][i] - uy[0][i];
    }

    // Q Criterion
    m_varList += ",Q";
    id = m_phys.size();
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    for(int i=0; i<m_Np; ++i) {
        m_phys[id][i] = -0.5 * (
            ux[0][i]*ux[0][i] + ux[1][i]*uy[0][i] + ux[2][i]*uz[0][i] +
            uy[0][i]*ux[1][i] + uy[1][i]*uy[1][i] + uy[2][i]*uz[1][i] +
            uz[0][i]*ux[2][i] + uz[1][i]*uy[2][i] + uz[2][i]*uz[2][i]
        );
    }
}

std::pair<int, int> IncFlow::GetDirection(double sigma, std::vector<std::vector<double> > & cores) {
    std::pair<int, int> res;
    res.first = -1;
    res.second = 0;
    if(cores.size()<2) return res;
    std::vector<double> ave(3, 0.);
    int count = 0;
    int plast = cores.size() -1;
    int p = plast-1;
    for(; p>=0 && Distance(cores[plast], cores[p]) <= sigma; --p) ;
    if(p<0) p = 0;
    for(int i=p; i<plast; ++i) {
        double dist = Distance(cores[p], cores[plast]);
        if(dist<1E-9) continue;
        dist = 1./dist;
        std::vector<double> tmp(3, 0.);
        AddVect(dist, cores[plast], -dist, cores[p], tmp);
        AddVect(1., tmp, 1., ave, ave);
    }
    
    double tmp = std::fabs(ave[0]);
    int itmp = 0;
    if(tmp < std::fabs(ave[1])) {
        itmp = 1;
        tmp = std::fabs(ave[1]);
    }
    if(tmp < std::fabs(ave[2])) {
        itmp = 2;
    }
    res.first = itmp;
    int sign = 1;
    if(ave[itmp] < 0) sign = -1;
    res.second = sign;
    printf("direction %d, %d\n", res.first, res.second);
    return res;
}
int IncFlow::ExtractCore(int f, double sigma, std::vector<std::vector<double> > & cores, int dir) {
    cores.clear();
    std::vector<std::vector<double> > odata;
    odata.push_back(m_phys[f]);
    Smoothing(sigma, odata);

    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<double> range = m_range;
    std::vector<int> padding(3);
    double paddingsize = 3.*sigma;
    for(int i=0; i<3; ++i) {
        if(N[i]>1) {
            padding[i] = std::max(myRound<double>(paddingsize/dx[i]), 3);
        } else {
            padding[i] = 0;
        }
    }

    int Trymax = m_N[0] + m_N[1] + m_N[2];
    std::set<std::pair<int, int> > searched;
    std::vector<int> intcenter;
    std::pair<int, int> plane(dir, 0);
    std::vector<int> planeN;
    std::vector<double> planedata;
    int count = 0;
    while(count<Trymax) {
        if(plane.second<0 || plane.second>=N[plane.first]) {
            break;
        }
        dir = plane.first;
        printf("search plane %d, %d\n", dir, plane.second);
        ExtractPlane(odata[0], plane, planeN, planedata);
        ShiftArray<double>(dx, 2-dir);
        ShiftArray<double>(range, 2*(2-dir));
        ShiftArray<int>(N, 2-dir);
        ShiftArray<int>(padding, 2-dir);
        ShiftArray<int>(intcenter, 2-dir);
        std::vector<double> newcenter;
        std::vector<int> Nslice = {N[0], N[1]};
        ExtractCore2Dplane(Nslice, intcenter, padding, planedata, newcenter);
        std::vector<double> physcenter = newcenter;
        intcenter.resize(2);
        for(int k=0; k<2; ++k) {
            physcenter[k] = range[2*k] + dx[k] * physcenter[k];
            intcenter[k] = myRound<double>(newcenter[k]);
        }
        physcenter.push_back(range[4] + (plane.second)*dx[2]);
        intcenter.push_back(plane.second);
        ShiftArray<double>(dx, dir-2);
        ShiftArray<double>(physcenter, dir-2);
        ShiftArray<double>(range, 2*(dir-2));
        ShiftArray<int>(N, dir-2);
        ShiftArray<int>(padding, dir-2);
        ShiftArray<int>(intcenter, dir-2);
        ShiftArray<int>(padding, dir-2);
        cores.push_back(physcenter);
        searched.insert(plane);
        printf("location %d, %d, %d\n", intcenter[0], intcenter[1], intcenter[2]);
        if(count) {
            std::pair<int, int> incplane = GetDirection(paddingsize, cores);
            plane.first = incplane.first;
            plane.second = intcenter[incplane.first] + incplane.second;
        } else {
            plane.second += 1;
        }
        printf("next plane %d, %d\n", plane.first, plane.second);
        if(searched.find(plane)!=searched.end()) {
            break;
        }
        ++count;
    }
    return count;
}

int IncFlow::ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
    const std::vector<int> &padding, std::vector<double> &data, std::vector<double> &core, bool ismin) {
    if(!ismin) {
        for(int i=0; i<N[0]*N[1]; ++i) {
            data[i] = -data[i];
        }
    }
    if(initial.size()>=2) {
        DoMask<double>(N, data.data(), initial, padding, false);
    }
    int Np = N[0] * N[1];
    int imin = FindMin<double>(Np, data.data());
    double pmin = data[imin];
    double thresh = 0.95 * pmin;
    std::vector<double> center;
    DoMaskShift<double>(Np, thresh, -1, data.data());
    core.clear();
    WeightedCenter<double>(N, data.data(), core);
    return core.size();
}
