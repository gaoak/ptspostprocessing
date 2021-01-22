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

void IncFlow::SetBody(std::string bodyname, std::vector<double> param) {
    m_body = Body(bodyname, param[0]);
}

std::pair<int, int> IncFlow::GetProceedDirection(const std::vector<double> &vor, double sign) {
    std::pair<int, int> res;
    res.first = FindAbsMax(3, vor.data());
    if(sign * vor[res.first] > 0) {
        res.second = 1;
    } else {
        res.second = -1;
    }
    printf("direction %d, %d\n", res.first, res.second);
    return res;
}

int IncFlow::ExtractCore(double sigma, std::vector<std::vector<double> > & cores,
    std::vector<double> & radius, std::vector<double> &circulation, int dir, int f) {
    cores.clear();
    std::vector<std::vector<double> > odata;
    //u, v, w, p, W_x, W_y, W_z, Q
    //0, 1, 2, 3, 4,   5,   6,   7
    odata.push_back(m_phys[4]);
    odata.push_back(m_phys[5]);
    odata.push_back(m_phys[6]);
    odata.push_back(m_phys[std::abs(f)]);
    if(f<0) {
        for(int i=0; i<m_Np; ++i) {
            odata[3][i] = -odata[3][i];
        }
    }
    Smoothing(sigma, odata);
    double vorticitysign = 1.;
    if(dir<0) {
        dir = -dir;
        vorticitysign = -1.;
    }

    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<double> range = m_range;
    std::vector<int> padding(3);
    double paddingsize = 6.*sigma;
    for(int i=0; i<3; ++i) {
        if(N[i]>1) {
            padding[i] = myRound<double>(paddingsize/dx[i]);
        } else {
            padding[i] = 0;
        }
    }

    int Trymax = m_N[0] + m_N[1] + m_N[2];
    std::set<int> searched;
    std::vector<int> intcenter;
    std::pair<int, int> plane(dir, 0);
    std::vector<int> planeN;
    std::vector<double> planedata;
    int count = 0;
    double tmpradius, tmpcirculation;
    std::vector<double> planevorticity;
    std::pair<int, int> incplane;
    incplane.first = dir;
    incplane.second = 1;
    while(count<Trymax) {
        if(plane.second<0 || plane.second>=N[plane.first]) {
            break;
        }
        dir = plane.first;
        printf("search plane %d, %d\n", dir, plane.second);
        ExtractPlane(odata[3], plane, planeN, planedata);
        ExtractPlane(odata[dir], plane, planeN, planevorticity);

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
        ExtractVortexParam2Dplane(Nslice, dx, intcenter, planevorticity, tmpradius, tmpcirculation);

        ShiftArray<double>(dx, dir-2);
        ShiftArray<double>(physcenter, dir-2);
        ShiftArray<double>(range, 2*(dir-2));
        ShiftArray<int>(N, dir-2);
        ShiftArray<int>(padding, dir-2);
        ShiftArray<int>(intcenter, dir-2);
        ShiftArray<int>(padding, dir-2);

        cores.push_back(physcenter);
        radius.push_back(tmpradius);
        circulation.push_back(std::fabs(tmpcirculation));
        int centerindex = Index(m_N, intcenter);
        if(searched.find(centerindex)!=searched.end()) {
            searched.insert(centerindex);
            break;
        }
        searched.insert(centerindex);
        printf("location %d, %d, %d\n", intcenter[0], intcenter[1], intcenter[2]);
        if(m_body.IsInBody(physcenter, m_dx[1])) {
            break;
        }

        std::vector<double> vor = {odata[0][centerindex], odata[1][centerindex], odata[2][centerindex]};
        incplane = GetProceedDirection(vor, vorticitysign);
        plane.first = incplane.first;
        plane.second = intcenter[incplane.first] + incplane.second;
        printf("next plane %d, %d\n", plane.first, plane.second);
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
    double thresh = 0.98 * pmin;
    std::vector<double> center;
    DoMaskShift<double>(Np, thresh, -1, data.data());
    core.clear();
    WeightedCenter<double>(N, data.data(), core);
    return core.size();
}

int IncFlow::ExtractVortexParam2Dplane(const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
    std::vector<double> &v, double &radius, double &circulation) {
        int indcore = Index(N, core);
        double sign = v[indcore] / std::fabs(v[indcore]);
        double thresh = 0.01 * sign * v[indcore];
        int ixmax, ixmin, iymax, iymin;
        for(ixmax=core[0]; ixmax<N[0]; ++ixmax) {
            int ind1 = Index(N, {ixmax, core[1]});
            if(sign * v[ind1] < thresh) {
                break;
            }
        }
        for(ixmin=core[0]; ixmin>=0; --ixmin) {
            int ind1 = Index(N, {ixmin, core[1]});
            if(sign * v[ind1] < thresh) {
                break;
            }
        }
        for(iymax = core[1]; iymax<N[1]; ++iymax) {
            int ind1 = Index(N, {core[0], iymax});
            if(sign * v[ind1] < thresh) {
                break;
            }
        }
        for(iymin = core[1]; iymin>=0; --iymin) {
            int ind1 = Index(N, {core[0], iymin});
            if(sign * v[ind1] < thresh) {
                break;
            }
        }
        ixmax = std::max(0, std::min(ixmax, N[0]-1));
        ixmin = std::max(0, std::min(ixmin, N[0]-1));
        iymax = std::max(0, std::min(iymax, N[1]-1));
        iymin = std::max(0, std::min(iymin, N[1]-1));
        int tmpmax = -1;
        double tmpvalue = -1E8;
        for(int j=iymin; j<=iymax; ++j) {
            for(int i=ixmin; i<=ixmax; ++i) {
                int n = Index(N, {i,j});
                double value = v[n]*sign;
                if(value < thresh) {
                    v[n] = 0;
                } else {
                    if(value > tmpvalue) {
                        tmpmax = n;
                        tmpvalue = value;
                    }
                }
            }
        }
        if(tmpmax>0) {
            invIndex(N, tmpmax, core);
        }

        double sum0 = 0., sum2 = 0.;
        for(int j=iymin; j<=iymax; ++j) {
            for(int i=ixmin; i<=ixmax; ++i) {
                int n = Index(N, {i,j});
                double x = (i - core[0])*dx[0];
                double y = (j - core[1])*dx[1];
                sum0 += v[n];
                sum2 += v[n] * (x*x + y*y);
            }
        }
        radius = std::sqrt(sum2/sum0);
        circulation = sum0 * dx[0] * dx[1];
        return 2;
    }