#include<tuple>
#include<limits>
#include "IncFlow.h"
#include "Util.h"
#include "UtilGraph.h"

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range,
        std::string bodyname, std::vector<double> param)
        :StructuredData(N, range), m_body(bodyname, param) {
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
    m_vars.push_back("W_x");
    m_vars.push_back("W_y");
    m_vars.push_back("W_z");
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
    m_vars.push_back("Q");
    id = m_phys.size();
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    for(int i=0; i<m_Np; ++i) {
        m_phys[id][i] = -0.5 * (
            ux[0][i]*ux[0][i] + ux[1][i]*uy[0][i] + ux[2][i]*uz[0][i] +
            uy[0][i]*ux[1][i] + uy[1][i]*uy[1][i] + uy[2][i]*uz[1][i] +
            uz[0][i]*ux[2][i] + uz[1][i]*uy[2][i] + uz[2][i]*uz[2][i]
        );
    }
    return m_Np;
}

int IncFlow::TransformCoord(const std::vector<double> &x0) {
    for(int k=0; k<3; ++k) {
        m_axis.m_o[k] += x0[k];
        for(int i=0; i<m_Np; ++i) {
            m_x[k][i] += x0[k];
        }
    }
    return m_Np;
}

std::pair<int, std::vector<int> > IncFlow::GetProceedDirection(const std::vector<double> &vor, double sign) {
    std::pair<int, std::vector<int>> res;
    res.first = FindAbsMax(3, vor.data());
    std::vector<int> inc(3);
    double ds = 0.;
    if(sign * vor[res.first] > 0) {
        ds = m_dx[res.first]/vor[res.first];
    } else {
        ds = -m_dx[res.first]/vor[res.first];
    }
    for(int i=0; i<3; ++i) {
        inc[i] = myRound<double>(vor[i] * ds / m_dx[i]);
    }
    res.second = inc;
    //printf("direction %d, (%d,%d,%d)\n", res.first, res.second[0], res.second[1], res.second[2]);
    return res;
}

int IncFlow::GetSubdomainRange(const std::vector<int> &center, double radius,
                               std::vector<int> &range) {
    range.resize(m_N.size() * 2);
    if(center.size()<m_N.size()) {
        for(int i=0; i<m_N.size(); ++i) {
            range[2*i  ] = 0;
            range[2*i+1] = m_N[i]-1;
        }
    } else {
        for(int i=0; i<m_N.size(); ++i) {
            int p = radius / m_dx[i];
            range[2*i  ] = std::min(m_N[i]-1, std::max(0, center[i] - p));
            range[2*i+1] = std::min(m_N[i]-1, std::max(0, center[i] + p));
        }
    }
    return m_N.size();
}

int IncFlow::ExtractCore(const double sigma, std::vector<std::vector<double> > & cores,
        std::set<int> &searched,
        std::vector<double> &inputcenter, const std::vector<int> &vf,
        const int field, const bool stoponwall, const double threshold) {
    std::vector<std::vector<double> > cores1, cores2;
    std::vector<double> inputc = inputcenter;
    ExtractCore(sigma, cores1, searched, inputcenter, vf, field,  1, stoponwall, threshold);
    ExtractCore(sigma, cores2, searched, inputc     , vf, field, -1, stoponwall, threshold);
    cores.clear();
    for(int i=(int)cores2.size()-1; i>0; --i) {
        cores.push_back(cores2[i]);
    }
    for(int i=0; i<(int)cores1.size()-1; ++i) {
        cores.push_back(cores1[i]);
    }
    return cores.size();
}

int IncFlow::SearchOneCoreXYZplane(
        std::vector<int> &intcenter, std::vector<double> &physcenter, std::vector<double> &info,
        const std::vector<int> &v, const double range, const bool ismax) {
    //printf("center is %d,%d,%d\n", intcenter[0], intcenter[1], intcenter[2]);
    int centerindex = Index(m_N, intcenter);
    std::vector<double> vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
    int dir = FindAbsMax(3, vor.data());
    std::pair<int, int> plane;
    plane.first = dir;
    plane.second = intcenter[dir];
    //printf("search plane %d, %d\n", dir, plane.second);
    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<int> planeN, subrange;
    std::vector<double> planedata, tmpradius, planevorticity;
    double tmpcirculation;
    GetSubdomainRange(intcenter, range, subrange);
    ExtractPlane(m_phys[v[3]], plane, subrange, planeN, planedata);
    ExtractPlane(m_phys[v[dir]], plane, subrange, planeN, planevorticity);

    ShiftArray<double>(dx, 2-dir);
    ShiftArray<int>(subrange, 2*(2-dir));
    ShiftArray<int>(N, 2-dir);
    ShiftArray<int>(intcenter, 2-dir);
    intcenter[0] -= subrange[0];
    intcenter[1] -= subrange[2];

    double tmpsign = planevorticity[Index(planeN, intcenter)];
    PurgeDifferentSign(planeN, planevorticity, planedata, tmpsign);
    FindLocMaxIn2DGraph(planeN, intcenter, planedata, intcenter, ismax);
    physcenter.resize(2);
    for(int k=0; k<2; ++k) {
        physcenter[k] = dx[k] * (intcenter[k] + subrange[2*k]);
    }
    physcenter.push_back(plane.second*dx[2]);
    intcenter.push_back(plane.second);

    ExtractVortexParam2Dplane(planeN, dx, intcenter, planevorticity, tmpradius, tmpcirculation);


    intcenter[0] += subrange[0];
    intcenter[1] += subrange[2];
    ShiftArray<double>(physcenter, dir-2);
    ShiftArray<int>(intcenter, dir-2);
    m_axis.ToPhysCoord(physcenter);
    info = {physcenter[0], physcenter[1], physcenter[2], tmpradius[0], tmpradius[1],
        std::fabs(tmpcirculation)};
    return 3;
}

int IncFlow::ExtractCore(const double sigma, std::vector<std::vector<double> > & cores, std::set<int> &searched,
    std::vector<double> &inputcenter, const std::vector<int> &vf, const int field, const int direction,
    const bool stoponwall, const double threshold) {
    if(vf.size()<3 || field==0 || direction==0 ||
       vf[0]<=0 || vf[1]<=0 || vf[2]<=0 || inputcenter.size()<3) {
        //printf("error incorrect parameters for ExtractCore\n");
        return -1;
    }

    // get parameter
    std::vector<std::vector<double> > odata;
    std::vector<int> v = {vf[0]-1, vf[1]-1, vf[2]-1, std::abs(field)-1};
    bool ismax = true;
    if(field<0) {
        ismax = false;
    }
    double vorticityproceedsign = 1.;
    if(direction<0) {
        vorticityproceedsign = -1.;
    }

    // set initial position and direction
    std::vector<double> physcenter = inputcenter;
    m_axis.ToCompCoord(physcenter);
    std::vector<int> intcenter(3);
    for(int i=0; i<3; ++i) {
        intcenter[i] = myRound<double>(physcenter[i]/m_dx[i]);
    }
    double radiusofsubrange = 0;
    for(int i=0; i<m_N.size(); ++i) {
        radiusofsubrange += m_dx[i] * (m_N[i] - 1);
    }

    //start extraction
    cores.clear();
    int Trymax = 2*(m_N[0] + m_N[1] + m_N[2]), count = 0;
    while(count<Trymax) {
        std::vector<double> coreinfo;
        SearchOneCoreXYZplane(intcenter, physcenter, coreinfo, v, radiusofsubrange, ismax);
        int centerindex = Index(m_N, intcenter);
        std::vector<double> vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
        std::pair<int, std::vector<int> > incplane = GetProceedDirection(vor, vorticityproceedsign);
        if((ismax && m_phys[v[3]][centerindex]<threshold) ||
           (!ismax && m_phys[v[3]][centerindex]>threshold) ) {
            break; //reach threshold
        }
        cores.push_back(coreinfo);
        if(count==0) {
            inputcenter = physcenter;
        }
        if(searched.find(centerindex)!=searched.end() && count) {
            break; //reaching previous point
        } else {
            searched.insert(centerindex);
        }
        if(stoponwall && m_body.IsInBody(physcenter, sigma)) {
            break; //touch wall
        }
        for(int i=0; i<3; ++i) {
            intcenter[i] += incplane.second[i];
        }
        if(intcenter[0]<0 || intcenter[0]>=m_N[0] ||
           intcenter[1]<0 || intcenter[1]>=m_N[1] ||
           intcenter[2]<0 || intcenter[2]>=m_N[2]) {
            break; //out of domain
        }
        ++count;
    }
    return 0;
}

IncFlow::IncFlow() {
    ;
}

int IncFlow::ExtractVortexParam2Dplane(const std::vector<int> &N, const std::vector<double> &dx,
        std::vector<int> core, std::vector<double> &planevorticity, std::vector<double> &radius,
        double &circulation) {
    std::vector<double> mu;
    ExtractPatchStat2DGraph(N, dx, core, planevorticity, mu, circulation);
    // l^2 - (xx + yy) l + xx yy -xy xy = 0
    double b = - mu[0] - mu[2];
    double c = mu[0]*mu[2] - mu[1]*mu[1];
    double det = std::sqrt(b*b - 4.*c);
    radius.clear();
    radius.push_back(std::sqrt(0.5*(-b - det)));
    radius.push_back(std::sqrt(0.5*(-b + det)));
    return 2;
}

int IncFlow::CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                             const std::vector<int> &skip, std::map<int, double> &field,
                             const IncFlow & big) {
    m_body = big.m_body;
    return StructuredData::CopyAsSubDomain(Ns, Ne, skip, field, big);
}

int IncFlow::OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot, const std::vector<double> &omega) {
    for(int i=0; i<m_Np; ++i) {
        std::vector<double> x = {m_x[0][i], m_x[1][i], m_x[2][i]};
        if(m_body.IsInBody(x, 10. * std::numeric_limits<double>::epsilon())) {
            std::vector<double> vel = AddVect(1., CrossVect(omega, AddVect(1., x, -1., pivot)), 1., u0);
            m_phys[0][i] = vel[0];
            m_phys[1][i] = vel[1];
            m_phys[2][i] = vel[2];
        }
    }
    return m_Np;
}