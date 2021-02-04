#include<tuple>
#include<set>
#include<limits>
#include "IncFlow.h"
#include "Util.h"

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
        m_range[2*k] += x0[k];
        m_range[2*k+1] += x0[k];
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

int IncFlow::ExtractCore(double sigma, std::vector<std::vector<double> > & cores,
    std::vector<std::vector<double> > & radius, std::vector<double> &circulation,
    std::vector<double> &inputcenter, std::vector<int> &vf, int field, int direction,
    bool stoponwall, double threshold) {
    if(vf.size()<3 || field==0 || direction==0 ||
       vf[0]<=0 || vf[1]<=0 || vf[2]<=0 || inputcenter.size()<3) {
        printf("error incorrect parameters for ExtractCore\n");
        return -1;
    }
    cores.clear();
    // get parameter
    std::vector<std::vector<double> > odata;
    std::vector<int> v = {vf[0]-1, vf[1]-1, vf[2]-1};
    int Qval = std::abs(field)-1;
    bool ismax = true;
    if(field<0) {
        ismax = false;
    }
    double vorticityproceedsign = 1.;
    if(direction<0) {
        direction = -direction - 1;
        vorticityproceedsign = -1.;
    } else {
        direction -= 1;
    }
    // set initial direction
    std::vector<int> intcenter(3);
    for(int i=0; i<3; ++i) {
        intcenter[i] = myRound<double>((inputcenter[i] - m_range[2*i]) / m_dx[i]);
    }
    int centerindex = Index(m_N, intcenter);
    std::vector<double> vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
    std::pair<int, std::vector<int> > incplane = GetProceedDirection(vor, vorticityproceedsign);
    std::pair<int, int> plane(incplane.first, intcenter[incplane.first]);

    int Trymax = 2*(m_N[0] + m_N[1] + m_N[2]), count = 0;
    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<double> range = m_range;
    std::set<int> searched;
    std::vector<int> planeN, subrange;
    std::vector<double> planedata, tmpradius, planevorticity;
    double radiusofsubrange = 0, tmpcirculation;
    for(int i=0; i<m_N.size(); ++i) {
        radiusofsubrange += m_range[2*i+1] - m_range[2*i];
    }
    while(count<Trymax) {
        if(plane.second<0 || plane.second>=N[plane.first]) {
            break;
        }
        int dir = plane.first;
        //printf("search plane %d, %d\n", dir, plane.second);
        GetSubdomainRange(intcenter, radiusofsubrange, subrange);
        ExtractPlane(m_phys[Qval], plane, subrange, planeN, planedata);
        ExtractPlane(m_phys[v[dir]], plane, subrange, planeN, planevorticity);

        ShiftArray<double>(dx, 2-dir);
        ShiftArray<double>(range, 2*(2-dir));
        ShiftArray<int>(subrange, 2*(2-dir));
        ShiftArray<int>(N, 2-dir);
        ShiftArray<int>(intcenter, 2-dir);
        intcenter[0] -= subrange[0];
        intcenter[1] -= subrange[2];

        double tmpsign = vorticityproceedsign;
        tmpsign = planevorticity[Index(planeN, intcenter)];
        PurgeDifferentSign(planeN, planevorticity, planedata, tmpsign);
        std::vector<double> physcenter;
        ExtractCore2Dplane(planeN, intcenter, planedata, physcenter, ismax);
        intcenter.resize(2);
        for(int k=0; k<2; ++k) {
            intcenter[k] = myRound<double>(physcenter[k]);
            physcenter[k] = range[2*k] + dx[k] * (physcenter[k] + subrange[2*k]);
        }
        physcenter.push_back(range[4] + (plane.second)*dx[2]);
        intcenter.push_back(plane.second);
        ExtractVortexParam2Dplane(planeN, dx, intcenter, planevorticity, tmpradius, tmpcirculation);
        tmpradius[0] = std::sqrt(tmpradius[0] * tmpradius[0] - 0. * sigma * sigma);
        tmpradius[1] = std::sqrt(tmpradius[1] * tmpradius[1] - 0. * sigma * sigma);
        radiusofsubrange = tmpradius[1] * 3.;

        intcenter[0] += subrange[0];
        intcenter[1] += subrange[2];
        ShiftArray<double>(dx, dir-2);
        ShiftArray<double>(physcenter, dir-2);
        ShiftArray<double>(range, 2*(dir-2));
        ShiftArray<int>(N, dir-2);
        ShiftArray<int>(intcenter, dir-2);

        centerindex = Index(m_N, intcenter);
        if((ismax && m_phys[Qval][centerindex]<threshold) || (!ismax && m_phys[Qval][centerindex]>threshold) ) {
            break;
        }
        cores.push_back(physcenter);
        radius.push_back(tmpradius);
        circulation.push_back(std::fabs(tmpcirculation));
        if(count==0) {
            inputcenter = physcenter;
        }
        if(searched.find(centerindex)!=searched.end()) {
            searched.insert(centerindex);
            break;
        }
        searched.insert(centerindex);
        //printf("location %f, %f, %f\n", physcenter[0], physcenter[1], physcenter[2]);
        if(stoponwall && m_body.IsInBody(physcenter, sigma)) {
            break;
        }

        vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
        incplane = GetProceedDirection(vor, vorticityproceedsign);
        plane.first = incplane.first;
        plane.second = intcenter[incplane.first] + incplane.second[incplane.first];
        for(int i=0; i<3; ++i) {
            intcenter[i] += incplane.second[i];
            if(intcenter[i]<0) {
                intcenter[i] = 0;
            }
            if(intcenter[i]>=m_N[i]) {
                intcenter[i] = m_N[i] - 1;
            }
        }
        //printf("next plane %d, %d\n", plane.first, plane.second);
        ++count;
    }
    return count;
}

IncFlow::IncFlow() {
    ;
}

int IncFlow::CopyAsSubDomain(const std::vector<int> &Ns, const std::vector<int> &Ne,
                             const std::vector<int> &skip, const IncFlow & big) {
    m_body = big.m_body;
    return StructuredData::CopyAsSubDomain(Ns, Ne, skip, big);
}

int IncFlow::ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
    std::vector<double> &data, std::vector<double> &core, bool ismax) {
    return FindLocMaxIn2DGraph(N, initial, data, core, ismax);
}

int IncFlow::PurgeDifferentSign(const std::vector<int> &N, const std::vector<double> &v,
                                std::vector<double> &data, double sign) {
    int Np = 1;
    for(int i=0; i<N.size(); ++i) Np *= N[i];
    if(Np>data.size()) Np = data.size();
    for(int i=0; i<Np; ++i) {
        if(v[i]*sign < 0.) {
            data[i] = 0.;
        }
    }
    return Np;
}

int IncFlow::ExtractVortexParam2Dplane(const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
    std::vector<double> &v, std::vector<double> &radius, double &circulation) {
    Fill2DGraph(N, v, core, 0.05, true);
    double sum0 = 0., sum1x = 0., sum1y = 0.;
    for(int j=0; j<N[1]; ++j) {
        int tmp = j * N[0];
        for(int i=0; i<N[0]; ++i) {
            int ind = i + tmp;
            sum0  += v[ind];
            sum1x += v[ind] * i;
            sum1y += v[ind] * j;
        }
    }
    if(std::fabs(sum0) < std::numeric_limits<double>::epsilon() ) {
        radius.clear();
        radius.push_back(0.);
        radius.push_back(0.);
        circulation = 0.;
        return 0;
    }

    sum1x /= sum0;
    sum1y /= sum0;
    double sumxx = 0., sumxy = 0., sumyy = 0.;
    for(int j=0; j<N[1]; ++j) {
        int tmp = j * N[0];
        for(int i=0; i<N[0]; ++i) {
            int ind = i + tmp;
            sumxx += v[ind] * (i - sum1x) * (i - sum1x) * dx[0] * dx[0];
            sumxy += v[ind] * (i - sum1x) * (j - sum1y) * dx[0] * dx[1];
            sumyy += v[ind] * (j - sum1y) * (j - sum1y) * dx[1] * dx[1];
        }
    }
    sumxx /= sum0;
    sumxy /= sum0;
    sumyy /= sum0;
    // l^2 - (xx + yy) l + xx yy -xy xy = 0
    double b = - sumxx - sumyy;
    double c = sumxx*sumyy - sumxy*sumxy;
    double det = std::sqrt(b*b - 4.*c);
    radius.clear();
    radius.push_back(std::sqrt(0.5*(-b - det)));
    radius.push_back(std::sqrt(0.5*(-b + det)));
    circulation = sum0 * dx[0] * dx[1];
    return 2;
}

int IncFlow::OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot, const std::vector<double> &omega) {
    for(int i=0; i<m_Np; ++i) {
        std::vector<double> x = {m_x[0][i], m_x[1][i], m_x[2][i]};
        if(m_body.IsInBody(x, 10. * std::numeric_limits<double>::epsilon())) {
            std::vector<double> vel = AddVect<double>(1., crossproduct<double>(omega, AddVect<double>(1., x, -1., pivot)), 1., u0);
            m_phys[0][i] = vel[0];
            m_phys[1][i] = vel[1];
            m_phys[2][i] = vel[2];
        }
    }
    return m_Np;
}