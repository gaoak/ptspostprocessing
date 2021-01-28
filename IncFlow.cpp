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
}

int IncFlow::TransformCoord(const std::vector<double> &x0) {
    for(int k=0; k<3; ++k) {
        m_range[2*k] += x0[k];
        m_range[2*k+1] += x0[k];
        for(int i=0; i<m_Np; ++i) {
            m_x[k][i] += x0[k];
        }
    }
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

int IncFlow::ExtractCore(double sigma, std::vector<std::vector<double> > & cores,
    std::vector<std::vector<double> > & radius, std::vector<double> &circulation,
    std::vector<double> &inputcenter, std::vector<int> &vf, int field, int direction,
    bool stoponwall, double threshold) {
    if(vf.size()<3 || field==0 || direction==0 ||
       vf[0]<=0 || vf[1]<=0 || vf[2]<=0) {
        printf("error incorrect parameters for ExtractCore\n");
        return -1;
    }
    cores.clear();
    std::vector<std::vector<double> > odata;
    odata.push_back(m_phys[vf[0]-1]);
    odata.push_back(m_phys[vf[1]-1]);
    odata.push_back(m_phys[vf[2]-1]);
    odata.push_back(m_phys[std::abs(field)-1]);
    bool ismax = true;
    if(field<0) {
        ismax = false;
    }
    //Smoothing(sigma, odata);
    double vorticityproceedsign = 1.;
    if(direction<0) {
        direction = -direction - 1;
        vorticityproceedsign = -1.;
    } else {
        direction -= 1;
    }

    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<double> range = m_range;

    int Trymax = 2*(m_N[0] + m_N[1] + m_N[2]);
    std::set<int> searched;
    std::vector<int> intcenter;
    if(inputcenter.size()>2) {
        intcenter.resize(3);
        for(int i=0; i<3; ++i) {
            intcenter[i] = myRound<double>((inputcenter[i] - m_range[2*i]) / m_dx[i]);
        }
    }
    std::pair<int, int> plane(direction, intcenter[direction]);
    std::vector<int> planeN;
    std::vector<double> planedata;
    int count = 0;
    double tmpcirculation;
    std::vector<double> tmpradius;
    std::vector<double> planevorticity;
    std::pair<int, std::vector<int> > incplane;
    incplane.first = direction;
    incplane.second = {0, 0, 0};
    incplane.second[direction] = 1;
    while(count<Trymax) {
        if(plane.second<0 || plane.second>=N[plane.first]) {
            break;
        }
        int dir = plane.first;
        //printf("search plane %d, %d\n", dir, plane.second);
        ExtractPlane(odata[3], plane, planeN, planedata);
        ExtractPlane(odata[dir], plane, planeN, planevorticity);

        ShiftArray<double>(dx, 2-dir);
        ShiftArray<double>(range, 2*(2-dir));
        ShiftArray<int>(N, 2-dir);
        ShiftArray<int>(intcenter, 2-dir);

        std::vector<double> newcenter;
        std::vector<int> Nslice = {N[0], N[1]};
        double tmpsign = vorticityproceedsign;
        if(intcenter.size()) {
            tmpsign = planevorticity[Index(Nslice, intcenter)];
        }
        PurgeDifferentSign(Nslice, planevorticity, planedata, tmpsign);
        ExtractCore2Dplane(Nslice, intcenter, planedata, newcenter, ismax);
        std::vector<double> physcenter = newcenter;
        intcenter.resize(2);
        for(int k=0; k<2; ++k) {
            physcenter[k] = range[2*k] + dx[k] * physcenter[k];
            intcenter[k] = myRound<double>(newcenter[k]);
        }
        physcenter.push_back(range[4] + (plane.second)*dx[2]);
        intcenter.push_back(plane.second);
        ExtractVortexParam2Dplane(Nslice, dx, intcenter, planevorticity, tmpradius, tmpcirculation);
        tmpradius[0] = std::sqrt(tmpradius[0] * tmpradius[0] - 0. * sigma * sigma);
        tmpradius[1] = std::sqrt(tmpradius[1] * tmpradius[1] - 0. * sigma * sigma);
        
        {
            //StructuredData tmp(Nslice, range);
            //tmp.AddPhysics("vor", planevorticity);
            //tmp.OutputData("vor.plt");
            //exit(0);
        }

        ShiftArray<double>(dx, dir-2);
        ShiftArray<double>(physcenter, dir-2);
        ShiftArray<double>(range, 2*(dir-2));
        ShiftArray<int>(N, dir-2);
        ShiftArray<int>(intcenter, dir-2);

        int centerindex = Index(m_N, intcenter);
        if((ismax && odata[3][centerindex]<threshold) || (!ismax && odata[3][centerindex]>threshold) ) {
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

        std::vector<double> vor = {odata[0][centerindex], odata[1][centerindex], odata[2][centerindex]};
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

int IncFlow::ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
    std::vector<double> &data, std::vector<double> &core, bool ismax) {
    return FindLocMaxIn2DGraph(N, initial, data, core, ismax);
}

int IncFlow::PurgeDifferentSign(const std::vector<int> &N, const std::vector<double> &v, std::vector<double> &data, double sign) {
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
        for(int i=0; i<N[0]; ++i) {
            int ind = Index(N, {i, j});
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
        for(int i=0; i<N[0]; ++i) {
            int ind = Index(N, {i, j});
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
}