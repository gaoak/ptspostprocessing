#include<tuple>
#include<limits>
#include "IncFlow.h"
#include "Util.h"
#include "UtilGraph.h"

IncFlow::IncFlow() {
    ;
}

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range,
        std::string bodyname, std::vector<double> param)
        :StructuredData(N, range), m_body(bodyname, param) {
}

IncFlow::IncFlow(std::string bodyname, std::vector<double> param)
        :m_body(bodyname, param) {
}

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range,
        const std::vector<std::vector<double> > & axis,
        std::string bodyname, std::vector<double> param)
        : StructuredData(N, range, axis), m_body(bodyname, param) {
}

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range,
        const std::vector<std::vector<double> > & axis)
        : StructuredData(N, range, axis) {
}

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range)
        :StructuredData(N, range) {
    ;
}

IncFlow::IncFlow(const IncFlow & flow)
    : StructuredData(flow),m_body(flow.m_body) {
}

IncFlow& IncFlow::operator=(const IncFlow & flow) {
    StructuredData::operator=(flow);
    m_body = flow.m_body;
    return *this;
}

int IncFlow::CalculateVorticity(int order) {
    if(GetVelocityDimension()==2) {
        return CalculateVorticity2D(order);
    } else if(GetNumCoords()==3) {
        return CalculateVorticity3D(order);
    } else {
        return -1;
    }
}

int IncFlow::CalculateVorticity3D(int order) {
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

int IncFlow::CalculateVorticity2D(int order) {
    std::vector<std::vector<double> > u, ux(2), uy(2);
    u.push_back(m_phys[0]);
    u.push_back(m_phys[1]);
    for(int i=0; i<2; ++i) {
        ux[i].resize(m_Np);
        uy[i].resize(m_Np);
    }
    Diff(u, ux, 0, order);
    Diff(u, uy, 1, order);
    int id;
    //W_x
    m_vars.push_back("W_z");
    id = m_phys.size();
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    for(int i=0; i<m_Np; ++i) {
        m_phys[id][i] = ux[1][i] - uy[0][i];
    }

    // Q Criterion
    m_vars.push_back("Q");
    id = m_phys.size();
    m_phys.push_back(std::vector<double>(m_Np, 0.));
    for(int i=0; i<m_Np; ++i) {
        m_phys[id][i] = -0.5 * (
            ux[0][i]*ux[0][i] + ux[1][i]*uy[0][i] +
            uy[0][i]*ux[1][i] + uy[1][i]*uy[1][i]
        );
    }
    return m_Np;
}

int IncFlow::TransformCoord(const std::vector<double> &x0) {
    for(int k=0; k<GetNumCoords(); ++k) {
        m_axis.m_o[k] += x0[k];
        for(int i=0; i<m_Np; ++i) {
            m_x[k][i] += x0[k];
        }
    }
    return m_Np;
}

std::pair<int, std::vector<int> > IncFlow::GetProceedDirectionInt(const std::vector<double> &vor, double sign) {
    //need check vor!=0 before calling this function
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

std::pair<int, std::vector<double> > IncFlow::GetProceedDirection(const std::vector<double> &vor, double sign) {
    //need check vor!=0 before calling this function
    std::pair<int, std::vector<double>> res;
    res.first = FindAbsMax(3, vor.data());
    std::vector<double> inc = {vor[0], vor[1], vor[2]};
    NormalizeVect(inc);
    double dx = std::min(m_dx[0], m_dx[1]);
    dx = std::min(dx, m_dx[2]);
    if(sign<0) {
        dx = -dx;
    }
    for(int i=0; i<3; ++i) {
        inc[i] = inc[i] * dx;
    }
    res.second = inc;
    //printf("direction %d, (%d,%d,%d)\n", res.first, res.second[0], res.second[1], res.second[2]);
    return res;
}

int IncFlow::GetSubdomainRange(const std::vector<int> &center, double radius,
                               std::vector<int> &range) {
    range.resize(m_N.size() * 2);
    if(center.size()<m_N.size()) {
        for(int i=0; i<(int)m_N.size(); ++i) {
            range[2*i  ] = 0;
            range[2*i+1] = m_N[i]-1;
        }
    } else {
        for(int i=0; i<(int)m_N.size(); ++i) {
            int p = radius / m_dx[i];
            range[2*i  ] = std::min(m_N[i]-1, std::max(0, center[i] - p));
            range[2*i+1] = std::min(m_N[i]-1, std::max(0, center[i] + p));
        }
    }
    return m_N.size();
}

int IncFlow::ExtractCoreByPoint(
        std::vector<std::vector<double> > & cores, std::set<int> &searched,
        std::vector<double> &inputcenter, const std::vector<int> vf,
        const int field, const bool stoponwall, const double threshold, const double walltol,
        VortexMethod vm) {
    std::vector<std::vector<double> > cores1, cores2;
    std::vector<double> inputc = inputcenter;
    VortexExtractionStopReason r1 = ExtractCoreByPointDirection(
        cores1, searched, inputc, vf, field,  1, stoponwall, threshold, walltol, vm);
    if(cores1.size()) {
        inputcenter = inputc;
    }
    VortexExtractionStopReason r2 = ExtractCoreByPointDirection(
        cores2, searched, inputc, vf, field, -1, stoponwall, threshold, walltol, vm);
    if(cores2.size()) {
        inputcenter = inputc;
    }
    cores.clear();
    for(int i=(int)cores2.size()-1; i>0; --i) {
        cores.push_back(cores2[i]);
    }
    for(int i=0; i<(int)cores1.size()-1; ++i) {
        cores.push_back(cores1[i]);
    }
    printf("ExtractCoreByPoint stop because of %d,%d\n", r1, r2);
    return cores.size();
}

int IncFlow::Extract2DVortex(std::vector<std::vector<int>> &intcenters,
        std::vector<std::vector<double>> &physcenters,
        std::vector<std::vector<double>> &info, const std::vector<int> &v,
        const std::pair<int, int> &plane, const double threshold) {

    bool ismax = false;
    SearchAllCoreXYZplane(intcenters, physcenters, info, v, ismax, plane, threshold);

    for(size_t i=0; i<intcenters.size(); ++i) {
        std::vector<double> tempinfo;
        SearchOneCoreXYZplane(intcenters[i], physcenters[i], tempinfo, v, 3., ismax, 2);
        info[i].insert(info[i].end(), tempinfo.begin(), tempinfo.end());
    }
    return info.size();
}

int IncFlow::SearchAllCoreXYZplane(std::vector<std::vector<int>> &intcenters,
    std::vector<std::vector<double>> &physcenters, std::vector<std::vector<double>> &info,
    const std::vector<int> &v, const bool ismax, const std::pair<int, int> plane, const double threshold) {
    if(plane.first<0 || plane.first>2) {
        printf("error: a correct direction [%d] must be assigned in SearchAllCoreXYZplane\n", plane.first);
        return -1;
    }
    intcenters.clear();
    physcenters.clear();
    info.clear();
    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<int> planeN, subrange = {0, m_N[0]-1, 0, m_N[1]-1, 0, m_N[2]-1};
    std::vector<double> planedata;
    ExtractPlane(m_phys[v[3]], plane, subrange, planeN, planedata);

    ShiftArray<double>(dx, 2-plane.first);
    ShiftArray<int>(subrange, 2*(2-plane.first));
    ShiftArray<int>(N, 2-plane.first);

    FindAllLocMaxIn2DGraph(planeN, planedata, intcenters, threshold, ismax);

    physcenters.resize(intcenters.size());
    for(size_t i=0; i<intcenters.size(); ++i) {
        double tempvalue = planedata[intcenters[i][0]+intcenters[i][1]*planeN[0]];
        intcenters[i][0] += subrange[0];
        intcenters[i][1] += subrange[2];
        intcenters[i].push_back(plane.second);
        physcenters[i].resize(2);
        for(int k=0; k<2; ++k) {
            physcenters[i][k] = dx[k] * (intcenters[i][k] + subrange[2*k]);
        }
        physcenters[i].push_back(plane.second*dx[2]);
        ShiftArray<int>(intcenters[i], plane.first-2);
        ShiftArray<double>(physcenters[i], plane.first-2);
        m_axis.ToPhysCoord(physcenters[i]);
        info.push_back({physcenters[i][0], physcenters[i][1], physcenters[i][2], tempvalue});
    }
    return info.size();
}

int IncFlow::SearchOneCoreXYZplane(
        std::vector<int> &intcenter, std::vector<double> &physcenter, std::vector<double> &info,
        const std::vector<int> &v, const double ran, const bool ismax, int dir) {
    //printf("center is %d,%d,%d\n", intcenter[0], intcenter[1], intcenter[2]);
    int centerindex = Index(m_N, intcenter);
    if(dir<0 || dir>2) {
        std::vector<double> vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
        dir = FindAbsMax(3, vor.data());
    }
    std::pair<int, int> plane;
    plane.first = dir;
    plane.second = intcenter[dir];
    //printf("search plane %d, %d\n", dir, plane.second);
    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<int> planeN, subrange;
    std::vector<double> planedata, tmpradius, planevorticity;
    double tmpcirculation;
    GetSubdomainRange(intcenter, ran, subrange);
    ExtractPlane(m_phys[v[3]], plane, subrange, planeN, planedata);
    ExtractPlane(m_phys[v[dir]], plane, subrange, planeN, planevorticity);

    ShiftArray<double>(dx, 2-dir);
    ShiftArray<int>(subrange, 2*(2-dir));
    ShiftArray<int>(N, 2-dir);
    ShiftArray<int>(intcenter, 2-dir);
    intcenter[0] -= subrange[0];
    intcenter[1] -= subrange[2];

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
    centerindex = Index(m_N, intcenter);
    for(size_t i=0; i<v.size(); ++i) {
        info.push_back(m_phys[v[i]][centerindex]);
    }
    return 3;
}

int IncFlow::InterpolateFrom(const IncFlow & origin, std::map<int,double> field) {
    m_body = origin.m_body;
    return StructuredData::InterpolateFrom(origin, field);
}

int IncFlow::SearchOneCorePerpendicular(
        std::vector<double> &physcenter, std::vector<double> &info,
        const std::vector<int> &v, const double ran, const bool ismax) {
    //[-ran, ran]
    double dx = std::min(m_dx[0], m_dx[1]);
    dx = std::min(dx, m_dx[2]);
    int Nr = std::min(int(ran/dx*2. + 0.5), SLICERESOLUTION);
    if(Nr<1) {
        return -1; //search fail
    }
    std::map<int, double> field;
    std::map<int, double> value;
    for(size_t i=0; i<m_phys.size(); ++i) {
        field[i] = 0.;
    }
    if(ismax) {
        field[v[3]] = -1e10;
    } else {
        field[v[3]] = 1e10;
    }
    physcenter.resize(3);
    InterpolatePoint(physcenter, field, value);
    std::vector<double> vor = {value[v[0]], value[v[1]], value[v[2]]};
    int dir = FindAbsMax(3, vor.data());
    if(std::fabs(vor[dir])<std::numeric_limits<double>::epsilon()) {
        return -2; //zero vorticity point
    }
    NormalizeVect(vor);
    std::vector<double> e0, e1;
    e0 = {1.,1.,1.};
    e0[dir] = 0.;
    e0 = CrossVect(vor, e0);
    NormalizeVect(e0);
    e1 = CrossVect(vor, e0);
    std::vector<std::vector<double> > axis;
    axis.push_back(e0); axis.push_back(e1); axis.push_back(vor);
    std::vector<double> range(6);
    for(int i=0; i<3; ++i) {
        range[2*i  ] = physcenter[i] - e0[i] * ran - e1[i] * ran;
        range[2*i+1] = 2. * ran;
    }
    std::vector<int> planeN = {Nr*2 + 1, Nr*2 + 1, 1};
    std::vector<int> tmpintcenter = {Nr, Nr, 0};
    IncFlow slice(planeN, range, axis);
    slice.InterpolateFrom(*this, field);
    std::vector<double> planevorticity(slice.GetTotPoints());
    std::vector<double> planedata(slice.GetTotPoints());
    for(size_t i=0; i<planevorticity.size(); ++i) {
        planevorticity[i] = vor[0] * slice.m_phys[v[0]][i] + vor[1] * slice.m_phys[v[1]][i]
                          + vor[2] * slice.m_phys[v[2]][i];
        planedata[i] = slice.m_phys[v[3]][i];
    }

    planeN.resize(2);
    FindLocMaxIn2DGraph(planeN, tmpintcenter, planedata, tmpintcenter, ismax);

    std::vector<double> tmpradius;
    double tmpcirculation;
    ExtractVortexParam2Dplane(planeN, slice.m_dx, tmpintcenter, planevorticity, tmpradius, tmpcirculation);
    physcenter = {tmpintcenter[0]*slice.m_dx[0], tmpintcenter[1]*slice.m_dx[1], 0.};
    slice.m_axis.ToPhysCoord(physcenter);
    info = {physcenter[0], physcenter[1], physcenter[2], tmpradius[0], tmpradius[1],
        std::fabs(tmpcirculation)};
    tmpintcenter.push_back(0);
    int centerindex = Index(slice.m_N, tmpintcenter);
    for(size_t i=0; i<v.size(); ++i) {
        info.push_back(slice.m_phys[v[i]][centerindex]);
    }
    //output test
    //slice.OutputData("perpendicular.plt");
    //printf("phys point %f,%f,%f\n", physcenter[0], physcenter[1], physcenter[2]);
    //exit(0);
    return 3;
}

VortexExtractionStopReason IncFlow::ExtractCoreByPointDirection(
        std::vector<std::vector<double> > & cores, std::set<int> &searched,
        std::vector<double> &inputcenter, const std::vector<int> vf, const int field, const int direction,
        const bool stoponwall, const double threshold, const double walltol, VortexMethod vm) {
    if(vf.size()<3 || field==0 || direction==0 ||
       vf[0]<=0 || vf[1]<=0 || vf[2]<=0 || inputcenter.size()<3) {
        printf("error incorrect parameters for ExtractCore\n");
        return StopParamsError;
    }
    if(vm == VorticityLine) {
        ExtractCoreByPointDirectionXYZ(cores, searched,
            inputcenter, vf, field, direction,
            stoponwall, threshold, walltol, true); //initial search
        return ExtractCoreByPointDirectionVorticityLine(cores, searched,
            inputcenter, vf, field, direction,
            stoponwall, threshold, walltol);
    } else {
        return ExtractCoreByPointDirectionXYZ(cores, searched,
            inputcenter, vf, field, direction,
            stoponwall, threshold, walltol);
    }
}

VortexExtractionStopReason IncFlow::ExtractCoreByPointDirectionXYZ(
        std::vector<std::vector<double> > & cores, std::set<int> &searched,
        std::vector<double> &inputcenter, const std::vector<int> vf, const int field, const int direction,
        const bool stoponwall, const double threshold, const double walltol,
        const bool onestep) {
    if(vf.size()<3 || field==0 || direction==0 ||
       vf[0]<=0 || vf[1]<=0 || vf[2]<=0 || inputcenter.size()<3) {
        printf("error incorrect parameters for ExtractCore\n");
        return StopParamsError;
    }

    // get parameter
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
    for(int i=0; i<(int)m_N.size(); ++i) {
        radiusofsubrange += m_dx[i] * (m_N[i] - 1);
    }

    //start extraction
    cores.clear();
    int Trymax = m_Np, count = 0;
    int dir = -1;
    while(count<Trymax) {
        if(intcenter[0]<0 || intcenter[0]>=m_N[0] ||
           intcenter[1]<0 || intcenter[1]>=m_N[1] ||
           intcenter[2]<0 || intcenter[2]>=m_N[2]) {
            return StopOutDomain;
        }
        std::vector<double> coreinfo;
        SearchOneCoreXYZplane(intcenter, physcenter, coreinfo, v, radiusofsubrange, ismax, dir);
        int centerindex = Index(m_N, intcenter);
        std::vector<double> vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
        if((ismax && m_phys[v[3]][centerindex]<threshold) ||
           (!ismax && m_phys[v[3]][centerindex]>threshold) ||
           NormVect(vor, 1) < std::numeric_limits<double>::epsilon()) {
            return StopThreshold;
        }
        if(count==0) {
            inputcenter = physcenter;
        }
        cores.push_back(coreinfo);
        if(onestep) {
            return StopSuccess;
        }
        std::pair<int, std::vector<int> > incplane = GetProceedDirectionInt(vor, vorticityproceedsign);
        dir = incplane.first;
        radiusofsubrange = 3. * coreinfo[CoreR2];
        if(searched.find(centerindex)!=searched.end() && count) {
            return StopRepeatPoint; //reaching previous point
        } else {
            searched.insert(centerindex);
        }
        if(stoponwall && m_body.IsInBody(physcenter, walltol)) {
            return StopReachWall;
        }
        for(int i=0; i<3; ++i) {
            intcenter[i] += incplane.second[i];
        }
        ++count;
    }
    return StopMaxTry;
}

VortexExtractionStopReason IncFlow::ExtractCoreByPointDirectionVorticityLine(
        std::vector<std::vector<double> > & cores, std::set<int> &searched,
        std::vector<double> &inputcenter, const std::vector<int> vf, const int field, const int direction,
        const bool stoponwall, const double threshold, const double walltol) {
    if(vf.size()<3 || field==0 || direction==0 ||
       vf[0]<=0 || vf[1]<=0 || vf[2]<=0 || inputcenter.size()<3) {
        printf("error incorrect parameters for ExtractCore\n");
        return StopParamsError;
    }

    // get parameter
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
    std::vector<double> tmp = inputcenter;
    m_axis.ToCompCoord(tmp);
    std::vector<int> intcenter(3);
    for(int i=0; i<3; ++i) {
        intcenter[i] = myRound<double>(tmp[i]/m_dx[i]);
        if(intcenter[i]<0 || intcenter[i]>=m_N[i]) {
            return StopOutDomain;
        }
    }
    double radiusofsubrange = 0;
    if(cores.size()) {
        radiusofsubrange = 3. * cores[0][CoreR2];
    } else {
        radiusofsubrange = 0.3 * SLICERESOLUTION * (m_dx[0] + m_dx[1] + m_dx[2]);
    }

    //start extraction
    cores.clear();
    int Trymax = m_Np, count = 0;
    while(count<Trymax) {
        std::vector<double> coreinfo;
        int serr = SearchOneCorePerpendicular(
            physcenter, coreinfo, v, radiusofsubrange, ismax);
        if(0 > serr) {
            printf("SearchOneCorePerpendicular error %d, r=%f\n", serr, radiusofsubrange);
            return StopPerpendicularError;
        }
        tmp = physcenter;
        m_axis.ToCompCoord(tmp);
        for(int i=0; i<3; ++i) {
            intcenter[i] = myRound<double>(tmp[i]/m_dx[i]);
            if(intcenter[i]<0 || intcenter[i]>=m_N[i]) {
                return StopOutDomain;
            }
        }
        int centerindex = Index(m_N, intcenter);
        std::vector<double> vor = {m_phys[v[0]][centerindex], m_phys[v[1]][centerindex], m_phys[v[2]][centerindex]};
        if((ismax && m_phys[v[3]][centerindex]<threshold) ||
           (!ismax && m_phys[v[3]][centerindex]>threshold) ||
           NormVect(vor, 1) < std::numeric_limits<double>::epsilon()) {
            return StopThreshold;
        }
        cores.push_back(coreinfo);
        std::pair<int, std::vector<double> > incplane = GetProceedDirection(vor, vorticityproceedsign);
        radiusofsubrange = 3. * coreinfo[CoreR2];
        if(count==0) {
            inputcenter = physcenter;
        }
        if(searched.find(centerindex)!=searched.end() && count) {
            ;//return StopRepeatPoint; //reaching previous point
        } else {
            searched.insert(centerindex);
        }
        if(stoponwall && m_body.IsInBody(physcenter, walltol)) {
            return StopReachWall;
        }
        for(int i=0; i<3; ++i) {
            physcenter[i] += incplane.second[i];
        }
        ++count;
    }
    return StopMaxTry;
}

VortexExtractionStopReason IncFlow::RefineCore(
        std::vector<std::vector<double> > & cores, const std::vector<int> vf, const int field) {
    //cores should be vortex cores found by ExtractCoreByPoint
    if(vf.size()<3 || field==0 || vf[0]<=0 || vf[1]<=0 || vf[2]<=0) {
        printf("error incorrect parameters for ExtractCore\n");
        return StopParamsError;
    }

    // get parameter
    std::vector<int> v = {vf[0]-1, vf[1]-1, vf[2]-1, std::abs(field)-1};
    bool ismax = true;
    if(field<0) {
        ismax = false;
    }

    for(size_t c=0; c<cores.size(); ++c) {
        std::vector<double> physcenter(3, 0.);
        for(int i=0; i<3; ++i) {
            physcenter[i] = cores[c][i];
        }
        double radiusofsubrange = 3. * cores[c][CoreR2];
        std::vector<double> coreinfo;
        if(SearchOneCorePerpendicular(physcenter, coreinfo, v, radiusofsubrange, ismax)) {
            for(size_t i=0; i<coreinfo.size(); ++i) {
                cores[c][i] = coreinfo[i];
            }
        }
    }
    return StopSuccess;
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

bool IncFlow::IsInBody(std::vector<double> p, double tol) {
    return m_body.IsInBody(p, tol);
}

int IncFlow::OverWriteBodyPoint(const std::vector<double> &u0, const std::vector<double> &pivot, const std::vector<double> &omega) {
    for(int i=0; i<m_Np; ++i) {
        std::vector<double> x = {m_x[0][i], m_x[1][i], GetCoordValue(2, i)};
        if(m_body.IsInBody(x, 10. * std::numeric_limits<double>::epsilon())) {
            std::vector<double> vel = AddVect(1., CrossVect(omega, AddVect(1., x, -1., pivot)), 1., u0);
            for (size_t d=0; d<vel.size(); ++d) {
                m_phys[d][i] = vel[d];
            }
        }
    }
    return m_Np;
}