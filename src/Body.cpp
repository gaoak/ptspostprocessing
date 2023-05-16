#include<cmath>
#include<limits>
#include "Body.h"
#include "Util.h"

Body::Body(std::string airfoil, std::vector<double> param)
    : m_airfoil(airfoil) {
    m_AoA = param[0];
    if(param.size()>=2) {
        m_span.push_back(param[1]);
    } else {
        m_span.push_back(std::numeric_limits<double>::min());
    }
    if(param.size()>=3) {
        m_span.push_back(param[2]);
    } else {
        m_span.push_back(std::numeric_limits<double>::max());
    }
    if(param.size()>=4) {
        m_roundtip = param[3];
    } else {
        m_roundtip = -10000.;
    }
}

Body::Body() {
    m_AoA = 0.;
}

Body& Body::operator=(const Body& b2) {
     m_airfoil = b2.m_airfoil;
     m_AoA = b2.m_AoA;
     m_span = b2.m_span;
     return *this;
}

 Body::Body(const Body & b2){
     *this = b2;
 }

bool Body::IsInBody(std::vector<double> p, double tol) {
    if(p[2]+tol<m_span[0] || p[2]-tol>m_span[1]) {
        return false;
    }
    p = transform(p, -m_AoA);
    if(p[0]<=0. || p[0]>=1.) return false;
    if(m_roundtip>0.) {
        double roundcurve = m_roundtip + sqrt(1 - (p[0]-1.)*(p[0]-1.) );
        if(p[1] > roundcurve) return false;
    }
    std::vector<double> foilup =  m_airfoil.up(p[0]);
    if(std::fabs(foilup[0] - p[0]) > 1E-6) {
        printf("error: asymmetric airfoil not supported\n");
    }
    if(std::fabs(p[1]) < tol+foilup[1]) {
        return true;
    } else {
        return false;
    }
}

bool Body::IsBelowBody(std::vector<double> p, double tol) {
    if(p[2]+tol<m_span[0] || p[2]-tol>m_span[1]) {
        return false;
    }
    if(p[0]<0 && p[1]<-0.1) return true;
    if(p[0]<1.1 && p[1]<-0.5) return true;
    if(p[0]<=0. || p[0]>=1.) return false;
    double y = p[0]*tan(-m_AoA);
    return p[1] <= y + tol;
}