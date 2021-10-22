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