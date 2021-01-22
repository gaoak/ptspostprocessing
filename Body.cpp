#include<cmath>
#include "Body.h"
#include "Util.h"

Body::Body(std::string airfoil, double AoA)
    : m_airfoil(airfoil) {
    m_AoA = AoA;
}

Body::Body()
    : m_airfoil("0012") {
    m_AoA = 0.;
}

bool Body::IsInBody(std::vector<double> p, double tol) {
    transform(p, -m_AoA);
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