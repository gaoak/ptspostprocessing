#ifndef HEADERFILEBODY_H
#define HEADERFILEBODY_H
#include "airfoil.h"

class Body {
public:
    Body();
    Body(std::string airfoil, std::vector<double> params);
    Body(const Body & b2);
    Body& operator=(const Body& b2);
    bool IsInBody(std::vector<double> p, double tol);
    bool IsBelowBody(std::vector<double> p, double tol);
protected:
    NACAmpxx m_airfoil;
    double m_AoA;
    std::vector<double> m_span;
    double m_roundtip;
};
#endif