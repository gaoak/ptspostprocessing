#include "airfoil.h"

class Body {
public:
    Body();
    Body(std::string airfoil, std::vector<double> params);
    bool IsInBody(std::vector<double> p, double tol);
protected:
    NACAmpxx m_airfoil;
    double m_AoA;
    std::vector<double> m_span;
};