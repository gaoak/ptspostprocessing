#include "CAD2D/airfoil.h"

class Body {
public:
    Body();
    Body(std::string airfoil, double AoA);
    bool IsInBody(std::vector<double> p, double tol);
protected:
    NACAmpxx m_airfoil;
    double m_AoA;
};