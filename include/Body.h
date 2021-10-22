#include "airfoil.h"

class Body {
public:
    ~Body();
    Body();
    Body(std::string airfoil, std::vector<double> params);
    Body(const Body & b2);
    bool IsInBody(std::vector<double> p, double tol);
protected:
    NACAmpxx m_airfoil;
    double m_AoA;
    std::vector<double> m_span;
};