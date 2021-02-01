#ifndef AIRFOIL_H
#define AIRFOIL_H
#include<vector>
#include<string>
class NACAmpxx {
public:
    NACAmpxx(double m, double p, double t);
    NACAmpxx(std::string name);
    std::vector<double> up(double x);
    std::vector<double> down(double x);
	double findx(double s, int up = 1);
	double finds(double x, int up = 1);
    std::vector<double> roundTrailingEdge(std::vector<double>&p0, double eps = 1.E-8);
    double roundTrailingSize();
protected:
    double findx(double s, std::vector<std::vector<double>> & arc);
    double finds(double x, std::vector<std::vector<double>> & arc);
    double halft(double x);
    double halfdt(double x);
    double chamber(double x);
    double theta(double x);
    double calculateTrailingRadius(double &xtmp);
    double testRadius(double x, double &r);
    void calculateArcTable();
    std::vector<std::vector<double>> m_arcu;
    std::vector<std::vector<double>> m_arcd;
    double m_m;
    double m_p;
    double m_t;
    double m_rRoundTrailing;
    double m_xRoundTrailing;
};

#endif
