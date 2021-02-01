#include<cmath>
#include<vector>
#include<iostream>
#include"airfoil.h"

NACAmpxx::NACAmpxx(double m, double p, double t) {
    m_m = m;
    m_p = p;
    if(m_m<0.001) {
        m_m = 0.;
        m_p = 100.;
    }
    m_t = t;
	calculateArcTable();
    m_rRoundTrailing = -1.;
}

NACAmpxx::NACAmpxx(std::string name) {
    if(name.size()<4) {
        std::cout << "incorrect NACA 4-digit format" << std::endl;
        exit(-1);
    }
    m_m = (name[0] - '0')*0.01;
    m_p = (name[1] - '0')*0.1;
    m_t = (name[2] - '0')*0.1 + (name[3] - '0')*0.01;
    if(m_m<0.001) {
        m_m = 0.;
        m_p = 100.;
    }
    calculateArcTable();
    m_rRoundTrailing = -1.;
}

double NACAmpxx::findx(double s, int surf) {
	if(surf>0) {
        return findx(s, m_arcu);
    } else {
        return findx(s, m_arcd);
    }
}

double NACAmpxx::findx(double s, std::vector<std::vector<double>> & arc) {
    s = fabs(s);
	if(s>= arc[arc.size()-1][1]) return 1.;
    int i = 1;
	while(arc[i][1]<s) ++i;
	return arc[i-1][0] + (s-arc[i-1][1])*(arc[i][0]-arc[i-1][0])/(arc[i][1]-arc[i-1][1]);
}

double NACAmpxx::finds(double x, int surf) {
	if(surf>0) {
        return finds(x, m_arcu);
    } else {
        return finds(x, m_arcd);
    }
}

double NACAmpxx::finds(double x, std::vector<std::vector<double>> & arc) {
    x = fabs(x);
	if(x>=1.) return arc[arc.size()-1][1];
    int i = 1;
	while(arc[i][0]<x) ++i;
	return arc[i-1][1] + (x-arc[i-1][0])*(arc[i][1]-arc[i-1][1])/(arc[i][0]-arc[i-1][0]);
}

double NACAmpxx::halft(double x) {
    x = fabs(x);
	std::vector<double> xs(5, 1.);
	for(int i=1; i<5; ++i) xs[i] = xs[i-1]*x;
    return 5.*m_t*(0.2969*sqrt(x) -0.1260*xs[1] - 0.3516*xs[2] + 0.2843*xs[3] - 0.1015*xs[4]);
}

double NACAmpxx::halfdt(double x) {
    x = fabs(x);
	std::vector<double> xs(5, 1.);
	for(int i=1; i<5; ++i) xs[i] = xs[i-1]*x;
    return 5.*m_t*(0.2969*0.5/sqrt(x) -0.1260 - 0.3516*2.*xs[1] + 0.2843*3.*xs[2] - 0.1015*4.*xs[3]);
}

double NACAmpxx::roundTrailingSize() {
    double x;
    calculateTrailingRadius(x);
    return 1. - x;
}

std::vector<double> NACAmpxx::roundTrailingEdge(std::vector<double>&p0, double eps) {
    if(m_rRoundTrailing < 0.) {
        m_rRoundTrailing = calculateTrailingRadius(m_xRoundTrailing);
    }
    std::vector<double> p1;
    bool notexit =  fabs(halft(p0[0]) - fabs(p0[1])) > eps;
    if(fabs(1.-p0[0])<eps && fabs(p0[1]) <= halft(1) + eps) notexit = false;
    if(m_rRoundTrailing < 0. || p0[0] <= m_xRoundTrailing || notexit) {
        p1.push_back(p0[0]);
        p1.push_back(p0[1]);
    } else {
        double r0 = p0[0] -1. + m_rRoundTrailing;
        double r1 = p0[1];
        double r = sqrt(r0*r0 + r1*r1);
        p1.push_back(1. + m_rRoundTrailing*(r0/r-1.));
        p1.push_back(     m_rRoundTrailing* r1/r);
    }
    return p1;
}

double NACAmpxx::calculateTrailingRadius(double &xtmp) {
    if(m_m > 1.E-6) {
        return -1.;
    }
    double eps = 1.E-14, rtmp, ftmp;
    double x[2];
    double f[2];
    x[0] = 0.9;
    x[1] = 0.9999; xtmp = x[1];
    f[0] = testRadius(x[0], rtmp);
    f[1] = testRadius(xtmp, rtmp);
    while(fabs(f[1])>eps) {
        xtmp = ( x[0]*f[1]-x[1]*f[0] )/( f[1]-f[0] );
        ftmp = testRadius(xtmp, rtmp);
        x[0] = x[1];
        f[0] = f[1];
        x[1] = xtmp;
        f[1] = ftmp;
    }
    return rtmp;
}

double NACAmpxx::testRadius(double x, double &r) {
    double y = halft(x);
    double t = halfdt(x);
    double xc = x + t*y;
    r = 1.-xc;
    return sqrt((x-xc)*(x-xc) + y*y) - r;
}

double NACAmpxx::chamber(double x) {
    double yc;
    if(x<=m_p) {
        yc = m_m/(m_p*m_p)*( 2.*m_p*x -x*x );
    } else {
        yc = m_m/(1.-m_p)/(1.-m_p)*( (1.-2.*m_p)+2.*m_p*x -x*x );
    }
    return yc;
}

std::vector<double> NACAmpxx::up(double x) {
    std::vector<double> res(2);
    if(x>1.) x = 1.;
    if(x<0.) x = 0.;
    double yt = halft(x);
    double yc = chamber(x);
    double the = theta(x);
    res[0] = x - yt*sin(the);
    res[1] = yc + yt*cos(the);
    return res;
}

std::vector<double> NACAmpxx::down(double x) {
    std::vector<double> res(2);
    if(x>1.) x = 1.;
    if(x<0.) x = 0.;
    double yt = halft(x);
    double yc = chamber(x);
    double the = theta(x);
    res[0] = x + yt*sin(the);
    res[1] = yc - yt*cos(the);
    return res;
}

double NACAmpxx::theta(double x) {
    double dyc;
    if(x<=m_p) {
        dyc = 2.*m_m/(m_p*m_p)*( m_p - x );
    } else {
        dyc = 2.*m_m/(1.-m_p)/(1.-m_p)*( m_p - x );
    }
    return atan(dyc);
}

static double distance(std::vector<double> p0, std::vector<double> p1) {
    return sqrt((p0[0] - p1[0])*(p0[0] - p1[0]) + (p0[1] - p1[1])*(p0[1] - p1[1]));
}
void NACAmpxx::calculateArcTable() {
	double xmid = 0.1;
	int N0 = 1000;
	int N1 = 900;
    std::vector<double> tmp(2, 0.);
    std::vector<double> p0(2, 0.);
    m_arcu.push_back(tmp);
    m_arcd.push_back(tmp);
    //////upper surface 
    double x, dx = xmid/N0/10., arcl = 0.;
    for(int i=1; i<=10*N0; ++i) {
        x = dx*i;
        std::vector<double> p1 = up(x);
		arcl += distance(p0, p1);
		if(i%10 == 0) {
			std::vector<double> tmp1(2);
			tmp1[0] = x;
			tmp1[1] = arcl;
            m_arcu.push_back(tmp1);
		}
		p0 = p1;
    }
	dx = (1.-xmid)/N1;
    for(int i=1; i<=N1; ++i) {
        x = dx*i + xmid;
        std::vector<double> p1 = up(x);
		arcl += distance(p0, p1);
		if(1) {
			std::vector<double> tmp1(2);
			tmp1[0] = x;
			tmp1[1] = arcl;
            m_arcu.push_back(tmp1);
		}
		p0 = p1;
    }
    //////lower surface
    dx = xmid/N0/10.;
    arcl = 0.;
    p0[0] = 0.; p0[1] = 0.;
    for(int i=1; i<=10*N0; ++i) {
        x = dx*i;
        std::vector<double> p1 = down(x);
		arcl += distance(p0, p1);
		if(i%10 == 0) {
			std::vector<double> tmp1(2);
			tmp1[0] = x;
			tmp1[1] = arcl;
            m_arcd.push_back(tmp1);
		}
		p0 = p1;
    }
	dx = (1.-xmid)/N1;
    for(int i=1; i<=N1; ++i) {
        x = dx*i + xmid;
        std::vector<double> p1 = down(x);
		arcl += distance(p0, p1);
		if(1) {
			std::vector<double> tmp1(2);
			tmp1[0] = x;
			tmp1[1] = arcl;
            m_arcd.push_back(tmp1);
		}
		p0 = p1;
    }
}

/*int main() {
    NACAmpxx airf("0012");
    double x;
    printf("%20.12f, %20.12f, %20.12f, %20.12f\n", x, airf.calculateTrailingRadius(x), 1.-airf.calculateTrailingRadius(x), x);
}*/