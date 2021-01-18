#ifndef DATAPROCESS_H
#define DATAPROCESS_H
#include<vector>
#define GAUSS 0
class KernelSmooth {
public:
    KernelSmooth(double sigma, double dx, double cutoff = 5., int kernel = GAUSS);
    int GetIMax();
    double GW(int i, int padding);
    std::vector<std::vector<double> > m_weight;
	// 0
	// 0 1
	// 0 1 2
	int DoSmooth(int N, const double * data, double * sdata);
    // multi dimensional smoothing
    int DoSmooth(const std::vector<int> &N, const double *data, double *sdata);
    int m_maxOffset;
private:
	double Sum(const double *data, int padding);
    int BoundPadding(int padding);
    double Kernel(double sigma, double x, int kernel = GAUSS);
};

class Derivative {
public:
    Derivative(int direction, double dx, double cutoff = 5., int order = 2);
};

#endif