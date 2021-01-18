#ifndef DATAPROCESS_H
#define DATAPROCESS_H
#include<vector>

class GaussKernel {
public:
    GaussKernel(double sigma, double dx, double cutoff = 5.);
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
};

#endif