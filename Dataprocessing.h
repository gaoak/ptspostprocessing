#ifndef DATAPROCESS_H
#define DATAPROCESS_H
#include<vector>

class GaussKernel {
public:
    GaussKernel(double sigma, double dx, double cutoff = 5.);
    int GetIMax();
    std::vector<std::vector<double> > m_weight;
	int DoSmooth(std::vector<double> &data, std::vector<double> & sdata);
	// 0
	// 0 1
	// 0 1 2
    int m_maxOffset;
private:
	double Sum(double *data, int padding);
};

#endif