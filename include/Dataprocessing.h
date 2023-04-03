#ifndef DATAPROCESS_H
#define DATAPROCESS_H
#include<vector>
#define GAUSS 0
class KernelSmooth {
public:
    KernelSmooth();
    int GetIMax();
    double GW(int i, int padding);
    std::vector<std::vector<double> > m_weight;
    // 0
    // 0 1
    // 0 1 2
    // multi dimensional smoothing
    int DoSmooth(const std::vector<double> &sigmaintegerwidth, const std::vector<int> &N, std::vector<std::vector<double> > &sdata);
    int m_maxOffset;
private:
    double Sum(const double *data, int padding);
    int BoundPadding(int padding);
    double Kernel(double sigma, double x, int kernel = GAUSS);
    void UpdateWeight(double sigmaintegerwidth = 10., double cutoff = 5., int kernel = GAUSS);
    int DoSmooth(int N, double * sdata);
};

class Derivative {
public:
    Derivative();
    int Diff(const std::vector<int> &N, const std::vector<std::vector<double> > &data, std::vector<std::vector<double> > &ddata, double dx, int order);
private:
    int Diff(int N, const double * data, double * ddata, double dx, int order);
    std::map<int, std::vector<double> > m_centralcoeff;
};

class Interpolation {
public:
    static void CalcWeight1D(const std::vector<double> &xi, std::vector<double> &w);
    static void CalcWeight2D(const std::vector<double> &xi, std::vector<double> &w);
    static void CalcWeight3D(const std::vector<double> &xi, std::vector<double> &w);
    static void CalcWeight(const std::vector<double> &x1, const int dim,
                           std::vector<double> &w);
};

class Interpolation1DNonuniform {
public:
    static void CalcWeight1D(const std::vector<double> &xi, std::vector<double> &x, std::vector<int> &index, std::vector<std::vector<double>> &weight);
};

#endif