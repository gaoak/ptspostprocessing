#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <map>
#include "Dataprocessing.h"
#include "Util.h"

KernelSmooth::KernelSmooth(double sigma, double dx, double cutoff, int kernel) {
    m_maxOffset = std::floor(cutoff*sigma/dx);
    std::vector<double> coeff(m_maxOffset + 1, 0.);
    std::vector<double> sum(m_maxOffset + 1, 1.);
    for(int i=0; i<=m_maxOffset; ++i) {
        coeff[i] = Kernel(sigma, i*dx, kernel);
        if(i) {
            sum[i] = sum[i-1] + 2. * coeff[i];
        }
    }
    for(int i=-0; i<=m_maxOffset; ++i) {
        std::vector<double> tmpw(i+1);
        for(int j=0; j<=i; ++j) {
            tmpw[j] = coeff[j]/sum[i];
        }
        m_weight.push_back(tmpw);
    }
}

double KernelSmooth::Kernel(double sigma, double x, int kernel) {
    if(kernel == GAUSS) {
        double tmp = x/sigma;
        return std::exp(-0.5*tmp*tmp);
    } else {
        return 1.;
    }
}

int KernelSmooth::BoundPadding(int padding) {
    if(padding<0) {
        return 0;
    } else if(padding>m_maxOffset) {
        return m_maxOffset;
    } else {
        return padding;
    }
}

double KernelSmooth::GW(int i, int padding) {
    if(i>padding) {
        printf("error: array exceed bound in GW\n");
    }
    padding = BoundPadding(padding);
    return m_weight[padding][std::abs(i)];
}

double KernelSmooth::Sum(const double *data, int padding) {
    double sum = *data;
    padding = BoundPadding(padding);
    sum = (*data) * m_weight[padding][0];
    for(int i=1; i<=padding; ++i) {
        sum += m_weight[padding][i] * (*(data-i) + *(data+i));
    }
    return sum;
}

int KernelSmooth::DoSmooth(int N, const double * data, double * sdata) {
    for(int i=0; i<N; ++i) {
        int padding = std::min(i, N-i);
        padding = BoundPadding(padding);
        sdata[i] = Sum(data+i, padding);
    }
    return N;
}

int KernelSmooth::DoSmooth(const std::vector<int> &Nraw, const double *data, double *sdata) {
    std::vector<int> N;
    for(int i=0; i<Nraw.size(); ++i) {
        if(Nraw[i]>1) {
            N.push_back(Nraw[i]);
        }
    }
    std::vector<int> Np(N.size(), 1);
    std::vector<int> oN(N.size());
    std::vector<int> sN(N.size());
    int Npt = 1;
    for(int dir=0; dir<N.size(); ++dir) {
        for(int i=0; i<N.size(); ++i) {
            if(i==dir) continue;
            Np[dir] *= N[i];
        }
        Npt *= N[dir];
        sN[dir] = N[dir];
    }
    std::vector<double> swap(Npt);
    for(int i=0; i<Npt; ++i) {
        swap[i] = data[i];
    }
    for(int dir=0; dir<N.size(); ++dir) {
        for(int p=0; p<Np[dir]; ++p) {
            DoSmooth(N[dir], swap.data()+p*N[dir], sdata+p*N[dir]);
        }
        LeftShiftIndex(sN, oN, sdata, swap.data());
        sN = oN;
    }
    for(int i=0; i<Npt; ++i) {
        sdata[i] = swap[i];
    }
}

int KernelSmooth::GetIMax() {
    return m_maxOffset;
}