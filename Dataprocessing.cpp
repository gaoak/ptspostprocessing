#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <map>
#include "Dataprocessing.h"
#include "Util.h"

KernelSmooth::KernelSmooth() {
    UpdateWeight();
}

void KernelSmooth::UpdateWeight(double sigmaintegerwidth, double cutoff, int kernel) {
    m_maxOffset = std::floor(cutoff*sigmaintegerwidth);
    std::vector<double> coeff(m_maxOffset + 1, 0.);
    std::vector<double> sum(m_maxOffset + 1, 1.);
    for(int i=0; i<=m_maxOffset; ++i) {
        coeff[i] = Kernel(sigmaintegerwidth, (double)i, kernel);
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

int KernelSmooth::DoSmooth(int N, double * sdata) {
    std::vector<double> data(N);
    for(int i=0; i<N; ++i) {
        data[i] = sdata[i];
    }
    for(int i=0; i<N; ++i) {
        int padding = std::min(i, N-i);
        padding = BoundPadding(padding);
        sdata[i] = Sum(data.data()+i, padding);
    }
    return N;
}

int KernelSmooth::DoSmooth(const std::vector<double> &sigmaintegerwidth, const std::vector<int> &Nraw, std::vector<std::vector<double> > &sdata) {
    std::vector<int> N;
    std::vector<int> sigma;
    for(int i=0; i<Nraw.size(); ++i) {
        if(Nraw[i]>1) {
            N.push_back(Nraw[i]);
            sigma.push_back(sigmaintegerwidth[i]);
        }
    }
    std::vector<int> Np(N.size(), 1);
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
    std::vector<std::vector<double> > swap;
    for(int i=0; i<sdata.size(); ++i) {
        std::vector<double> tmpv(Npt, 0.);
        swap.push_back(tmpv);
    }
    for(int d=0; d<sdata.size(); ++d) {
        for(int i=0; i<Npt; ++i) {
            swap[d][i] = sdata[d][i];
        }
    }
    for(int dir=0; dir<N.size(); ++dir) {
        UpdateWeight(sigma[dir]);
        for(int d=0; d<sdata.size(); ++d) {
            for(int p=0; p<Np[dir]; ++p) {
                DoSmooth(N[dir], swap[d].data()+p*N[dir]);
            }
        }
        ShiftIndex(sN, swap, -1);
    }
    for(int d=0; d<sdata.size(); ++d) {
        for(int i=0; i<Npt; ++i) {
            sdata[d][i] = swap[d][i];
        }
    }
}

int KernelSmooth::GetIMax() {
    return m_maxOffset;
}

Derivative::Derivative() {
    m_centralcoeff = {
        {1, {0., 1./2.}},
        {2, {0., 2./3., -1./12.}},
        {3, {0., 3./4., -3./20., 1./60.0}},
        {4, {0., 4./5., -1./5.0, 4./105., -1/280.}}
    };
}

int Derivative::Diff(int N, const double * uu, double * du, double dx, int order) {
    if(N<2) {
        printf("error: array length is 1 when calculating derivative\n");
        return 0;
    }
    dx = 1./dx;
    du[0] = dx * (uu[1] - uu[0]);
    du[N-1] = dx * (uu[N-1] - uu[N-2]);
    for(int i=1; i<N-1; ++i) {
        int padding = std::min(order/2, i);
        padding = std::min(padding, N-1-i);
        du[i] = 0.;
        for(int k=1;k<padding;++k) {
            du[i] += m_centralcoeff[padding][k] * (uu[i+k] - uu[i-k]);
        }
        du[i] *= dx;
    }
    return N;
}

int Derivative::Diff(const std::vector<int> &N, const std::vector<std::vector<double> > &data, std::vector<std::vector<double> > &ddata, double dx, int order) {
    int Np = 1;
    for(int i=1; i<N.size(); ++i) {
        Np *= N[i];
    }
    for(int d=0; d<data.size(); ++d) {
        for(int p=0; p<Np; ++p) {
            Diff(N[0], data[d].data() + p*N[0], ddata[d].data() + p*N[0], dx, order);
        }
    }
}