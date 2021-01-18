#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <map>
#include "Dataprocessing.h"
#include "Util.h"

GaussKernel::GaussKernel(double sigma, double dx, double cutoff) {
    m_maxOffset = std::floor(cutoff*sigma/dx);
    std::vector<double> gauss(m_maxOffset + 1, 0.);
    std::vector<double> sum(m_maxOffset + 1, 1.);
    for(int i=0; i<=m_maxOffset; ++i) {
        double tmp = i*dx/sigma;
        tmp = std::exp(-0.5*tmp*tmp);
        gauss[i] = tmp;
        if(i) {
            sum[i] = sum[i-1] + 2. * gauss[i];
        }
    }
    for(int i=-0; i<=m_maxOffset; ++i) {
        std::vector<double> tmpw(i+1);
        for(int j=0; j<=i; ++j) {
            tmpw[j] = gauss[j]/sum[i];
        }
        m_weight.push_back(tmpw);
    }
}

int GaussKernel::BoundPadding(int padding) {
    if(padding<0) {
        return 0;
    } else if(padding>m_maxOffset) {
        return m_maxOffset;
    } else {
        return padding;
    }
}

double GaussKernel::GW(int i, int padding) {
    if(i>padding) {
        printf("error: array exceed bound in GW\n");
    }
    padding = BoundPadding(padding);
    return m_weight[padding][std::abs(i)];
}

double GaussKernel::Sum(const double *data, int padding) {
    double sum = *data;
    padding = BoundPadding(padding);
    sum = (*data) * m_weight[padding][0];
    for(int i=1; i<=padding; ++i) {
        sum += m_weight[padding][i] * (*(data-i) + *(data+i));
    }
    return sum;
}

int GaussKernel::DoSmooth(int N, const double * data, double * sdata) {
    for(int i=0; i<N; ++i) {
        int padding = std::min(i, N-i);
        padding = BoundPadding(padding);
        sdata[i] = Sum(data+i, padding);
    }
    return N;
}

int GaussKernel::DoSmooth(const std::vector<int> &N, const double *data, double *sdata) {
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
        if(dir<N.size()-1) {
            LeftShiftIndex(sN, oN, sdata, swap.data());
            sN = oN;
        }
    }
}



int GaussKernel::GetIMax() {
    return m_maxOffset;
}

std::map<bool, std::string> testresults = {{true, "pass"},{false, "fail"}};

int test_weight(int padding) {
    GaussKernel ker(1., 0.3);
    double sum = 0.;
    for(int i=-padding; i<=padding; ++i) {
        double tmp = ker.GW(i, padding);
        printf("%g,",tmp);
        sum+= tmp;
    }
    printf("\ntest weight: sum %g, %s\n", sum-1., testresults[fabs(sum-1.)<1E-14].c_str());
}

int test_index(std::vector<int> &N) {
    int Np = N[0];
    for(int i=1; i<N.size(); ++i) {
        Np *= N[i];
    }
    std::vector<int> ind(N.size());
    int failcount = 0;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        int tmp = Index(N, ind);
        printf("%d, (", i);
        for(int j=0; j<ind.size(); ++j) {
            printf("%d,",ind[j]);
        }
        printf("), %d\n", tmp);
        failcount += i != tmp;
    }
    printf("test index: %d, %s\n", failcount, testresults[!failcount].c_str());
}

int main() {
    test_weight(6);
    std::vector<int> N={4,2,3};
    test_index(N);
    N={2,1,3};
    test_index(N);
    N={1,2,3};
    test_index(N);
    N={1,2,1};
    test_index(N);
    N={4,2,1};
    test_index(N);
    N={2,3};
    test_index(N);
    N={1,2};
    test_index(N);
    N={4,1};
    test_index(N);
    N={4,2,3,5};
    test_index(N);
}