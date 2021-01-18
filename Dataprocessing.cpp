#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include "Dataprocessing.h"

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
    double sum = 0.;
    for(int i=-0; i<=m_maxOffset; ++i) {
        std::vector<double> tmpw(i+1);
        for(int j=0; j<=i; ++j) {
            tmpw[j] = gauss[j]/sum[i];
        }
        m_weight.push_back(tmpw);
    }
}

double GaussKernel::Sum(double *data, int padding) {
    double sum = *data;
    if(padding<0 || padding>m_maxOffset) {
        printf("error: illegal padding %d, max[%d]\n", padding, m_maxOffset);
        return sum;
    }
    sum = (*data) * m_weight[padding][0];
    for(int i=1; i<=padding; ++i) {
        sum += m_weight[padding][i] * (*(data-i) + *(data+i));
    }
    return sum;
}

int GaussKernel::DoSmooth(std::vector<double> &data, std::vector<double> & sdata) {
    ;
}

int GaussKernel::GetIMax() {
    return m_maxOffset;
}

int main() {
    GaussKernel ker(1., 0.3);
    double sum = 0.;
    for(int i=-ker.GetIMax(); i<=ker.GetIMax(); ++i) {
        double tmp = ker.GW(i);
        printf("%g,",tmp);
        sum+= tmp;
    }
    printf("sum %g\n", sum-1.);
}