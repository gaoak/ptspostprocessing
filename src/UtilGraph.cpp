#include "UtilGraph.h"
#include "Util.h"
#include<limits>
#include<cmath>

int PurgeDifferentSign(const std::vector<int> &N, const std::vector<double> &v,
                       std::vector<double> &data, double sign) {
    int Np = 1;
    for(size_t i=0; i<N.size(); ++i) Np *= N[i];
    if(Np>(int)data.size()) Np = data.size();
    for(int i=0; i<Np; ++i) {
        if(v[i]*sign < 0.) {
            data[i] = 0.;
        }
    }
    return Np;
}

int ExtractPatchStat2DGraph(const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
    std::vector<double> &v, std::vector<double> &mu, double &sum) {
    Fill2DGraph(N, v, core, 0.05, true);
    double sum0 = 0., sum1x = 0., sum1y = 0.;
    for(int j=0; j<N[1]; ++j) {
        int tmp = j * N[0];
        for(int i=0; i<N[0]; ++i) {
            int ind = i + tmp;
            sum0  += v[ind];
            sum1x += v[ind] * i;
            sum1y += v[ind] * j;
        }
    }
    if(std::fabs(sum0) < std::numeric_limits<double>::epsilon() ) {
        mu.clear();
        mu.push_back(0.);
        mu.push_back(0.);
        sum = 0.;
        return 0;
    }

    sum1x /= sum0;
    sum1y /= sum0;
    double sumxx = 0., sumxy = 0., sumyy = 0.;
    for(int j=0; j<N[1]; ++j) {
        int tmp = j * N[0];
        for(int i=0; i<N[0]; ++i) {
            int ind = i + tmp;
            sumxx += v[ind] * (i - sum1x) * (i - sum1x);
            sumxy += v[ind] * (i - sum1x) * (j - sum1y);
            sumyy += v[ind] * (j - sum1y) * (j - sum1y);
        }
    }
    sumxx *= dx[0] * dx[0];
    sumxy *= dx[0] * dx[1];
    sumyy *= dx[1] * dx[1];
    mu.clear();
    mu.push_back(sumxx / sum0);
    mu.push_back(sumxy / sum0);
    mu.push_back(sumyy / sum0);
    sum = sum0 * dx[0] * dx[1];
    return 4;
}

int Fill2DGraph(const std::vector<int> &rawN, std::vector<double> &value, const std::vector<int> &init, const double &eps, bool monotone)
{
    if(eps>=1. || eps <= 0.) {
        printf("error: threshold value %f should be in (0, 1)\n", eps);
        return -1;
    }
    std::vector<int> N(2);
    N[0] = rawN[0];
    N[1] = rawN[1];
    double sign = 1.;
    std::vector<int> intcore(2);
    if(value[Index(N, init)] < 0.) {
        sign = -1.;
        FindLocMaxIn2DGraph(N, init, value, intcore, false);
    } else {
        FindLocMaxIn2DGraph(N, init, value, intcore, true);
    }
    double threshold = value[Index(N, intcore)] * eps * sign;
    //prepare memory
	std::vector<double> val(N[0]*N[1], 0);
	std::vector<int> res(N[0]*N[1], 0);
    for(int j=0; j<N[1]; ++j) {
        int tmp = j* N[0];
        for(int i=0; i<N[0]; ++i) {
            int ind = i + tmp;
            if(value[ind] * sign > threshold) {
                val[ind] = value[ind] * sign;
            } else {
                val[ind] = -1.;
            }
        }
    }
    //start fill
    res[Index(N, intcore)] = 1;
    bool proceed = true;
    while(proceed) {
        proceed = false;
        for(int j=0; j<N[1]; ++j) {
            int tmp = j* N[0];
            for(int i=0; i<N[0]; ++i) {
                int ind = i + tmp;
                if(res[ind]) continue;
                std::vector<int> p;
                if(i>0) {
                    p.push_back(ind - 1);
                }
                if(i<N[0]-1) {
                    p.push_back(ind + 1);
                }
                if(j>0) {
                    p.push_back(ind - N[0]);
                }
                if(j<N[1]-1) {
                    p.push_back(ind + N[0]);
                }
                for(int n=0; n<(int)p.size(); ++n) {
                    int indt = p[n];
                    bool checkmono = !monotone || val[ind] <= val[indt];
                    if( res[indt] && val[ind]>0. && checkmono) {
                        res[ind] = 1;
                        proceed = true;
                        break;
                    }
                }
            }
        }
    }

    int count = 0;
    for(int j=0; j<N[1]; ++j) {
        int tmp = j * N[0];
        for(int i=0; i<N[0]; ++i) {
            int ind = i + tmp;
            if(res[ind] == 0) {
                value[ind] = 0.;
            } else {
                ++count;
            }
        }
    }
	return count;
}

int FindLocMaxIn2DGraph(const std::vector<int> &N, const std::vector<int> &initial,
    std::vector<double> &data, std::vector<int> &core, bool ismax) {
    if(!ismax) {
        for(int i=0; i<N[0]*N[1]; ++i) {
            data[i] = -data[i];
        }
    }
    int Np = N[0] * N[1];
    if(initial.size()<2) {
        int imax = FindMax<double>(Np, data.data());
        std::vector<int> p;
        invIndex(N, imax, p);
        core = p;
    } else {
        std::vector<int> p = initial;
        double maxv = data[Index(N, p)];
        bool proceed = true;
        int tmp;
        while(proceed) {
            proceed = false;
            if(p[0]>0 && data[tmp = Index(N, {p[0]-1, p[1]})] > maxv) {
                maxv = data[tmp];
                p[0] -= 1;
                proceed = true;
            } else if(p[0]<N[0]-1 && data[tmp = Index(N, {p[0]+1, p[1]})] > maxv) {
                maxv = data[tmp];
                p[0] += 1;
                proceed = true;
            }
            if(p[1]>0 && data[tmp = Index(N, {p[0], p[1]-1})] > maxv) {
                maxv = data[tmp];
                p[1] -= 1;
                proceed = true;
            } else if(p[1]<N[1]-1 && data[tmp = Index(N, {p[0], p[1]+1})] > maxv) {
                maxv = data[tmp];
                p[1] += 1;
                proceed = true;
            }
        }
        core = p;
    }
    core.resize(2);
    if(!ismax) {
        for(int i=0; i<N[0]*N[1]; ++i) {
            data[i] = -data[i];
        }
    }
    return core.size();
}