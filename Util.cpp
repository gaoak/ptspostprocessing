#include<vector>
#include<string>
#include<cmath>
#include "Util.h"

double Distance(const std::vector<double> &a, const std::vector<double> &b) {
    double res = 0.;
    for(int i=0; i<a.size(); ++i) {
        res += (b[i]-a[i])*(b[i]-a[i]);
    }
    return std::sqrt(res);
}

int FindAbsMax(int N, const double* data) {
    int ind=0;
    double tmp = std::fabs(data[0]);
    for(int i=1; i<N; ++i) {
        if(tmp < std::fabs(data[i])) {
            ind = i;
            tmp = std::fabs(data[i]);
        }
    }
    return ind;
}

std::vector<double> transform(const std::vector<double> &p, double AoA) {
    std::vector<double> res = p;
    res[0] = p[0]*cos(AoA) + p[1]*sin(AoA);
    res[1] =-p[0]*sin(AoA) + p[1]*cos(AoA);
    return res;
}

void parserUInt(const char * cstr, std::vector<int> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if(cstr[i]>='0' && cstr[i]<='9') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    int k;
    for(int i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);// data  in [s e-1]
        if(sscanf(cuts.c_str(), "%d", &k)<1) {
            printf("error: parser int %s \n",  cuts.c_str());
        }
        if(i>0 && (digs[i] - dige[i-1])==1 && cstr[digs[i]-1]=='-') {
            for(int j=value[value.size()-1]+1; j<k; ++j) {
                value.push_back(j);
            }
        }
        value.push_back(k);
    }
}

void parserString(const char * cstr, std::vector<std::string> & value, char sep) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if(cstr[i]!=sep && cstr[i]!=0) {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    double k;
    for(int i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);
        value.push_back(cuts);
    }
}

void parserDouble(const char * cstr, std::vector<double> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if((cstr[i]>='0' && cstr[i]<='9') ||
            cstr[i]=='.' ||
            cstr[i]=='e' || cstr[i]=='E' ||
            cstr[i]=='+' || cstr[i]=='-') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    double k;
    for(int i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);
        if(sscanf(cuts.c_str(), "%lf", &k)<1) {
            printf("error: parser double %s, in \n%s\n", cuts.c_str(), cstr);
        }
        value.push_back(k);
    }
}

double StringToDouble(std::string str) {
    double k;
    if(sscanf(str.c_str(), "%lf", &k)<1) {
        printf("error: parser double %s\n", str.c_str());
    }
    return k;
}

bool CompactArray(std::vector<int> &N) {
    std::vector<int> sN = N;
    N.clear();
    for(int i=0; i<sN.size(); ++i) {
        if(1!=sN[i]) {
            N.push_back(sN[i]);
        }
    }
    return N.size()!=sN.size();
}

bool NeedShift(std::vector<int> N, int dir) {
    for(int i=0; i<N.size(); ++i) {
        if(N[i]!=1) {
            N[i] = 10 + i;
        }
    }
    std::vector<int> sN = N;
    ShiftArray<int>(sN, dir);
    CompactArray(N);
    CompactArray(sN);
    for(int i=0; i<N.size(); ++i) {
        if(N[i]!=sN[i]) {
            return true;
        }
    }
    return false;
}

int Index(const std::vector<int> &N, const std::vector<int> & index) {
    int res = index[0];
    int Np = N[0];
    for(int i=1; i<N.size(); ++i) {
        if(index[i]<0 || index[i]>=N[i]) {
            printf("error: illegal index %d in Index %d\n", index[i], N[i]);
        }
        res += Np*index[i];
        Np *= N[i];
    }
    return res;
}

void invIndex(const std::vector<int> &N, int index, std::vector<int> & res) {
    res.resize(N.size());
    res[0] = N[0];
    for(int i=1; i<N.size()-1; ++i) {
        res[i] = res[i-1] * N[i];
    }
    for(int i=res.size()-1; i>0; --i) {
        res[i] = index / res[i-1];
        index = index % res[i-1];
    }
    res[0] = index;
}

int Fill2DGraph(const std::vector<int> &rawN, std::vector<double> &value, const std::vector<int> &init, const double &eps, bool monotone)
{
    if(eps>=1. || eps <= 0.) {
        printf("error: threshold value %f should be in (0, 1)\n", eps);
        return -1;
    }
    std::vector<int> N(2);
    std::vector<double> core;
    N[0] = rawN[0];
    N[1] = rawN[1];
    double sign = 1.;
    if(value[Index(N, init)] < 0.) {
        sign = -1.;
        FindLocMaxIn2DGraph(N, init, value, core, false);
    } else {
        FindLocMaxIn2DGraph(N, init, value, core, true);
    }
    std::vector<int> intcore(2);
    intcore[0] = myRound(core[0]);
    intcore[1] = myRound(core[1]);
    double threshold = value[Index(N, intcore)] * eps * sign;
    //prepare memory
	std::vector<double> val(N[0]*N[1], 0);
	std::vector<int> res(N[0]*N[1], 0);
    for(int j=0; j<N[1]; ++j) {
        for(int i=0; i<N[0]; ++i) {
            int ind = Index(N, {i, j});
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
            for(int i=0; i<N[0]; ++i) {
                int ind = Index(N, {i, j});
                if(res[ind]) continue;
                std::vector<int> p;
                if(i>0) {
                    p.push_back(Index(N, {i-1, j}));
                }
                if(i<N[0]-1) {
                    p.push_back(Index(N, {i+1, j}));
                }
                if(j>0) {
                    p.push_back(Index(N, {i, j-1}));
                }
                if(j<N[1]-1) {
                    p.push_back(Index(N, {i, j+1}));
                }
                for(int n=0; n<p.size(); ++n) {
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
        for(int i=0; i<N[0]; ++i) {
            int ind = Index(N, {i, j});
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
    std::vector<double> &data, std::vector<double> &core, bool ismax) {
    if(!ismax) {
        for(int i=0; i<N[0]*N[1]; ++i) {
            data[i] = -data[i];
        }
    }
    int Np = N[0] * N[1];
    if(initial.size()<2) {
        int imin = FindMin<double>(Np, data.data());
        double pmin = data[imin];
        double thresh = 0.98 * pmin;
        std::vector<double> center;
        DoMaskShift<double>(Np, thresh, -1, data.data());
        core.clear();
        WeightedCenter<double>(N, data.data(), core);
    } else {
        std::vector<int> p = initial;
        double maxv = data[Index(N, p)];
        bool proceed = true;
        while(proceed) {
            proceed = false;
            if(p[0]>0 && data[Index(N, {p[0]-1, p[1]})] > maxv) {
                maxv = data[Index(N, {p[0]-1, p[1]})];
                p[0] -= 1;
                proceed = true;
            }
            if(p[0]<N[0]-1 && data[Index(N, {p[0]+1, p[1]})] > maxv) {
                maxv = data[Index(N, {p[0]+1, p[1]})];
                p[0] += 1;
                proceed = true;
            }
            if(p[1]>0 && data[Index(N, {p[0], p[1]-1})] > maxv) {
                maxv = data[Index(N, {p[0], p[1]-1})];
                p[1] -= 1;
                proceed = true;
            }
            if(p[1]<N[1]-1 && data[Index(N, {p[0], p[1]+1})] > maxv) {
                maxv = data[Index(N, {p[0], p[1]+1})];
                p[1] += 1;
                proceed = true;
            }
        }
        core.resize(2);
        core[0] = p[0];
        core[1] = p[1];
    }
    if(!ismax) {
        for(int i=0; i<N[0]*N[1]; ++i) {
            data[i] = -data[i];
        }
    }
    return core.size();
}