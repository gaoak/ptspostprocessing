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
            printf("error: parser double %s\n", cuts.c_str());
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