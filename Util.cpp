#include<vector>
#include<string>
#include "Util.h"
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

void ShiftArray(std::vector<int> &a, int dir) {
    if(a.size()==0) {
        return;
    }
    if(dir<0) {
        dir = (a.size() + dir % a.size()) % a.size();
    } else {
        dir = dir % a.size();
    }
    std::vector<int> tmp = a;
    for(int i=0; i<a.size(); ++i) {
        a[(i+dir)%a.size()] = tmp[i];
    }
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
    ShiftArray(sN, dir);
    CompactArray(N);
    CompactArray(sN);
    for(int i=0; i<N.size(); ++i) {
        if(N[i]!=sN[i]) {
            return true;
        }
    }
    return false;
}

int ShiftIndex(std::vector<int> &N, std::vector<std::vector<double> > &odata, int dir) {
    // i,j,k to j,k,i
    if(odata.size()==0) {
        return 0;
    }
    std::vector<int> sN = N;
    ShiftArray(N, dir);
    if(!NeedShift(sN, dir)) {
        return 0;
    }

    size_t dim = sN.size();
    int Np = 1;
    for(int i=0; i<dim; ++i) {
        Np *= sN[i];
    }
    std::vector<std::vector<double> > data;
    for(int i=0; i<odata.size(); ++i) {
        data.push_back(odata[i]);
    }
    
    std::vector<int> ind(dim);
    std::vector<int> indo(dim);
    for(int i=0; i<Np; ++i) {
        invIndex(sN, i, ind);
        indo = ind;
        ShiftArray(indo, dir);
        int it = Index(N, indo);
        for(int n=0; n<data.size(); ++n) {
            odata[n][it] = data[n][i];
        }
    }
    return Np;
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