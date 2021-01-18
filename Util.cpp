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

void LeftShiftArray(const std::vector<int> &a, std::vector<int> &res) {
    for(int i=0; i<res.size(); ++i) {
        res[i] = a[(i+1)%res.size()];
    }
}

int LeftShiftIndex(std::vector<int> &N, std::vector<int> &oN, const double *data, double * odata) {
    // i,j,k to j,k,i
    size_t dim = N.size();
    LeftShiftArray(N, oN);
    int Np = 1;
    for(int i=0; i<dim; ++i) {
        Np *= N[i];
    }
    std::vector<int> ind(dim);
    std::vector<int> indo(dim);
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        LeftShiftArray(ind, indo);
        odata[Index(oN, indo)] = data[i];
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