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

int SwitchCoord(std::vector<int> &N, std::vector<double> & data) {
    for(int i=N.size(); i<3; ++i) {
        N.push_back(1);
    }
    int Np = N[0] * N[1] * N[2];
    if(data.size() < Np) {
        printf("error: data (%d) is not enought for SwithCoord (%d)\n", (int)data.size(), Np);
    }
    // x, y, z -> y, z, x
}