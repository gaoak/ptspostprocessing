#include<string>
#include<vector>
#include "PlungingMotion.h"
using namespace std;

int main(int argc, char *argv[]) {
    string dataconfigue("dataconfigue");
    int mode = 0;
    for(int c=1; c<argc; ++c) {
        if(0==string("data").compare(argv[c])) {
            dataconfigue = argv[c+1];
        }
        if(0==string("exp").compare(argv[c])) {
            mode = 1;
        }
        if(0==string("airfoil").compare(argv[c])) {
            mode = 2;
        }
    }
    PlungingMotion plungdata(dataconfigue);
    int dir = 1;
    for(int c=1; c<argc; ++c) {
        if(0==string("processinverse").compare(argv[c])) {
            dir = -1;
        }
    }
    if(mode==1) {
        plungdata.ProcessEXPWingData(dir);
    } else if(mode==0) {
        plungdata.ProcessCFDWingData(dir);
    } else if(mode==2) {
        plungdata.ProcessCFDAirfoilData(dir);
    }
    printf("finished\n");
    return 0;
}