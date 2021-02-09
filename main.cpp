#include<string>
#include<vector>
#include "PlungingMotion.h"
using namespace std;

int main(int argc, char *argv[]) {
    string dataconfigue("dataconfigue");
    bool expdata = false;
    for(int c=1; c<argc; ++c) {
        if(0==string("data").compare(argv[c])) {
            dataconfigue = argv[c+1];
        }
        if(0==string("exp").compare(argv[c])) {
            expdata = true;
        }
    }
    PlungingMotion plungdata(dataconfigue);
    int dir = 1;
    for(int c=1; c<argc; ++c) {
        if(0==string("processinverse").compare(argv[c])) {
            dir = -1;
        }
    }
    if(expdata) {
        plungdata.ProcessEXPWingData(dir);
    } else {
        plungdata.ProcessCFDWingData(dir);
    }
    printf("finished\n");
    return 0;
}