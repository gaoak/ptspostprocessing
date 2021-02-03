#include<string>
#include<vector>
#include "PlungingMotion.h"
using namespace std;

int main(int argc, char *argv[]) {
    string dataconfigue("dataconfigue");
    for(int c=1; c<argc; ++c) {
        if(0==string("data").compare(argv[c])) {
            dataconfigue = argv[c+1];
        }
    }
    PlungingMotion plungdata(dataconfigue);
    int dir = 1;
    for(int c=1; c<argc; ++c) {
        if(0==string("processinverse").compare(argv[c])) {
            dir = -1;
        }
    }
    plungdata.ProcessFlowData(dir);
    printf("finished\n");
    return 0;
}