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
    std::vector<double> sigma, value;
    for(int c=1; c<argc; ++c) {
        if(0==string("dump").compare(argv[c])) {
            plungdata.Dumppoints();
        }
        if(0==string("process").compare(argv[c])) {
            plungdata.ProcessFlowData(1);
        }
        if(0==string("processinverse").compare(argv[c])) {
            plungdata.ProcessFlowData(-1);
        }
    }
    printf("finished\n");
    return 0;
}