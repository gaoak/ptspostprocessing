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
    plungdata.Dumppoints();
    printf("finished\n");
    return 0;
}