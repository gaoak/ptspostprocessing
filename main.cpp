#include<iostream>
#include<string>
#include<vector>
#include<map>
#include "Dataprocessing.h"
#include "Util.h"
#include "StructuredData.h"
#include "PlungingMotion.h"
using namespace std;
std::map<bool, std::string> testresults = {{true, "pass"},{false, "fail"}};

int test_weight(int padding) {
    KernelSmooth ker;
    double sum = 0.;
    for(int i=-padding; i<=padding; ++i) {
        double tmp = ker.GW(i, padding);
        printf("%g,",tmp);
        sum += tmp;
    }
    printf("\ntest weight: sum %g, %s\n", sum-1., testresults[fabs(sum-1.)<1E-14].c_str());
}

int test_index(std::vector<int> &N) {
    int Np = N[0];
    for(int i=1; i<N.size(); ++i) {
        Np *= N[i];
    }
    std::vector<int> ind(N.size());
    int failcount = 0;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        int tmp = Index(N, ind);
        printf("%d, (", i);
        for(int j=0; j<ind.size(); ++j) {
            printf("%d,",ind[j]);
        }
        printf("), %d\n", tmp);
        failcount += i != tmp;
    }
    printf("test index: %d, %s\n", failcount, testresults[!failcount].c_str());
}
double highfreq(std::vector<double> p) {
    return 0.;
    return sin(200.*M_PI*p[2]) + sin(200.*M_PI*p[0]) + sin(200.*M_PI*p[1]);
}
double lowfreq(std::vector<double> p) {
    return sin(2.*M_PI*p[1]) + sin(2.*M_PI*p[2]);
}
double deriv(std::vector<double> p) {
    return 2.*M_PI*cos(2.*M_PI*p[1]) - p[4];
}
double sinfunc(std::vector<double> p) {
     return highfreq(p) + lowfreq(p);
}
double exact(std::vector<double> p) {
    return lowfreq(p);
}
int test_structuredData() {
    std::vector<int> N = {12, 128 + 1, 2};
    std::vector<double> range = {0., 1., 0., 1., 0., 1.};
    StructuredData sdata(N, range);
    //sdata.OutputCSV("coor.csv");
    //sdata.OutputTec360("coor.plt");
    sdata.AddPhysics("sin", (void*) sinfunc);
    //sdata.OutputTec360("addphys.plt");
    std::vector<int> field = {0};
    //sdata.Smoothing(0.02, field, false);

    field[0] = 0;
    sdata.Diff(field, 1, 6);
    //field = {2,3};
    //sdata.Smoothing(0.02, field, true);

    sdata.AddPhysics("deriv", (void*) deriv);
    std::map<int, double> def = {{2, 2.*M_PI}, {3, 2.*M_PI}};
    field = {2};
    def = {{2,0.}};
    sdata.MaskBoundary(0.1, field, def);
    sdata.OutputTec360("smoothphys.plt");
    printf("sum %g\n", sdata.GetPhysNorm(2, -1));
}

int TestSummary() {
    test_weight(6);
    std::vector<int> N = {4,2,3};
    test_index(N);
    N={2,1,3};
    test_index(N);
    N={1,2,3};
    test_index(N);
    N={1,2,1};
    test_index(N);
    N={4,2,1};
    test_index(N);
    N={2,3};
    test_index(N);
    N={1,2};
    test_index(N);
    N={4,1};
    test_index(N);
    N={4,2,3,5};
    test_index(N);
}

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
            plungdata.ProcessFlowData();
        }
    }
    printf("finished\n");
    return 0;
}