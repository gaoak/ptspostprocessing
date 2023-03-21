#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<tuple>
#include "Dataprocessing.h"
#include "Util.h"
#include "StructuredData.h"
#include "IncFlow.h"

using namespace std;

std::map<bool, std::string> testresults = {{true, "pass"},{false, "+++++++fail!!!!!!!!!!!!!!!!!!!!!!"}};

int testfftw() {
    std::vector<double> data0(8, 0.);
    for (size_t i=0; i<data0.size(); ++i)
    {
        double t = i * 2. * M_PI / data0.size();
        data0[i] = sin(t) + 0.5* cos(3 * t)  + sin(3 * t) + 2.;
    }
    std::vector<double> spectral(data0.size());
    std::vector<double> beta(data0.size());
    doRealFFT(data0.data(), spectral.data(), beta.data(), 2.*M_PI, data0.size());
    for (size_t i=0; i<data0.size(); ++i)
    {
        printf("%3lu: %f, %f, %g\n", i, data0[i], beta[i], spectral[i]);
    }
    return 0;
}

int test_weight(int padding) {
    KernelSmooth ker;
    double sum = 0.;
    for(int i=-padding; i<=padding; ++i) {
        double tmp = ker.GW(i, padding);
        sum += tmp;
    }
    printf("\ntest weight: sum %g, %s\n", sum-1., testresults[fabs(sum-1.)<1E-14].c_str());
    return 0;
}

int test_index(std::vector<int> &N) {
    int Np = N[0];
    for(int i=1; i<(int)N.size(); ++i) {
        Np *= N[i];
    }
    std::vector<int> ind(N.size());
    int failcount = 0;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        int tmp = Index(N, ind);
        //printf("%d, (", i);
        for(int j=0; j<(int)ind.size(); ++j) {
            ;//printf("%d,",ind[j]);
        }
        //printf("), %d\n", tmp);
        failcount += i != tmp;
    }
    printf("test index: %d, %s\n", failcount, testresults[!failcount].c_str());
    return 0;
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
    return 0;
}

double highfreq(std::vector<double> p) {
    return 0.*(sin(200.*M_PI*p[2]) + sin(200.*M_PI*p[0]) + sin(200.*M_PI*p[1]));
}
double lowfreq(std::vector<double> p) {
    return sin(2.*M_PI*p[1]) * sin(2.*M_PI*p[2]);
}
double deriv(std::vector<double> p) {
    return 2.*M_PI*cos(2.*M_PI*p[1]) * sin(2.*M_PI*p[2]);
}
double sinfunc(std::vector<double> p) {
     return highfreq(p) + lowfreq(p);
}
double sinfunc2d(std::vector<double> p) {
     return sin(2.*M_PI*p[0]) * sin(2.*M_PI*p[1]);
}
double exact(std::vector<double> p) {
    return lowfreq(p);
}
double linear(std::vector<double> p) {
    return p[2]*p[1]*p[0];
}
int test_structuredData() {
    std::vector<int> N = {33,65,33};
    std::vector<double> range = {-0.5,2.,-0.5, 2.,0,5.5};
    //origin0, length0, origin1, length1, origin2, length2
    StructuredData sdata(N, range);
    sdata.OutputData("coor.csv");
    sdata.OutputData("coor.plt");
    sdata.AddPhysics("sin", (void(*)()) sinfunc);
    double sum = sdata.GetPhysNorm(0,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-0.238695)<1E-6].c_str());
    std::vector<int> field = {0};
    sdata.Smoothing(0.02, field, true);
    sum = sdata.GetPhysNorm(0,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-0.238695)<1E-6].c_str());

    field[0] = 0;
    sdata.Diff(field, 1, 6);
    sum = sdata.GetPhysNorm(1,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-9.70294)<1E-5].c_str());

    sdata.AddPhysics("deriv", (void(*)()) deriv);
    sum = sdata.GetPhysNorm(2,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-9.71776)<1E-5].c_str());
    std::map<int, double> def = {{2, 2.*M_PI}, {3, 2.*M_PI}};
    field = {2};
    def = {{2,0.}};
    sdata.MaskBoundary(0.1, field, def);
    sdata.OutputData("testfile.plt");
    sum = sdata.GetPhysNorm(2,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-7.08265)<1E-5].c_str());
    sum = sdata.GetPhysNorm(2, -1);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-1.)<1E-10].c_str());
    //extract plane
    range = {0.0,0.5,-0.5,0.5,0,2.75};
    std::vector<double> planedata;
    std::pair<int, int> plane = std::make_pair(0, 16);
    ShiftArray<double>(range, -2);
    std::vector<int> planeN;
    std::vector<double> phys2(sdata.GetTotPoints());
    for(size_t i=0; i<phys2.size(); ++i) {
        phys2[i] = sdata.GetPhysValue(2, i);
    }
    sdata.ExtractPlane(phys2, plane, {8,16,0,32,0,16}, planeN, planedata);
    StructuredData tmp(planeN, range);
    tmp.AddPhysics("deriv", planedata);
    tmp.OutputData("plane0.plt");
    sum = tmp.GetPhysNorm(0,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-8.48457)<1E-5].c_str());
    plane = std::make_pair(1, 32);
    ShiftArray<double>(range, -2);
    sdata.ExtractPlane(phys2, plane, {8,16,0,32,0,16}, planeN, planedata);
    StructuredData tmp1(planeN, range);
    tmp1.AddPhysics("deriv", planedata);
    tmp1.OutputData("plane1.plt");
    sum = tmp1.GetPhysNorm(0,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-19.7392)<1E-4].c_str());
    plane = std::make_pair(2, 16);
    ShiftArray<double>(range, -2);
    sdata.ExtractPlane(phys2, plane, {8,16,0,32,0,16}, planeN, planedata);
    StructuredData tmp2(planeN, range);
    tmp2.AddPhysics("deriv", planedata);
    tmp2.OutputData("plane2.plt");
    sum = tmp2.GetPhysNorm(0,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-16.9691)<1E-4].c_str());
    //3D rotation data
    sdata = StructuredData(N, range, {{1,1,0.},{-1,1,0},{1,1,1}});
    sdata.AddPhysics("sin", (void(*)()) sinfunc2d);
    field = {0};
    sdata.Smoothing(0.02, field, true);
    sdata.OutputData("3drot.plt");
    ///////////test 2D
    std::vector<int> N2d = {33,65};
    std::vector<double> range2d = {-0.5,2.,-0.5,2.};
    StructuredData sdata2d(N2d, range2d);
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(0, Index(N2d, {32,64}))-1.5)<1E-4].c_str());
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(1,Index(N2d, {32,64}))-1.5)<1E-4].c_str());
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(0,0)+0.5)<1E-4].c_str());
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(1, 0)+0.5)<1E-4].c_str());
    sdata2d.OutputData("2d.plt");
    sdata2d = StructuredData(N2d, range2d, {{1,1},{-1,1}});
    sdata2d.AddPhysics("sin", (void(*)()) sinfunc2d);
    field = {0};
    sdata2d.Smoothing(0.02, field, true);
    sdata2d.OutputData("2drot.plt");
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(0,Index(N2d, {32,64}))+0.5)<1E-4].c_str());
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(1,Index(N2d, {32,64}))+0.5-2*sqrt(2.0))<1E-4].c_str());
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(0,0)+0.5)<1E-4].c_str());
    printf("test structuredData 2d %s\n", testresults[fabs(sdata2d.GetCoordValue(1,0)+0.5)<1E-4].c_str());
    sum = sdata2d.GetPhysNorm(0,2);
    printf("test structuredData %g, %s\n", sum, testresults[fabs(sum-0.245466)<1E-6].c_str());
    return 0;
}

int test_fileio() {
    StructuredData sdata({33,65,33},{-0.5,2.,-0.5,2.,0,5.5});
    sdata.InputData("testfile.plt");
    double sum = sdata.GetPhysNorm(2,2);
    printf("test fileio %g, %s\n", sum, testresults[fabs(sum-7.08265)<1E-5].c_str());
    sdata.OutputData("plt.dat");
    sdata.OutputData("plt.csv");
    sdata.OutputData("plt.plt");
    //
    sdata.InputData("plt.csv");
    sum = sdata.GetPhysNorm(2,2);
    printf("test fileio %g, %s\n", sum, testresults[fabs(sum-7.08265)<1E-5].c_str());
    sdata.OutputData("plt_csv.dat");
    //
    sdata.InputData("plt.plt");
    sum = sdata.GetPhysNorm(2,2);
    printf("test fileio %g, %s\n", sum, testresults[fabs(sum-7.08265)<1E-5].c_str());
    sdata.OutputData("plt_plt.dat");
    //2d io
    sdata.InputData("2drot.plt");
    sum = sdata.GetPhysNorm(0,2);
    printf("test fileio %g, %s\n", sum, testresults[fabs(sum-0.245466)<1E-6].c_str());
    sdata.OutputData("2drot_reload.dat");
    sdata.OutputData("2drot_reload.plt");
    //3d rot
    sdata.InputData("3drot.plt");
    sdata.OutputData("3drot_reload.plt");
    sdata.OutputData("3drot_reload.dat");
    return 0;
}

int test_subdomain() {
    StructuredData sdata({33,65,33},{-0.5,2.,-0.5,2.,0,5.5});
    std::map<int, double> field = {{0,0.}, {1,0.}, {2,0.}};
    sdata.InputData("testfile.plt");
    double sum = sdata.GetPhysNorm(2,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-7.08265)<1E-5].c_str());
    StructuredData sub1;
    sub1.CopyAsSubDomain({8,16,8},{17,33,17},{1,1,1}, field, sdata);
    sub1.OutputData("sub1.plt");
    sum = sub1.GetPhysNorm(2,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-10.7204)<1E-4].c_str());
    sub1.CopyAsSubDomain({8,16,8},{17,33,17},{2,2,1}, field, sdata);
    sub1.OutputData("sub2.plt");
    sum = sub1.GetPhysNorm(2,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-11.2498)<1E-4].c_str());
    //
    IncFlow sflow({33,65,33},{-0.5,2.,-0.5,2.,0,5.5}, "0012", {15,0.,5.});
    sflow.InputData("testfile.plt");
    sum = sflow.GetPhysNorm(2,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-7.08265)<1E-5].c_str());
    IncFlow sub2;
    sub2.CopyAsSubDomain({8,16,8},{17,33,17},{1,1,1}, field, sflow);
    sum = sub2.GetPhysNorm(2,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-10.7204)<1E-4].c_str());
    sub2.OutputData("ssub1.plt");
    sub2.CopyAsSubDomain({8,16,8},{17,33,17},{2,2,1}, field, sflow);
    sum = sub2.GetPhysNorm(2,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-11.2498)<1E-4].c_str());
    sub2.OutputData("ssub2.plt");
    //test coarsen
    StructuredData coar;
    field = {{1,0.}, {2,0.}, {3,0},{4,0},{9,0}};
    coar.CopyAsSubDomain({0,0,0},{33,65,33},{2,2,2}, field, sdata);
    sum = coar.GetPhysNorm(1,2);
    printf("test subdomain %g, %s\n", sum, testresults[fabs(sum-6.60311)<1E-5].c_str());
    coar.OutputData("coarsen.plt");
    return 0;
}

int test_interpolation() {
    std::map<int, double> field = {{0,0.},{1,0},{2,0},{3,0},{4,0}};
    StructuredData sdata({17,17,17},{-0.5,2.,-0.5,2.,0,5.5}), data2;
    sdata.AddPhysics("linear", (void(*)())linear);
    sdata.AddPhysics("sin", (void(*)())sinfunc);
    double sum = sdata.GetPhysNorm(0,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-4.06189)<1E-5].c_str());
    sum = sdata.GetPhysNorm(1,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-0.221453)<1E-6].c_str());
    sdata.OutputData("interporigin.plt");
    data2 = sdata;
    data2.InterpolateFrom(sdata, field);
    sum = data2.GetPhysNorm(0,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-4.06189)<1E-5].c_str());
    sum = data2.GetPhysNorm(1,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-0.221453)<1E-6].c_str());
    data2.OutputData("interp3dsame.plt");
    data2 = StructuredData({65,65,129},{-0.5,2.,-0.5,2.,0,5.5}, {{1,1,0.},{-1,1,1},{1,0,1}});
    data2.InterpolateFrom(sdata, field);
    sum = data2.GetPhysNorm(0,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-0.895969)<1E-6].c_str());
    sum = data2.GetPhysNorm(1,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-0.0436265)<1E-7].c_str());
    data2.OutputData("interp3d.plt");
    data2 = StructuredData({65,65,1},{-0.5,2.,-0.5,2.,1}, {{1,1,0.},{-1,1,1},{1,0,1}});
    data2.InterpolateFrom(sdata, field);
    sum = data2.GetPhysNorm(0,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-0.0938391)<1E-7].c_str());
    sum = data2.GetPhysNorm(1,2);
    printf("test interpolation %g, %s\n", sum, testresults[fabs(sum-0.0512146)<1E-7].c_str());
    data2.OutputData("interp2d.plt");
    //2D interp
    StructuredData subdata2({22,22,1},{-0.5,2.,-0.5,2.,1}, {{1,1,0.},{-1,1,1},{1,0,1}});
    subdata2.InterpolateFrom(data2, field);
    subdata2.OutputData("interpsub2d.plt");
    //point
    std::map<int, double> value;
    sdata.InterpolatePoint ({-0.595777387263203506, 0.919303221756027966, 1.75754030295902952},
                            field, value);
    printf("test interpolation %s\n", testresults[fabs(value[0]-0.)<1E-5].c_str());
    printf("test interpolation %s\n", testresults[fabs(value[1]-0.)<1E-5].c_str());
    return 0;
}

int main() {
    TestSummary();
    test_structuredData();
    test_fileio();
    test_subdomain();
    test_interpolation();
    return 0;
}