#include<fstream>
#include<map>
#include<cmath>
#include<algorithm>
#include<set>
#include "IncFlow.h"
#include "Util.h"
using namespace std;

double func(std::vector<double> p) {
    p[0] = 0.;
    return p[0];
}

int main(int argc, char *argv[]) {
    if (argc<4) {
        printf("argc is %d [<4]\n", argc);
    }
    string baseflowname(argv[1]);
    double amplitude = 0.05;
    sscanf(argv[2], "%lf", &amplitude);
    string modename(argv[3]);

    IncFlow baseflow;
    IncFlow modeflow;
    baseflow.InputData(baseflowname, true);
    modeflow.InputData(modename, true);
    if (modeflow.GetNumPhys()==2) {
        modeflow.AddPhysics("w", (void*)func);
    }
    for (int f=0; f<baseflow.GetNumPhys() - 1; ++f) { //no p
        for (int p=0; p < baseflow.GetTotPoints(); ++p) {
            double v = baseflow.GetPhysValue(f, p) + amplitude * modeflow.GetPhysValue(f, p);
            modeflow.SetPhysValue(v, f, p);
        }
    }
    int fw = modeflow.GetNumPhys() - 1;
    for (int p=0; p < baseflow.GetTotPoints(); ++p) {
        double v = amplitude * modeflow.GetPhysValue(fw, p);
        modeflow.SetPhysValue(v, fw, p);
    }
    modeflow.CalculateVorticity();
    modeflow.OutputData("summaryflow.plt");
    return 0;
}