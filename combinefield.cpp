#include<fstream>
#include<map>
#include<cmath>
#include<algorithm>
#include<set>
#include "IncFlow.h"
#include "Util.h"
using namespace std;
/************
 * description
 * combine two fields by appending new physics fields
 * Usage:
 * combinefield field1.plt field2.plt
 ***********/

int main(int argc, char *argv[]) {
    if (argc<3) {
        printf("argc is %d [<3]\n", argc);
        exit(-1);
    }
    string baseflowname(argv[1]);
    string appdflowname(argv[2]);
    bool addall = true;
    if(argc>3) {
        string strflag(argv[3]);
        addall = (strflag=="addall");
    }

    StructuredData baseflow;
    StructuredData appdflow;
    baseflow.InputData(baseflowname, true);
    appdflow.InputData(appdflowname, true);
    std::map<int, double> field;
    for(int i=0; i<appdflow.GetNumPhys(); ++i) {
        field[i] = 0.;
    }
    appdflow.InterpolateFrom(baseflow, field);
    for(int i=0; i<appdflow.GetNumPhys(); ++i) {
        if (addall || !baseflow.HasField(appdflow.GetPhysVarName(i))) {
            baseflow.AddPhysics(appdflow.GetPhysVarName(i), appdflow.GetPhys(i));
        }
    }
    baseflow.OutputData("comb_" + baseflowname);
    return 0;
}