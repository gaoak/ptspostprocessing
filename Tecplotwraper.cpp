#include<TECIO.h>
#include<MASTER.h>
#include<string>
#include<vector>
#include "Tecplotwraper.h"
int OutputTec360_calllib(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<void*> data,
                 int isdouble,
                 int debug,
                 int filetype,
                 int fileformat)
{
    std::string varlist = variables[0];
    for(int i=1; i<variables.size(); ++i) {
        varlist += "," + variables[i];
    }
    INTEGER4 Debug = debug;
    INTEGER4 VIsDouble = isdouble;
    INTEGER4 FileType = filetype;
    INTEGER4 FileFormat = fileformat; // 0 == PLT, 1 == SZPLT
    INTEGER4 I = 0; /* Used to track return codes */
    /*
    * Open the file and write the tecplot datafile
    * header information
    */
    I = TECINI142((char*)"OutputTec360", // dataset title, seems not important
                (char*)varlist.c_str(),  // variables list
                (char*)filename.c_str(), // output filename
                (char*)".", // Scratch Directory
                &FileFormat,
                &FileType,
                &Debug,
                &VIsDouble);
    if(I) return -1;

    /*Ordered Zone Parameters*/
    INTEGER4 IMax = N[0];
    INTEGER4 JMax = N[1];
    INTEGER4 KMax = N[2];
    INTEGER4 ZoneType = 0;
    INTEGER4 ICellMax = 0;
    INTEGER4 JCellMax = 0;
    INTEGER4 KCellMax = 0;
    double   SolTime  = 0.;
    INTEGER4 StrandID = 0;
    INTEGER4 ParentZn = 0;
    INTEGER4 IsBlock = 1;
    INTEGER4 NFConns = 0;
    INTEGER4 FNMode = 0;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 ShrConn = 0;
    I = TECZNE142((char*)"Ordered Zone", // zone name, seems not important
                &ZoneType, // 0 is ordered zone
                &IMax,
                &JMax,
                &KMax,
                &ICellMax,
                &JCellMax,
                &KCellMax,
                &SolTime,
                &StrandID,
                &ParentZn,
                &IsBlock,
                &NFConns,
                &FNMode,
                &TotalNumFaceNodes,
                &TotalNumBndryFaces,
                &TotalNumBndryConnections,
                NULL,
                NULL,
                NULL,
                &ShrConn);
    if(I) return -1;

    INTEGER4 III = IMax * JMax * KMax;
    for(int i=0; i<data.size(); ++i) {
        I = TECDAT142(&III, data[i], &VIsDouble);
        if(I) return -1;
    }
    I = TECEND142();
    if(I) return -1;
    
    return 0;
}