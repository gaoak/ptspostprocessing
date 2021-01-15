#include "TECIO.h"
#include "MASTER.h"
#include<string>
#include<vector>
#include<cmath>
#include "StructuredData.h"

StructuredData::StructuredData(int i, int j, int k) {
    m_i = i;
    m_j = j;
    m_k = k;
    m_Np = m_i * m_j * m_k;
    if(m_Np==0) {
        printf("error: 0 points\n");
    }
    for(int i=0; i<3; ++i) {
        m_x.push_back(std::vector<double>(m_Np, 0.));
    }
    m_varList = "x, y, z";
}

int StructuredData::OutputCSV(std::string filename) {
    std::ofstream file(filename.c_str());
    file << "# " << m_varList << "\n";
    file << std::scientific << std::setprecision(16);
    for(int index=0; index<m_Np; ++index) {
        file << m_x[0][index] << "," << m_x[1][index] << "," << m_x[2][index];
        for(int v=0; v<m_phys.size(); ++v) {
            file << "," << m_phys[v][index];
        }
        file << "\n";
    }
    return m_Np + 1;
}

int StructuredData::ParserCSVHeader(const char * header) {
    int i = 0;
    for(; header[i]=='#' || header[i]==' '; ++i);
    m_varList = header + i;
    int vcount = -2;
    for(;header[i];++i) {
        vcount += header[i]==',';
    }
    return vcount;
}

int StructuredData::LoadCSV(std::string filename) {
    std::ifstream file(filename.c_str());
    char buffer[1000];
    std::vector<double> value;
    file.getline(buffer, sizeof(buffer));
    int vcount = ParserCSVHeader(buffer);
    m_phys.clear();
    for(int i=0; i<vcount; ++i) {
        m_phys.push_back(std::vector<double>(m_Np, 0.));
    }

    file.getline(buffer, sizeof(buffer));
    int index = 0;
    while(!file.eof()) {
        parserDouble(buffer, value);
        for(int i=0; i<vcount; ++i) {
            m_phys[i][index] = value[3+i];
        }
        ++index;
        file.getline(buffer, sizeof(buffer));
    }
    return index;
}

int StructuredData::GenPoints(double x0, double x1, double y0, double y1, double z0, double z1) {
    for(int k=0; k<m_k; ++k) {
        for(int j=0; j<m_j; ++j) {
            for(int i=0; i<m_i; ++i) {
                int index = ArrayIndex(i, j, k);
                if(m_i>1) {
                    m_x[0][index] = x0 + i*(x1 - x0)/(m_i-1);
                } else {
                    m_x[0][index] = x0;
                }
                if(m_j>1) {
                    m_x[1][index] = y0 + j*(y1 - y0)/(m_j-1);
                } else {
                    m_x[1][index] = y0;
                }
                if(m_k>1) {
                    m_x[2][index] = z0 + k*(z1 - z0)/(m_k-1);
                } else {
                    m_x[2][index] = z0;
                }
            }
        }
    }
    return m_Np;
}

int StructuredData::ArrayIndex(int i, int j, int k) {
    if(i<0 || i>=m_i || j<0 || j>=m_j || k<0 || k>=m_k) {
        printf("error: ilegal index (%d,%d,%d) for array size (%d,%d,%d)\n", i, j, k, m_i, m_j, m_k);
    }
    return i + j*m_i + k*m_i*m_j;
}

void parserDouble(const char * cstr, std::vector<double> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if((cstr[i]>='0' && cstr[i]<='9') ||
            cstr[i]=='.' ||
            cstr[i]=='e' || cstr[i]=='E' ||
            cstr[i]=='+' || cstr[i]=='-') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    double k;
    for(int i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);
        if(sscanf(cuts.c_str(), "%lf", &k)<1) {
            printf("error: parser double %s\n", cuts.c_str());
        }
        value.push_back(k);
    }
}

int OutputTec360(std::string filename, std::string variables,
                 int i, int j, int k, std::vector<void*> data,
                 int isdouble,
                 int debug,
                 int filetype,
                 int fileformat)
{
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
                (char*)variables.c_str(),  // variables list
                (char*)filename.c_str(), // output filename
                (char*)".", // Scratch Directory
                &FileFormat,
                &FileType,
                &Debug,
                &VIsDouble);
    if(I) return -1;

    /*Ordered Zone Parameters*/
    INTEGER4 IMax = i;
    INTEGER4 JMax = j;
    INTEGER4 KMax = k;
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

int StructuredData::OutputTec360(std::string filename)
{
    int isdouble = 1;
    std::vector<void*> data;
    for(int i=0; i<m_x.size(); ++i) {
        data.push_back((void*) (m_x[i].data()) );
    }
    for(int i=0; i<m_phys.size(); ++i) {
        data.push_back((void*) (m_phys[i].data()) );
    }
    ::OutputTec360(filename, m_varList, m_i, m_j, m_k, data, isdouble);
    return m_x.size() + m_phys.size();
}