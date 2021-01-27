#include<iostream>
#include<iomanip>
#include "FileIO.h"
#include "Util.h"

int OutputCSV(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data) {
    std::ofstream file(filename.c_str());
    file << "# ";
    for(int i=0; i<variables.size(); ++i) {
        file << variables[i];
        if(i!=variables.size()-1) {
            file << ",";
        } else {
            file << "\n";
        }
    }
    file << std::scientific << std::setprecision(16);
    int Np = 1;
    for(int i=0; i<N.size(); ++i) {
        Np *= N[i];
    }
    for(int index=0; index<Np; ++index) {
        for(int v=0; v<data.size(); ++v) {
            file << data[v][index];
            if(v!=data.size()-1) {
                file << ",";
            } else {
                file << "\n";
            }
        }
    }
    file.close();
    return Np;
}

int ParserCSVHeader(const char * header, std::vector<std::string> &vars) {
    int i = 0;
    for(; header[i]=='#' || header[i]==' '; ++i);
    std::string varList = header + i;
    parserString(varList.c_str(), vars, ',');
    return vars.size();
}

int InputCSV(const std::string filename, std::vector<std::string> &variables,
                 std::vector<int> &N, std::vector<std::vector<double> > &data,
                 int &isdouble) {
    std::ifstream file(filename.c_str());
    if(!file.is_open()) {
        printf("error: unable to open file %s\n", filename.c_str());
    }
    int Np = 1;
    for(int i=0; i<N.size(); ++i) {
        Np *= N[i];
    }
    char buffer[1000];
    std::vector<double> value;
    file.getline(buffer, sizeof(buffer));
    variables.clear();
    int vcount = ParserCSVHeader(buffer, variables);
    data.clear();
    for(int i=0; i<vcount; ++i) {
        data.push_back(std::vector<double>(Np, 0.));
    }

    file.getline(buffer, sizeof(buffer));
    int index = 0;
    while(!file.eof()) {
        parserDouble(buffer, value);
        for(int i=0; i<vcount; ++i) {
            data[i][index] = value[i];
        }
        ++index;
        file.getline(buffer, sizeof(buffer));
    }
    file.close();
    return index;
}

int OutputTec360_ascii(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data,
                 int isdouble)
{
    std::string varlist = variables[0];
    for(int i=1; i<variables.size(); ++i) {
        varlist += "," + variables[i];
    }
    std::ofstream ofile(filename.c_str());
    ofile << "variables = " << varlist << std::endl;
    ofile << "zone I = " << N[0] << ", J = " << N[1] << ", K = " << N[2] << std::endl;
    int Np = N[0] * N[1] * N[2];
    for(int i=0; i<Np; ++i) {
        for(int j=0; j<data.size(); ++j) {
            ofile << data[j][i] << " ";
        }
        ofile << "\n";
    }
    ofile.close();
    return 0;
}

int BinaryWrite(std::ofstream &ofile, std::string str) {
    int tmp = 0;
    for(int i=0; i<str.size(); ++i) {
        tmp = str[i];
        ofile.write((char*)&tmp, 4);
    }
    tmp = 0;
    ofile.write((char*)&tmp, 4);
}

int OutputTec360_binary(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data,
                 int isdouble)
{
    std::ofstream odata;
    odata.open(filename, std::ios::binary);
	if(!odata.is_open())
	{
        printf("error unable to open file %s\n", filename.c_str());
		return -1;
	}
    char tecplotversion[] = "#!TDV112";
	odata.write((char*)tecplotversion, 8);
    int value1 = 1;
	odata.write((char*)&value1, 4);
    int filetype = 0;
	odata.write((char*)&filetype, 4);
	//read file title and variable names
	int tempi = 0;
    std::string filetitle = "";
    BinaryWrite(odata, filetitle);
    int nvar = variables.size();
    odata.write((char*)&nvar, 4);//number of variables
    std::vector<std::string> vartitle;
	for(int i=0; i<nvar; ++i)
	{
        BinaryWrite(odata, variables[i]);
	}
    float marker299I = 299.0f;
	odata.write((char*)&marker299I, 4);
	//zone title
    std::string zonetitle("ZONE 0");
	BinaryWrite(odata, zonetitle);
    int parentzone = -1;
	odata.write((char*)&parentzone, 4);
    int strandid = -1;
	odata.write((char*)&strandid, 4);
    double soltime = 0.0;
    odata.write((char*)&soltime, 8);
    int unused = -1;
    odata.write((char*)&unused, 4);
    int zonetype = 0;
    odata.write((char*)&zonetype, 4);
    int zero = 0;
    odata.write((char*)&zero, 4);
    odata.write((char*)&zero, 4);
    odata.write((char*)&zero, 4);
	for(int i=0; i<3; ++i) {
        int tmp = N[i];
        odata.write((char*)&tmp, 4);
    }

    odata.write((char*)&zero, 4);
    float marker357 = 357.0f;
	odata.write((char*)&marker357, 4);
    float marker299II = 299.0f;
	odata.write((char*)&marker299II, 4);
	std::vector<int> binarydatatype(nvar, 1 + (isdouble>0));
	odata.write((char*)binarydatatype.data(), 4*nvar);
    odata.write((char*)&zero, 4);
    odata.write((char*)&zero, 4);
    int minus1 = -1;
    odata.write((char*)&minus1, 4);
    
    int datanumber, datasize;
	datanumber = N[0] * N[1] * N[2];
	datasize = N[0] * N[1] * N[2] * 8;
    for(int i=0; i<nvar; ++i) {
        double minv = 0., maxv=1.;
        maxv = data[i][FindMax<double>(datanumber, data[i].data())];
        minv = data[i][FindMin<double>(datanumber, data[i].data())];
        odata.write((char*)&minv, 8);
        odata.write((char*)&maxv, 8);
    }

    std::vector<float> vardata(datanumber);
    for(int i=0; i<nvar; ++i) {
        if(isdouble) {
            odata.write((char*)data[i].data(), datasize);
        } else {
            std::vector<float> fdata(datanumber);
            for(int j=0; j<datanumber; ++j) {
                fdata[j] = data[i][j];
            }
            odata.write((char*)fdata.data(), datasize);
        }
    }
    odata.close();
    return 0;
}

int InputTec360_binary(const std::string filename, std::vector<std::string> &variables,
                 std::vector<int> &N, std::vector<std::vector<double> > &data,
                 int &isdouble) {
	std::ifstream indata;
    indata.open(filename, std::ios::binary);
	if(!indata.is_open())
	{
        printf("error unable to open file %s\n", filename.c_str());
		return -1;
	}
    char tecplotversion[8];
	indata.read((char*)tecplotversion, 8);
    int value1;
	indata.read((char*)&value1, 4);
    int filetype;
	indata.read((char*)&filetype, 4);
	//read file title and variable names
	int tempi;
    std::string filetitle;
	while(1)
	{
        indata.read((char*)&tempi, 4);
        if(tempi==0) break;
        char c = tempi;
        filetitle.push_back(c);
	}
    int nvar;
    indata.read((char*)&nvar, 4);//number of variables
    variables.clear();
	for(int i=0; i<nvar; ++i)
	{
        std::string vname;
        while(1)
        {
            indata.read((char*)&tempi, 4);
            if(tempi==0) break;
            char c = tempi;
            vname.push_back(c);
        }
        variables.push_back(vname);
	}
    float marker299I;
	indata.read((char*)&marker299I, 4);
	if(marker299I!=299.)
	{
		printf("error in reading file %s\n", filename.c_str());
		return -1;
	}
	//zone title
    std::string zonentitle;
	while(1)
	{
        indata.read((char*)&tempi, 4);
        if(tempi==0) break;
        char c = tempi;
        zonentitle.push_back(c);
	}
    int parentzone;
	indata.read((char*)&parentzone, 4);
    int strandid;
	indata.read((char*)&strandid, 4);
    double soltime;
    indata.read((char*)&soltime, 8);
    int unused;
    indata.read((char*)&unused, 4);
    int zonetype;
    indata.read((char*)&zonetype, 4);
    int zero;
    indata.read((char*)&zero, 4);
    indata.read((char*)&zero, 4);
    indata.read((char*)&zero, 4);
    N.resize(3);
	indata.read((char*)N.data(), 3*4);
    indata.read((char*)&zero, 4);
    float marker357;
	indata.read((char*)&marker357, 4);
    float marker299II;
	indata.read((char*)&marker299II, 4);
	if(marker357!=357.||marker299II!=299.)
	{
		printf("error in reading file %s\n", filename.c_str());
		return -1;
	}
	std::vector<int> binarydatatype(nvar);
	indata.read((char*)binarydatatype.data(), 4*nvar);
    isdouble = binarydatatype[0] == 2;
    indata.read((char*)&zero, 4);
    indata.read((char*)&zero, 4);
    int minus1;
    indata.read((char*)&minus1, 4);

    for(int i=0; i<nvar; ++i) {
        double minv, maxv;
        indata.read((char*)&minv, 8);
        indata.read((char*)&maxv, 8);
    }

    int datanumber, datasize;
	datanumber = N[0] * N[1] * N[2];
	datasize = N[0] * N[1] * N[2] * 4;
    for(int i=0; i<nvar; ++i) {
        if(isdouble) {
            int datasize = data.size();
            data.resize(datasize+1);
            data[datasize].resize(datanumber);
            indata.read((char*)data[datasize].data(), datasize * 2);
        } else {
            std::vector<float> vardata(datanumber);
            indata.read((char*)vardata.data(), datasize);
            int datasize = data.size();
            data.resize(datasize+1);
            data[datasize].resize(datanumber);
            for(int j=0; j<datanumber; ++j) {
                data[datasize][j] = vardata[j];
            }
        }
    }
    indata.close();
    return 0;
}