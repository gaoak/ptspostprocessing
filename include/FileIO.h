#ifndef FILEIO_H
#define FILEIO_H
#include<string>
#include<vector>
#include<fstream>
#include<map>
#ifdef TEC360USEDOUBLE
#define TEC360USEDOUBLE 1
#else
#define TEC360USEDOUBLE 0
#endif

struct DataPack;

struct DataPack {
    std::vector<int> N;
    std::vector<std::vector<double> > data;
};

int BinaryWrite(std::ofstream &ofile, std::string str);

int OutputTec360_ascii(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data,
                 int isdouble = 1);

int OutputTec360_binary(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data,
                 int isdouble);

int OutputCSV(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data);

int ParserCSVHeader(const char * header, std::vector<std::string> &vars);

int InputTec360_FSILBM2D(const std::string filename, std::vector<DataPack> &zones);

int InputTec360_binary(const std::string filename, std::vector<std::string> &variables,
                 std::vector<int> &N, std::vector<std::vector<double> > &data,
                 int &isdouble, std::map<int, int> &vm);

int InputCSV(const std::string filename, std::vector<std::string> &variables,
                 std::vector<int> &N, std::vector<std::vector<double> > &data,
                 int &isdouble, std::map<int, int> &vm);
#endif
