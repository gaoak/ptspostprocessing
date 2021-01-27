#ifndef FILEIO_H
#define FILEIO_H
#include<string>
#include<vector>
#include<fstream>

int BinaryWrite(std::ofstream &ofile, std::string str);

int OutputTec360_ascii(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data,
                 int isdouble = 1);

int OutputTec360_binary(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data,
                 int isdouble);

int InputTec360_binary(const std::string filename, std::vector<std::string> &variables,
                 std::vector<int> &N, std::vector<std::vector<double> > &data,
                 int &isdouble);

int OutputCSV(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<std::vector<double> > &data);

int ParserCSVHeader(const char * header, std::vector<std::string> &vars);

int InputCSV(const std::string filename, std::vector<std::string> &variables,
                 std::vector<int> &N, std::vector<std::vector<double> > &data,
                 int &isdouble);
#endif
