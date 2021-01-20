#ifndef UTIL_H
#define UTIL_H
#include<vector>
void parserUInt(const char * cstr, std::vector<int> & value);
void parserDouble(const char * cstr, std::vector<double> & value);
int Index(const std::vector<int> &N, const std::vector<int> & index);
void invIndex(const std::vector<int> &N, int index, std::vector<int> & res);
void ShiftArray(std::vector<int> &a, int dir);
int ShiftIndex(std::vector<int> &N, std::vector<std::vector<double> > &odata, int dir);
#endif