#ifndef UTIL_H
#define UTIL_H
#include<vector>
void parserDouble(const char * cstr, std::vector<double> & value);
int Index(const std::vector<int> &N, const std::vector<int> & index);
void invIndex(const std::vector<int> &N, int index, std::vector<int> & res);
void LeftShiftArray(const std::vector<int> &a, std::vector<int> &res);
int LeftShiftIndex(std::vector<int> &N, std::vector<int> &oN, const double *data, double * odata);
#endif