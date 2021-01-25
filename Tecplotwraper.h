#ifndef TECPLOTWRAPER_H
#define TECPLOTWRAPER_H
#include<string>
#include<vector>

int OutputTec360_calllib(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<void*> data,
                 int isdouble,
                 int debug = 0,
                 int filetype = 0,
                 int fileformat = 0);
#endif
