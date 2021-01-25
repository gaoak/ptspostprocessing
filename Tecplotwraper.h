#ifndef TECPLOTWRAPER_H
#define TECPLOTWRAPER_H
#include<string>
#include<vector>
#ifdef __linux__
#define TECPLOTEXT ".plt"
#else
#define TECPLOTEXT ".dat"
#endif
int OutputTec360(const std::string filename, const std::vector<std::string> &variables,
                 const std::vector<int> &N, const std::vector<void*> data,
                 int isdouble,
                 int debug = 0,
                 int filetype = 0,
                 int fileformat = 0);
#endif
