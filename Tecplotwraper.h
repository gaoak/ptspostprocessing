#ifndef TECPLOTWRAPER_H
#define TECPLOTWRAPER_H
int OutputTec360(std::string filename, std::string variables,
                 int i, int j, int k, std::vector<void*> data,
                 int isdouble = 1,
                 int debug = 1,
                 int filetype = 0,
                 int fileformat = 0);
#endif