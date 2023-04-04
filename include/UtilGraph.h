#ifndef HEADERFILEUTILGRAPH_H
#define HEADERFILEUTILGRAPH_H
#include<vector>
int Fill2DGraph(const std::vector<int> &rawN, std::vector<double> &value, const std::vector<int> &init,
                const double &eps, bool monotone);
int FindLocMaxIn2DGraph(const std::vector<int> &N, const std::vector<int> &initial,
    std::vector<double> &data, std::vector<int> &core, bool ismax);
int PurgeDifferentSign(const std::vector<int> &N, const std::vector<double> &v,
                       std::vector<double> &data, double sign);
int ExtractPatchStat2DGraph(const std::vector<int> &N, const std::vector<double> &dx, std::vector<int> core,
    std::vector<double> &v, std::vector<double> &mu, double &sum);
int FindAllLocMaxIn2DGraph(const std::vector<int> &N,
    std::vector<double> &data, std::vector<std::vector<int> > &cores, double threshold, bool ismax);
#endif