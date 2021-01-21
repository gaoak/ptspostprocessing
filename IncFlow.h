#ifndef INCFLOW_H
#define INCFLOW_H
#include "StructuredData.h"
class IncFlow : public StructuredData {
public:
    IncFlow(const std::vector<int> &N, const std::vector<double> &range);
    int ExtractCore(int f, double sigma, std::vector<std::vector<double> > & cores, int dir = 2);
    int ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
        const std::vector<int> &padding, std::vector<double> &data,
        std::vector<double> &core, bool ismin = true);
};
#endif