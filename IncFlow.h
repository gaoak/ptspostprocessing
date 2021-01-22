#ifndef INCFLOW_H
#define INCFLOW_H
#include "StructuredData.h"
#include "Body.h"
class IncFlow : public StructuredData {
public:
    IncFlow(const std::vector<int> &N, const std::vector<double> &range);
    void SetBody(std::string bodyname, std::vector<double> param);
    int ExtractCore(int f, double sigma, std::vector<std::vector<double> > & cores, int dir = 2);
    int ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
        const std::vector<int> &padding, std::vector<double> &data,
        std::vector<double> &core, bool ismin = true);
    int CalculateVorticity(int order = 2);
    std::pair<int, int> GetProceedDirection(double sigma, std::vector<std::vector<double> > & cores);
protected:
    Body m_body;
};
#endif