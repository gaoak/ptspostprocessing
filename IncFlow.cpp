#include<tuple>
#include<set>
#include "IncFlow.h"
#include "Util.h"

IncFlow::IncFlow(const std::vector<int> &N, const std::vector<double> &range) 
        :StructuredData(N, range) {
}

int IncFlow::ExtractCore(int f, double sigma, std::vector<std::vector<double> > & cores, int dir) {
    cores.clear();
    std::vector<std::vector<double> > odata;
    odata.push_back(m_phys[f]);
    Smoothing(sigma, odata);

    std::vector<int> N = m_N;
    std::vector<double> dx = m_dx;
    std::vector<double> range = m_range;
    std::vector<int> padding(3);
    for(int i=0; i<3; ++i) {
        if(N[i]>1) {
            padding[i] = std::round(3.*sigma/dx[i]);
        } else {
            padding[i] = 0;
        }
    }

    int Trymax = m_N[0] + m_N[1] + m_N[2];
    std::set<std::pair<int, int> > searched;
    std::vector<int> intcenter;
    std::pair<int, int> plane(dir, 0);
    std::vector<int> planeN;
    std::vector<double> planedata;
    int count = 0;
    while(count<Trymax) {
        if(plane.second<0 || plane.second>=N[plane.first]) {
            break;
        }
        dir = plane.first;
        printf("search plane %d, %d\n", dir, plane.second);
        ExtractPlane(odata[0], plane, planeN, planedata);
        ShiftArray<double>(dx, 2-dir);
        ShiftArray<double>(range, 2-dir);
        ShiftArray<int>(N, 2-dir);
        ShiftArray<int>(padding, 2-dir);
        ShiftArray<int>(intcenter, 2-dir);
        std::vector<double> newcenter;
        std::vector<int> Nslice = {N[0], N[1]};
        ExtractCore2Dplane(Nslice, intcenter, padding, planedata, newcenter);
        std::vector<double> physcenter = newcenter;
        intcenter.resize(2);
        for(int k=0; k<2; ++k) {
            physcenter[k] = range[2*k] + dx[k] * physcenter[k];
            intcenter[k] = std::round(newcenter[k]);
        }
        physcenter.push_back(range[4] + (plane.second)*dx[2]);
        intcenter.push_back(plane.second);
        ShiftArray<double>(dx, dir-2);
        ShiftArray<double>(physcenter, dir-2);
        ShiftArray<double>(range, dir-2);
        ShiftArray<int>(N, dir-2);
        ShiftArray<int>(padding, dir-2);
        ShiftArray<int>(intcenter, dir-2);
        ShiftArray<int>(padding, dir-2);
        cores.push_back(physcenter);
        searched.insert(plane);
        if(count) {
            int num = cores.size() - 1;
            double tmp = std::fabs(cores[num][0] - cores[num-1][0]);
            int itmp = 0;
            if(tmp < std::fabs(cores[num][1] - cores[num-1][1])) {
                itmp = 1;
                tmp = std::fabs(cores[num][1] - cores[num-1][1]);
            }
            if(tmp < std::fabs(cores[num][2] - cores[num-1][2])) {
                itmp = 2;
            }
            plane.first = itmp;
            int sign = 1;
            if(cores[num][itmp] - cores[num-1][itmp] < 0) sign = -1;
            plane.second = intcenter[itmp] + sign;

        } else {
            plane.second += 1;
        }
        if(searched.find(plane)!=searched.end()) {
            break;
        }
        ++count;
    }
    return count;
}

int IncFlow::ExtractCore2Dplane(const std::vector<int> &N, const std::vector<int> &initial,
    const std::vector<int> &padding, std::vector<double> &data, std::vector<double> &core, bool ismin) {
    if(!ismin) {
        for(int i=0; i<N[0]*N[1]; ++i) {
            data[i] = -data[i];
        }
    }
    if(initial.size()>=2) {
        DoMask<double>(N, data.data(), initial, padding, false);
    }
    int Np = N[0] * N[1];
    int imin = FindMin<double>(Np, data.data());
    double pmin = data[imin];
    double thresh = 0.9 * pmin;
    std::vector<double> center;
    DoMaskShift<double>(Np, thresh, -1, data.data());
    core.clear();
    WeightedCenter<double>(N, data.data(), core);
    return core.size();
}