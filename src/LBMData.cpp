#include "LBMData.h"
#include "FileIO.h"
#include "Dataprocessing.h"
#include "StructuredData.h"
int LBMData::Extractxy() {
    std::vector<int> N = m_zones[0].N;
    std::vector<int> Nl(N.size(), 1);
    m_x.resize(N.size());
    for(size_t i=0; i<N.size(); ++i) {
        m_x[i].resize(N[i], 1.);
        if(i>0) {
            Nl[i] *= Nl[i-1] * N[i-1];
        }
    }
    for(size_t dim=0; dim<N.size(); ++dim) {
        for(int i=0; i<N[dim]; ++i) {
            m_x[dim][i] = m_zones[0].data[dim][i*Nl[dim]];
        }
    }
    return 0;
}

LBMData::LBMData(const std::string filename, const std::vector<std::vector<int>> &N){
  m_zones.resize(N.size());
  for(size_t i=0; i<N.size(); ++i) {
    m_zones[i].N = N[i];
  }
  InputTec360_FSILBM2D(filename, m_zones);
  ShiftIndex<double>(m_zones[0].N, m_zones[0].data, 1);
  Extractxy();
}

int LBMData::Interpolation(StructuredData &data, std::vector<std::vector<double>> &u1) {
    std::vector<double> o = data.GetOrigion();
    std::vector<double> dx = data.GetDetx();
    std::vector<int> N = data.GetN();
    std::vector<std::vector<double>> x1(N.size());
    for(size_t dim = 0; dim<N.size(); ++dim) {
        x1[dim].resize(N[dim]);
        for(int i=0; i<N[dim]; ++i) {
            x1[dim][i] = o[dim] + dx[dim] * i;
        }
    }
    return Interpolation(m_x, m_zones[0].data, x1, u1);
}

int LBMData::Interpolation(std::vector<std::vector<double>> &x1, std::vector<std::vector<double>> &u1) {
    return Interpolation(m_x, m_zones[0].data, x1, u1);
}

void GenerateStencil(std::vector<int> &Np, std::vector<int> &n1,
    std::vector<std::vector<int>> &index,
    std::vector<std::vector<std::vector<double>>>& weight,
    std::vector<int> &stencil, std::vector<double> &w) {
    int dim = Np.size();
    if(dim == 3 && Np[2]==1) {
        dim = 2;
    }
    if(dim==2 && Np[1]==1) {
        dim = 1;
    }
    if(dim==1) {
        stencil.push_back(index[0][n1[0]]-1);
        stencil.push_back(index[0][n1[0]]);
        w = weight[0][n1[0]];
    }else if(dim==2) {
        std::vector<int> id{index[0][n1[0]], index[1][n1[1]], 0};
        stencil.resize(4);
        w.resize(4);
        stencil[0] = Index(Np, id);
        w[0] = weight[0][n1[0]][1] * weight[1][n1[1]][1];
        stencil[1] = stencil[0] - 1;
        w[1] = weight[0][n1[0]][0] * weight[1][n1[1]][1];
        stencil[2] = stencil[0] - Np[0];
        w[2] = weight[0][n1[0]][1] * weight[1][n1[1]][0];
        stencil[3] = stencil[2] - 1;
        w[3] = weight[0][n1[0]][0] * weight[1][n1[1]][0];
    }else if(dim==3) {
        std::vector<int> id{index[0][n1[0]], index[1][n1[1]], index[2][n1[2]]};
        stencil.resize(8);
        w.resize(8);
        stencil[0] = Index(Np, id);
        w[0] = weight[0][n1[0]][1] * weight[1][n1[1]][1] * weight[2][n1[2]][1];
        stencil[1] = stencil[0] - 1;
        w[1] = weight[0][n1[0]][0] * weight[1][n1[1]][1] * weight[2][n1[2]][1];
        stencil[2] = stencil[0] - Np[0];
        w[2] = weight[0][n1[0]][1] * weight[1][n1[1]][0] * weight[2][n1[2]][1];
        stencil[3] = stencil[2] - 1;
        w[3] = weight[0][n1[0]][0] * weight[1][n1[1]][0] * weight[2][n1[2]][1];

        stencil[4] = stencil[0] - Np[0] * Np[1];
        w[4] = weight[0][n1[0]][1] * weight[1][n1[1]][1] * weight[2][n1[2]][0];
        stencil[5] = stencil[4] - 1;
        w[5] = weight[0][n1[0]][0] * weight[1][n1[1]][1] * weight[2][n1[2]][0];
        stencil[6] = stencil[4] - Np[0];
        w[6] = weight[0][n1[0]][1] * weight[1][n1[1]][0] * weight[2][n1[2]][0];
        stencil[7] = stencil[6] - 1;
        w[7] = weight[0][n1[0]][0] * weight[1][n1[1]][0] * weight[2][n1[2]][0];
    }
}

/***
 * input x = {{xmin:dx:xmax}, {ymin:dy:ymax}, ...}, only 1-D arrays
 * input u field data on the spatial grid, x goes first, then y, finally z
 * input x1, targeting grid 
 * input u1, targeting fields
***/
int LBMData::Interpolation(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &u,
    std::vector<std::vector<double>> &x1, std::vector<std::vector<double>> &u1) {
    std::vector<int> Np;
    std::vector<int> Np1;
    int Ntot1 = 1;
    for(size_t i = 0; i < x1.size(); ++i) {
        Np.push_back(x[i].size());
        Np1.push_back(x1[i].size());
        Ntot1 *= Np1[i];
    }
    u1.resize(u.size());
    for(size_t v=0; v<u1.size(); ++v) {
        u1[v].resize(Ntot1, 0.);
    }
    std::vector<std::vector<int>> index(x.size());
    std::vector<std::vector<std::vector<double>>> weight(x.size());
    for(size_t dim = 0; dim < x1.size(); ++dim) {
        Interpolation1DNonuniform::CalcWeight1D(x[dim], x1[dim], index[dim], weight[dim]);
    }
    for(int n=0; n<Ntot1; ++n) {
        std::vector<int> stencil;
        std::vector<double> w;
        std::vector<int> ind1;
        invIndex(Np1, n, ind1);
        GenerateStencil(Np, ind1, index, weight, stencil, w);
        for(int v=0; v<(int)u.size(); ++v) {
            for(size_t k=0; k<stencil.size(); ++k) {
                u1[v][n] += u[v][stencil[k]]*w[k];
            }
        }
    }
    return (int)u1.size();
}