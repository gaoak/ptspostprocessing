#ifndef HEADERFILEUTIL_H
#define HEADERFILEUTIL_H
#include<string>
#include<vector>
//overloading functions
int myMod(int x, int N);
//parser data from string
double StringToDouble(std::string str);
void parserUInt(const char * cstr, std::vector<int> & value);
void parserInt(const char * cstr, std::vector<int> & value);
void parserDouble(const char * cstr, std::vector<double> & value);
void parserString(const char * cstr, std::vector<std::string> & value, char sep = ' ');
//high dimensional index
int Index(const std::vector<int> &N, const std::vector<int> & index);
void invIndex(const std::vector<int> &N, int index, std::vector<int> & res);
bool NeedShift(std::vector<int> N, int dir);
template<typename T>
void ShiftArray(std::vector<T> &a, int dir);
template<typename T>
int ShiftIndex(std::vector<int> &N, std::vector<std::vector<T> > &odata, int dir);
//vector operation
int FindAbsMax(int N, const double* data);
template<typename T>
int FindMax(int N, const T* data);
template<typename T>
int FindMin(int N, const T* data);
std::vector<double> transform(const std::vector<double> &p, double AoA);
double DotVect(const std::vector<double> &a, const std::vector<double> &b);
double NormVect(const std::vector<double> &data, int p);
void NormalizeVect(std::vector<double> &data);
void AddVect(const double a1, const std::vector<double> &a,
             const double b1, const std::vector<double> &b,
             std::vector<double> &res);
std::vector<double> AddVect(const double a1, const std::vector<double> &a,
             const double b1, const std::vector<double> &b);
std::vector<double> CrossVect(std::vector<double> v1, std::vector<double> v2);
//////////////////////template function definition
template<typename T>
int FindMax(int N, const T* data) {
    int ind=0;
    T tmp = data[0];
    for(int i=1; i<N; ++i) {
        if(tmp < data[i]) {
            ind = i;
            tmp = data[i];
        }
    }
    return ind;
}

template<typename T>
int FindMin(int N, const T* data) {
    int ind=0;
    T tmp = data[0];
    for(int i=1; i<N; ++i) {
        if(tmp > data[i]) {
            ind = i;
            tmp = data[i];
        }
    }
    return ind;
}

template <typename T>
int DoMask(const std::vector<int> &N, T* data, const std::vector<int> &center, const std::vector<int> &padding, const bool relation) {
    int Np = 1;
    for(int i=0; i<(int)N.size(); ++i) {
        Np *= N[i];
    }
    std::vector<int> ind;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        int k=0;
        for(; k<(int)N.size(); ++k) {
            if (std::abs(ind[k] - center[k]) > padding[k]) {
                break;
            }
        }
        if(k==(int)N.size() && relation) {
            data[i] = 0;
        }
        if(k<(int)N.size() && !relation) {
            data[i] = 0;
        }
    }
    return Np;
}

template <typename T>
int WeightedCenter(const std::vector<int> &N, T* data, std::vector<T> &center) {
    int Np = 1;
    for(int i=0; i<(int)N.size(); ++i) {
        Np *= N[i];
    }
    center = std::vector<T>(N.size(), 0);
    std::vector<int> ind;
    T sum = 0;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        sum += data[i];
        for(int k=0; k<(int)N.size(); ++k) {
            center[k] += data[i] * ind[k];
        }
    }
    for(int k=0; k<(int)N.size(); ++k) {
        center[k] /= sum;
    }
    return Np;
}

template <typename T>
int DoMask(const int N, const T threshold, const int comp, T* data) {
    for(int i=0; i<N; ++i) {
        if(comp > 0 && threshold > data[i]) {
            data[i] = 0;
        }
        if(comp == 0 && threshold == data[i]) {
            data[i] = 0;
        }
        if(comp < 0 && threshold < data[i]) {
            data[i] = 0;
        }
    }
    return N;
}

template <typename T>
int DoMaskShift(const int N, const T threshold, const int comp, T* data) {
    for(int i=0; i<N; ++i) {
        data[i] -= threshold;
        if(comp > 0 && 0 > data[i]) {
            data[i] = 0;
        }
        if(comp == 0 && 0 == data[i]) {
            data[i] = 0;
        }
        if(comp < 0 && 0 < data[i]) {
            data[i] = 0;
        }
    }
    return N;
}

template<typename T>
int ShiftIndex(std::vector<int> &N, std::vector<std::vector<T> > &odata, int dir) {
    // i,j,k to j,k,i
    if(odata.size()==0) {
        return 0;
    }
    std::vector<int> sN = N;
    ShiftArray<int>(N, dir);
    if(!NeedShift(sN, dir)) {
        return 0;
    }

    int dim = sN.size();
    dir = myMod(dir, dim);
    std::vector<std::vector<double> > data(odata.size());
    for(int i=0; i<(int)odata.size(); ++i) {
        data[i] = odata[i];
    }
    if(dim==2) {
        if(dir!=1) {
            return 0;
        }
        int Np = sN[0] * sN[1];
        for(int d=0; d<(int)odata.size(); ++d) {
            int count = 0;
            for(int i=0; i<sN[0]; ++i) {
                for(int j=i; j<Np; j+=sN[0]) {
                    odata[d][count++] = data[d][j];
                }
            }
        }
        return Np;
    } else if(dim==3) {
        int Np = sN[0] * sN[1] * sN[2];
        int N01 = sN[0] * sN[1];
        if(dir==2) {
            for(int d=0; d<(int)odata.size(); ++d) {
                int count = 0;
                for(int i=0; i<sN[0]; ++i) {
                    for(int k=i; k<Np; k+=N01) {
                        int jmax = k + N01;
                        for(int j= k; j<jmax; j+=sN[0]) {
                            odata[d][count++] = data[d][j];
                        }
                    }
                }
            }
        } else if(dir==1) {
            for(int d=0; d<(int)odata.size(); ++d) {
                int count = 0;
                for(int j=0; j<N01; j+=sN[0]) {
                    for(int i=0; i<sN[0]; ++i) {
                        for(int k=i+j; k<Np; k+=N01) {
                            odata[d][count++] = data[d][k];
                        }
                    }
                }
            }
        }
        return Np;
    } else {
        printf("error: unsupported space dimension %d\n", dim);
        return 0;
    }
}

template<typename T>
void ShiftArray(std::vector<T> &a, int dir) {
    int N = a.size();
    if(N==0) {
        return;
    }
    dir = myMod(dir, N);
    std::vector<T> tmp = a;
    for(int i=0; i<N; ++i) {
        a[(i+dir)%N] = tmp[i];
    }
}

template<typename T>
int myRound(T x) {
    if(x>0) {
        return (int) (x + 0.5);
    } else {
        return -(int) (-x + 0.5);
    }
}

int getRealPowerSpectral(const double *data, double *spectral,
    double *beta,double Tlen, int Nfft);

int doRealFFT(const double *data, double *spectral,
    double *beta,double Tlen, int Nfft);

bool startWithNumber(const std::string &str);

bool startWithNumber(const char *str, size_t len);
#endif