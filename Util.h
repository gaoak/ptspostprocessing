#ifndef UTIL_H
#define UTIL_H
#include<string>
#include<vector>
double StringToDouble(std::string str);
void parserUInt(const char * cstr, std::vector<int> & value);
void parserDouble(const char * cstr, std::vector<double> & value);
void parserString(const char * cstr, std::vector<std::string> & value, char sep = ' ');
int Index(const std::vector<int> &N, const std::vector<int> & index);
void invIndex(const std::vector<int> &N, int index, std::vector<int> & res);
bool NeedShift(std::vector<int> N, int dir);
template<typename T>
int FindMax(int N, const T* data);
template<typename T>
int FindMin(int N, const T* data);
template<typename T>
void ShiftArray(std::vector<T> &a, int dir);
template<typename T>
int ShiftIndex(std::vector<int> &N, std::vector<std::vector<T> > &odata, int dir);
double Distance(const std::vector<double> &a, const std::vector<double> &b);
template<typename T>
std::vector<T> AddVect(const T a1, const std::vector<T> &a, const T b1, const std::vector<T> &b);
std::vector<double> transform(const std::vector<double> &p, double AoA);

template<typename T>
std::vector<T> AddVect(const T a1, const std::vector<T> &a, const T b1, const std::vector<T> &b) {
    std::vector<T> res(a.size());
    for(int i=0; i<a.size(); ++i) {
        res[i] = a1*a[i] + b1*b[i];
    }
    return res;
}

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

int FindAbsMax(int N, const double* data);

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

template<typename T>
std::vector<T> crossproduct(std::vector<T> v1, std::vector<T> v2)
{
    std::vector<T> res(3, 0.);
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return res;
}

template <typename T>
int DoMask(const std::vector<int> &N, T* data, const std::vector<int> &center, const std::vector<int> &padding, const bool relation) {
    int Np = 1;
    for(int i=0; i<N.size(); ++i) {
        Np *= N[i];
    }
    std::vector<int> ind;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        int k=0;
        for(; k<N.size(); ++k) {
            if (std::abs(ind[k] - center[k]) > padding[k]) {
                break;
            }
        }
        if(k==N.size() && relation) {
            data[i] = 0;
        }
        if(k<N.size() && !relation) {
            data[i] = 0;
        }
    }
    return Np;
}

template <typename T>
int WeightedCenter(const std::vector<int> &N, T* data, std::vector<T> &center) {
    int Np = 1;
    for(int i=0; i<N.size(); ++i) {
        Np *= N[i];
    }
    center = std::vector<T>(N.size(), 0);
    std::vector<int> ind;
    T sum = 0;
    for(int i=0; i<Np; ++i) {
        invIndex(N, i, ind);
        sum += data[i];
        for(int k=0; k<N.size(); ++k) {
            center[k] += data[i] * ind[k];
        }
    }
    for(int k=0; k<N.size(); ++k) {
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

    size_t dim = sN.size();
    int Np = 1;
    for(int i=0; i<dim; ++i) {
        Np *= sN[i];
    }
    std::vector<std::vector<double> > data;
    for(int i=0; i<odata.size(); ++i) {
        data.push_back(odata[i]);
    }
    
    std::vector<int> ind(dim);
    std::vector<int> indo(dim);
    for(int i=0; i<Np; ++i) {
        invIndex(sN, i, ind);
        indo = ind;
        ShiftArray<int>(indo, dir);
        int it = Index(N, indo);
        for(int n=0; n<data.size(); ++n) {
            odata[n][it] = data[n][i];
        }
    }
    return Np;
}

template<typename T>
void ShiftArray(std::vector<T> &a, int dir) {
    int N = a.size();
    if(N==0) {
        return;
    }
    if(dir<0) {
        dir = (N + dir % N ) % N;
    } else {
        dir = dir % N;
    }
    std::vector<T> tmp = a;
    for(int i=0; i<N; ++i) {
        a[(i+dir)%N] = tmp[i];
    }
}

template<typename T>
int myRound(T x) {
    return (int) (x + 0.5);
}
#endif