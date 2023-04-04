#include<vector>
#include<string>
#include<cmath>
#include <fftw3.h>
#include "Util.h"
#include "Dataprocessing.h"

int getRealPowerSpectral(const double *data, double *spectral,
    double *beta, double Tlen, int Nfft)
{
    // prepare wavenumber, lz = 2*pi/beta
    double temp = 2.*M_PI/Tlen;
    for(int i=0;i<Nfft;++i)
    {
        if(i<=Nfft-i)
        {
            beta[i] = temp*double(i);
        }
    }
    //fftw initiates
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    fftw_plan pf;
    pf = fftw_plan_dft_1d(Nfft, in, out, FFTW_FORWARD, FFTW_MEASURE);
    //copy data
    for(int i=0;i<Nfft;++i)
    {
        in[i][0] = data[i];
        in[i][1] = 0.;
    }
    fftw_execute(pf); 
    // power spectral, integral of energy
    temp = 1. / (double(Nfft)*double(Nfft));
    spectral[0] = Tlen * (out[0][0]*out[0][0] + out[0][1]*out[0][1]) * temp;
    for(int i=1;i<Nfft;++i)
    {
        if(i<Nfft-i)
        {
            spectral[i] = 2. * Tlen * (out[i][0]*out[i][0] + out[i][1]*out[i][1]) * temp;
        }
        else if (i==Nfft-i)
        {
            spectral[i] = 0.5 * Tlen * out[i][0]*out[i][0] * temp;
        }
    }
    return 0;
}

int doRealFFT(const double *data, double *spectral,
    double *beta,double Tlen, int Nfft)
{
    // prepare wavenumber, lz = 2*pi/beta, beta is circular frequency, i.e. omega
    double temp = 2.*M_PI/Tlen;
    beta[0] = 0.;
    beta[1] = temp * (Nfft / 2);
    for(int i=1;i<Nfft/2;++i)
    {
        int k = i * 2;
        beta[k  ] = temp*double(i); //cos mode
        beta[k+1] =-temp*double(i); //sin mode
    }
    //fftw initiates
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    fftw_plan pf;
    pf = fftw_plan_dft_1d(Nfft, in, out, FFTW_FORWARD, FFTW_MEASURE);
    //copy data
    for(int i=0;i<Nfft;++i)
    {
        in[i][0] = data[i];
        in[i][1] = 0.;
    }
    fftw_execute(pf); 
    // Fourier coefficients, nektar++ format
    // [0], [cos N/2] not 0 in Nektar++, [cos 1], [sin 1] not -1 in Nektar++, [cos 2], [sin 2], ...
    temp = 1. / double(Nfft);
    spectral[0] = out[0][0] * temp;
    spectral[1] = out[Nfft/2][0] * temp;
    for(int i=1;i<Nfft/2;++i)
    {
        int k = i * 2;
        spectral[k  ] = 2. * out[i][0] * temp;
        spectral[k+1] =-2. * out[i][1] * temp;
    }
    return 0;
}

double DotVect(const std::vector<double> &a, const std::vector<double> &b) {
    double res = 0.;
    for(size_t i=0; i<a.size(); ++i) {
        res += a[i]*b[i];
    }
    return res;
}

double NormVect(const std::vector<double> &data, int p) {
    if(p<=0) {
        return (int)data.size();
    }
    double sum = 0.;
    if(p==1) {
        for(size_t i=0; i<data.size(); ++i) {
            sum += std::fabs(data[i]);
        }
    } else if(p==2) {
        for(size_t i=0; i<data.size(); ++i) {
            sum += data[i] * data[i];
        }
    } else {
        for(size_t i=0; i<data.size(); ++i) {
            sum += std::pow(std::fabs(data[i]), p);
        }
    }
    return sum;
}

void NormalizeVect(std::vector<double> &data) {
    double norm = std::sqrt(NormVect(data, 2));
    for(size_t i=0; i<data.size(); ++i) {
        data[i] /= norm;
    }
}

void AddVect(const double a1, const std::vector<double> &a,
             const double b1, const std::vector<double> &b,
             std::vector<double> &res) {
    size_t n = std::min(a.size(), b.size());
    if(res.size()<n) {
        res.resize(n);
    }
    for(size_t i=0; i<n; ++i) {
        res[i] = a1*a[i] + b1*b[i];
    }
}

std::vector<double> AddVect(const double a1, const std::vector<double> &a,
             const double b1, const std::vector<double> &b) {
    size_t n = std::min(a.size(), b.size());
    std::vector<double> res(n);
    for(size_t i=0; i<n; ++i) {
        res[i] = a1*a[i] + b1*b[i];
    }
    return res;
}

std::vector<double> CrossVect(std::vector<double> v1, std::vector<double> v2) {
    std::vector<double> res(3, 0.);
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return res;
}

int FindAbsMax(int N, const double* data) {
    int ind=0;
    double tmp = std::fabs(data[0]);
    for(int i=1; i<N; ++i) {
        if(tmp < std::fabs(data[i])) {
            ind = i;
            tmp = std::fabs(data[i]);
        }
    }
    return ind;
}

std::vector<double> transform(const std::vector<double> &p, double AoA) {
    std::vector<double> res = p;
    res[0] = p[0]*cos(AoA) + p[1]*sin(AoA);
    res[1] =-p[0]*sin(AoA) + p[1]*cos(AoA);
    return res;
}

void parserUInt(const char * cstr, std::vector<int> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if(cstr[i]>='0' && cstr[i]<='9') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    int k;
    for(int i=0; i<(int) digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);// data  in [s e-1]
        if(sscanf(cuts.c_str(), "%d", &k)<1) {
            printf("error: parser int %s \n",  cuts.c_str());
        }
        if(i>0 && (digs[i] - dige[i-1])==1 && cstr[digs[i]-1]=='-') {
            for(int j=value[(int)value.size()-1]+1; j<k; ++j) {
                value.push_back(j);
            }
        }
        value.push_back(k);
    }
}

void parserInt(const char * cstr, std::vector<int> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if((cstr[i]>='0' && cstr[i]<='9') ||
            cstr[i]=='-' || cstr[i]=='+') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    int k;
    for(size_t i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);// data  in [s e-1]
        if(sscanf(cuts.c_str(), "%d", &k)<1) {
            printf("error: parser int %s \n",  cuts.c_str());
        }
        value.push_back(k);
    }
}

void parserString(const char * cstr, std::vector<std::string> & value, char sep) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if(cstr[i]!=sep && cstr[i]!=0) {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    for(size_t i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);
        value.push_back(cuts);
    }
}

void parserDouble(const char * cstr, std::vector<double> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if((cstr[i]>='0' && cstr[i]<='9') ||
            cstr[i]=='.' ||
            cstr[i]=='e' || cstr[i]=='E' ||
            cstr[i]=='+' || cstr[i]=='-') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    double k;
    for(size_t i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);
        if(sscanf(cuts.c_str(), "%lf", &k)<1) {
            printf("error: parser double %s, in \n%s\n", cuts.c_str(), cstr);
        }
        value.push_back(k);
    }
}

double StringToDouble(std::string str) {
    double k;
    if(sscanf(str.c_str(), "%lf", &k)<1) {
        printf("error: parser double %s\n", str.c_str());
    }
    return k;
}

bool CompactArray(std::vector<int> &N) {
    std::vector<int> sN = N;
    N.clear();
    for(size_t i=0; i<sN.size(); ++i) {
        if(1!=sN[i]) {
            N.push_back(sN[i]);
        }
    }
    return N.size()!=sN.size();
}

bool NeedShift(std::vector<int> N, int dir) {
    for(size_t i=0; i<N.size(); ++i) {
        if(N[i]!=1) {
            N[i] = 10 + i;
        }
    }
    std::vector<int> sN = N;
    ShiftArray<int>(sN, dir);
    CompactArray(N);
    CompactArray(sN);
    for(size_t i=0; i<N.size(); ++i) {
        if(N[i]!=sN[i]) {
            return true;
        }
    }
    return false;
}

int myMod(int x, int N) {
    if(N<=0 || x==0) {
        return x;
    }
    if(x<0) {
        return (N + x % N ) % N;
    } else {
        return x % N;
    }
}

int Index(const std::vector<int> &N, const std::vector<int> & index) {
    int res = index[0];
    int Np = N[0];
    for(size_t i=1; i<N.size(); ++i) {
        if(index[i]<0 || index[i]>=N[i]) {
            printf("error: illegal index %d in Index %d\n", index[i], N[i]);
        }
        res += Np*index[i];
        Np *= N[i];
    }
    return res;
}

void invIndex(const std::vector<int> &N, int index, std::vector<int> & res) {
    res.resize(N.size());
    res[0] = N[0];
    for(int i=1; i<(int)N.size()-1; ++i) {
        res[i] = res[i-1] * N[i];
    }
    for(int i=(int)res.size()-1; i>0; --i) {
        res[i] = index / res[i-1];
        index = index % res[i-1];
    }
    res[0] = index;
}

bool startWithNumber(const char *str, int len) {
    size_t i=0;
    for(; i<len; ++i) {
        if(str[i]!=' ') break;
    }
    if(str[i]=='+' || str[i]=='-' || str[i]=='.' || (str[i]>=0 && str[i]<=9)) {
        return true;
    }
    return false;
}

bool startWithNumber(std::string &str) {
    return startWithNumber(str.c_str(), str.size());
}