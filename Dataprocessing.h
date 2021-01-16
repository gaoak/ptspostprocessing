#ifndef DATAPROCESS_H
#define DATAPROCESS_H
#include <fftw3.h>
#define PROC_SUCCESS 0
#define PROC_CANNOT_OPENFILE -1001
#define PROC_BAD_MEMORY -1002
#define PROC_MIS_CALLED -1003
#define PROC_SINGULAR -1004
#ifndef STATUS
	typedef int STATUS; 
#endif
enum FilterType
{
	eroundoff,
	ebox,
	eGaussian,
	esharp,
	eCauchy,
	ePao,
	eNone
};
STATUS loadfile(const char *filename, int &num, double * & time, double * & data);
STATUS loadforcemomentumfile(const char * filename, int &num, double *&time, double *force_data[3], double *momentum_data[3]);
STATUS findacmes(const double *data, const int num, int & numacmes, int * & acme_index);
STATUS parser_arguments(int argc, char * argv[], char * filename, double & nperiods);
STATUS find_maxspectral(fftw_complex * spectral, int Nfft, int firstn, int * & index);
STATUS filter(fftw_complex * spectral, double * omega, int Nfft, FilterType type);
STATUS derivatives12(double * data, double * d_data, double * dd_data, double Tlen, int Nfft, FilterType ftype = eGaussian);


class GaussKernel {
public:
    GaussKernel(double sigma, double dx, double cutoff = 5.);
    double GW(int offset);
    int GetIMax();
    std::map<int, double> m_weight;
    int m_maxOffset;
};

#endif