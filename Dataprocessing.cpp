#include "Dataprocessing.h"
#include<cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

int filter_Gaussian(fftw_complex * spectral, double * omega, int Nfft, double sigma)
{
	//Gaussian filter
	double Delta, temp, ftmp;
	temp = sigma*sigma/2.;
	for(int i=1;i<Nfft;i++)
	{
		ftmp = exp(-temp*omega[i]*omega[i]);
		spectral[i][0] *= ftmp;
		spectral[i][1] *= ftmp;
	}
	return 0;
}

STATUS derivatives1D(double * data, double * d_data, double Tlen, int Nfft, FilterType ftype, double sigma)
{
	//using FFT to get derivatives
	//fftw initiates
    fftw_complex *in, *out, *tilt_data;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    tilt_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
    if(NULL==in||NULL==out||NULL==tilt_data)
    {
    	return PROC_BAD_MEMORY;
    }
    fftw_plan pf, pb;
    pf = fftw_plan_dft_1d(Nfft, in, out, FFTW_FORWARD, FFTW_MEASURE);
	pb = fftw_plan_dft_1d(Nfft, in, out, FFTW_BACKWARD, FFTW_MEASURE);
	if(Tlen<1.E-16) return PROC_SINGULAR;
    double temp = 2.*M_PI/Tlen;
    int i;
    double *omega;
    try
	{
		omega = new double[Nfft];
	}
	catch(const std::bad_alloc &ex)
	{
		return PROC_BAD_MEMORY;
	}
    for(i=0;i<Nfft;i++)
    {
       	if(i<=Nfft-i)
    	{
    		omega[i] = temp*double(i);
    	}
    	else
    	{
    		omega[i] = temp*double(i-Nfft);
    	}
    }
    //execute FFT
    //get spectral of data
    for(i=0;i<Nfft;i++)
    {
    	in[i][0] = data[i];
    	in[i][1] = 0.;
    }
    fftw_execute(pf); 
    STATUS return_status;
    return_status = filter(out, omega, Nfft, ftype);
    if(PROC_SUCCESS != return_status)//filter
    {
    	return return_status;
    }
    for(i=0;i<Nfft;i++)
    {
    	tilt_data[i][0] = out[i][0]/double(Nfft);
    	tilt_data[i][1] = out[i][1]/double(Nfft);
    }
    //get derivative of data
    for(i=0;i<Nfft;i++)
    {
    	in[i][0] = -omega[i]*tilt_data[i][1];
    	in[i][1] =  omega[i]*tilt_data[i][0];
    }
    fftw_execute(pb); 
    for(i=0;i<Nfft;i++)
    {
    	d_data[i] = out[i][0];
    }
    //get double derivative of data
    for(i=0;i<Nfft;i++)
    {
    	temp = -omega[i]*omega[i];
    	in[i][0] = temp*tilt_data[i][0];
    	in[i][1] = temp*tilt_data[i][1];
    }
    fftw_execute(pb); 
    for(i=0;i<Nfft;i++)
    {
    	dd_data[i] = out[i][0];
    }
    //free allocated spaces
	fftw_destroy_plan(pf);
	fftw_destroy_plan(pb);
    fftw_free(in); fftw_free(out);fftw_free(tilt_data);
	return PROC_SUCCESS;
}


GaussKernel::GaussKernel(double sigma, double dx, double cutoff) {
    m_maxOffset = std::floor(cutoff*sigma/dx);
    double sum = 0.;
    for(int i=-m_maxOffset; i<=m_maxOffset; ++i) {
        double tmp = i*dx/sigma;
        tmp = std::exp(-0.5*tmp*tmp);
        m_weight[i] = tmp;
        sum += tmp;
    }
    for(int i=-m_maxOffset; i<=m_maxOffset; ++i) {
        m_weight[i] /= sum;
    }
}

double GaussKernel::GW(int offset) {
    if(m_weight.count(offset)) {
        return m_weight[offset];
    } else {
        return 0.;
    }
}

int GaussKernel::GetIMax() {
    return m_maxOffset;
}

int main() {
    GaussKernel ker(1., 0.3);
    double sum = 0.;
    for(int i=-ker.GetIMax(); i<=ker.GetIMax(); ++i) {
        double tmp = ker.GW(i);
        printf("%g,",tmp);
        sum+= tmp;
    }
    printf("sum %g\n", sum-1.);
}