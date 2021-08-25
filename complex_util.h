#include "fftw3.h"

fftw_complex* c_multiply(fftw_complex a,fftw_complex b,fftw_complex *res);
fftw_complex* c_divide(fftw_complex a,fftw_complex b,fftw_complex *res);
fftw_complex* c_add(fftw_complex a,fftw_complex b,fftw_complex *res);
fftw_complex* c_sub(fftw_complex a,fftw_complex b,fftw_complex *res);
double c_phase(fftw_complex a);
double c_mag(fftw_complex a);
fftw_complex* c_pow(fftw_complex a, float pw,fftw_complex *res);

