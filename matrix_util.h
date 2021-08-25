#include"fftw3.h"
float invert_matrix(float in[3][3], float out[3][3]);
float c_invert_matrix(fftw_complex in[3][3], fftw_complex out[3][3]);
void c_matrix_times_vector(fftw_complex mat[3][3],fftw_complex vec[3],fftw_complex *result);

