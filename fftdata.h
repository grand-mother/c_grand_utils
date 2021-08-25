#include<math.h>

void fft_init(int len);
void fft_forward(float *input, float *output);
void fft_backward(float *input, float *output);
void mag_and_phase(float *in, float *out_mag,float * out_phase);
void trace_from_mag_phase(float *in_mag,float * in_phase,float *out_trace);
void envelope(float *in,float *out);
void filter_data(float *in, float *out, float freqsample, float freqmin, float freqmax);
