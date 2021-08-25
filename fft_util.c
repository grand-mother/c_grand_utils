/*! \file fft_util.c
This file contains routines that aid implementation of fftw functionality in GRAND
*/
#include<string.h>
#include "math.h"
#include "fftw3.h"

fftw_plan fftpf,fftpb;

fftw_complex *fftin=NULL,*fftout=NULL;

int fft_len=0;


int iswap=0;
short *datbuf=NULL;

/**
 * Initialize the fft. If fftw was already initialized, first memory is released.
 * @param[in] len the length of the trace(s) to be Fourier transformed
 */
void fft_init(int len)
{
  if(fftin != NULL) fftw_free(fftin);
  if(fftout != NULL) fftw_free(fftout);
  if(fftpf != NULL) fftw_destroy_plan(fftpf);
  if(fftpb != NULL) fftw_destroy_plan(fftpb);
  fftin=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*len);
  fftout=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*len);
  fftpf=fftw_plan_dft_1d(len,fftin,fftout,FFTW_FORWARD,FFTW_MEASURE);
  fftpb=fftw_plan_dft_1d(len,fftout,fftin,FFTW_BACKWARD,FFTW_MEASURE);
  fft_len = len;
}

/**
 * perform a forward FFT (time to frequency)
 * @param[in] input the time series (complex number as 2 floats)
 * @param[out] output the frequency series (complex number as 2 floats)
 */
void fft_forward(float *input, float *output)
{
  int i;
  for(i=0;i<fft_len;i++){
    fftin[i][0] = input[2*i];
    fftin[i][1] = input[2*i+1];
  }
  fftw_execute(fftpf);
  for(i=0;i<fft_len;i++){
    output[2*i] = fftout[i][0];
    output[2*i+1] = fftout[i][1];
  }
}

/**
 * perform a backward FFT (frequency to time)
 * @param[in] input the frequency series (complex number as 2 floats)
 * @param[out] output the time series (complex number as 2 floats)
 */
void fft_backward(float *input,float *output)
{
  int i;

  for(i=0;i<fft_len;i++){
    fftout[i][0] = input[2*i];
    fftout[i][1] = input[2*i+1];
  }
  fftw_execute(fftpb);
  for(i=0;i<fft_len;i++){
    output[2*i] = fftin[i][0]/fft_len;
    output[2*i+1] = fftin[i][1]/fft_len;
  }
}

/**
 * perform a forward FFT (time to frequency)
 * @param[in] in the time series (only real values!!)
 * @param[out] out_mag magnitudes of the frequency series
 * @param[out] out_phase phases of the frequency series
 */
void mag_and_phase(float *in, float *out_mag,float * out_phase)
{
  int i;
  for(i=0;i<fft_len;i++){
    fftin[i][0] = in[i];
    fftin[i][1] = 0.;
  }
  fftw_execute(fftpf);
  for(i=0;i<fft_len;i++){
    out_mag[i] = sqrt(fftout[i][0]*fftout[i][0]+fftout[i][1]*fftout[i][1]);
    out_phase[i] = atan2(fftout[i][1],fftout[i][0]);
  }
}

/**
 * perform a backward FFT (frequency to time)
 * @param[in] in_mag magnitudes of the frequency series
 * @param[in] in_phase phases of the frequency series
 * @param[out] out_trace the time series (only real part)
 */
void trace_from_mag_phase(float *in_mag,float * in_phase,float *out_trace)
{
  int i;
  for(i=0;i<fft_len;i++){
    fftout[i][0] = in_mag[i]*cos(in_phase[i]);
    fftout[i][1] = in_mag[i]*sin(in_phase[i]);
  }
  fftw_execute(fftpb);
  for(i=0;i<fft_len;i++){
    out_trace[i] = fftin[i][0]/fft_len;
  }
}

/**
 * create a Hilbert Envelope of a time series
 * @param[in] in input time series (float)
 * @param[out] out Hilbert Envelope (float)
 */
void envelope(float *in, float *out)
{
  int i;
  float tmp;
  for(i=0;i<fft_len;i++){
    fftin[i][0] = in[i];
    fftin[i][1] = 0.;
  }
  fftw_execute(fftpf);
  for(i=0;i<fft_len;i++){
    if(i<(fft_len/2)) {
      tmp = fftout[i][1];
      fftout[i][1] = -fftout[i][0];
    }
    else {
      tmp = -fftout[i][1];
      fftout[i][1] = fftout[i][0];
    }
    fftout[i][0] = tmp;
  }
  fftw_execute(fftpb);
  for(i=0;i<fft_len;i++){
    out[i] = fftin[i][0]*fftin[i][0]+fftin[i][1]*fftin[i][1];
    out[i] = sqrt(out[i]/(fft_len*fft_len)+in[i]*in[i]);
  }
}

/**
 * block filter the data
 * @param[in] in input time series (float)
 * @param[in] freqsample sampling frequency
 * @param[in] freqmin all frequencies below freqmin are blocked
 * @param[in] freqmax all frequencies above freqmax are blocked
 * @param[out] out filtered time series (float)
 */
void filter_data(float *in, float *out, float freqsample, float freqmin, float freqmax)
{
  int i;
  //printf("freqs: %g %g %g %d\n",freqsample,freqmin,freqmax,fft_len);
  for(i=0;i<fft_len;i++){
    fftin[i][0] = in[i];
    fftin[i][1] = 0;
  }
  fftw_execute(fftpf);
  if(freqmin>0){
    fftout[0][0] = 0;
    fftout[0][1] = 0;
  }
  for(i=1;i<fft_len/2;i++){
    if(freqsample*i/fft_len<freqmin ||freqsample*i/fft_len>freqmax){
      fftout[i][0] = 0;
      fftout[i][1] = 0;
      fftout[fft_len-i][0] = 0;
      fftout[fft_len-i][1] = 0;
    }
  }
  fftw_execute(fftpb);
  for(i=0;i<fft_len;i++){
    out[i] = fftin[i][0]/fft_len;
  }

}
