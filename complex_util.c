/*! \file complex_util.c
 * some simple calculations using fftw_complex representation of complex numbers
 */
#include "complex_util.h"
#include "math.h"

/**
 * Multiplying complex numbers
 *  @param[in] a
 *  @param[in] b
 *  @param[out] res=a*b
 */
fftw_complex* c_multiply(fftw_complex a,fftw_complex b,fftw_complex *res)
{
  res[0][0] = a[0]*b[0] - a[1]*b[1];
  res[0][1] = a[1]*b[0] + a[0]*b[1];
  return(res);
}

/**
 * Dividing complex numbers
 *  @param[in] a
 *  @param[in] b
 *  @param[out] res=a/b
 */
fftw_complex* c_divide(fftw_complex a,fftw_complex b,fftw_complex *res)
{
  double bnorm;
  
  bnorm = b[0]*b[0]+b[1]*b[1];
  res[0][0] = (a[0]*b[0] + a[1]*b[1])/bnorm;
  res[0][1] = (a[1]*b[0] - a[0]*b[1])/bnorm;
  return(res);
}

/**
 * Adding complex numbers
 *  @param[in] a
 *  @param[in] b
 *  @param[out] res=a+b
 */
fftw_complex* c_add(fftw_complex a,fftw_complex b,fftw_complex *res)
{
  res[0][0] = a[0]+b[0];
  res[0][1] = a[1]+b[1];
  return(res);
}

/**
 * Subtracting complex numbers
 *  @param[in] a
 *  @param[in] b
 *  @param[out] res=a-b
 */
fftw_complex* c_sub(fftw_complex a,fftw_complex b,fftw_complex *res)
{
  res[0][0] = a[0]-b[0];
  res[0][1] = a[1]-b[1];
  return(res);
}

/**
 * Getting the angle of a complex number in the complex plane
 *  @param[in] a
 *  @result angle in radians
 */
double c_phase(fftw_complex a)
{
  double res;
  
  res = atan2(a[1],a[0]);
  return(res);
}

/**
 * Getting the magnitude of a complex number
 *  @param[in] a
 *  @result |a|
 */
double c_mag(fftw_complex a)
{
  double res;
  
  res = sqrt(a[0]*a[0]+a[1]*a[1]);
  return(res);
}

/**
 * raising a complex number to an arbitrary power
 *  @param[in] a the complex number
 *  @param[in] pw the power
 *  @param[out] res=a^pw
 */
fftw_complex* c_pow(fftw_complex a, float pw,fftw_complex *res)
{
  double phase = c_phase(a);
  double mag = c_mag(a);
  phase*=pw;
  mag = pow(mag,pw);
  res[0][0] = mag*cos(phase);
  res[0][1] = mag*sin(phase);
  return(res);
}
