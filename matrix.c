/*! \file matrix.c
This file contains routines that aid matrix manipulation
*/
#include "complex_util.h"

/**
 * Invert a real 3x3 matrix
 * @param[in] in original matrix
 * @param[out] out inverted matrix
 */
float invert_matrix(float in[3][3], float out[3][3])
{
  int i, j;
  float determinant = 0;
  float *orig_mat[3],*inv_mat[3];
  
  orig_mat[0] = (float *)in;
  orig_mat[1] = orig_mat[0]+3;
  orig_mat[2] = orig_mat[1]+3;
  inv_mat[0] = (float *)out;
  inv_mat[1] = inv_mat[0]+3;
  inv_mat[2] = inv_mat[1]+3;
  //finding determinant
  for(i = 0; i < 3; i++){
    determinant = determinant + (orig_mat[0][i] * (orig_mat[1][(i+1)%3] * orig_mat[2][(i+2)%3] - orig_mat[1][(i+2)%3] * orig_mat[2][(i+1)%3]));
  }
  printf("Det = %g\n",determinant);
  if(determinant == 0) return(-1);
  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++)
      inv_mat[i][j] = ((orig_mat[(j+1)%3][(i+1)%3] * orig_mat[(j+2)%3][(i+2)%3]) - (orig_mat[(j+1)%3][(i+2)%3] * orig_mat[(j+2)%3][(i+1)%3]))/ determinant;
  }
  return(determinant);
}

/**
 * Invert a complex 3x3 matrix
 * @param[in] in original matrix
 * @param[out] out inverted matrix
 */
float c_invert_matrix(fftw_complex in[3][3], fftw_complex out[3][3])
{
  int i, j;
  fftw_complex determinant;
  fftw_complex c_temp[5];
  
  determinant[0] = 0;
  determinant[1] = 0;
  //finding determinant
  for(i = 0; i < 3; i++){
    c_multiply(in[1][(i+1)%3],in[2][(i+2)%3],&c_temp[0]);
    c_multiply(in[1][(i+2)%3],in[2][(i+1)%3],&c_temp[1]);
    c_sub(c_temp[0],c_temp[1],&c_temp[2]);
    c_multiply(in[0][i],c_temp[2],&c_temp[3]);
    c_add(determinant,c_temp[3],&determinant);
   }
  //printf("Det = (%g,%g)\n",determinant[0],determinant[1]);
  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++){
      c_multiply(in[(j+1)%3][(i+1)%3],in[(j+2)%3][(i+2)%3],&c_temp[0]);
      c_multiply(in[(j+1)%3][(i+2)%3],in[(j+2)%3][(i+1)%3],&c_temp[1]);
      c_sub(c_temp[0],c_temp[1],&c_temp[2]);
      c_divide(c_temp[2],determinant,&out[i][j]);
    }
  }
  
  return(c_mag(determinant));
}

/**
 * Multiply a complex 3x3 matrix with a 3D complex vector
 * @param[in] mat  matrix
 * @param[in] vec input vector
 * @param[out] result resulting vector
 */
void c_matrix_times_vector(fftw_complex mat[3][3],fftw_complex vec[3],fftw_complex *result)
{
  
  int i,j;
  fftw_complex c_temp;
  for(i=0;i<3;i++){
    result[i][0] = 0;
    result[i][1] = 0;
    for(j=0;j<3;j++){
      c_multiply(mat[i][j],vec[j],&c_temp);
      c_add(result[i],c_temp,&result[i]);
    }
  }
}
