/*! \file antenna_util.c
* In this file we have routines to convert an electric field to voltage "apply antenna"
 * and to go from voltage to electric field "invert_antenna"
 * it requires the event structure from GRAND and the antenna definition
 * External package fftw is required
 * In addition, some calculations with complex numbers require the complex_util part of the GRAND library
 */
#include<string.h>
#include "fftw3.h"
#include "fftdata.h"
#include "antenna.h"
#include "GRANDevent.h"
#include "matrix_util.h"
#include "complex_util.h"
#define RADDEG 57.2957795131

extern int fft_len;
extern fftw_complex *fftin,*fftout;
extern fftw_plan fftpf,fftpb;

/**
 * Applying an antenna model to a GRAND detector in an event. The result of this routine will be that:
 *
 *    In the detectorinfo the antennamodel and electronics model are set (now fixed to GP300 and GRANDProto_V 2)\n
 *    For all 3 antenna arms an FFT of the simulated E-field is stored in SimE_fftmag and SimE_fftPhase\n
 *    This is multiplied with the effective length of the antenna, using the simulated direction of the incoming primary particle\n
 *    Next an inverse FFT is performed and the result is stored in the SimVoltage arrays for each antenna arm
 *  @param[in] evnt The GRAND event
 *  @param[in] iant The index of the detectordata  on which the antenna model will be applied
 *  @param[in] antenna The antenna model
 *
 */
void apply_antenna(GRANDEvent *evnt, int iant, struct grand_antenna *antenna)
{
  int i,iarm;
  GRANDDetectorData *tr = &(evnt->DetectorData[iant]);
  double Effl_mag[3]={0,0,0},Effl_phase[3]={0,0,0};
  fftw_complex Effl_mat[3][3];
  fftw_complex Efield[tr->SimPoints][3],Vval[tr->SimPoints][3];
  int fft_mid = tr->SimPoints/2+tr->SimPoints%2;

  sprintf(evnt->DetectorInfo[iant].DetectorModel,"GP300"); //fixed for now
  sprintf(evnt->DetectorInfo[iant].ElectronicsModel,"GRANDProto_V2"); //fixed for now
  memset((void *)Effl_mat,0,sizeof(Effl_mat));
  if(tr->SimPoints != fft_len) fft_init(tr->SimPoints);
  for(iarm=0;iarm<3;iarm++){
    for(i=0;i<fft_len;i++){
      fftin[i][0] = tr->SimEfield[iarm][i];
      fftin[i][1] = 0;
    }
    fftw_execute(fftpf);
    for(i=0;i<fft_len;i++){
      Efield[i][iarm][0] = fftout[i][0];
      Efield[i][iarm][1] = fftout[i][1];
      tr->SimE_fftMag[iarm][i] = c_mag(fftout[i]);
      tr->SimE_fftPhase[iarm][i] = c_phase(fftout[i]);
    }
  }
  for(i=0;i<fft_mid;i++){
    for(iarm=0;iarm<3;iarm++){
      antenna->effective_length(antenna,(enum grand_antenna_arm) iarm,(tr->SimMSPS/fft_len)*i, (evnt->SimEventInfo.Azimuth),(evnt->SimEventInfo.Zenith),Effl_mag,Effl_phase);
      for(int ior = 0;ior<3;ior++){ // create the mult. matrix
        Effl_mat[iarm][ior][0] = Effl_mag[ior]*cos(Effl_phase[ior]);
        Effl_mat[iarm][ior][1] = Effl_mag[ior]*sin(Effl_phase[ior]);
      }
    }
    c_matrix_times_vector(Effl_mat,Efield[i],Vval[i]);
  }
  for(iarm=0;iarm<3;iarm++){
    for(i=0;i<fft_mid;i++){
      fftout[i][0] = Vval[i][iarm][0];
      fftout[i][1] = Vval[i][iarm][1];
      if(i>0){
        fftout[fft_len-i][0] = Vval[i][iarm][0];
        fftout[fft_len-i][1] = -Vval[i][iarm][1];
      }
    }
    fftw_execute(fftpb);
    for(i=0;i<fft_len;i++){
      tr->SimVoltage[iarm][i] = fftin[i][0]/fft_len;
    }
  }
}

/**
 * inverting an antenna model using voltage traces in a GRAND detector in an event. The result of this routine will be that:
 *
 *    In the detectorinfo the antennamodel and electronics model are set (now fixed to GP300 and GRANDProto_V 2)\n
 *    For all 3 antenna arms an FFT of the raw voltages is created\n
 *    For each frequency the effective length of all antenna arms is obtained, using the reconstructed zenith and azimuth angles.\n
 *    For each frequency Matrix inversion is attempted to obtain the full 3D E-field from the 3 voltage values\n
 *    If the matrix inversion fails, the 3D e-field is obtained from 2 voltage values combined with the demand that the E-field is perpendicular to the incoming shower\n
 *    The result is stored in the RecE_fftPhase and RecE_fftMag arrays\n
 *    Next an inverse FFT is performed and the result is stored in the RecEfield arrays for each antenna arm
 *  @param[in] evnt The GRAND event
 *  @param[in] iant The index of the detectordata  on which the antenna model will be applied
 *  @param[in] antenna The antenna model
 *
 */
void inverse_antenna(GRANDEvent *evnt, int iant, struct grand_antenna *antenna)
{
  int i,iarm,a[3];
  GRANDDetectorData *tr = &(evnt->DetectorData[iant]);
  double Effl_mag[3],Effl_phase[3];
  fftw_complex shdir[3];
  fftw_complex det,ctemp[5];
  fftw_complex Effl_mat[3][3],Effl_inv[3][3];
  fftw_complex Effl_2mat[2][2],Effl_2inv[2][2];
  float Effl_rmat[3][3],Effl_rinv[3][3];
  fftw_complex Efield[tr->RawPoints][3],Vval[tr->RawPoints][3];
  int fft_mid = tr->RawPoints/2+tr->RawPoints%2;

  sprintf(evnt->DetectorInfo[iant].DetectorModel,"GP300"); //fixed for now
  sprintf(evnt->DetectorInfo[iant].ElectronicsModel,"GRANDProto_V2"); //fixed for now
  float rzen = evnt->RecEventInfo.Zenith/RADDEG;
  float razim =evnt->RecEventInfo.Azimuth/RADDEG;
  shdir[0][0] = cos(razim)*sin(rzen);
  shdir[0][1] = 0;
  shdir[1][0] = sin(razim)*sin(rzen);
  shdir[1][1] = 0;
  shdir[2][0] = cos(rzen);
  shdir[2][1] = 0;
  if(tr->RawPoints != fft_len) fft_init(tr->RawPoints);
  for(iarm = 0;iarm<3;iarm++){
    for(i=0;i<fft_len;i++){
      fftin[i][0] = tr->RawVoltage[iarm][i];
      fftin[i][1] = 0;
    }
    fftw_execute(fftpf);
    for(i=0;i<fft_len;i++){
      Vval[i][iarm][0] = fftout[i][0];
      Vval[i][iarm][1] = fftout[i][1];
    }
  }
  for(i=0;i<fft_mid;i++){
    for(iarm=0;iarm<3;iarm++){
      Efield[i][iarm][0] = 0; // initialize to zero
      Efield[i][iarm][1] = 0;
      antenna->effective_length(antenna,(enum grand_antenna_arm) iarm,(tr->RawMSPS/fft_len)*i, (evnt->RecEventInfo.Azimuth),(evnt->RecEventInfo.Zenith),Effl_mag,Effl_phase);
      if((tr->RawMSPS/fft_len)*i<35 || (tr->RawMSPS/fft_len)*i>200){ // get rid of instabilities
        Vval[i][iarm][0] = 0;
        Vval[i][iarm][1] = 0;
        if(i>0){
          Vval[fft_len-i][iarm][0] = 0;
          Vval[fft_len-i][iarm][1] = 0;
        }
      }
      for(int ior = 0;ior<3;ior++){ // create the mult. matrix
        Effl_mat[iarm][ior][0] = Effl_mag[ior]*cos(Effl_phase[ior]);
        Effl_rmat[iarm][ior] = Effl_mag[ior];
        Effl_mat[iarm][ior][1] = Effl_mag[ior]*sin(Effl_phase[ior]);
      }
    }
    if(Effl_rmat[0][0]<1E-10 && Effl_rmat[1][1]<1E-10 && Effl_rmat[2][2]<1E-10)continue; //gets rid of small and large frequencies in simulation
    if((c_invert_matrix(Effl_mat, Effl_inv)/(Effl_rmat[0][0]+Effl_rmat[1][1]+Effl_rmat[2][2]))< 1.E-10){
      //printf("Use 2D Inverted matrix\n");
      if((Effl_rmat[0][0]< Effl_rmat[1][1]) && (Effl_rmat[0][0]< Effl_rmat[2][2])){
        a[2] = 0;
        a[0] = 1;
        a[1] = 2;
      }
      if((Effl_rmat[1][1]< Effl_rmat[0][0]) && (Effl_rmat[1][1]< Effl_rmat[2][2])){
        a[2] = 1;
        a[0] = 0;
        a[1] = 2;
      }
      if((Effl_rmat[2][2]< Effl_rmat[1][1]) && (Effl_rmat[2][2]< Effl_rmat[0][0])){
        a[2] = 2;
        a[0] = 0;
        a[1] = 1;
      }
      c_sub(Effl_mat[a[0]][a[0]],*c_multiply(*c_divide(shdir[a[0]],shdir[a[2]],&ctemp[0]),Effl_mat[a[0]][a[2]],&ctemp[1]),&Effl_2mat[0][0]);
      c_sub(Effl_mat[a[0]][a[1]],*c_multiply(*c_divide(shdir[a[1]],shdir[a[2]],&ctemp[0]),Effl_mat[a[0]][a[2]],&ctemp[1]),&Effl_2mat[0][1]);
      c_sub(Effl_mat[a[1]][a[0]],*c_multiply(*c_divide(shdir[a[0]],shdir[a[2]],&ctemp[0]),Effl_mat[a[1]][a[2]],&ctemp[1]),&Effl_2mat[1][0]);
      c_sub(Effl_mat[a[1]][a[1]],*c_multiply(*c_divide(shdir[a[1]],shdir[a[2]],&ctemp[0]),Effl_mat[a[1]][a[2]],&ctemp[1]),&Effl_2mat[1][1]);
      c_sub(*c_multiply(Effl_2mat[0][0],Effl_2mat[1][1],&ctemp[0]),*c_multiply(Effl_2mat[1][0],Effl_2mat[0][1],&ctemp[1]),&det);
      c_divide(*c_sub(*c_multiply(Effl_2mat[1][1],Vval[i][a[0]],&ctemp[0]),*c_multiply(Effl_2mat[0][1],Vval[i][a[1]],&ctemp[1]),&ctemp[2]),det,&Efield[i][a[0]]);
      c_divide(*c_sub(*c_multiply(Effl_2mat[0][0],Vval[i][a[1]],&ctemp[0]),*c_multiply(Effl_2mat[1][0],Vval[i][a[0]],&ctemp[1]),&ctemp[2]),det,&Efield[i][a[1]]);
      c_divide(*c_add(*c_multiply(shdir[a[0]],Efield[i][a[0]],&ctemp[0]),*c_multiply(shdir[a[1]],Efield[i][a[1]],&ctemp[1]),&ctemp[2]),shdir[a[2]],
               &Efield[i][a[2]]);
      Efield[i][a[2]][0] = -Efield[i][a[2]][0];
      Efield[i][a[2]][1] = -Efield[i][a[2]][1];
      continue;
    }
    c_matrix_times_vector(Effl_inv,Vval[i],Efield[i]);
  }
  for(iarm=0;iarm<3;iarm++){
    for(i=0;i<fft_mid;i++){
      fftout[i][0] = Efield[i][iarm][0];
      fftout[i][1] = Efield[i][iarm][1];
      tr->RecE_fftMag[iarm][i] = c_mag(Efield[i][iarm]);
      tr->RecE_fftPhase[iarm][i] = c_phase(Efield[i][iarm]);
      if(i>0){
        fftout[fft_len-i][0] = Efield[i][iarm][0];
        fftout[fft_len-i][1] = -Efield[i][iarm][1];
        tr->RecE_fftMag[iarm][fft_len-i] = c_mag(Efield[fft_len-i][iarm]);
        tr->RecE_fftPhase[iarm][fft_len-i] = c_phase(Efield[fft_len-i][iarm]);
      }
    }
    fftw_execute(fftpb);
    //printf("Write Antenna data %s\n",evnt->antennainfo[iant].ID);
    for(i=0;i<fft_len;i++){
      tr->RecEfield[iarm][i] = fftin[i][0]/fft_len;
    }
  }
}
