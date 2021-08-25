/*! \file reco_util.c
This file contains a geometrical reconstruction to be used in GRAND. Most routines are not to be used by the user. Only reconstruct_event is really for external use. The ROOT external package is required
*/
#include <stdlib.h>
#include <stdio.h>
#include<string.h>
#include <math.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TClass.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"
#include "GRANDevent.h"

#define SPLIGHT 0.299792458 // m/ns = km/musec
#define PI 3.141592653589793
#define RAD_TO_DEG 57.295779513082321
#define GIGA 1000000000

GRANDEvent *event;


/**
 * Calculate the distance between 2 point in 3D space. Routine intended for internal use only
 * @param[in] p1 point 1
 * @param[in] p2 point 2
 * @result distance between the points
 */
double distance(double *p1,double *p2)
{
  return(sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2])));
}

/**
 * Calculate the opening angle between 2 vectors,each given as 2 angles. Routine intended for internal use only
 * @param[in] phi1 azimuth first vector (radians)
 * @param[in] theta1 zenith first vector (radians)
 * @param[in] phi2 azimuth second vector (radians)
 * @param[in] theta2 zenith second vector (radians)
 * @result opening angle in radians
 */
double calcopang(double phi1,double theta1,double phi2,double theta2) //all inputs in radians!
{
  double co = (cos(phi1)*cos(phi2)+sin(phi1)*sin(phi2))*sin(theta1)*sin(theta2)
  +cos(theta1)*cos(theta2);
  return(acos(co));
}

/**
 * calculate the time difference (in nanoseconds) between 2 GPS time stamps. Routine intended for internal use only
 * @param[in] gps1
 * @param[in] gps2
 * @result gps1-gps2 (nanoseconds)
 */
double calcdt(GPS gps1, GPS gps2)
{
  double result;
  double dsec, dnsec;
  
  if(gps1.Second>gps2.Second){
    dsec = gps1.Second-gps2.Second;
    result = 1.E9*(double)dsec;
  }
  else {
    dsec = gps2.Second-gps1.Second;
    result = -1.E9*(double)dsec;
  }
  if(gps1.NanoSec>gps2.NanoSec){
    dnsec = gps1.NanoSec-gps2.NanoSec;
    result += (double)dnsec;
  }
  else{
    dnsec = gps2.NanoSec-gps1.NanoSec;
    result -= (double)dnsec;
  }
  return(result);
}

/**
 * provides an analytical least squares solution for the zenith and azimuth angle for an event, assuming all detectors are at the same altitude. Routine intended for internal use only. It uses the positions and pulse times of the individual triggered detectors as input
 * @param[out] phi azimuth (radians)
 * @param[out] theta zenith (radians)
 */
void calc_uv(double *phi, double *theta)
{
  // we assume here that all the stations tanks are at the same altitude
  // times are in ns
  // coordinates are in meters
  int i;
  double a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0, a8 = 0, a9 = 0;
  double u,v,t0;
  double rmstimesq;
  double coords[2];
  GRANDDetectorData *det,*det0;
  GRANDDetectorInfo *detI,*detI0;

  *theta = -1;
  *phi = -1;
  det0 = &(event->DetectorData[0]);
  detI0 = &(event->DetectorInfo[0]);
  for(i=0;i<event->n_det;i++ ) {
    det = &(event->DetectorData[i]);
    if(!det->IsTriggered) continue;
    detI = &(event->DetectorInfo[i]);
    rmstimesq =det->e_DetTime*det->e_DetTime + det0->e_DetTime*det0->e_DetTime;
    a1 += (detI->Position[0]-detI0->Position[0]) * (detI->Position[0]-detI0->Position[0])/rmstimesq;
    a2 += (detI->Position[0]-detI0->Position[0]) * (detI->Position[1]-detI0->Position[1])/rmstimesq;
    a3 += (detI->Position[0]-detI0->Position[0])/rmstimesq;
    a4 += (detI->Position[1]-detI0->Position[1]) * (detI->Position[1]-detI0->Position[1]) /rmstimesq;
    a5 += (detI->Position[1]-detI0->Position[1]) /rmstimesq;
    a6 += 1./rmstimesq;
    a7 += (detI->Position[0]-detI0->Position[0]) * (det->DetTime-det0->DetTime)/rmstimesq;
    a8 += (detI->Position[1]-detI0->Position[1]) * (det->DetTime-det0->DetTime)/rmstimesq;
    a9 += (det->DetTime-det0->DetTime)/rmstimesq;
    //printf("%d: %g %g %g %g %g %g %g %g %g\n",i,a1,a2,a3,a4,a5,a6,a7,a8,a9);
  }
  a1 /= SPLIGHT;
  a2 /= SPLIGHT;
  a4 /= SPLIGHT;
  
  double denominateur = a1*a4*a6*SPLIGHT
  -a1*a5*a5-a2*a2*a6*SPLIGHT
  +2*a2*a3*a5-a3*a3*a4;
  
  u = (-a7*(a4*a6*SPLIGHT-a5*a5)
       +a8*(a2*a6*SPLIGHT-a3*a5)
       -a9*(a2*a5-a3*a4)*SPLIGHT) / denominateur;
  
  v = (a7*(a2*a6*SPLIGHT-a3*a5)
       -a8*(a1*a6*SPLIGHT-a3*a3)
       +a9*(a1*a5-a2*a3)*SPLIGHT) / denominateur;
  
  t0 = (a7*(a2*a5-a3*a4)
        -a8*(a1*a5-a2*a3)
        +a9*(a1*a4-a2*a2)*SPLIGHT) / denominateur;
  
  double sq = sqrt(u*u+v*v);
  if(sq>1. & sq<=1.01) sq = 1; //new for GP300
  if( sq <= 1 ){
    *theta = asin(sq);
    *phi = fmod(atan2(v,u)+2.*PI,2.*PI);
  }else{
    printf("Calc_uv: U=%g V=%g T=%g sq = %g denom = %g\n",u,v,t0,sq,denominateur);
    //for(i=0;i<nstat;i++ ) {
    //  printf("%g %g %g %g\n",x[i],y[i],t[i]-t[0],sqrt((x[i]-x[0])*(x[i]-x[0])+(y[i]-y[0])*(y[i]-y[0]))/SPLIGHT);  //rmstimesq;
    //}
  }
}

/**
 * function used in the ROOT implementation of Minuit to fit a plane wave using the peak times and the locations of all triggered detectors
 */
void fplane(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Int_t i;
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  float inprod;
  GRANDDetectorData *det;
  GRANDDetectorInfo *detI;

  for (i=0;i<event->n_det; i++) {
    det = &(event->DetectorData[i]);
    if(!det->IsTriggered) continue;
    detI = &(event->DetectorInfo[i]);
    inprod = (par[3]-detI->Position[0])*cos(par[0])*sin(par[1]);
    inprod += (par[4]-detI->Position[1])*sin(par[0])*sin(par[1]);
    inprod += (par[5]-detI->Position[2])*cos(par[1]);
    delta  = (det->DetTime-par[2]-inprod/SPLIGHT);
    chisq += delta*delta/(det->e_DetTime*det->e_DetTime);
  }
  f = chisq;
}


/**
 * function used in the ROOT implementation of Minuit to fit an expanding sphere using the peak times and the locations of all triggered detectors
 */
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Int_t i;
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  Double_t dxmax = pow(10.,par[2]);
  double center[3],core[3];
  GRANDDetectorData *det;
  GRANDDetectorInfo *detI;

  core[0] = par[3];
  core[1] = par[4];
  core[2] = par[5];
  center[0] = par[3]+cos(par[0])*sin(par[1])*dxmax;
  center[1] = par[4]+sin(par[0])*sin(par[1])*dxmax;
  center[2] = par[5]+cos(par[1])*dxmax;
  for (i=0;i<event->n_det; i++) {
    det = &(event->DetectorData[i]);
    if(!det->IsTriggered) continue;
    detI = &(event->DetectorInfo[i]);
    delta  = (det->DetTime-par[6]-(distance(detI->Position,center) - distance(core,center))/SPLIGHT)
    /det->e_DetTime;
    chisq += delta*delta;
  }
  //printf("FCN: %g %g %g %g\n",par[0]*RAD_TO_DEG,par[1]*RAD_TO_DEG,1./curv, chisq);
  f = chisq;
}

/**
 * calculate the core position through the weighted average of the raw voltages\n
 * Next obtain the nearest detector to the core and use its timing as the core time\n
 * Finally calculate the time difference wrt the core time for each detector
 */
int reconstruct_core()
{
  double Emag,Wsum=0;
  GRANDDetectorData *det;
  GRANDDetectorInfo *detI;
  int iclose;
  double dclose,dist;
  double Vsum;

  for(int i=0;i<3;i++) event->RecEventInfo.Core[i] = 0;
  for(int i=0;i<event->n_det;i++){
    det = &(event->DetectorData[i]);
    if(!det->IsTriggered) continue;
    detI = &(event->DetectorInfo[i]);
    Vsum = 0;
    for(int it=0;it<det->RawPoints;it++){
      Emag = sqrt(det->RawVoltage[0][it]*det->RawVoltage[0][it]+
                  det->RawVoltage[1][it]*det->RawVoltage[1][it]+
                  det->RawVoltage[2][it]*det->RawVoltage[2][it]);
      if(Emag>10) Vsum+=Emag; //avoid adding noise...
    }
    Wsum+=pow(Vsum,2);
    for(int ic=0;ic<3;ic++) event->RecEventInfo.Core[ic] += pow(Vsum,2)*detI->Position[ic];
  }
  if(Wsum == 0) return(-1);
  for(int ic=0;ic<3;ic++)event->RecEventInfo.Core[ic] /=Wsum;
  dclose = -1;
  for(int i=0;i<event->n_det;i++){
    detI = &(event->DetectorInfo[i]);
    det = &(event->DetectorData[i]);
    if(!det->IsTriggered) continue;
    dist = distance(event->RecEventInfo.Core,detI->Position);
    if(dist<dclose || dclose <0){
      dclose = dist;
      iclose = i;
    }
  }
  event->RecEventInfo.CoreTime.Second =event->DetectorData[iclose].RawGPS.Second;
  event->RecEventInfo.CoreTime.NanoSec =
    event->DetectorData[iclose].RawGPS.NanoSec+(int)event->DetectorData[iclose].TPulse;
  for(int i=0;i<event->n_det;i++){
    det = &(event->DetectorData[i]);
    det->DetTime = calcdt(det->RawGPS,event->RecEventInfo.CoreTime)+det->TPulse;
  }
  return(0);
}
/**
 Perform a geometrical reconstruction of the event as follows:\n
 1. reconstruct the core position; this also sets the time of each station wrt the core. This is used in:
 2. analytical plane wave calculation to get initial guesses of the zenith and azimuth. Start values for
 3. minuit plane wave fit, which uses the  core and core timing as a reference for all. The result is a start for
 4. minuit expanding sphere fit as follows:\n
  a. fit distance to xmax (leaving all others fixed)\n
  b. fit the zenith and azimuth, fixing all others\n
  c. fit the core position, fixing all others\n
  d. fit core,direction and distance
 5. Record all reconstructed parameters in the event; set the proper core time and time differences
 */
int reconstruct_event(GRANDEvent *evt)
{
  Double_t arglist[10];
  Double_t phi,theta,tcore;
  GRANDDetectorData *det;
  int itry;
  Int_t ierflg = 0;
  //make sure the event is known
  event = evt;
  //
  //analytical plane wave approach, assuming Z locations of antennas are all the same
  //
  reconstruct_core();
  calc_uv(&phi,&theta);
  if(phi <0 && theta <0) return(-1); //cannot initiate!
  evt->RecEventInfo.Azimuth = phi*RAD_TO_DEG;
  evt->RecEventInfo.Zenith = theta*RAD_TO_DEG;
  //
  // Fit direction assuming a plane wave
  //
  TMinuit *gMplane = new TMinuit(6);  //initialize TMinuit with a maximum of 6 params
  gMplane->SetFCN(fplane);
  gMplane->Command("SET PRINT -1");
  gMplane->mnexcm("SET NOWarnings",0,0,ierflg);
  arglist[0] = 1;
  gMplane->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  Double_t pstart[6] = {phi, theta , 0 , event->RecEventInfo.Core[0],event->RecEventInfo.Core[1],event->RecEventInfo.Core[2]};
//    GP300[90].location[0],GP300[90].location[1],GP300[90].location[2]};
  Double_t pstep[6] = {0.05 , 0.05 , 200. , 100,100,100};
  Double_t pbmin[6] = {0.,     0.,-2000.,event->RecEventInfo.Core[0]-10000., event->RecEventInfo.Core[1]-10000.,event->RecEventInfo.Core[2]-500.};
  Double_t pbmax[6] = {2*PI,PI/2.,2000, event->RecEventInfo.Core[0]+10000,    event->RecEventInfo.Core[1]+10000.,event->RecEventInfo.Core[2]+500.};
  gMplane->mnparm(0, "Phi", pstart[0], pstep[0], pbmin[0],pbmax[0],ierflg);
  gMplane->mnparm(1, "Theta", pstart[1], pstep[1], pbmin[1],pbmax[1],ierflg);
  gMplane->mnparm(2, "TCore", pstart[2], pstep[2], pbmin[2],pbmax[2],ierflg);
  gMplane->mnparm(3, "CoreX", pstart[3], pstep[3], pbmin[3],pbmax[3],ierflg);
  gMplane->mnparm(4, "CoreY", pstart[4], pstep[4], pbmin[4],pbmax[4],ierflg);
  gMplane->mnparm(5, "CoreZ", pstart[5], pstep[5], pbmin[5],pbmax[5],ierflg);
  // Minuit parameters: Core position (assumed to be fixed)
  gMplane->FixParameter(3);
  gMplane->FixParameter(4);
  gMplane->FixParameter(5);
  arglist[0] = 1000;
  arglist[1] = 2*evt->n_det; //tolerance to minimum
  ierflg = -1;
  itry = 0;
  while(ierflg != 0 && itry<3){
    gMplane->mnexcm("MIGRAD", arglist ,2,ierflg);
    itry++;
  }
  gMplane->GetParameter(0,phi,pstep[0]);
  gMplane->GetParameter(1,theta,pstep[1]);
  gMplane->GetParameter(2,tcore,pstep[2]);
  evt->RecEventInfo.Azimuth = phi*RAD_TO_DEG;
  evt->RecEventInfo.Zenith = theta*RAD_TO_DEG;
  delete gMplane;
  //
  // Fit expanding sphere
  //
  TMinuit *gMinuit = new TMinuit(7);  //initialize TMinuit with a maximum of 6 params
  gMinuit->SetFCN(fcn);
  arglist[0] = 1;
  gMinuit->Command("SET PRINT -1");
  gMinuit->mnexcm("SET NOWarnings",0,0,ierflg);
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  Double_t vstart[7] = {phi, theta , 5 ,event->RecEventInfo.Core[0],event->RecEventInfo.Core[1],event->RecEventInfo.Core[2],0};
  //    GP300[90].location[0],GP300[90].location[1],GP300[90].location[2]};
  Double_t step[7] = {0.05 , 0.05 , 0.1 , 100,    100,   10,100};
  Double_t bmin[7] = {0.,     0.,-10,event->RecEventInfo.Core[0]-10000., event->RecEventInfo.Core[1]-10000.,event->RecEventInfo.Core[2]-10.,-3000};
  Double_t bmax[7] = {2*PI,1.1*PI/2.,20, event->RecEventInfo.Core[0]+10000,    event->RecEventInfo.Core[1]+10000.,event->RecEventInfo.Core[2]+10.,3000};
  gMinuit->mnparm(0, "Phi", vstart[0], step[0], bmin[0],bmax[0],ierflg);
  gMinuit->mnparm(1, "Theta", vstart[1], step[1], bmin[1],bmax[1],ierflg);
  gMinuit->mnparm(2, "10log(dxmax)", vstart[2], step[2], bmin[2],bmax[2],ierflg);
  gMinuit->mnparm(3, "CoreX", vstart[3], step[3], bmin[3],bmax[3],ierflg);
  gMinuit->mnparm(4, "CoreY", vstart[4], step[4], bmin[4],bmax[4],ierflg);
  gMinuit->mnparm(5, "CoreZ", vstart[5], step[5], bmin[5],bmax[5],ierflg);
  gMinuit->mnparm(6, "CoreTime", vstart[6], step[6], bmin[6],bmax[6],ierflg);
  // Minuit parameters: Center of shower and Core position (assumed to be fixed)
  gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->FixParameter(3);
  gMinuit->FixParameter(4);
  gMinuit->FixParameter(5);
  
  // Now ready for minimization step
  //step 1 reconstruct Distance
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = -1;
  itry = 0;
  //printf("Step 1\n");
  while(itry<3&&ierflg !=0){
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    itry++;
  }
  //step 2 reconstruct angle
  gMinuit->Release(0);
  gMinuit->Release(1);
  gMinuit->FixParameter(2);
  gMinuit->FixParameter(6);
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = -1;
  itry = 0;
  //printf("Step 2\n");
  while(itry<3&&ierflg !=0){
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    itry++;
  }
  //step 3 reconstruct core
  //gMinuit->Release(2);
  gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->Release(3);
  gMinuit->Release(4);
  gMinuit->Release(6);
  //gMinuit->Release(5);
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = -1;
  itry = 0;
  //printf("Step 3\n");
  while(itry<3&&ierflg !=0){
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    itry++;
  }
  //step 4 reconstruct both
  gMinuit->Release(0);
  gMinuit->Release(1);
  gMinuit->Release(2);
  //gMinuit->FixParameter(3);
  //gMinuit->FixParameter(4);
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = -1;
  itry = 0;
  //printf("Step 3\n");
  while(itry<3&&ierflg !=0){
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    itry++;
  }
  //if(ierflg != 4) printf("Error Flag... %d\n",ierflg);
  gMinuit->GetParameter(0,evt->RecEventInfo.Azimuth,evt->RecEventInfo.e_Azimuth);
  gMinuit->GetParameter(1,evt->RecEventInfo.Zenith,evt->RecEventInfo.e_Zenith);
  gMinuit->GetParameter(2,evt->RecEventInfo.XmaxDistance,evt->RecEventInfo.e_XmaxDistance);
  gMinuit->GetParameter(3,evt->RecEventInfo.Core[0],evt->RecEventInfo.e_Core[0]);
  gMinuit->GetParameter(4,evt->RecEventInfo.Core[1],evt->RecEventInfo.e_Core[1]);
  gMinuit->GetParameter(5,evt->RecEventInfo.Core[2],evt->RecEventInfo.e_Core[2]);
  gMinuit->GetParameter(6,tcore,step[6]);
  event->RecEventInfo.CoreTime.NanoSec += tcore;
  if(event->RecEventInfo.CoreTime.NanoSec>GIGA){
    event->RecEventInfo.CoreTime.NanoSec -=GIGA;
    event->RecEventInfo.CoreTime.Second +=1;
  }
  if(event->RecEventInfo.CoreTime.NanoSec<0){
    event->RecEventInfo.CoreTime.NanoSec +=GIGA;
    event->RecEventInfo.CoreTime.Second -=1;
  }
  for(int i=0;i<event->n_det;i++){
    det = &(event->DetectorData[i]);
    det->DetTime = calcdt(det->RawGPS,event->RecEventInfo.CoreTime)+det->TPulse;
  }
  //gMinuit->DeleteArrays();
  delete gMinuit;
  evt->RecEventInfo.Azimuth =RAD_TO_DEG*evt->RecEventInfo.Azimuth;
  evt->RecEventInfo.e_Azimuth =RAD_TO_DEG*evt->RecEventInfo.e_Azimuth;
  evt->RecEventInfo.Zenith =RAD_TO_DEG*evt->RecEventInfo.Zenith;
  evt->RecEventInfo.e_Zenith =RAD_TO_DEG*evt->RecEventInfo.e_Zenith;
  evt->RecEventInfo.XmaxDistance = pow(10,evt->RecEventInfo.XmaxDistance);
  evt->RecEventInfo.e_XmaxDistance = evt->RecEventInfo.XmaxDistance*log(10.)*evt->RecEventInfo.e_XmaxDistance;
  return(0);
}
