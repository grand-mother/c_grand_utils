/*! \file hardware_util.c
This file contains an implementation of the hardware (electronics except antenna) The ROOT external package is required for random number generation
*/
#include "GRANDevent.h"
#include "TRandom2.h"

#define GRANDMSPS 500
#define THRESADC 100
#define GPSRES 10
#define GIGA 1000000000

/**
 * go from simulated voltages to raw ADC and voltage values. At this moment there is no proper trigger simulation, nor does it read the hardware configuration. It is a very basic "digitization" based upon the sampling rate and bit depth of the ADC. The trigger is met when the maximal value is larger than a hardcoded threshold.
 * @param[in] evnt the event
 * @param[in] iant the detector number
 */
void apply_hardware(GRANDEvent *evnt, int iant)
{
  int imax,amax;
  float dax;
  int i,iraw;
  GRANDDetectorData *det = &(evnt->DetectorData[iant]);
  double ADC_muV = (1.E-6)*10*16384/1.8;
  static TRandom2 RandTime;

  det->IsTriggered = false;
  imax = 0;
  amax = 0;
  for(int iarm=0;iarm<3;iarm++){
    i = 0;
    iraw = 0;
    dax = det->SimMSPS/det->RawMSPS;
    while(i<det->SimPoints){
      det->RawADC[iarm][iraw] = (int)(ADC_muV*det->SimVoltage[iarm][i]);
      //printf("Filling RawADC  %g %g\n",det->RawADC[iarm][iraw],ADC_muV*det->SimVoltage[iarm][i]);
      if(det->RawADC[iarm][iraw]>amax){
        imax = iraw;
        amax = det->RawADC[iarm][iraw];
      }
      if(det->RawADC[iarm][iraw]<-amax){
        imax = iraw;
        amax = -det->RawADC[iarm][iraw];
      }
      det->RawVoltage[iarm][iraw] = det->RawADC[iarm][iraw]/(ADC_muV);
      iraw++;
      i = iraw*dax;
    }
  }
  if(amax>THRESADC) {
    //printf("A triggered detector %d\n",iant);
    det->IsTriggered = true;
  } 
  det->RawGPS.Second = det->SimGPS.Second;
  det->RawGPS.NanoSec = det->SimGPS.NanoSec + RandTime.Gaus(0,GPSRES);
  if( det->RawGPS.NanoSec>=GIGA){
    det->RawGPS.NanoSec -=GIGA;
    det->RawGPS.Second +=1;
  }
  if( det->RawGPS.NanoSec<0){
    det->RawGPS.NanoSec +=GIGA;
    det->RawGPS.Second -=1;
  }
  det->TPulse = (imax*(1000/det->RawMSPS));
  det->e_DetTime = GPSRES;
  det->RawPoints = iraw;
}
