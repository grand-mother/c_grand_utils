/*! \file aires_util.c
This file contains the input of AIRES events and the conversion to GRAND. It requires hdf5. For use outside the library are the routines read_aires_event, clear_aires_event and convert_aires_GRAND
*/
#include<stdlib.h>
#include<stdio.h>
#include <string.h>
#include"hdf5.h"
#include "AIRESevent.h"
#include "GRANDevent.h"
#define GRANDMSPS 500

/**
 * Convert an ZHAires event (as given in the files from Matias) to a GRAND event
 * @param[in] Aires the ZHaires event
 * @param[in] GRAND the GRAND event in which the AIRES information is copied
 */
void convert_aires_GRAND(AIRESEvent *Aires, GRANDEvent *GRAND)
{
  char *c;
  int i,id,iarm;
  int Sim_tot_data = 0,Raw_tot_data = 0;
  // Fill general event info
  c = Aires->eventinfo.EventName;
  for(i=0; &c[i] !=Aires->eventinfo.EventID && i<79; i++)GRAND->GenEventInfo.EventName[i] = c[i];
  GRAND->GenEventInfo.EventName[i] = 0;
  c = Aires->eventinfo.EventID;
  for(i=0; &c[i] !=Aires->eventinfo.Primary && i<79; i++)GRAND->GenEventInfo.EventID[i] = c[i];
  GRAND->GenEventInfo.EventID[i] = 0;
  c = Aires->eventinfo.Site;
  for(i=0; &c[i] !=Aires->eventinfo.Date && i<79; i++)GRAND->GenEventInfo.Site[i] = c[i];
  GRAND->GenEventInfo.Site[i] = 0;
  c = Aires->eventinfo.Date;
  for(i=0; &c[i] !=Aires->eventinfo.Latitude && i<79; i++)GRAND->GenEventInfo.Date[i] = c[i];
  GRAND->GenEventInfo.Date[i] = 0;
  sscanf(Aires->eventinfo.Latitude,"%lg",&GRAND->GenEventInfo.Latitude);
  sscanf(Aires->eventinfo.Longitude,"%lg",&GRAND->GenEventInfo.Longitude);
  GRAND->GenEventInfo.GroundAltitude = *Aires->eventinfo.GroundAltitude;
  GRAND->GenEventInfo.BField =*Aires->eventinfo.BField;
  GRAND->GenEventInfo.BFieldIncl =*Aires->eventinfo.BFieldIncl;
  GRAND->GenEventInfo.BFieldDecl =*Aires->eventinfo.BFieldDecl;
  c = Aires->eventinfo.AtmosphericModel;
  for(i=0; &c[i] !=Aires->eventinfo.AtmosphericModelParameters && i<79; i++)
    GRAND->GenEventInfo.AtmosphericModel[i] = c[i];
  GRAND->GenEventInfo.AtmosphericModel[i] = 0;
  c = Aires->eventinfo.AtmosphericModelParameters;
  for(i=0; &c[i] !=(char *)Aires->eventinfo.EnergyInNeutrinos && i<99; i++)
    GRAND->GenEventInfo.AtmosphericModelParameters[i] = c[i];
  GRAND->GenEventInfo.AtmosphericModel[i] = 0;
  // Fill simulation event info
  c = Aires->eventinfo.Primary;
  for(i=0; &c[i] !=(char *)Aires->eventinfo.Energy && i<79; i++)
    GRAND->SimEventInfo.Primary[i] = c[i];
  GRAND->SimEventInfo.Primary[i] = 0;
  GRAND->SimEventInfo.Energy = *Aires->eventinfo.Energy;
  GRAND->SimEventInfo.Zenith = 180 - (*Aires->eventinfo.Zenith);
  GRAND->SimEventInfo.Azimuth = 180+ (*Aires->eventinfo.Azimuth);
  if(GRAND->SimEventInfo.Azimuth>360) GRAND->SimEventInfo.Azimuth -=360;
  GRAND->SimEventInfo.XmaxDistance = *Aires->eventinfo.XmaxDistance;
  for(i=0;i<3;i++) GRAND->SimEventInfo.XmaxPosition[i] =Aires->eventinfo.XmaxPosition[i];
  GRAND->SimEventInfo.XmaxAltitude = *Aires->eventinfo.XmaxAltitude;
  GRAND->SimEventInfo.SlantXmax = *Aires->eventinfo.SlantXmax;
  GRAND->SimEventInfo.InjectionAltitude = *Aires->eventinfo.InjectionAltitude;
  GRAND->SimEventInfo.EnergyInNeutrinos = *Aires->eventinfo.EnergyInNeutrinos;
  //The detector positions
  GRAND->n_det = Aires->n_ant;
  if(GRAND->DetectorInfo != NULL) free((void *)GRAND->DetectorInfo);
  GRAND->DetectorInfo = (GRANDDetectorInfo *)malloc(GRAND->n_det*sizeof(GRANDDetectorInfo));
  if(GRAND->DetectorData != NULL) free((void *)GRAND->DetectorData);
  GRAND->DetectorData = (GRANDDetectorData *)malloc(GRAND->n_det*sizeof(GRANDDetectorData));
  Sim_tot_data = 0;
  Raw_tot_data = 0;
  for(id=0;id<GRAND->n_det;id++){
    c =Aires->antennainfo[id].ID;
    for(i=0; &c[i] !=(char *)Aires->antennainfo[id].pos && i<79; i++)
      GRAND->DetectorInfo[id].ID[i] = c[i];
    GRAND->DetectorInfo[id].ID[i] = 0;
    for(i=0;i<3;i++) GRAND->DetectorInfo[id].Position[i] = Aires->antennainfo[id].pos[i];
    GRAND->DetectorData[id].SimPoints =Aires->antennatrace[id].n_point;
    i =Aires->antennatrace[id].n_point;
    GRAND->DetectorData[id].SimMSPS =1000*(i-1)/(Aires->antennatrace[id].time[i-1]-Aires->antennatrace[id].time[0]);
    GRAND->DetectorData[id].RawMSPS = GRANDMSPS;
    GRAND->DetectorData[id].RawPoints =1+(int)(float)(Aires->antennatrace[id].n_point*GRAND->DetectorData[id].RawMSPS/GRAND->DetectorData[id].SimMSPS);
    GRAND->DetectorData[id].SimGPS.Second = 100;
    GRAND->DetectorData[id].SimGPS.NanoSec =(unsigned int)(500000000 + Aires->antennatrace[id].time[0]);
    Sim_tot_data += GRAND->DetectorData[id].SimPoints;
    Raw_tot_data += GRAND->DetectorData[id].RawPoints;
  }
  if(GRAND->TraceBuffer != NULL) free((void *)GRAND->TraceBuffer);
  GRAND->TraceBuffer = (float *)malloc(Sim_tot_data*12*sizeof(float)+Raw_tot_data*15*sizeof(float));
  // fill traces
  i = 0;
  for(id=0;id<GRAND->n_det;id++){
    GRAND->DetectorData[id].SimEfield[0] = &GRAND->TraceBuffer[i];
    memcpy((void *)GRAND->DetectorData[id].SimEfield[0],(void *)Aires->antennatrace[id].Ex,
           GRAND->DetectorData[id].SimPoints*sizeof(float));
    i+= GRAND->DetectorData[id].SimPoints;
    GRAND->DetectorData[id].SimEfield[1] = &GRAND->TraceBuffer[i];
    memcpy((void *)GRAND->DetectorData[id].SimEfield[1],(void *)Aires->antennatrace[id].Ey,
           GRAND->DetectorData[id].SimPoints*sizeof(float));
    i+= GRAND->DetectorData[id].SimPoints;
    GRAND->DetectorData[id].SimEfield[2] = &GRAND->TraceBuffer[i];
    memcpy((void *)GRAND->DetectorData[id].SimEfield[2],(void *)Aires->antennatrace[id].Ez,
           GRAND->DetectorData[id].SimPoints*sizeof(float));
    i+= GRAND->DetectorData[id].SimPoints;
    for(iarm=0;iarm<3;iarm++){
      GRAND->DetectorData[id].SimE_fftMag[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].SimPoints;
      GRAND->DetectorData[id].SimE_fftPhase[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].SimPoints;
      GRAND->DetectorData[id].SimVoltage[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].SimPoints;
    }
    for(iarm=0;iarm<3;iarm++){
      GRAND->DetectorData[id].RawADC[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].RawPoints;
      GRAND->DetectorData[id].RawVoltage[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].RawPoints;
      GRAND->DetectorData[id].RecEfield[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].RawPoints;
      GRAND->DetectorData[id].RecE_fftMag[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].RawPoints;
      GRAND->DetectorData[id].RecE_fftPhase[iarm] = &GRAND->TraceBuffer[i];
      i+= GRAND->DetectorData[id].RawPoints;
    }
  }
}
/**
 * release all the allocated memories connected to an Aires event and set all elements to zero
 * @param[in] event the AIRES event to be cleared
 */
void clear_aires_event(AIRESEvent *event)
{
  if(event->runinfo.buffer != NULL)
    free(event->runinfo.buffer);
  if(event->eventinfo.buffer != NULL)
    free(event->eventinfo.buffer);
  if(event->showersiminfo.buffer != NULL)
    free(event->showersiminfo.buffer);
  if(event->signalsiminfo.buffer != NULL)
    free(event->signalsiminfo.buffer);
  if(event->showertable.lateralbuffer != NULL)
    free(event->showertable.lateralbuffer);
  if(event->showertable.lateralprofile != NULL)
    free(event->showertable.lateralprofile);
  if(event->showertable.longbuffer != NULL)
    free(event->showertable.longbuffer);
  if(event->showertable.longprofile != NULL)
    free(event->showertable.longprofile);
  if(event->antennabuffer != NULL)
    free(event->antennabuffer);
  if(event->antennainfo != NULL)
    free(event->antennainfo);
  if(event->tracebuffer != NULL)
    free(event->tracebuffer);
  if(event->antennatrace != NULL)
    free(event->antennatrace);
  if(event->antennaP2Pbuffer != NULL)
    free(event->antennaP2Pbuffer);
  if(event->antennap2p != NULL)
    free(event->antennap2p);
  memset(event,0,sizeof(AIRESEvent));
}
/**
 * read in the runinfo table from an Aires event file
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the runinfo is stored
 */
int read_runinfo(hid_t file, AIRESEvent *event)
{
  hid_t group;
  hid_t dataset,dataspace,datatype; //freq,phi,theta
  int rank,Class,sz;
  herr_t status;
  char *name;
  
  
  if((group = H5Gopen(file, "/", H5P_DEFAULT))<0) return(-1);
  //1. read Run info
  if((dataset = H5Dopen(group, "RunInfo", H5P_DEFAULT))<0) {
    clear_aires_event(event);
    memset(event,0,sizeof(AIRESEvent));
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->runinfo.size < sz){
    if(event->runinfo.size>0) free(event->runinfo.buffer);
    event->runinfo.buffer = (char *)malloc(sz);
    if(event->runinfo.buffer ==NULL)return(-3);
    event->runinfo.size = sz;
  }
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,event->runinfo.buffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"EventName") == 0)
      event->runinfo.EventName =
      (char *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"EventID") == 0)
      event->runinfo.EventID =
      (char *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Primary") == 0)
      event->runinfo.Primary =
      (char *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Energy") == 0)
      event->runinfo.Energy =
      (double *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Zenith") == 0)
      event->runinfo.Zenith =
      (double *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Azimuth") == 0)
      event->runinfo.Azimuth =
      (double *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"XmaxDistance") == 0)
      event->runinfo.XmaxDistance =
      (double *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"SlantXmax") == 0)
      event->runinfo.SlantXmax =
      (double *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"HadronicModel") == 0)
      event->runinfo.HadronicModel =
      (char *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"InjectionAltitude") == 0)
      event->runinfo.InjectionAltitude =
      (double *)&(event->runinfo.buffer[H5Tget_member_offset(datatype,i)]);
    H5free_memory(name);
  }
  event->runinfo.EventID[0] = 0; //for now a quick fix!
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  return(0);
}

/**
 * read in the antennainfo table from an Aires event file
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_antennainfo(hid_t file, AIRESEvent *event)
{
  hid_t group;
  hid_t dataset,dataspace,datatype; //freq,phi,theta
  int rank,Class,sz;
  herr_t status;
  hsize_t dims,max_dims;
  char *name;
  int iant;
  
  
  if((group = H5Gopen(file, event->runinfo.EventName, H5P_DEFAULT))<0) return(-1);
  //1. read Run info
  if((dataset = H5Dopen(group, "AntennaInfo", H5P_DEFAULT))<0) {
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->antennabuffer == NULL){
    event->antennabuffer = (char *)malloc(sz*dims);
    if(event->antennabuffer == NULL)return(-3);
    event->antennainfo = (Antennainfo *)malloc(sizeof(Antennainfo)*dims);
    if(event->antennainfo ==NULL)return(-3);
  }
  event->n_ant = dims;
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,event->antennabuffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"ID") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennainfo[iant].ID =
        (char *)&event->antennabuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"X") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennainfo[iant].pos =
        (float *)&event->antennabuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"T0") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennainfo[iant].T0 =
        (float *)&event->antennabuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"SlopeA") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennainfo[iant].SlopeA =
        (double *)&event->antennabuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"SlopeB") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennainfo[iant].SlopeB =
        (double *)&event->antennabuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  return(0);
}

/**
 * read in the antenna traces group from an Aires event file. This contains the e-field, voltages and filtered voltages for all simulated antennas
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_antennatraces(hid_t file, AIRESEvent *event)
{
  hid_t group,agroup;
  hid_t dataset,dataspace,datatype; //freq,phi,theta
  int rank,Class,sz;
  herr_t status;
  hsize_t dims,max_dims;
  char *name,*ch;
  char groupname[200];
  int iant,i;
  int npoint,ipoint;
  float *tbuf;
  
  
  sprintf(groupname,"%s/AntennaTraces",event->runinfo.EventName);
  if((group = H5Gopen(file, groupname, H5P_DEFAULT))<0) return(-1);
  npoint = 0;
  for(iant = 0;iant<event->n_ant;iant++){
    i = 0;
    for(ch = event->antennainfo[iant].ID;ch!=(char *)event->antennainfo[iant].pos;ch++)
      groupname[i++] =*ch;
    groupname[i] = 0;
    if((agroup = H5Gopen(group,groupname, H5P_DEFAULT))<0) return(-1);
  //1. Determine the size of the e-field vectors
    if((dataset = H5Dopen(agroup, "efield", H5P_DEFAULT))<0) {
      return(-1);
    }
    dataspace = H5Dget_space(dataset);
    datatype = H5Dget_type (dataset);
    Class = H5Tget_class(datatype);
    if(Class != H5T_COMPOUND) return(-2);
    rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
    npoint+=dims;
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Gclose(agroup);
  }
  event->tracebuffer = (char *)malloc((1+3*3)*npoint*sizeof(float)); //T, (E,V,Vf)*(x,y,z)
  event->antennatrace = (AntennaTrace *)malloc(event->n_ant*sizeof(AntennaTrace));
  tbuf = (float *)malloc(4*npoint*sizeof(float)); //tempory buffer, only for this routine!
  ipoint = 0;
  for(iant = 0;iant<event->n_ant;iant++){
    i = 0;
    for(ch = event->antennainfo[iant].ID;ch!=(char *)event->antennainfo[iant].pos;ch++)
      groupname[i++] =*ch;
    groupname[i] = 0;
    if((agroup = H5Gopen(group,groupname, H5P_DEFAULT))<0) return(-1);
  //1. Determine the size of the e-field vectors
    if((dataset = H5Dopen(agroup, "efield", H5P_DEFAULT))<0) {
      return(-1);
    }
    dataspace = H5Dget_space(dataset);
    datatype = H5Dget_type (dataset);
    Class = H5Tget_class(datatype);
    if(Class != H5T_COMPOUND) return(-2);
    rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
    event->antennatrace[iant].n_point = dims;
    event->antennatrace[iant].time = (float *)&event->tracebuffer[ipoint];
    event->antennatrace[iant].Ex = &event->antennatrace[iant].time[dims];
    event->antennatrace[iant].Ey = &event->antennatrace[iant].Ex[dims];
    event->antennatrace[iant].Ez = &event->antennatrace[iant].Ey[dims];
    event->antennatrace[iant].Vx = &event->antennatrace[iant].Ez[dims];
    event->antennatrace[iant].Vy = &event->antennatrace[iant].Vx[dims];
    event->antennatrace[iant].Vz = &event->antennatrace[iant].Vy[dims];
    event->antennatrace[iant].Vfx = &event->antennatrace[iant].Vz[dims];
    event->antennatrace[iant].Vfy = &event->antennatrace[iant].Vfx[dims];
    event->antennatrace[iant].Vfz = &event->antennatrace[iant].Vfy[dims];
    ipoint+=dims*10*sizeof(float);
    status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,tbuf);
    for(int i=0;i<dims;i++){
      event->antennatrace[iant].time[i] = tbuf[4*i];
      event->antennatrace[iant].Ex[i] = tbuf[4*i+1];
      event->antennatrace[iant].Ey[i] = tbuf[4*i+2];
      event->antennatrace[iant].Ez[i] = tbuf[4*i+3];
    }
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    if((dataset = H5Dopen(agroup, "voltage", H5P_DEFAULT))<0) {
      return(-1);
    }
    dataspace = H5Dget_space(dataset);
    datatype = H5Dget_type (dataset);
    Class = H5Tget_class(datatype);
    if(Class != H5T_COMPOUND) return(-2);
    status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,tbuf);
    rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
    for(int i=0;i<dims;i++){
      event->antennatrace[iant].Vx[i] = tbuf[4*i+1];
      event->antennatrace[iant].Vy[i] = tbuf[4*i+2];
      event->antennatrace[iant].Vz[i] = tbuf[4*i+3];
    }
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    if((dataset = H5Dopen(agroup, "filteredvoltage", H5P_DEFAULT))<0) {
      return(-1);
    }
    dataspace = H5Dget_space(dataset);
    datatype = H5Dget_type (dataset);
    Class = H5Tget_class(datatype);
    if(Class != H5T_COMPOUND) return(-2);
    status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,tbuf);
    rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
    for(int i=0;i<dims;i++){
      event->antennatrace[iant].Vfx[i] = tbuf[4*i+1];
      event->antennatrace[iant].Vfy[i] = tbuf[4*i+2];
      event->antennatrace[iant].Vfz[i] = tbuf[4*i+3];
    }
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Gclose(agroup);
  }
  free(tbuf);
  H5Gclose(group);
  return(0);
}
/**
 * read in the p2p (peak-to-peak) info table from an Aires event file
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_p2p(hid_t file, AIRESEvent *event)
{
  hid_t group;
  hid_t dataset,dataspace,datatype;
  int rank,Class,sz;
  herr_t status;
  hsize_t dims,max_dims;
  char *name;
  int iant;
  
  if((group = H5Gopen(file, event->runinfo.EventName, H5P_DEFAULT))<0) return(-1);
  if((dataset = H5Dopen(group, "AntennaP2PInfo", H5P_DEFAULT))<0) {
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->antennaP2Pbuffer == NULL){
    event->antennaP2Pbuffer = (char *)malloc(sz*dims);
    if(event->antennaP2Pbuffer == NULL)return(-3);
    event->antennap2p = (AntennaP2P *)malloc(sizeof(AntennaP2P)*dims);
    if(event->antennap2p ==NULL)return(-3);
  }
  event->n_ant = dims;
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,event->antennaP2Pbuffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"ID") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].ID =
        (char *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"P2P_efield") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].P2P_efield =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"P2P_voltage") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].P2P_voltage =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"P2P_filtered") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].P2P_filtered =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"HilbertPeakE") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].HilbertPeakE =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"HilbertPeakTimeE") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].HilbertPeakTimeE =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"HilbertPeakV") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].HilbertPeakV =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"HilbertPeakTimeV") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].HilbertPeakTimeV =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"HilbertPeakFV") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].HilbertPeakVf =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    if(strcmp(name,"HilbertPeakTimeFV") == 0){
      for(iant=0;iant<event->n_ant;iant++)event->antennap2p[iant].HilbertPeakTimeVf =
        (float *)&event->antennaP2Pbuffer[H5Tget_member_offset(datatype,i)+iant*sz];
    }
    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  return(0);
}

/**
 * readthe showertables group from an Aires event file. This contains the lateral and longitudinal profile of a simulated shower
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_showertables(hid_t file, AIRESEvent *event)
{
  hid_t group,groupa;
  hid_t dataset,dataspace,datatype;
  int rank,Class,sz;
  herr_t status;
  hsize_t dims,max_dims;
  char *name;
  int ib;
  
  if((group = H5Gopen(file, event->runinfo.EventName, H5P_DEFAULT))<0) return(-1);
  if((groupa = H5Gopen(group,"ShowerTables", H5P_DEFAULT))<0) return(-1);

  if((dataset = H5Dopen(groupa, "NLateralProfile", H5P_DEFAULT))<0) {
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->showertable.lateralbuffer == NULL){
    event->showertable.lateralbuffer = (char *)malloc(sz*dims);
    if(event->showertable.lateralbuffer == NULL)return(-3);
    event->showertable.lateralprofile = (LateralProfile *)malloc(sizeof(LateralProfile)*dims);
    if(event->showertable.lateralprofile ==NULL)return(-3);
  }
  event->showertable.n_lateral = dims;
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,
                   event->showertable.lateralbuffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"Distance") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Distance =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Ngamma") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Ngamma =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Ne_plus_minus") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Ne_plus_minus =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Ne_plus") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Ne_plus =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Nmu_plus_minus") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Nmu_plus_minus =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Nmu_plus") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Nmu_plus =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Nall_charged") == 0){
      for(ib=0;ib<event->showertable.n_lateral;ib++)
        event->showertable.lateralprofile[ib].Nall_charged =(float *)&event->showertable.
        lateralbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  if((dataset = H5Dopen(groupa, "NLongitudinalProfile", H5P_DEFAULT))<0) {
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Sget_simple_extent_dims (dataspace, &dims, &max_dims);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->showertable.longbuffer == NULL){
    event->showertable.longbuffer = (char *)malloc(sz*dims);
    if(event->showertable.longbuffer == NULL)return(-3);
    event->showertable.longprofile = (LongProfile *)malloc(sizeof(LongProfile)*dims);
    if(event->showertable.longprofile ==NULL)return(-3);
  }
  event->showertable.n_long = dims;
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,
                   event->showertable.longbuffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"SlantDepth") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].SlantDepth =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"VerticalDepth") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].VerticalDepth =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Ngamma") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Ngamma =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Ne_plus_minus") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Ne_plus_minus =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Ne_plus") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Ne_plus =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Nmu_plus_minus") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Nmu_plus_minus =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Nmu_plus") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Nmu_plus =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Npi_plus_minus") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Npi_plus_minus =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Npi_plus") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Npi_plus =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    if(strcmp(name,"Nall_charged") == 0){
      for(ib=0;ib<event->showertable.n_long;ib++)
        event->showertable.longprofile[ib].Nall_charged =(float *)&event->showertable.
        longbuffer[H5Tget_member_offset(datatype,i)+ib*sz];
    }
    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(groupa);
  H5Gclose(group);
  return(0);
}

/**
 * read in the signalsiminfo table from an Aires event file. This provides parameters used in the simulation
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_signalsiminfo(hid_t file, AIRESEvent *event)
{
  hid_t group,groupa;
  hid_t dataset,dataspace,datatype; //freq,phi,theta
  int rank,Class,sz;
  herr_t status;
  char *name;
  
  
  if((group = H5Gopen(file, event->runinfo.EventName, H5P_DEFAULT))<0) return(-1);
  if((groupa = H5Gopen(group,"ShowerTables", H5P_DEFAULT))<0) return(-1);

  //1. read Run info
  if((dataset = H5Dopen(group, "SignalSimInfo", H5P_DEFAULT))<0) return(-1);
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->signalsiminfo.size < sz){
    if(event->signalsiminfo.size > 0) free(event->signalsiminfo.buffer);
    event->signalsiminfo.buffer = (char *)malloc(sz);
    if(event->signalsiminfo.buffer ==NULL)return(-3);
    event->signalsiminfo.size = sz;
  }
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,event->signalsiminfo.buffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"FieldSimulator") == 0)
      event->signalsiminfo.FieldSimulator =
      (char *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"RefractionIndexModel") == 0)
      event->signalsiminfo.RefractionIndexModel =
      (char *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"RefractionIndexModelParameters") == 0)
      event->signalsiminfo.RefractionIndexModelParameters =
      (char *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"TimeBinSize") == 0)
      event->signalsiminfo.TimeBinSize =
      (double *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"TimeWindowMin") == 0)
      event->signalsiminfo.TimeWindowMin =
      (double *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"TimeWindowMax") == 0)
      event->signalsiminfo.TimeWindowMax =
      (double *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"OtherParameters") == 0)
      event->signalsiminfo.OtherParameters =
      (char *)&(event->signalsiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(groupa);
  H5Gclose(group);
  return(0);
}
/**
 * read in the eventinfo table from an Aires event file. General event information (angles, particle type, xmax, etc)
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_eventinfo(hid_t file, AIRESEvent *event)
{
  hid_t group;
  hid_t dataset,dataspace,datatype; //freq,phi,theta
  int rank,Class,sz;
  herr_t status;
  char *name;
  
  
  if((group = H5Gopen(file, event->runinfo.EventName, H5P_DEFAULT))<0) return(-1);
  if((dataset = H5Dopen(group, "EventInfo", H5P_DEFAULT))<0) {
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->eventinfo.size < sz){
    if(event->eventinfo.size > 0) free(event->eventinfo.buffer);
    event->eventinfo.buffer = (char *)malloc(sz);
    if(event->eventinfo.buffer ==NULL)return(-3);
    event->eventinfo.size = sz;
  }
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,event->eventinfo.buffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"EventName") == 0)
      event->eventinfo.EventName =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"EventID") == 0)
      event->eventinfo.EventID =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Primary") == 0)
      event->eventinfo.Primary =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Energy") == 0)
      event->eventinfo.Energy =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Zenith") == 0)
      event->eventinfo.Zenith =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Azimuth") == 0)
      event->eventinfo.Azimuth =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"XmaxDistance") == 0)
      event->eventinfo.XmaxDistance =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"XmaxPosition") == 0)
      event->eventinfo.XmaxPosition =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"XmaxAltitude") == 0)
      event->eventinfo.XmaxAltitude =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"SlantXmax") == 0)
      event->eventinfo.SlantXmax =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"InjectionAltitude") == 0)
      event->eventinfo.InjectionAltitude =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"GroundAltitude") == 0)
      event->eventinfo.GroundAltitude =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Site") == 0)
      event->eventinfo.Site =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Date") == 0)
      event->eventinfo.Date =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Latitude") == 0)
      event->eventinfo.Latitude =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"Longitude") == 0)
      event->eventinfo.Longitude =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"BField") == 0)
      event->eventinfo.BField =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"BFieldIncl") == 0)
      event->eventinfo.BFieldIncl =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"BFieldDecl") == 0)
      event->eventinfo.BFieldDecl =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"AtmosphericModel") == 0)
      event->eventinfo.AtmosphericModel =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"AtmosphericModelParameters") == 0)
      event->eventinfo.AtmosphericModelParameters =
      (char *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"EnergyInNeutrinos") == 0)
      event->eventinfo.EnergyInNeutrinos =
      (double *)&(event->eventinfo.buffer[H5Tget_member_offset(datatype,i)]);

    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  return(0);
}

/**
 * read in the showersiminfo table from an Aires event file. Parameters dictating the shower development
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */int read_showersiminfo(hid_t file, AIRESEvent *event)
{
  hid_t group;
  hid_t dataset,dataspace,datatype; //freq,phi,theta
  int rank,Class,sz;
  herr_t status;
  char *name;
  
  
  if((group = H5Gopen(file, event->runinfo.EventName, H5P_DEFAULT))<0) return(-1);
  if((dataset = H5Dopen(group, "ShowerSimInfo", H5P_DEFAULT))<0) {
    return(-1);
  }
  dataspace = H5Dget_space(dataset);
  datatype = H5Dget_type (dataset);
  Class = H5Tget_class(datatype);
  if(Class != H5T_COMPOUND) return(-2);
  rank = H5Tget_nmembers(datatype);
  sz = H5Tget_size(datatype);
  if(event->showersiminfo.size < sz){
    if(event->showersiminfo.size > 0) free(event->showersiminfo.buffer);
    event->showersiminfo.buffer = (char *)malloc(sz);
    if(event->showersiminfo.buffer ==NULL)return(-3);
    event->showersiminfo.size = sz;
  }
  status = H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,event->showersiminfo.buffer);
  for(int i=0;i<rank;i++){
    name =H5Tget_member_name(datatype,i);
    if(strcmp(name,"ShowerSimulator") == 0)
      event->showersiminfo.ShowerSimulator =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"HadronicModel") == 0)
      event->showersiminfo.HadronicModel =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"RandomSeed") == 0)
      event->showersiminfo.RandomSeed =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"RelativeThining") == 0) //watch spelling-error!
      event->showersiminfo.RelativeThinning =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"WeightFactor") == 0)
      event->showersiminfo.WeightFactor =
      (double *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"GammaEnergyCut") == 0)
      event->showersiminfo.GammaEnergyCut =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"ElectronEnergyCut") == 0)
      event->showersiminfo.ElectronEnergyCut =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"MuonEnergyCut") == 0)
      event->showersiminfo.MuonEnergyCut =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"MesonEnergyCut") == 0)
      event->showersiminfo.MesonEnergyCut =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"NucleonEnergyCut") == 0)
      event->showersiminfo.NucleonEnergyCut =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"CPUTime") == 0)
      event->showersiminfo.CPUTime =
      (double *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);
    if(strcmp(name,"OtherParameters") == 0)
      event->showersiminfo.OtherParameters =
      (char *)&(event->showersiminfo.buffer[H5Tget_member_offset(datatype,i)]);

    H5free_memory(name);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  return(0);
}

/**
 * read in the complete Aires event. This is the only read-in routine you need as it calls the othes
 * @param[in] file the hdf5 file
 * @param[in] event the AIRES event in which the info is stored
 */
int read_aires_event(char *fname, AIRESEvent *event)
{
  hid_t file;
  int rcode;
 
  
  if((file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT))<0) return(-1);
  if((rcode = read_runinfo(file,event))< 0) return(rcode);
  if((rcode = read_eventinfo(file,event))< 0) return(rcode);
  if((rcode = read_showersiminfo(file,event))<0) return(rcode);
  if((rcode = read_signalsiminfo(file,event))<0) return(rcode);
  if((rcode = read_antennainfo(file,event))<0) return(rcode);
  if((rcode = read_antennatraces(file,event))<0) return(rcode);
  if((rcode = read_p2p(file,event))<0) return(rcode);
  if((rcode = read_showertables(file,event))<0) return(rcode);

  //printf("Read all info!\n");
  H5Fclose(file);
  return(0);
}


