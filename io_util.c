/*! \file io_util.c
In this file we should have input and output routines for GRAND events. Right now, only an output routine to hdf5 is written. The development of this file depends on the agreed file-format in GRAND.
*/
#include <stdlib.h>
#include <string.h>
#include "GRANDevent.h"
#include "hdf5.h"
/**
 * writing a GRAND event into a hdf5 file. The name of the hdf5 file is determined from the event name.
 * @param[in] event GRAND event
 */
int write_GRANDevent(GRANDEvent *event)
{
  hid_t file_id,evt_id,tr_id;
  hid_t gentype,genspace,dataset,strtype;
  hid_t simtype,simspace;
  hid_t rectype,recspace;
  hid_t dettype,detspace;
  hid_t gpstype,dattype;
  hid_t datspace,datchunk;
  hid_t postype,tracetype;
  int offset;
  hsize_t dim[3];
  char name[80],*buffer;
  float *Efield,*E_fftMag,*E_fftPhase,*Volt,*ADC;
  
  sprintf(name,"%s.hdf5",event->GenEventInfo.EventName);
  if((file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT))<0){
    return(-2);
  }
  if((evt_id = H5Gcreate(file_id,event->GenEventInfo.EventName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT))<0){
    H5Fclose(file_id);
    return(-2);
  }
  dim[0] = 3;
  postype = H5Tarray_create(H5T_NATIVE_DOUBLE,1,dim);//3-d array
  strtype =H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype,80);
  gentype = H5Tcreate (H5T_COMPOUND, sizeof(GRANDGenEventInfo));
  H5Tinsert(gentype,"EventName",HOFFSET(GRANDGenEventInfo,EventName),strtype);
  H5Tinsert(gentype,"EventID",HOFFSET(GRANDGenEventInfo,EventID),strtype);
  H5Tinsert(gentype,"Site",HOFFSET(GRANDGenEventInfo,Site),strtype);
  H5Tinsert(gentype,"Date",HOFFSET(GRANDGenEventInfo,Date),strtype);
  H5Tinsert(gentype,"Latitude",HOFFSET(GRANDGenEventInfo,Latitude),H5T_NATIVE_DOUBLE);
  H5Tinsert(gentype,"Longitude",HOFFSET(GRANDGenEventInfo,Longitude),H5T_NATIVE_DOUBLE);
  H5Tinsert(gentype,"GroundAltitude",HOFFSET(GRANDGenEventInfo,GroundAltitude),H5T_NATIVE_DOUBLE);
  H5Tinsert(gentype,"BField",HOFFSET(GRANDGenEventInfo,BField),H5T_NATIVE_DOUBLE);
  H5Tinsert(gentype,"BFieldIncl",HOFFSET(GRANDGenEventInfo,BFieldIncl),H5T_NATIVE_DOUBLE);
  H5Tinsert(gentype,"BFieldDecl",HOFFSET(GRANDGenEventInfo,BFieldDecl),H5T_NATIVE_DOUBLE);
  H5Tinsert(gentype,"AtmosphericModel",HOFFSET(GRANDGenEventInfo,AtmosphericModel),strtype);
  H5Tclose(strtype);
  strtype =H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype,100);
  H5Tinsert(gentype,"AtmosphericModelParameters",HOFFSET(GRANDGenEventInfo,AtmosphericModelParameters),strtype);
  H5Tclose(strtype);
  dim[0] = 1;
  genspace = H5Screate_simple(1,dim,NULL);
  dataset = H5Dcreate(evt_id,"GeneralInfo",gentype,genspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, gentype, H5S_ALL, H5S_ALL,H5P_DEFAULT,&event->GenEventInfo );
  H5Dclose(dataset);
  H5Sclose(genspace);
  H5Tclose(gentype);
  simtype = H5Tcreate (H5T_COMPOUND, sizeof(GRANDSimEventInfo));
  strtype =H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype,80);
  H5Tinsert(simtype,"Primary",HOFFSET(GRANDSimEventInfo,Primary),strtype);
  H5Tclose(strtype);
  H5Tinsert(simtype,"Energy",HOFFSET(GRANDSimEventInfo,Energy),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"Zenith",HOFFSET(GRANDSimEventInfo,Zenith),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"Azimuth",HOFFSET(GRANDSimEventInfo,Azimuth),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"XmaxDistance",HOFFSET(GRANDSimEventInfo,XmaxDistance),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"XmaxPosition",HOFFSET(GRANDSimEventInfo,XmaxPosition),postype);
  H5Tinsert(simtype,"XmaxAltitude",HOFFSET(GRANDSimEventInfo,XmaxAltitude),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"SlantXmax",HOFFSET(GRANDSimEventInfo,SlantXmax),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"InjectionAltitude",HOFFSET(GRANDSimEventInfo,InjectionAltitude),H5T_NATIVE_DOUBLE);
  H5Tinsert(simtype,"EnergyInNeutrinos",HOFFSET(GRANDSimEventInfo,EnergyInNeutrinos),H5T_NATIVE_DOUBLE);
  dim[0] = 1;
  simspace = H5Screate_simple(1,dim,NULL);
  dataset = H5Dcreate(evt_id,"SimulationInfo",simtype,simspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, simtype, H5S_ALL, H5S_ALL,H5P_DEFAULT,&event->SimEventInfo );
  H5Dclose(dataset);
  H5Sclose(simspace);
  H5Tclose(simtype);
  gpstype =H5Tcreate(H5T_COMPOUND, sizeof(GPS));
  H5Tinsert(gpstype,"Second",HOFFSET(GPS,Second),H5T_NATIVE_LONG);
  H5Tinsert(gpstype,"NanoSec",HOFFSET(GPS,NanoSec),H5T_NATIVE_UINT);
  rectype = H5Tcreate (H5T_COMPOUND, sizeof(GRANDRecEventInfo));
  H5Tinsert(rectype,"Energy",HOFFSET(GRANDRecEventInfo,Energy),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"e_Energy",HOFFSET(GRANDRecEventInfo,e_Energy),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"Zenith",HOFFSET(GRANDRecEventInfo,Zenith),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"e_Zenith",HOFFSET(GRANDRecEventInfo,e_Zenith),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"Azimuth",HOFFSET(GRANDRecEventInfo,Azimuth),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"e_Azimuth",HOFFSET(GRANDRecEventInfo,e_Azimuth),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"Core",HOFFSET(GRANDRecEventInfo,Core),postype);
  H5Tinsert(rectype,"e_Core",HOFFSET(GRANDRecEventInfo,e_Core),postype);
  H5Tinsert(rectype,"XmaxDistance",HOFFSET(GRANDRecEventInfo,XmaxDistance),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"e_XmaxDistance",HOFFSET(GRANDRecEventInfo,e_XmaxDistance),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"SlantXmax",HOFFSET(GRANDRecEventInfo,SlantXmax),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"e_SlantXmax",HOFFSET(GRANDRecEventInfo,e_SlantXmax),H5T_NATIVE_DOUBLE);
  H5Tinsert(rectype,"CoreTime",HOFFSET(GRANDRecEventInfo,CoreTime),gpstype);
  dim[0] = 1;
  recspace = H5Screate_simple(1,dim,NULL);
  dataset = H5Dcreate(evt_id,"ReconstructionInfo",rectype,recspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, rectype, H5S_ALL, H5S_ALL,H5P_DEFAULT,&event->RecEventInfo );
  H5Dclose(dataset);
  H5Sclose(recspace);
  H5Tclose(rectype);
  dettype = H5Tcreate (H5T_COMPOUND, sizeof(GRANDDetectorInfo));
  strtype =H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype,30);
  H5Tinsert(dettype,"ID",HOFFSET(GRANDDetectorInfo,ID),strtype);
  H5Tclose(strtype);
  strtype =H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype,80);
  H5Tinsert(dettype,"DetectorModel",HOFFSET(GRANDDetectorInfo,DetectorModel),strtype);
  H5Tinsert(dettype,"ElectronicsModel",HOFFSET(GRANDDetectorInfo,ElectronicsModel),strtype);
  H5Tclose(strtype);
  H5Tinsert(dettype,"Position",HOFFSET(GRANDDetectorInfo,Position),postype);
  dim[0] = event->n_det;
  detspace = H5Screate_simple(1,dim,NULL);
  dataset = H5Dcreate(evt_id,"DetectorInfo",dettype,detspace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, dettype, H5S_ALL, H5S_ALL,H5P_DEFAULT,event->DetectorInfo );
  H5Dclose(dataset);
  H5Sclose(detspace);
  H5Tclose(dettype);
  tr_id = H5Gcreate(evt_id,"Traces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(int iant=0;iant<event->n_det;iant++){
    dattype = H5Tcreate (H5T_COMPOUND, 2*sizeof(GPS)+3*sizeof(float)+2*sizeof(double)+sizeof(bool)+
                         12*event->DetectorData[iant].SimPoints*sizeof(float)+
                         15*event->DetectorData[iant].RawPoints*sizeof(float));
    offset = 0;
    H5Tinsert(dattype,"SimGPS",offset,gpstype);
    offset += sizeof(GPS);
    H5Tinsert(dattype,"SimMSPS",offset,H5T_NATIVE_FLOAT);
    offset+=sizeof(float);
    dim[0] = 3;
    dim[1] = event->DetectorData[iant].SimPoints;
    tracetype = H5Tarray_create(H5T_NATIVE_FLOAT,2,dim);//3-d array
    H5Tinsert(dattype,"SimEfield",offset,tracetype);
    offset +=3*event->DetectorData[iant].SimPoints*sizeof(float);
    H5Tinsert(dattype,"SimE_fftMag",offset,tracetype);
    offset +=3*event->DetectorData[iant].SimPoints*sizeof(float);
    H5Tinsert(dattype,"SimE_fftPhase",offset,tracetype);
    offset +=3*event->DetectorData[iant].SimPoints*sizeof(float);
    H5Tinsert(dattype,"SimVoltage",offset,tracetype);
    offset +=3*event->DetectorData[iant].SimPoints*sizeof(float);
    H5Tclose(tracetype);
    H5Tinsert(dattype,"RawGPS",offset,gpstype);
    offset += sizeof(GPS);
    H5Tinsert(dattype,"RawMSPS",offset,H5T_NATIVE_FLOAT);
    offset+=sizeof(float);
    H5Tinsert(dattype,"TPulse",offset,H5T_NATIVE_FLOAT);
    offset+=sizeof(float);
    H5Tinsert(dattype,"DetTime",offset,H5T_NATIVE_DOUBLE);
    offset+=sizeof(double);
    H5Tinsert(dattype,"e_DetTime",offset,H5T_NATIVE_DOUBLE);
    offset+=sizeof(double);
    H5Tinsert(dattype,"IsTriggered",offset,H5T_NATIVE_HBOOL);
    offset+=sizeof(bool);
    dim[0] = 3;
    dim[1] = event->DetectorData[iant].RawPoints;
    tracetype = H5Tarray_create(H5T_NATIVE_FLOAT,2,dim);//3-d array
    H5Tinsert(dattype,"RawADC",offset,tracetype);
    offset +=3*event->DetectorData[iant].RawPoints*sizeof(float);
    H5Tinsert(dattype,"RawVoltage",offset,tracetype);
    offset +=3*event->DetectorData[iant].RawPoints*sizeof(float);
    H5Tinsert(dattype,"RecEfield",offset,tracetype);
    offset +=3*event->DetectorData[iant].RawPoints*sizeof(float);
    H5Tinsert(dattype,"RecE_fftMag",offset,tracetype);
    offset +=3*event->DetectorData[iant].RawPoints*sizeof(float);
    H5Tinsert(dattype,"RecE_fftPhase",offset,tracetype);
    offset +=3*event->DetectorData[iant].RawPoints*sizeof(float);
    H5Tclose(tracetype);
    buffer = (char *)malloc(offset);
    dim[0] =1;
    datspace =  H5Screate_simple(1,dim,NULL);
    datchunk = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(datchunk, H5D_CHUNKED);
    H5Pset_chunk(datchunk, 1,dim);
    H5Pset_deflate(datchunk,6); //compression level 6
    dataset = H5Dcreate(tr_id,event->DetectorInfo[iant].ID,dattype,datspace,
                        H5P_DEFAULT, datchunk, H5P_DEFAULT);
    offset = 0;
    memcpy(&buffer[offset],&(event->DetectorData[iant].SimGPS),sizeof(GPS));
    offset += sizeof(GPS);
    memcpy(&buffer[offset],&(event->DetectorData[iant].SimMSPS),sizeof(float));
    offset += sizeof(float);
    Efield =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].SimPoints*sizeof(float);
    E_fftMag =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].SimPoints*sizeof(float);
    E_fftPhase =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].SimPoints*sizeof(float);
    Volt =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].SimPoints*sizeof(float);
    for(int iarm=0;iarm<3;iarm++){
      for(int i=0;i<event->DetectorData[iant].SimPoints;i++){
        Efield[iarm*event->DetectorData[iant].SimPoints+i]=
        event->DetectorData[iant].SimEfield[iarm][i];
        E_fftMag[iarm*event->DetectorData[iant].SimPoints+i]=
        event->DetectorData[iant].SimE_fftMag[iarm][i];
        E_fftPhase[iarm*event->DetectorData[iant].SimPoints+i]=
        event->DetectorData[iant].SimE_fftPhase[iarm][i];
        Volt[iarm*event->DetectorData[iant].SimPoints+i]=
        event->DetectorData[iant].SimVoltage[iarm][i];
      }
    }
    memcpy(&buffer[offset],&(event->DetectorData[iant].RawGPS),sizeof(GPS));
    offset += sizeof(GPS);
    memcpy(&buffer[offset],&(event->DetectorData[iant].RawMSPS),sizeof(float));
    offset += sizeof(float);
    memcpy(&buffer[offset],&(event->DetectorData[iant].TPulse),sizeof(float));
    offset += sizeof(float);
    memcpy(&buffer[offset],&(event->DetectorData[iant].DetTime),sizeof(double));
    offset += sizeof(double);
    memcpy(&buffer[offset],&(event->DetectorData[iant].e_DetTime),sizeof(double));
    offset += sizeof(double);
    memcpy(&buffer[offset],&(event->DetectorData[iant].IsTriggered),sizeof(bool));
    offset += sizeof(bool);
    ADC =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].RawPoints*sizeof(float);
    Volt =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].RawPoints*sizeof(float);
    Efield =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].RawPoints*sizeof(float);
    E_fftMag =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].RawPoints*sizeof(float);
    E_fftPhase =(float *)&buffer[offset];
    offset += 3*event->DetectorData[iant].RawPoints*sizeof(float);
    for(int iarm=0;iarm<3;iarm++){
      for(int i=0;i<event->DetectorData[iant].RawPoints;i++){
        ADC[iarm*event->DetectorData[iant].RawPoints+i]=
        event->DetectorData[iant].RawADC[iarm][i];
        Volt[iarm*event->DetectorData[iant].RawPoints+i]=
        event->DetectorData[iant].RawVoltage[iarm][i];
        Efield[iarm*event->DetectorData[iant].RawPoints+i]=
        event->DetectorData[iant].RecEfield[iarm][i];
        E_fftMag[iarm*event->DetectorData[iant].RawPoints+i]=
        event->DetectorData[iant].RecE_fftMag[iarm][i];
        E_fftPhase[iarm*event->DetectorData[iant].RawPoints+i]=
        event->DetectorData[iant].RecE_fftPhase[iarm][i];
      }
    }
    H5Dwrite(dataset,dattype,H5S_ALL, H5S_ALL,H5P_DEFAULT,buffer);
    H5Dclose(dataset);
    free(buffer);
    H5Pclose(datchunk);
    H5Tclose(dattype);
    H5Sclose(datspace);
  }
  H5Tclose(gpstype);
  H5Tclose(postype);
  H5Gclose(tr_id);
  H5Gclose(evt_id);
  H5Fclose(file_id);
  return(0);
}
