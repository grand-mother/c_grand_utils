/*! \file io_util.c
In this file we should have input and output routines for GRAND events. Right now, only an output routine to hdf5 is written and an input from ROOT format. The development of this file depends on the agreed file-format in GRAND.
*/
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TClass.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "GRANDevent.h"

using namespace std;

#define PI 3.141592653589793

typedef struct{
  unsigned int run_number;
  unsigned int  run_mode;
  unsigned int first_event;
  unsigned int  first_event_time;
  unsigned int last_event;
  unsigned int  last_event_time;
  char data_source[80];
  char data_generator[80];
  char data_generator_version[80];
  char site[80];
  float site_long;
  float site_lat;
  float origin_geoid[3];
}runhdr;

typedef struct{
  unsigned int run_number;
  char refractivity_model[80];
  double refrac_param[2];
  float t_pre;
  float t_post;
  float t_bin_size;

}runsim;
runhdr run_hdr;
runsim run_simdata;


void grand_root_read_runhdr(TFile *g)
{
  memset(&run_hdr,0,sizeof(runhdr));
  TTree *trun = (TTree *)g->Get("trun");
  trun->SetBranchAddress("run_number",&run_hdr.run_number);
  trun->SetBranchAddress("run_mode",&run_hdr.run_mode);
  trun->SetBranchAddress("first_event",&run_hdr.first_event);
  trun->SetBranchAddress("first_event_time",&run_hdr.first_event_time);
  trun->SetBranchAddress("last_event",&run_hdr.last_event);
  trun->SetBranchAddress("last_event_time",&run_hdr.last_event_time);
  string *data_source = new string();
  trun->SetBranchAddress("data_source",&data_source);
  string *data_generator = new string();
  trun->SetBranchAddress("data_generator",&data_generator);
  string *data_generator_version = new string();
  trun->SetBranchAddress("data_generator_version",&data_generator_version);
  string *site = new string();
  trun->SetBranchAddress("site",&site);
  trun->SetBranchAddress("site_long",&run_hdr.site_long);
  trun->SetBranchAddress("site_lat",&run_hdr.site_lat);
  trun->SetBranchAddress("origin_geoid",run_hdr.origin_geoid);
  for(int i=0;i< trun->GetEntries();i++){
    if(trun->GetEntry(i)<= 0) break;
    strcpy(run_hdr.data_source,data_source->c_str());
    strcpy(run_hdr.data_generator,data_generator->c_str());
    strcpy(run_hdr.data_generator_version,data_generator_version->c_str());
    strcpy(run_hdr.site,site->c_str());
  }
  trun->Delete();
}

void grand_root_read_runsimdata(TFile *g)
{
  memset(&run_simdata,0,sizeof(runsim));
  TTree *trunefieldsimdata = (TTree *)g->Get("trunefieldsimdata");
  trunefieldsimdata->SetBranchAddress("run_number",&run_simdata.run_number);
  string *model = new string();
  trunefieldsimdata->SetBranchAddress("refractivity_model",&model);
  vector<double> *model_param = new vector<double>(2);
  trunefieldsimdata->SetBranchAddress("refractivity_model_parameters",&model_param);
  trunefieldsimdata->SetBranchAddress("t_pre",&run_simdata.t_pre);
  trunefieldsimdata->SetBranchAddress("t_post",&run_simdata.t_post);
  trunefieldsimdata->SetBranchAddress("t_bin_size",&run_simdata.t_bin_size);
  for(int i=0;i< trunefieldsimdata->GetEntries();i++){
    if(trunefieldsimdata->GetEntry(i)<=0) break;
    strcpy(run_simdata.refractivity_model,model->c_str());
    run_simdata.refrac_param[0] = *((double *)model_param->data());
    run_simdata.refrac_param[1] = *((double *)(1+model_param->data()));
  }
  trunefieldsimdata->Delete();
}

int grand_root_read_event(TFile *g,int ievt, GRANDEvent *event)
{
  unsigned int n_det,simtrlen,rawtrlen;
  int ibuf;

  simtrlen = ((run_simdata.t_post-run_simdata.t_pre)/run_simdata.t_bin_size)-1;
  rawtrlen =(simtrlen*run_simdata.t_bin_size+1)/2; // to be changed 2 = 2 ns or 500 MSPS

  TTree *teventefield = (TTree *)g->Get("teventefield");
  teventefield->SetBranchAddress("du_count",&n_det);
  if(ievt>=teventefield->GetEntries()) return(-1);
  if(teventefield->GetEntry(ievt)<=0) return(-1);
  printf("DU count %d\n",n_det);
  if(n_det != event->n_det && event->DetectorData != NULL){
    free(event->DetectorData);
    event->DetectorData = NULL;
    free(event->DetectorInfo);
    event->DetectorInfo = NULL;
    free(event->TraceBuffer);
    event->TraceBuffer = NULL;
 }
  if(event->DetectorData == NULL) {
    event->DetectorData = (GRANDDetectorData *)malloc(n_det*sizeof(GRANDDetectorData));
  }
  if(event->DetectorInfo == NULL) {
    event->DetectorInfo = (GRANDDetectorInfo *)malloc(n_det*sizeof(GRANDDetectorInfo));
  }
  if(event->TraceBuffer != NULL) free((void *)event->TraceBuffer);
  event->TraceBuffer = (float *)malloc(n_det*(simtrlen*12+rawtrlen*15)*sizeof(float));
  event->n_det = n_det;
  teventefield->SetBranchAddress("run_number",&event->GenEventInfo.RunNumber);
  teventefield->SetBranchAddress("event_number",&event->GenEventInfo.EventNumber);
  teventefield->SetBranchAddress("event_type",&event->GenEventInfo.EventType);
  teventefield->SetBranchAddress("time_seconds",&event->GenEventInfo.EventTime.Second);
  teventefield->SetBranchAddress("time_nanoseconds",&event->GenEventInfo.EventTime.NanoSec);
  vector<unsigned short> *du_id = new vector<unsigned short>(n_det);
  teventefield->SetBranchAddress("du_id",&du_id);
  vector<unsigned int> *du_seconds = new vector<unsigned int>(n_det);
  teventefield->SetBranchAddress("du_seconds",&du_seconds);
  vector<unsigned int> *du_nanoseconds = new vector<unsigned int>(n_det);
  teventefield->SetBranchAddress("du_nanoseconds",&du_nanoseconds);
  vector<unsigned short> *trigger_position = new vector<unsigned short>(n_det);
  teventefield->SetBranchAddress("trigger_position",&trigger_position);
  vector<unsigned short> *trigger_flag = new vector<unsigned short>(n_det);
  teventefield->SetBranchAddress("trigger_flag",&trigger_flag);
  vector<unsigned short> *trigger_pattern = new vector<unsigned short>(n_det);
  teventefield->SetBranchAddress("trigger_pattern",&trigger_pattern);
  vector<unsigned short> *trigger_rate = new vector<unsigned short>(n_det);
  teventefield->SetBranchAddress("trigger_rate",&trigger_rate);
  vector<float> *atm_pressure = new vector<float>(n_det);
  teventefield->SetBranchAddress("atm_pressure",&atm_pressure);
  vector<float> *atm_temperature = new vector<float>(n_det);
  teventefield->SetBranchAddress("atm_temperature",&atm_temperature);
  vector<float> *atm_humidity = new vector<float>(n_det);
  teventefield->SetBranchAddress("atm_humidity",&atm_humidity);
  vector<float> *gps_long = new vector<float>(n_det);
  teventefield->SetBranchAddress("gps_long",&gps_long);
  vector<float> *gps_lat = new vector<float>(n_det);
  teventefield->SetBranchAddress("gps_lat",&gps_lat);
  vector<float> *gps_alt = new vector<float>(n_det);
  teventefield->SetBranchAddress("gps_alt",&gps_alt);
  vector<float> *pos_x = new vector<float>(n_det);
  teventefield->SetBranchAddress("pos_x",&pos_x);
  vector<float> *pos_y = new vector<float>(n_det);
  teventefield->SetBranchAddress("pos_y",&pos_y);
  vector<float> *pos_z = new vector<float>(n_det);
  teventefield->SetBranchAddress("pos_z",&pos_z);
  vector<vector<float>> tracex(n_det,vector<float>(simtrlen));
  vector<vector<float>>*trace_x = &tracex;
  teventefield->SetBranchAddress("trace_x",&trace_x);
  vector<vector<float>> tracey(n_det,vector<float>(simtrlen));
  vector<vector<float>>*trace_y = &tracey;
  teventefield->SetBranchAddress("trace_y",&trace_y);
  vector<vector<float>> tracez(n_det,vector<float>(simtrlen));
  vector<vector<float>>*trace_z = &tracez;
  teventefield->SetBranchAddress("trace_z",&trace_z);
  // check on digi_prepost!
  if(teventefield->GetEntry(ievt)<=0) return(-1);
  ibuf = 0;
  for(int i=0;i<n_det;i++){
    sprintf(event->DetectorInfo[i].ID,"DU%d",du_id->at(i));
    event->DetectorInfo[i].Position[0] = pos_x->at(i);
    event->DetectorInfo[i].Position[1] = pos_y->at(i);
    event->DetectorInfo[i].Position[2] = pos_z->at(i);
    //GPS positions not in memory yet!
    //Atmosphere to be done
    event->DetectorData[i].SimGPS.Second = du_seconds->at(i);
    event->DetectorData[i].SimGPS.NanoSec = du_nanoseconds->at(i);
    // triggerinfo te be put in memory
    event->DetectorData[i].SimPoints = simtrlen;
    event->DetectorData[i].RawPoints = rawtrlen;
    event->DetectorData[i].SimMSPS = 1000/run_simdata.t_bin_size;
    event->DetectorData[i].RawMSPS = 500; // hardcoded 500 MSPS
    for(int iarm = 0;iarm<3;iarm++){
      event->DetectorData[i].SimEfield[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=simtrlen;
      event->DetectorData[i].SimE_fftMag[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=simtrlen;
      event->DetectorData[i].SimE_fftPhase[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=simtrlen;
      event->DetectorData[i].SimVoltage[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=simtrlen;
      event->DetectorData[i].RawADC[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=rawtrlen;
      event->DetectorData[i].RawVoltage[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=rawtrlen;
      event->DetectorData[i].RecEfield[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=rawtrlen;
      event->DetectorData[i].RecE_fftMag[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=rawtrlen;
      event->DetectorData[i].RecE_fftPhase[iarm] = &event->TraceBuffer[ibuf];
      ibuf+=rawtrlen;
    }
    for(int id=0;id<simtrlen;id++){
      event->DetectorData[i].SimEfield[0][id] = (trace_x->at(i)).at(id);
      event->DetectorData[i].SimEfield[1][id] = (trace_y->at(i)).at(id);
      event->DetectorData[i].SimEfield[2][id] = (trace_z->at(i)).at(id);
    }
  }
  TTree *teventshowersimdata = (TTree *)g->Get("teventshowersimdata");
  vector<float> *prim_energy = new vector<float>(1);
  teventshowersimdata->SetBranchAddress("prim_energy",&prim_energy);
  teventshowersimdata->SetBranchAddress("shower_azimuth",&event->SimEventInfo.Azimuth);
  teventshowersimdata->SetBranchAddress("shower_zenith",&event->SimEventInfo.Zenith);
  teventshowersimdata->SetBranchAddress("xmax_distance",&event->SimEventInfo.XmaxDistance);
  teventshowersimdata->SetBranchAddress("xmax_grams",&event->SimEventInfo.SlantXmax);
  teventshowersimdata->SetBranchAddress("xmax_pos_shc",event->SimEventInfo.XmaxPosition);
  teventshowersimdata->SetBranchAddress("xmax_alt",&event->SimEventInfo.XmaxAltitude);
  teventshowersimdata->SetBranchAddress("magnetic_field",event->GenEventInfo.BCart);
  if(teventshowersimdata->GetEntry(ievt)<=0) return(-1);
  event->SimEventInfo.Energy = prim_energy->at(0);
  float *bf  = event->GenEventInfo.BCart;
  event->GenEventInfo.BField = sqrt(bf[0]*bf[0]+bf[1]*bf[1]+bf[2]*bf[2]);
  event->GenEventInfo.BFieldDecl = acos(bf[2]/event->GenEventInfo.BField)*180/PI;
  event->GenEventInfo.BFieldIncl = atan2(bf[1],bf[0])*180/PI;
  event->SimEventInfo.Azimuth  =  180 + event->SimEventInfo.Azimuth; //to  or from....
  event->SimEventInfo.Zenith  =  180 - event->SimEventInfo.Zenith;
  return(1);
}



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
