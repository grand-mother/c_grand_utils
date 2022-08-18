/*! \file GRANDevent.h */


typedef struct{
  long Second;          //!< GPS seconds (hardware or simulation date)
  unsigned int NanoSec; //!< GPS nanoseconds (hardware or simulation)
}GPS;

typedef struct{
  char EventName[80]; //!< event name (simulation)
  char EventID[80]; //!< event ID simulation/raw data
  char Site[80];  //!< site for which the simulation was done/ at which data was taken
  char Date[80];  //!< date for which the simulation was done/at which data was taken
  double Latitude;  //!< latitude of site (degrees)
  double Longitude; //!< longitude of site (degrees)
  double GroundAltitude;  //!< altitude of site (m)
  double BField;  //!< magnitude of B-field (uT)
  double BFieldIncl;  //!< inclination B field (degrees)
  double BFieldDecl;  //!< declination B-field (degrees)
  char AtmosphericModel[80]; //!< atmospheric model used in simulation/reconstruction
  char AtmosphericModelParameters[100]; //!< optional: parameters of the atmospheric model
}GRANDGenEventInfo;

typedef struct{
  char Primary[80]; //!< primary particle in simulation
  double Energy;  //!< energy of primary particle (EeV)
  double Zenith; //!< zenith angle incoming particle (degrees, 0 = coming from above)
  double Azimuth; //!< azimuth angle primary particle (degrees, 0=coming from North)
  double XmaxDistance; //!< distance between Xmax and the core position
  double XmaxPosition[3]; //!<The Xmax coordinates (m) (x,y,z)
  double XmaxAltitude; //!< the altitude of Xmax (m) above the  Earth
  double SlantXmax; //!< amount of atmosphere encountered by the shower/primary particle  until Xmax (g/cm^2)
  double InjectionAltitude; //!< altitude above Earth (m) at which primary particle was injected
  double EnergyInNeutrinos; //!< amount of energy in Neutrinos (EeV)
}GRANDSimEventInfo;

typedef struct{
  double Energy, e_Energy; //!<Shower energy and uncertainty (EeV)
  double Zenith,e_Zenith; //!<Shower zenith angle and uncertainty (degrees, 0=from above)
  double Azimuth,e_Azimuth; //!<Shower azimuth and uncertainty (degrees, 0=from North)
  double Core[3],e_Core[3];//!<Shower core (x,y,z) at surface; defined by Z-coordinate (m)
  double XmaxDistance,e_XmaxDistance; //!<distance between core and Xmax (m)
  double SlantXmax,e_SlantXmax; //!<amount of atmosphere encountered by the shower/primary particle  until Xmax (g/cm^2)
  GPS CoreTime; //!< time when the shower was at the core position
}GRANDRecEventInfo;

typedef struct{
  char ID[30];  //!< identifier of the detector unit
  char DetectorModel[80]; //!< description of the hardware (eg antenna type, particle detector type, ..)
  char ElectronicsModel[80]; //!< description of the electronics used
  double Position[3];//!< (x,y,z) of the detector unit in the field/simulation
}GRANDDetectorInfo;

typedef struct{
  int SimPoints; //!< trace length in simulation
  int RawPoints; //!< trace length in raw data
  GPS SimGPS; //!< simulated time of the start of the detector data
  GPS RawGPS; //!< raw time of the start of the detector data
  float TPulse; //!< offset between start and the peak of the signal
  double DetTime; //!< time of the peak of the signal relative to the time that the shower hit the surface
  double e_DetTime; //!< uncertainty on this peak time
  bool IsTriggered; //!< flag if this unit was triggered
  float SimMSPS;  //!< simulated sampling speed (Mega Samples per Second)
  float RawMSPS; //!< raw sampling speed (Mega Samples per Second)
  float *SimEfield[3]; //!< simulated time traces of the electric field (x,y,z) in uV/m
  float *SimE_fftMag[3]; //!<magnitude of the  FFT of the simulated E-field
  float *SimE_fftPhase[3]; //!< phase of the FFT of the simulated E-field (radians)
  float *SimVoltage[3]; //!< Voltage after applying the antenna model/electronics to the E-field (uV)
  float *RawADC[3]; //!< time traces (x,y,z) raw ADC values (from electronics or digitizing and resampling the simulation)
  float *RawVoltage[3]; //!< Voltage time traces after applying electronics description to ADC (uV)
  float *RecEfield[3]; //!< time traces after applying inverting antenna description to voltage (x,y,z) (uV/m)
  float *RecE_fftMag[3]; //!< magnitude of the FFT of the reconstructed E field (x,y,z)
  float *RecE_fftPhase[3]; //!< phase of the FFT of the reconstructed E field (x,y,z) in radians
}GRANDDetectorData;

typedef struct{
  int n_det; //!< number of detectors in the event
  GRANDGenEventInfo GenEventInfo; //!< general information (site, atmosphere)
  GRANDSimEventInfo SimEventInfo; //!< simulation specific info
  GRANDRecEventInfo RecEventInfo; //!< reconstruction specific info
  GRANDDetectorInfo *DetectorInfo; //!< detector information (type, location)
  float *TraceBuffer; //!< buffer holding the data (length not a-priori known)
  GRANDDetectorData *DetectorData; //!< all traces and related detector readout parameters
}GRANDEvent;

