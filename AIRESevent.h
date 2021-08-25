/*! \file AIRESevent.h */
typedef struct{
  int size; //!< size of the RunInfo Table
  char *buffer; //!< Buffer in which the data is stored
  char *EventName; //!< Name of the Aires event (same as in EventInfo)
  char *EventID; //!< Identifier of the Aires event (same as in EventInfo)
  char *Primary; //!< Primary particle used in simulation  (same as in EventInfo)
  double *Energy; //!< Energy of the primary particle (EeV)  (same as in EventInfo)
  double *Zenith; //!< Zenith angle (degrees) in AIRES convention (in the direction the particle is moving)  (same as in EventInfo)
  double *Azimuth;//!< Azimuth angle (degrees) in AIRES convention (in the direction the particle is moving)  (same as in EventInfo)
  double *XmaxDistance; //!< Distance between the core position and Xmax (m)  (same as in EventInfo)
  double *SlantXmax; //!< amout of atmosphere traversed by incoming particle and shower until Xmax (g/cm^2)  (same as in EventInfo)
  char *HadronicModel; //!< hadronic interaction model used  (same as in ShowersimInfo)
  double *InjectionAltitude; //!< Altitude above sea level of the injection of the primary particle  (same as in EventInfo)
}Runinfo;

typedef struct{
  int size; //!< size of the Eventinfo table
  char *buffer; //!< buffer in which the data is stored
  char *EventName; //!< Name of the Aires event (same as in RunInfo)
  char *EventID; //!< identifier of the Aires event (same as in Runinfo)
  char *Primary; //!< Primary particle used in simulation  (same as in RunInfo)
  double *Energy;//!< Energy of the primary particle (EeV)  (same as in RunInfo)
  double *Zenith;//!< Zenith angle (degrees) in AIRES convention (in the direction the particle is moving)  (same as in RunInfo)
  double *Azimuth;//!< Azimuth angle (degrees) in AIRES convention (in the direction the particle is moving)  (same as in RunInfo)
  double *XmaxDistance; //!< Distance between the core position and Xmax (m)  (same as in RunInfo)
  double *XmaxPosition; //!< Coordinates of the Xmax position (x,y,z)
  double *XmaxAltitude;//!< Altitude of Xmax above sea level (m)
  double *SlantXmax; //!< amout of atmosphere traversed by incoming particle and shower until Xmax (g/cm^2)  (same as in RunInfo)
  double *InjectionAltitude; //!< Altitude above sea level of the injection of the primary particle  (same as in RunInfo)
  double *GroundAltitude; //!< Altitude of the Earth surface above sea level (m)
  char *Site; //!<name of the site
  char *Date; //!<date of the simulated event (day/month/year, eg 13/May/2021)
  char *Latitude; //!<Latitude of the site in degrees
  char *Longitude; //!<Longitude of the site in degrees
  double *BField; //!< Magnitude of the B-field in uT
  double *BFieldIncl; //!< Inclination B-field (degrees)
  double *BFieldDecl; //!<Declination B-field (degrees)
  char *AtmosphericModel; //!<Name of the atmospheric model used
  char *AtmosphericModelParameters; //!<parameters set in the atmospheric model
  double *EnergyInNeutrinos; //!< amount of energy in shower-neutrinos (EeV)
}Eventinfo;

typedef struct{
  int size; //!< size of the Showersiminfo table
  char *buffer; //!< buffer in which the data is stored
  char *ShowerSimulator; //!< name and version of Aires
  char *HadronicModel; //!<hadronic interaction model used  (same as in RunInfo)
  char *RandomSeed; //!< Seed used in simulation
  char *RelativeThinning; //!< Thinning factor
  double * WeightFactor; //!< Wight factor
  char *GammaEnergyCut; //!<  Minimal gamma energy (MeV)
  char *ElectronEnergyCut; //!< Minimal electron energy (MeV)
  char *MuonEnergyCut; //!< Minimal muon energy (MeV)
  char *MesonEnergyCut; //!< Minimal meson energy (MeV)
  char *NucleonEnergyCut; //!< Minimal nucleonenergy (MeV)
  double *CPUTime; //!< CPU time used for simulation (s)
  char *OtherParameters; //!< Potential other parameters
}Showersiminfo;

typedef struct{
  int size; //!<Size of Signalsiminfo table
  char *buffer; //!< buffer in which the data is stored
  char *FieldSimulator; //!< program used to simulate the electric field
  char *RefractionIndexModel; //!< model simulating the refractive index in air
  char *RefractionIndexModelParameters; //!<parameters of the model
  double *TimeBinSize; //!< size of the time steps (ns)
  double *TimeWindowMin; //!< minimal time  of the window (ns)
  double *TimeWindowMax; //!<maximal time of window (ns)
  char *OtherParameters; //!< placeholder
}Signalsiminfo;

typedef struct{
  char *ID; //!< identifier of antenna
  float *pos;//!< (x,y,z) position of the antenna (m)
  float *T0;//!< time offset (ns)
  double *SlopeA,*SlopeB; //!< I don't know what these slopes are (m)
}Antennainfo;

typedef struct{
  int n_point; //!< number of points in the simulation
  float *Ex,*Ey,*Ez; //!< time series E-field (uV/m)
  float *Vx,*Vy,*Vz; //!< time series Voltage values (uV)
  float *Vfx,*Vfy,*Vfz; //!< time series filtered voltages (uV)
  float *time; //!< the array of time values (ns)
}AntennaTrace;

typedef struct{
  char *ID; //!< identifier of the antenna
  float *P2P_efield; //!< The P2P of E_tot, E_x, E_y, E_z (uV/m)
  float *P2P_voltage; //!<The P2P of V_tot, V_x, V_y, V_z (uV)
  float *P2P_filtered; //!<The P2P of the filtered V_tot, V_x, V_y, V_z (uV)
  float *HilbertPeakE;//!< The peak of the Hilbert envelope of the E-field (uV/m)
  float *HilbertPeakTimeE; //!< The time at which the Hilbert envelope of the E-field peaks (ns)
  float *HilbertPeakV; //!< The peak of the Hilbert envelope of V_tot (uV)
  float *HilbertPeakTimeV;//!< The time at which the Hilbert envelope of V_tot peaks (ns)
  float *HilbertPeakVf; //!< The peak of the Hilbert envelope of the filtered V_tot (uV)
  float *HilbertPeakTimeVf; //!< The time at which the Hilbert envelope of the filtered V_tot peaks (ns)
}AntennaP2P;

typedef struct{
  float *Distance; //!< Distance from shower axis (m)
  float *Ngamma; //!< Number of photons
  float *Ne_plus_minus; //!< Number of electrons+positrons
  float *Ne_plus; //!< Number of positrons
  float *Nmu_plus_minus;//!< Number of positively and negatively charged muons
  float *Nmu_plus; //!< Number of positively charged muons
  float *Nall_charged; //!< Number of charged particles
}LateralProfile;

typedef struct{
  float *SlantDepth;//!< Amount of atmosphere traversed along the path of shower/primary (g/cm^2)
  float *VerticalDepth;//!< vetrical depth of Xmax from top of atmosphere (g/cm^2)
  float *Ngamma; //!< Number of photons
  float *Ne_plus_minus;//!< Number of electrons+positrons
  float *Ne_plus; //!< Number of positrons
  float *Nmu_plus_minus;//!< Number of positively and negatively charged muons
  float *Nmu_plus;//!< Number of positively charged muons
  float *Npi_plus_minus;//!< Number of positively and negatively charged pions
  float *Npi_plus;//!< Number of positively charged pions
  float *Nall_charged;//!< Number of charged particles
}LongProfile;

typedef struct{
  int n_lateral; //!< number of lateral samples
  char *lateralbuffer; //!< buffer containing the lateral profiles
  LateralProfile *lateralprofile; //!< Pointer to the n_lateral profiles
  int n_long; //!< number of longitudinal samples
  char *longbuffer; //!< buffer containing the longitudinal profiles
  LongProfile *longprofile; //!< Pointer to the n_long profiles
}ShowerTable;

typedef struct{
  Runinfo runinfo; //!< run information
  Eventinfo eventinfo; //!< event information
  Showersiminfo showersiminfo; //!< shower simulation information
  Signalsiminfo signalsiminfo; //!< radio signal simulation information
  int n_ant; //!< Number of antennas
  char *antennabuffer; //!< buffer holding the antenna information of each antenna
  Antennainfo *antennainfo; //!< Pointer to all antenna information structures
  char *tracebuffer; //!< buffer holding all traces
  AntennaTrace *antennatrace; //!< Pointer to all antenna trace structures
  char *antennaP2Pbuffer; //!< buffer holding all peak-to-peak data
  AntennaP2P *antennap2p; //<! pointer to the peak-to-peak structures
  ShowerTable showertable; //!< the shower table information
}AIRESEvent;

