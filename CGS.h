/***************************************************************************
                          CGS.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_CGS
#define H_CGS

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "vdsrclass.h"

//---------------------------- Class UnitsCGS -----------------------
class UnitsCGS
{
friend class Detector;
friend class Particle;
friend class Beam;
friend class Domain;
friend class ExternalFields;
private:
  Domain *domain(void) {return Domain::p_D;};
  ExternalFields *externalFields(void) {return ExternalFields::p_ExternalFields;};
  char c_Begin[4];                   // begin of variables
  double f_Qe;
  double f_Me;
  double f_Mp;
  double f_hbar;
  double f_QeSI;

  double f_OmegaLaser,f_Ncrit;
  double f_PhotonLaser_eV,f_PhotonPlasma_eV;
  double f_Wavelength;
  double f_Ts;
  double f_Hx;
  double f_Hy;
  double f_Hz;
  double f_Clight;
  double f_r_e;                       // Electron radius.        
  double f_r_B;                       // Bohr radius.             
  double f_OmegaAtomic;               // Atomic Unit. of freq.    
  double f_Electron_ComptonWavelength;
  char   c_End[4];                    // end of variables
public:
  static UnitsCGS *p_UnitsCGS;
  double GetQe() {return f_Qe;};
  double GetMe() {return f_Me;};
  double GetMp() {return f_Mp;};
  double GetOmegaLaser() {return f_OmegaLaser;};
  double GetCritDensity() {return f_Ncrit;};
  double GetWavelength() {return f_Wavelength;};
  double GetTs() {return f_Ts;};
  double GetHx() {return f_Hx;};
  double GetHy() {return f_Hy;};
  double GetHz() {return f_Hz;};
  double GetClight() {return f_Clight;};

  double Getre()   {return f_r_e;};  
  double GetOmegaAtomic() { return f_OmegaAtomic; }
  double GetrB() { return f_r_B; }
  double Get_Electron_ComptonWavelength() {return f_Electron_ComptonWavelength;};

  UnitsCGS(void);
};
#endif
