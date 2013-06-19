#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "vdsr.h"

//---------------------------- UnitsCGS::UnitsCGS -----------------------
UnitsCGS* UnitsCGS::p_UnitsCGS = NULL;
UnitsCGS::UnitsCGS(void)
{
  UnitsCGS::p_UnitsCGS = this;
  f_Qe = -4.8e-10;
  f_QeSI = -1.6e-19;
  f_hbar=6.626e-34/(2*PI);   // unit in J.s
  f_Me = 0.911e-27;          // unit in g
  f_Mp = 1836.*f_Me;

  f_Wavelength = externalFields()->GetWavelength()*1.0e-4;  // here f_Wavelength in the unit of cm, original externalFields()->GetWavelength() in the unit of um

  f_Clight = 3e10;   // unit in cm
  f_Ts = domain()->GetTs()*f_Wavelength/f_Clight;


  f_OmegaLaser = f_Clight*2.*PI/f_Wavelength;
  f_Ncrit = f_Me*f_OmegaLaser*f_OmegaLaser/(4.*PI*f_Qe*f_Qe);

  //For ADK Probability   from  i_vlpl 090216
  f_OmegaAtomic = 4.134E16; // Atomic Unit. of freq.
  f_r_B = 5.2918E-9; // Bohr radius. unit in cm
  f_r_e = 2.8179E-13; // Electron radius. unit in cm

  f_PhotonLaser_eV=f_OmegaLaser*f_hbar/(-f_QeSI);
  f_PhotonPlasma_eV=f_PhotonLaser_eV*sqrt(externalFields()->GetDensity_Normalize());
  f_Electron_ComptonWavelength=2.4263E-10; // unit in cm   h/m_e c
}
