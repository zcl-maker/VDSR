/***************************************************************************
                          Detector.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_DETECTOR
#define H_DETECTOR

#include "vdsrclass.h"
#include "vdsr.h"



//---------------------------- Detector class -----------------------
class Detector : public NList{
friend class Pixel;
friend class Domain;
friend class Particle;
public:
	UnitsCGS *p_UnitsCGS;
	static Detector *p_Detector;
	long   l_Pixel_number;
	int    i_my_rank,i_OmegaLog;
	double f_total_emission_energy,f_receive_emission_energy;
	Pixel  **Pixel_Array;

	int    GetmyPE(){return i_my_rank;};
	void   WriteRadiation();
	void   Run(int i_Radiation_Type);
	void   Global_radiation_reduction();
	void   Global_Intensity_reduction();
	void   Projector_Rotation();
    Detector(char *infile, int rank);
    ~Detector() {;};
private:
	int    i_RealDetector,i_InitialEndIntegration;
	long   l_nEbins,l_nPhibins,l_nThetabins;
    double *f_Omega,f_omega_max,f_omega_min,f_theta_max,f_theta_min,f_phi_max,f_phi_min,
		   f_domega,f_dphi,f_dtheta,f_PhotonPlasma_eV,f_R_dp,f_Xlength,f_Ylength,f_DCenterX,f_DCenterY;
	long l_Emission_number,l_Recorder_number,l_Emission_number_tot,l_Recorder_number_tot;
};
#endif