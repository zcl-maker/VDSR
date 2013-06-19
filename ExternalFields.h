/***************************************************************************
                          ExternalFields.h -
						  In the original VDSR code before May.17, 2012, this class was named LaserPlasma
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_ExternalFields
#define H_ExternalFields

#include "vdsrclass.h"
#include "vdsr.h"
#include "stdlib.h"


//---------------------------- Particle class -----------------------
class ExternalFields : public NList{
	friend class Particle;
	friend class Domain;
	friend class Beam;

public:
	static ExternalFields *p_ExternalFields;
	ExternalFields(char *infile, int rank);      //Pnumber represent the particle number
	double GetWavelength() {return f_Lambda_L0;};
	double GetDensity_Normalize() {return f_Density_Normalize;};
	double DampingCoefficient(double t);
	virtual ~ExternalFields() {;};
	

private:
	int i_my_rank,i_Particle_Field,i_Particle_Trajectory,i_DampingFlag,i_Radiation_Reaction,i_Laser_Polarization,i_AzCalculation,i_Laser_Polarization_2;
    double f_Lambda_L0,f_Length_Normalize,f_Density_Normalize,f_Density_Plasma,f_ExternalField_End,f_gamma_wake,f_r_bubble,f_ExternalField_Duration,f_ExternalField_Intensity,f_Laser_FocusX,
		   f_Laser_FocusY,f_Laser_FocusZ,f_Laser_Kx,f_Laser_Ky,f_Laser_Kz,f_Laser_FocusW,f_Laser_Mode,f_Laser_IniPhase,f_Laser_LtoFocus,
		   f_BeginDampingTime,f_EndDampingTime,f_Undu_B, f_Undu_Lambda, f_Undu_Length,f_Undulator_a,f_Undulator_b,
		   f_ExternalField_Duration_2,f_ExternalField_Intensity_2,f_Laser_FocusX_2,f_Lambda_L2,
		   f_Laser_FocusY_2,f_Laser_FocusZ_2,f_Laser_Kx_2,f_Laser_Ky_2,f_Laser_Kz_2,f_Laser_FocusW_2,f_Laser_Mode_2,f_Laser_IniPhase_2,f_Laser_LtoFocus_2;
	double f_Laser1_GaussianEnergy;
	double *Data_ArrayEx,*Data_ArrayEz,*Data_ArrayBy;
	void Set_Field(double t,double x,double y,double z,double &Ex,double &Ey,double &Ez,double &Bx, double &By, double &Bz);
	void Cal_Field(double t,double x,double y,double z,double &Ex,double &Ey,double &Ez,double &Bx, double &By, double &Bz, int i_FieldType);
	double Sn(double n,double omega0,double omega,double phi,double phiG){return (pow(omega0/omega,n)*sin(phi+n*phiG));};
	double Cn(double n,double omega0,double omega,double phi,double phiG){return (pow(omega0/omega,n)*cos(phi+n*phiG));};
	double ax(double x,double y,double z,double t);
	double az(double x,double y,double z,double t);
	void Write_Field(long l_StepOrder, double t, double x_start, double y_start, double z_start, double dx, double dy, double dz, long l_Field_Nx, long l_Field_Ny, long l_Field_Nz, long local_N_start, long local_N_end);
    void Write_Field_Parallel(long l_StepOrder, double t, double x_start, double y_start, double z_start, double dx, double dy, double dz, long l_Field_Nx, long l_Field_Ny, long l_Field_Nz, long local_N_start, long local_N_end);
	void Cal_OutField(long l_StepOrder, double t, double x_start, double y_start, double z_start, double dx, double dy, double dz, long l_Field_Nx, long l_Field_Ny, long l_Field_Nz, long local_N_start, long local_N_end);
	double CP_MC_LaserField_GetNp();
	double CP_MC_LaserField_GetPhi(double x,double y,double z,double t);
	double CP_MC_LaserField_GetlocalEnvelope_nplocal(double x,double y,double z,double t);
	double CP_MC_PhotonK_Selection();
};
#endif