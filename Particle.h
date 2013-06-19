#ifndef H_PARTICLE
#define H_PARTICLE


#include "vdsr.h"
#include "vdsrclass.h"

#include <complex>

/*
#define EPS 6.0e-8
#define MAXIT 100
#define FPMIN 1.0e-30
#define XMIN 1.5
#define PIBY2 (3.1415926/2.0)
#define TRUE 1
#define ONE complex<double>(1.0,0.0)
*/
using std::complex;

/***************************************************************************
                          Particle.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/

//---------------------------- Particle class -----------------------
class Particle{
	friend class Pixel;
	friend class Domain;
	friend class Beam;
	friend class Detector;
public:
	long l_Step,l_Recorder;
    SpecialFunc *p_SpecialFunc;
	void Set_XP(double x0, double y0, double z0, double vx0, double vy0, double vz0);   // read or calculate particle trajecotry and momenta
	void Set_XP0(double x0, double y0, double z0, double vx0, double vy0, double vz0);  // only set initial particle positions and momenta
	void Read_XP_VLPL(char *str_HistoryFname,long i);  // here i is the particle label
	void Read_XP_VORPAL(char *str_HistoryFname,long i);
	void Write_XP(long l_Porder);	
	void Interpolation(long k);
	Particle(long l_Porder);      //Pnumber represent the particle number
	virtual ~Particle() {;};
	Photon *Photon_String_Head, *Photon_String_Current;
	UnitsCGS *p_UnitsCGS;
    double tp_particle,f_averageemissionnumber_classical;
	double polarization[3],polarizationprime[3],kexi1lab,kexi2lab,kexi3lab;
private:
	FILE *Ftrajectory;
	int i_Particle_Field;
	long l_Porder,l_emission_number;
    double *x,*y,*z,*vx,*vy,*vz,*Field,weight,dt,ImIjk[3],ReIjk[3],LFold[3],b1[3],a1[3],a2[3],RK[6],Charge2Number_weight,f_Amplify_Factor;
	double gamma,px,py,pz,pxnew,pynew,pznew,gammanew;                      // temporal variables in the Lab Frame
	double gammaprime;
	union
    {
      double p3prime[1];
      double pxprime;
    };
	double pyprime,pzprime;  // temporal variables in the Transformed Frame
	double f_k2p;

	void InitialEndModification(double theta, double phi, double omega);
	void StepIntegration(double theta, double phi, double omega, long k);
	void StepIntegration(double theta, double phi, double omega, long k, int i); // here use the Taylor expansion to calculate this stepintegration
	void StepIntegration(double theta, double phi, double omega, long k, double i); // here use the Fresnel integration to calculate this stepintegration
	//void fresnel(double x, double *s, double *c);  // 111012, this function has been moved to the new class SpecialFunc. 
	double f_field(double t, double x,double y,double z, double px, double py, double pz);
	double k_photon_lab0[4],k_photon[4],k_photon_prime[4],TMatrix[9];


	void   CP_Motion(double px_new,double py_new, double pz_new, long k);  // k is the current time step, k+1 is the next time step
    void   CP_MC_PhotonEmissionRecorder(double X_lab,double Y_lab, double kx_photon, double ky_photon, double kz_photon, long k);  // finished
    void   CP_MC_ProjectSinglePhoton();
	double CP_MC_ScatteringPossibility(double kx,double ky,double kz,double npklocal, long kstep);  // input photon information, calculation scattering probability for current electron
	double CP_MC_CrossSectionTot(double px,double py,double pz,double kx,double ky, double kz);  // total scattering cross section
	double CP_MC_CrossSectionEgPhi(); // differential cross section  dsigema/dEgdphi
	double CP_MC_CrossSectionEg();    // differential cross section  dsigema/dEg
	double CP_MC_FormularEg2Theta(double Egprime,double Epprime);  // using Compton Formula from Eg' to calculate Theta' in the electron rest frame
	void   CP_MC_GetSCPhotonPolarization(double Pt,double Pc,double X_lab, double Y_lab,double phif);
	void   CP_MC_Lorentz_Transform(double beta1,double beta2,double beta3, double &x0, double &x1, double &x2, double &x3, double &x0p, double &x1p, double &x2p, double &x3p); // From K frame to K'frame, K' parallel to K, and K' has a velocity of beta in K frame
	void   CP_MC_Lorentz_Transform_EPolarization(double beta1,double beta2,double beta3,double k1,double k2,double k3); // From K frame to K'frame, K' parallel to K, and K' has a velocity of beta in K frame, photon polarization in K frame is (p1,p2,p3), in K' frame is (p1',p2',p3')
	void   CP_MC_Lorentz_Rotation(double *TMatrix,double *x3); //TMatrix is a 3X3 matrix defines the rotation matrix, x0, x1,x2 will be rotated,  
    double CP_MC_EgprimeSampling(double Epprime);
    double CP_MC_Fintegration(double Epprime,double x);
	double CP_MC_GetEgprime(double random_F,double Epprime,double Egprimemin, double Egprimemax);
	double CP_MC_PhiprimeSampling(double Egprime,double Epprime,double tao, double Pt);
	double CP_MC_GetPhiprime(double c1,double tao,double random_G);
	double CP_MC_Gintegration(double c1, double tao, double phi);
	void   CP_MC_GetLocalNp(double x,double y,double z, double t,double *k_photon_Lab);
	void   CP_MC_ComptonRadiation(long kstep);
	double Np_total;
	void   Matrix_Inverse(double *TMatrix);
	void   Set_TMatrix(double kx,double ky,double kz,double nx,double ny,double nz);
	double determinant_TMatrix3(double *TMatrix);
	void   Print_TMatrix3(double *TMatrix);


	void   CP_LBNL_ComptonRadiation(long kstep);
	void   CP_LBNL_PhotonEmissionRecorder(double kx,double ky,double kz,long kstep);
    double CP_LBNL_GetChi_e(long k);
	double CP_LBNL_ScatteringPossibility(double p0,double Chi_e);  // return the scattering possibility during a single time step
	void   CP_LBNL_ProjectSinglePhoton();
	double CP_LBNL_EgSampling(double Chi_e);



};

#endif