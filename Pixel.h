/***************************************************************************
                          Pixel.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_PIXEL
#define H_PIXEL

#include "vdsrclass.h"
#include "vdsr.h"


//---------------------------- Particle class -----------------------
class Pixel : public NList{
friend class Particle;
friend class Detector;

public:
	long l_NumberOmega;

	Detector *detector() {return Detector::p_Detector;};
	Pixel(double Pixel_theta, double Pixel_phi);
    void ProjectSinglePa(Particle *pa, long k);
	void GetIntensity();
	double ReturnIntensityX(int i_omega);
	double ReturnIntensityY(int i_omega);
	double ReturnIntensityZ(int i_omega);
	void InitialEndModification(Particle *pa);
    virtual ~Pixel() {;};

private:
    double theta,phi,omega_max,f_omega_min,f_domega,*IntensityX,*IntensityY,*IntensityZ,*Im_IntegrationX,*Re_IntegrationX,*Im_IntegrationY,*Re_IntegrationY,*Im_IntegrationZ,*Re_IntegrationZ;
	void CleanImRe_IntegrationXYZ();

};
#endif