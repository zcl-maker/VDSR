/***************************************************************************
                          Photon.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_PHOTON
#define H_PHOTON

#include "vdsrclass.h"
#include "vdsr.h"
#include "stdlib.h"


//---------------------------- Particle class -----------------------
class Photon {
	friend class Particle;
	friend class Domain;
	friend class Pixel;
public:
	Photon(double X_lab,double Y_lab, double kx, double ky, double kz);      // k is normalized by kp=2pi/lambda_p
	virtual ~Photon() {;};
private:
	double X_lab,Y_lab,k_vector[3];
	Photon *p_Photon_Next,*p_Photon_Precede; 
};

#endif