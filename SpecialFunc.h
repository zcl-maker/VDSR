/***************************************************************************
                          SpecialFunc.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_SPECIALFUNC
#define H_SPECIALFUNC

#include "vdsrclass.h"
#include "vdsr.h"
#include "stdlib.h"

#define PIBY2 (3.1415926/2.0)
#define TRUE 1
#define ONE complex<double>(1.0,0.0)
using std::complex;

//---------------------------- Particle class -----------------------
class SpecialFunc {
	friend class Particle;
public:
	static SpecialFunc *p_SpecialFunc;
	SpecialFunc();      
	virtual ~SpecialFunc() {;};

	void fresnel(double x, double *s, double *c); //Computes the Fresnel intergrals S(x) and C(x) for all real x.   From P 256,  <<Numercal Recipes in C>>   F(x)=\int_0^x {cos(PI*t^2/2) dt}

    void airy(double x, double *ai, double *bi, double *aip, double *bip);  
	//Returns Airy functions Ai(x), Bi(x), and their derivatives Ai'(x), Bi'(x).

	void bessik(double x,double xnu,double *ri,double *rk,double *rip,double *rkp);
	//Returns the modified Bessel functions ri = Ixnu, rk = Kxnu and their derivatives rip = I'xnu, rkp = K'xnu, for positive x and for xnu >= 0.

	void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);
	//Returns the Bessel functions rj = J¦Í, ry = Y¦Í and their derivatives rjp = J'xnu, ryp = Y'xnu, for positive x and for xnu>= 0.
	
	void beschb(double x,double *gam1,double *gam2,double *gampl,double *gammi);

	double chebev(double a, double b, double c[], int m, double x);

	double IntAiry(double z0, double z1);

private:
};

#endif