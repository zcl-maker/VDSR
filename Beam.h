/***************************************************************************
                          Beam.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_BEAM
#define H_BEAM

#include "vdsrclass.h"
#include "vdsr.h"
#include "stdlib.h"


//---------------------------- Particle class -----------------------
class Beam : public NList{
	friend class Particle;
	friend class Domain;
public:
	static Beam *p_Beam;
	Beam(char *infile, int rank, long Pnumber, int numprocs, int flag);      //Pnumber represents the particle number, flag=1 means setting beam initial positions from input file, =0 means reading trajectories from another file
	virtual ~Beam() {;};
private:
	int    i_ShapeX,i_ShapeY,i_ShapeZ,i_ShapePX,i_ShapePY,i_ShapePZ,i_ShapeR,i_ShapePR,i_my_rank,i_numprocs;
	long   l_Pnumber, l_myPnumber,l_start_position;
    double f_Qcharge,f_Energy,f_EnergySpread,f_RadiusX,f_RadiusY,f_RadiusZ,f_RadiusR,f_CenterX,f_CenterY,
		   f_CenterZ,f_CenterR,f_RadiusPX,f_RadiusPY,f_RadiusPZ,f_RadiusPR,f_CenterPX,f_CenterPY,f_CenterPZ,f_CenterPR,f_ExpTemp,f_Exp_min,f_Exp_max,*f_sampling_data,*f_using_data;
	void   Sampling(double f_Center,double f_Radius,long l_number, int i_shape);
	double F_rand();
};

#endif