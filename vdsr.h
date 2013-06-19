/***************************************************************************
                         vdsr.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/#ifndef H_VDSR
#define H_VDSR

#include <iostream>
using namespace std;

#include <vector>   
#include <cmath>    
#include <assert.h> 
#include <math.h>
#include "myhdfshell.h"
//#include <omp.h>


#define NIL 0
#define PI 3.1415926
#define TWOPI (2*PI)
#define H5Dcreate H5Dcreate1
#define H5Dopen H5Dopen1
#undef _DEBUG

#define V_MPI   //add by mchen for pc parallized

#ifdef V_MPI
#include <mpi.h>
#endif

//#define PARALLELHDF5

#include "buffers.h"
#include "namelist.h"

#include "vdsrclass.h"
#include "Beam.h"
#include "domain.h"
#include "Particle.h"
#include "Detector.h"
#include "Pixel.h"
#include "ExternalFields.h"
#include "CGS.h"
#include "Controls.h"
#include "Photon.h"
#include "SpecialFunc.h"


#endif

