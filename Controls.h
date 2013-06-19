/***************************************************************************
                          Controls.h -
    email                : mchen911@ustc.edu
 ***************************************************************************/
#ifndef H_CONTROLS
#define H_CONTROLS

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "vdsrclass.h"
#include "vdsr.h"

//---------------------------- Controls class -----------------------

class Controls : public NList {
   friend class Domain;
   friend class Detector;
   friend class ExternalFields;
   friend class Particle;
   friend class Pixel;

private:
   Domain *domain() {return Domain::p_D;};
   char c_Begin[4];                     // begin of variables
   double f_CPU;
   double f_CPU_start;
   double f_CPU_finish;
   time_t t_CPU_start;
   time_t t_CPU_finish;
   time_t t_WallClockStart;
   time_t t_WallClock;
   double f_WallClockElapsed;

   double f_CPUstop;


   char   c_End[4];                    // end of variables


public:

   static Controls *p_Controls;
   double GetCPU(void);                // return CPU time, unit in second
   double GetWallClockElapsed(void) { return f_WallClockElapsed;};   // return wall elapsed time, unit in second

   Controls (char *infile, int rank);
};

#endif

