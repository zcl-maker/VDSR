#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include <time.h>
clock_t clock(void);

#include "vdsr.h"
Controls* Controls::p_Controls = NULL;
//---------------------------- Controls::Controls -----------------------
Controls::Controls (char *infile, int rank) : NList ("Controls")
{
   Controls::p_Controls = this;
   AddEntry("f_CPUstop", &f_CPUstop);

   FILE *p_File = NULL;


   if (rank==NIL)
   {
      p_File = fopen(infile,"rt");     
      if (p_File == NULL)
      {
		  cout<<"Controls no input file: "<<infile<<" "<<endl;
		  // following shoule be noted,since ExternalFields class is defined earlier than domain class
		  //Domain::p_D->out_Flog<<"ExternalFields no input file: "<<infile<<"\n";
          exit (-1);
      }
   }

   if (p_File)
   {
      rewind(p_File); 
      read(p_File); 
   }
 
#ifdef V_MPI
   CBuffer *buf = new CBuffer;
   buf->reset();
   pack_nls(buf);
   Domain::p_D->BroadCast(buf);
   if (rank) unpack_nls(buf);
   delete buf;
#endif
  

   f_CPU = 0.;
   f_CPU_start = 0.;
   f_CPU_finish = 0.;

   f_CPU_start = clock();
   time( &t_WallClockStart );
   time( &t_WallClock );
   if(rank==0){
	   cout<<"Controls finished!"<<endl;
   }
}

//--- Controls:: ----------------------->
double Controls::GetCPU(void)
{
   clock_t cl;
   double ftime;

   f_CPU_finish = clock();
   /* printf("\n cl = %d, Divider = %d \n",cl,CLOCKS_PER_SEC);*/
   /*time = (double)cl;*/

   ftime = (f_CPU_finish - f_CPU_start);

   if (ftime < 0.) {
      ftime += 2147483647;
   }

   f_CPU_start = clock();
   f_CPU += ftime/CLOCKS_PER_SEC;

   /* printf("\n time = %g \n",time);*/

   time( &t_WallClock );

   f_WallClockElapsed = difftime( t_WallClock, t_WallClockStart );

   return f_CPU;
}

