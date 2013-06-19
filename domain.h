#ifndef H_DOMAIN
#define H_DOMAIN

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "vdsrclass.h"

//---------------------------- Domain class -----------------------
class Domain : public NList{
 friend class Particle;

 public:
  static Domain *p_D;
  char   *p_BufferMPI;
  int    i_Radiation_Calculation,i_my_rank,i_CP_ParticleEnergyChange;
  long   l_BufferMPIsize,l_Step,l_Step_Cut1,l_Step_Cut2,pNumber,l_Pnumber;
  double f_Amplify_Factor,f_MaxScattering_Prob;
  Particle **pa;
  ofstream out_Flog;

  int    GetmyPE(){return i_my_rank;};
  double GetTs(){return (f_dt);};

 private:
  int  i_Code,i_Read_Trajectory,i_Integration,i_numprocs,i_OutputField,i_OutputFieldEvolution;
  long l_OutputPTrajec,l_OutputBeamDis,l_OutputField,l_Field_Nx,l_Field_Ny,l_Field_Nz;
  char c_Begin[4];
  double f_dt,f_MovingWindowt,f_xLength,f_yLength,f_zLength,f_xPosition,f_yPosition,f_zPosition,f_x0Output,f_y0Output,f_z0Output;
  char c_End[4];

  FILE *p_File;
  FILE *p_SaveFile;
  FILE *p_MovieFile;
  FILE *FBeamdistribution;
  ofstream out_Fig8;

  char *str_File;
  char *str_SName;
  char *str_DName;
  char *str_LogName;
  char str_FileName[200];
  char str_HistoryFname[200];
  char *str_DataDirectory;
  char *str_LogDirectory;
  char *str_MovieDirectory;
  char *str_HistoryDirectory;

  void WriteBeamDistribution(long l_StepOrder);
  void WriteBeamDistribution_Parallel(long l_StepOrder);
  void OutputFieldDistribution(int i_OutputField, long l_StepOrder);
  void OutputFieldEvolution(int i_my_rank);
  void OutputBeamDiagnosis(long l_StepOrder);
  void Read_VLPL_HISTORY(char *infile, int rank, int numprocs);
  void Read_VORPAL_HISTORY(char *infile, int rank, int numprocs);

public:                      // functions

   FILE* GetIniFile() {return p_File;}
   ofstream& Getout_Flog() {return out_Flog;}  


   CBuffer*	 GetBufMPP();
   void		 ResetBufMPP();

   void BroadCast(CBuffer *b);  
   Domain(char *infile, int rank, int numprocs);
   ~Domain();
};

#endif
