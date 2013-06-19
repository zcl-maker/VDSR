#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "vdsr.h"

static char *MeshDataName = "History_Particles";

Domain* Domain::p_D = NULL;
int MakeNewDirectory(char* directory);

//---------------------------- Domain::Domain -----------------------
Domain::Domain (char *infile, int rank, int numprocs) : NList("Domain")
{
   int i = 0;
   i_my_rank=rank;
   p_BufferMPI = NULL;
   l_BufferMPIsize = 0;
   i_numprocs=numprocs;

   Domain::p_D = this;

   str_SName = ".save";
   str_DName = ".dat";
   str_LogName = ".log";

   str_DataDirectory = ".";
   str_LogDirectory = ".";          
   str_MovieDirectory = "."; 
   str_HistoryDirectory = ".";
 
   int ierr=0;  

   AddEntry("f_Ts", &f_dt, 0.);
   AddEntry("f_MovingWindowt", &f_MovingWindowt,100.0);
   AddEntry("l_Step", &l_Step, 10);
   AddEntry("l_Step_Cut1", &l_Step_Cut1, 0);
   AddEntry("l_Step_Cut2", &l_Step_Cut2, 0);
   AddEntry("l_Pnumber",&l_Pnumber,1);
   AddEntry("l_OutputPTrajec",&l_OutputPTrajec,1); 
   AddEntry("l_OutputBeamDis",&l_OutputBeamDis,10);
   AddEntry("i_Integration", &i_Integration,0); 
   AddEntry("i_Code",&i_Code,1);
   AddEntry("i_Read_Trajectory",&i_Read_Trajectory,0);
   AddEntry("i_Radiation_Calculation",&i_Radiation_Calculation,1);
   AddEntry("l_OutputField",&l_OutputField,100000);
   AddEntry("f_xLength",&f_xLength,100.0);
   AddEntry("f_yLength",&f_yLength,100.0);
   AddEntry("f_zLength",&f_zLength,100.0);
   AddEntry("f_x0Output",&f_x0Output,0.0);
   AddEntry("f_y0Output",&f_y0Output,0.0);
   AddEntry("f_z0Output",&f_z0Output,0.0);
   AddEntry("l_Field_Nx",&l_Field_Nx,100);
   AddEntry("l_Field_Ny",&l_Field_Ny,1);
   AddEntry("l_Field_Nz",&l_Field_Nz,100);
   AddEntry("i_OutputField",&i_OutputField,0);
   AddEntry("f_xPosition",&f_xPosition,1.0);
   AddEntry("f_yPosition",&f_yPosition,0.0);
   AddEntry("f_zPosition",&f_zPosition,1000.0);
   AddEntry("i_OutputFieldEvolution",&i_OutputFieldEvolution,0);
   AddEntry("f_Amplify_Factor",&f_Amplify_Factor,1.0);
   AddEntry("f_MaxScattering_Prob",&f_MaxScattering_Prob,1.0);
   AddEntry("i_CP_ParticleEnergyChange",&i_CP_ParticleEnergyChange,1);

   l_OutputPTrajec=l_Pnumber/l_OutputPTrajec;

   p_File = NULL;

   sprintf(str_FileName,"%s//v%d.log",str_LogDirectory,rank); 
   if(i_my_rank==0) out_Flog.open(str_FileName, ios::out);   //Here added if(i_my_rank==0) 

   if (rank==NIL)
   {
      p_File = fopen(infile,"rt");     
      if (p_File == NULL)
      {
         out_Flog << "Error. Domain: No such file " << infile << "\n";  
         exit (-1);
      }
	  else{
         out_Flog << "Domain: Input file " << infile << " is found!\n";  

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
   BroadCast(buf);
   if (rank) unpack_nls(buf);
   delete buf;

#endif


   if(i_my_rank==0){
	   if(i_Integration==0){ 
			cout<<"i_Integration=0, Taylor Expansion for integration will be used!"<<endl;
			out_Flog<<"i_Integration=0, Taylor Expansion for integration will be used!\n\n";
	   }
	   else{
            cout<<"i_Integration="<<i_Integration<<", Fresnel Integration will be used!"<<endl;
			out_Flog<<"i_Integration="<<i_Integration<<", Fresnel Integration will be used!\n\n";
	   }
   }
		
   if(i_Read_Trajectory==1){                        
	   if(i_Code==1){
		   if(rank==0) {
			   cout<<"VLPL's output file History.h5 will be read for the particles' trajectories."<<endl;
			   out_Flog << "VLPL's output file History.h5 will be read for the particles' trajectories.\n";  
		   }
		   Read_VLPL_HISTORY(infile, rank, numprocs);
	   }
	   else if(i_Code==2){
		   if(rank==0) {
			   cout<<"VORPAL's output file History.h5 will be read for the particles' trajectories."<<endl;
			   out_Flog<<"VORPAL's output file History.h5 will be read for the particles' trajectories.\n";
		   }
		   Read_VORPAL_HISTORY(infile, rank, numprocs);
	   }
	   else if(i_Code==3)
		   cout<<"Sorry we cannot process the input from i_Code==3."<<endl;
	   else 
		   cout<<"Sorry we cannot process the input."<<endl;

   }
   else if(i_Read_Trajectory==0){   // Set the particle initial positions according to the beam parameters and then calculate the particle trajectories.
	   
	   if(i_OutputFieldEvolution==1)   OutputFieldEvolution(rank);
	   if(rank==0){
		   cout<<"\n------------ BeforeFieldOutput CPU TIME= "<<Controls::p_Controls->GetCPU()<<endl;
		   out_Flog<<"\n------------ BeforeFieldOutput CPU TIME= "<<Controls::p_Controls->GetCPU()<<"\n";
	   }
	   if(i_OutputField>0){
		   for(long i=0;i<l_Step;i=i+l_OutputField){
			   OutputFieldDistribution(i_OutputField,i);
			   if(rank==0){
				   cout<<"\n------------ AfterField["<<i<<"] Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<endl;
				   out_Flog<<"\n------------ AfterField["<<i<<"] Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<"\n";
			   }
		   }
	   }

       
	   if(rank==0) {
			   cout<<"\n****** Beam Initialization......\n"<<endl;
			   out_Flog<<"\n****** Beam Initialization......\n";
		   }

	   pNumber=l_Pnumber/numprocs;                                          // this is the particle within one process
       if(numprocs-i_my_rank<l_Pnumber%numprocs+1)                          // Here means the last few processes will possibly take one particle more than the first few processes.
		   pNumber++; 
	   if(pNumber>0){
		   pa=new Particle*[pNumber];
		   Beam *beam = new Beam(infile, i_my_rank, l_Pnumber,numprocs,1);
	   

		   for(long i=0;i<pNumber;i++){                                      // initial particle positions from the sampling value
			   if(rank==0 && i%1000==0){
				   cout<<"Particle Trajectory Calculation begins "<<100.0*i/pNumber<<"% at CPU TIME:"<<Controls::p_Controls->GetCPU()<<endl;
				   out_Flog<<"Particle Trajectory Calculation begins "<<100.0*i/pNumber<<"% at CPU TIME:"<<Controls::p_Controls->GetCPU()<<"\n";
			   }
				pa[i]=new Particle(beam->l_start_position+i);
				if(i_Radiation_Calculation==0 || i_Radiation_Calculation==1)
					pa[i]->Set_XP(beam->f_using_data[i*6+0],beam->f_using_data[i*6+1],beam->f_using_data[i*6+2],beam->f_using_data[i*6+3],beam->f_using_data[i*6+4],beam->f_using_data[i*6+5]);            if(i_Radiation_Calculation==2 || i_Radiation_Calculation==3)
					pa[i]->Set_XP0(beam->f_using_data[i*6+0],beam->f_using_data[i*6+1],beam->f_using_data[i*6+2],beam->f_using_data[i*6+3],beam->f_using_data[i*6+4],beam->f_using_data[i*6+5]);
	       }
          
		   for(long i=0;i<l_Step;i=i+l_OutputBeamDis){
			   if(rank==0){
					cout<<"\n------------ Before BeamDis["<<i<<"] Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<endl;
                    out_Flog<<"\n------------ Before BeamDis["<<i<<"] Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<"\n";
			   }
#ifdef PARALLELHDF5
			   if(rank==0 && i==0){
					cout<<"\n********************************************\nPARALLELHDF5 has been defined!"<<endl;
					cout<<"\n********************************************\nPARALLELHDF5 has been defined!\n";
			   }
			    WriteBeamDistribution_Parallel(i);
#else
				WriteBeamDistribution(i);
#endif

				if(rank==0){
					cout<<"------------ After BeamDis["<<i<<"] Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<endl;
					out_Flog<<"------------ After BeamDis["<<i<<"] Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<"\n";
				}
		   }
	   }
		
   }
   else {
	   if(i_my_rank==0)   cout <<"Sorry I do not know how to get the particles' trajectory. Done and Exit!"<<endl;
	   exit(-1);
   }

   if(rank==0){
					cout<<"------------ Single Particle Trajectory Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<endl;
					out_Flog<<"------------ Single Particle Trajectory Output CPU TIME= "<<Controls::p_Controls->GetCPU()<<" Walltime Elapsed "<<Controls::p_Controls->GetWallClockElapsed()<<"\n";
				}
   for(long i=0;i<pNumber;i++){ 
	   if((Beam::p_Beam->l_start_position+i)%l_OutputPTrajec==0) pa[i]->Write_XP(Beam::p_Beam->l_start_position+i);  //every l_OutputPTrajec particle output one particle
   }

   if(i_my_rank==0){
		cout<<"dt="<<f_dt<<" l_Step="<<l_Step<<"  l_Pnumber="<<l_Pnumber<<endl;
		cout << "***********Domain is created***************"<<endl;
		cout<<endl;
	    out_Flog<<"dt="<<f_dt<<" l_Step="<<l_Step<<"  l_Pnumber="<<l_Pnumber<<"\n";
		out_Flog << "***********Domain is created***************\n";
		out_Flog<<"\n";
   }
}

void Domain::Read_VLPL_HISTORY(char *infile, int rank, int numprocs)
{
	 // following read History.h5 from VLPL-PIC simulation
	 f_dt=0.0;
     hid_t       F_HistoryFile, fdataset;         /* File and dataset            */
     hid_t       fdataspace;   /* Dataspace handles           */
     hsize_t     dimsf[3], dimsfi[3], maxdims[3];              /* Dataset dimensions          */
     hid_t rank1 = 1;
     hid_t rank2 = 2;
     hid_t rank3 = 3;
	 sprintf(str_HistoryFname,"%s//History.h5",str_HistoryDirectory);
	 F_HistoryFile = H5Fopen (str_HistoryFname, H5F_ACC_RDONLY, H5P_DEFAULT);
	 fdataset = H5Dopen(F_HistoryFile,MeshDataName);
	 assert (fdataset >= 0);
	 fdataspace= H5Dget_space(fdataset);
	 H5Sget_simple_extent_dims(fdataspace,dimsf,maxdims);
	 l_Step=dimsf[1]-l_Step_Cut1-l_Step_Cut2;
	 l_Pnumber=abs(l_Pnumber)>dimsf[0]?dimsf[0]:abs(l_Pnumber);
	 pNumber=l_Pnumber/numprocs;            // this is the particle within one process
     if(numprocs-i_my_rank<l_Pnumber%numprocs+1)                          // Here means the last few processes will possibly take one particle more than the first few processes.
		pNumber++; 
	 if(pNumber>0){
		   pa=new Particle*[pNumber];
		   Beam *beam = new Beam(infile, i_my_rank,l_Pnumber,numprocs,0);  // the last 0 means, beam will not initialize the particles' initial positions
		   for(long i=0;i<pNumber;i++){      // initial particle positions from the sampling value
				pa[i]=new Particle(beam->l_start_position+i);
				pa[i]->Read_XP_VLPL(str_HistoryFname,beam->l_start_position+i);
	       }

		   for(long i=0;i<l_Step;i=i+l_OutputBeamDis)
			    WriteBeamDistribution(i);
	 }
	 H5Fclose(F_HistoryFile);
};

void Domain::Read_VORPAL_HISTORY(char *infile, int rank, int numprocs)
{
	 // following read History.h5 from VORPAL-PIC simulation
	 //f_dt=0.0;
     hid_t       F_HistoryFile, fdataset, idataset;         /* File and dataset            */
     hid_t       fdataspace;   /* Dataspace handles           */
     hsize_t     dimsf[3], dimsfi[3], dimsp[3],maxdims[3];              /* Dataset dimensions          */
     hid_t rank1 = 1;
     hid_t rank2 = 2;
     hid_t rank3 = 3;
	 char *History_particle_name;
	 History_particle_name = "etraj";
	 sprintf(str_HistoryFname,"%s//test_History.h5",str_HistoryDirectory);
	 F_HistoryFile = H5Fopen (str_HistoryFname, H5F_ACC_RDONLY, H5P_DEFAULT);
	 fdataset = H5Dopen(F_HistoryFile,History_particle_name);
	 assert (fdataset >= 0);
	 fdataspace= H5Dget_space(fdataset);
	 H5Sget_simple_extent_dims(fdataspace,dimsf,maxdims);
	 l_Step=dimsf[0]-l_Step_Cut1-l_Step_Cut2;  
	 l_Pnumber=abs(l_Pnumber)>dimsf[1]?dimsf[1]:abs(l_Pnumber);
	 pNumber=l_Pnumber/numprocs;            // this is the particle within one process
     if(numprocs-i_my_rank<l_Pnumber%numprocs+1)                          // Here means the last few processes will possibly take one particle more than the first few processes.
		pNumber++; 
	 if(pNumber>0){
		   pa=new Particle*[pNumber];
		   Beam *beam = new Beam(infile, i_my_rank,l_Pnumber,numprocs,0);  // the last 0 means, beam will not initialize the particles' initial positions
		   for(long i=0;i<pNumber;i++){      // initial particle positions from the sampling value
				pa[i]=new Particle(beam->l_start_position+i);
				pa[i]->Read_XP_VORPAL(str_HistoryFname,beam->l_start_position+i);
	       }

		   for(long i=0;i<l_Step;i=i+l_OutputBeamDis)
			    WriteBeamDistribution(i);
	 }
	 H5Fclose(F_HistoryFile);
	 
};

//---------------------------- Domain::BroadCast -----------------------
void Domain::BroadCast(CBuffer *b) 
{
#ifdef V_MPI
   int root = 0;
#ifdef _DEBUG
   printf("Start BroadCast\nlen = %d\n\n", b->getlen());
#endif

   int ierr = MPI_Bcast(b->getbuf(), b->getlen(),
      MPI_BYTE, root, MPI_COMM_WORLD);
   switch (ierr)
   {
   case MPI_SUCCESS:
      // No error
      break;
   case MPI_ERR_COMM:
      out_Flog << "Invalid communicator. A common error is to use a null communicator in a call. ierr = " << ierr << ", PE = " << GetmyPE() << endl;
      break;
   case MPI_ERR_COUNT:
      out_Flog << "Invalid count argument. Count arguments must be non-negative; a count of zero is often valid. ierr = " << ierr << ", PE = " << GetmyPE() << endl;
      break;
   case MPI_ERR_TYPE:
      out_Flog << "Invalid datatype argument. May be an uncommitted MPI_Datatype (see MPI_Type_commit). ierr = " << ierr << ", PE = " << GetmyPE() << endl;  					
      break;
   case MPI_ERR_BUFFER:
      out_Flog << "Invalid buffer pointer. Usually a null buffer where one is not valid. ierr = " << ierr << ", PE = " << GetmyPE() << endl;  					  					
      break;
   case MPI_ERR_ROOT:
      out_Flog << "Invalid root. The root must be specified as a rank in the communicator. Ranks must be between zero and the size of the communicator minus one. ierr = " << ierr << ", PE = " << GetmyPE() << endl;  					  					 	
      break;
   default:
      out_Flog << "Unknown error Domain::BroadCast MPI_Bcast. ierr = " << ierr << ", PE = " << GetmyPE() << endl;				
   };

#ifdef _DEBUG
   printf("End BroadCast\n\n");
#endif
#endif
}


void Domain::WriteBeamDistribution(long l_StepOrder){
    double gamma;
	double *writedata;
	long l_TransferPNumber,i,j,k;
	writedata= new double[6*(pNumber+1)];
	MPI_Status status;

	if(i_my_rank==0){
		char s_filename[50];
		sprintf(s_filename,"BeamDistribution_at_%ld.dat",l_StepOrder);
		FBeamdistribution=fopen(s_filename,"w");
		for(k=0;k<pNumber;k++){
			    gamma=1.0/sqrt(1-pa[k]->vx[l_StepOrder]*pa[k]->vx[l_StepOrder]-pa[k]->vy[l_StepOrder]*pa[k]->vy[l_StepOrder]-pa[k]->vz[l_StepOrder]*pa[k]->vz[l_StepOrder]);
				fprintf(FBeamdistribution,"%15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n",pa[k]->x[l_StepOrder],
					    pa[k]->y[l_StepOrder],pa[k]->z[l_StepOrder],pa[k]->vx[l_StepOrder]*gamma,pa[k]->vy[l_StepOrder]*gamma,pa[k]->vz[l_StepOrder]*gamma);
		}
		for(i=1;i<i_numprocs;i++){
			MPI_Send(&pNumber,1,MPI_LONG,i,0,MPI_COMM_WORLD);                           // ask process i to transfer the data
			MPI_Recv(&l_TransferPNumber,1,MPI_LONG,i,i,MPI_COMM_WORLD,&status);         // Receive how many particles' information will be transfered   l_TransferPNumber
            MPI_Recv(writedata,l_TransferPNumber*6,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);         // Receive particles' information   writedata
			for(k=0;k<l_TransferPNumber;k++){
				for(j=0;j<6;j++)
					fprintf(FBeamdistribution,"%15.5e ",writedata[k*6+j]);
				fprintf(FBeamdistribution,"\n");
			}
		}
		fclose(FBeamdistribution);
	}
	else{
        MPI_Recv(&l_TransferPNumber,1,MPI_LONG,0,0,MPI_COMM_WORLD,&status);            // get the signal from process 0
		MPI_Send(&pNumber,1,MPI_LONG,0,i_my_rank,MPI_COMM_WORLD);              // Send how many particles' information will be transfered   pNumber
		for(k=0;k<pNumber;k++){
			    gamma=1.0/sqrt(1-pa[k]->vx[l_StepOrder]*pa[k]->vx[l_StepOrder]-pa[k]->vy[l_StepOrder]*pa[k]->vy[l_StepOrder]-pa[k]->vz[l_StepOrder]*pa[k]->vz[l_StepOrder]);
                writedata[k*6+0]=pa[k]->x[l_StepOrder];
				writedata[k*6+1]=pa[k]->y[l_StepOrder];
				writedata[k*6+2]=pa[k]->z[l_StepOrder];
				writedata[k*6+3]=pa[k]->vx[l_StepOrder]*gamma;
				writedata[k*6+4]=pa[k]->vy[l_StepOrder]*gamma;
				writedata[k*6+5]=pa[k]->vz[l_StepOrder]*gamma;
		}
		MPI_Send(writedata,pNumber*6,MPI_DOUBLE,0,i_my_rank,MPI_COMM_WORLD);              // Send particles' information     writedata
	}
}

void Domain::WriteBeamDistribution_Parallel(long l_StepOrder){
    double gamma, *writedata;
	char fname[200];
	long i,j,k;
	writedata= new double[6*pNumber];
    for(k=0;k<pNumber;k++){
		gamma=1.0/sqrt(1-pa[k]->vx[l_StepOrder]*pa[k]->vx[l_StepOrder]-pa[k]->vy[l_StepOrder]*pa[k]->vy[l_StepOrder]-pa[k]->vz[l_StepOrder]*pa[k]->vz[l_StepOrder]);
        writedata[k*6+0]=pa[k]->x[l_StepOrder];
		writedata[k*6+1]=pa[k]->y[l_StepOrder];
		writedata[k*6+2]=pa[k]->z[l_StepOrder];
		writedata[k*6+3]=pa[k]->vx[l_StepOrder]*gamma;
		writedata[k*6+4]=pa[k]->vy[l_StepOrder]*gamma;
		writedata[k*6+5]=pa[k]->vz[l_StepOrder]*gamma;
	}    

	hid_t       file, fdataset;                     /* File and dataset            */
	hid_t       fdataspace,memspace;                /* Dataspace handles           */
	hsize_t     dimsf[2], dimsfmem[2];              /* Dataset dimensions          */
	herr_t      status;                             /* Error checking              */
	hid_t rank1 = 1;
	hid_t rank2 = 2;
	hid_t rank3 = 3;
	hsize_t     count[2],offset[2],count_out[2],offset_out[2];   

	hid_t acc_tpl1;
	herr_t ret;


    dimsf[1] = 6;
    dimsfmem[1] = 6;

    dimsf[0] = l_Pnumber;
    dimsfmem[0] = pNumber;
   
    /* Creating the dataspace                                         */
    fdataspace = H5Screate_simple(rank2, dimsf, NULL); 
    assert (fdataspace >= 0);

    memspace = H5Screate_simple(rank2, dimsfmem, NULL); 
    assert (memspace >= 0);

	offset[1] = 0; 
    count[1] = 6;  
	offset_out[1]=0;
	count_out[1]=6;

    offset[0] = Beam::p_Beam->l_start_position; 
    count[0] = pNumber;  
	offset_out[0]=0;
	count_out[0]=pNumber;

	//Open the file and the dataset.

	acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
	assert(acc_tpl1>=0);
#ifdef PARALLELHDF5
    ret = H5Pset_fapl_mpio(acc_tpl1,MPI_COMM_WORLD,MPI_INFO_NULL);
    assert(ret >=0);	
#endif

    sprintf(fname,"BeamDis_at_%d.h5",l_StepOrder);
	file = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT,acc_tpl1);
    if (file <= 0 && i_my_rank==0) {
		 Domain::p_D->out_Flog << "BeamDis cannot be created: " << fname << "\n";
         Domain::p_D->out_Flog.flush();
		 cout << "BeamDis cannot be created: " << fname << "\n";
    }
	else{
		if(i_my_rank==0){
		    //cout<<"my_rank="<<i_my_rank<<"  local_N_Start="<<local_N_Start<<"  local_N_End="<<local_N_End<<" ExternalFieldsOutputFieldDis begin!"<<endl;
			Domain::p_D->out_Flog << "BeamDis has been created: " << fname << "\n";
			Domain::p_D->out_Flog.flush();
		}
	}
	assert (file >= 0);
#ifdef PARALLELHDF5
	ret = H5Pclose(acc_tpl1);
#endif

    fdataspace = H5Screate_simple(rank2, dimsf, NULL); 
    assert (fdataspace >= 0);

    memspace = H5Screate_simple(rank2, dimsfmem, NULL); 
    assert (memspace >= 0);


	// Creating the dataset within the dataspace =======================                      
	fdataset = H5Dcreate(file,"BeamDis", H5T_IEEE_F64LE,fdataspace,H5P_DEFAULT);
	assert (fdataset >= 0);
  

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, writedata);      
    assert (status >= 0);


	status = H5Dclose(fdataset);
	assert (status >= 0);
	status = H5Sclose(fdataspace);
	assert (status >= 0);
	status = H5Sclose(memspace);
	assert (status >= 0);
    status = H5Fclose(file);
	assert (status >= 0);
}

void Domain::OutputFieldDistribution(int i_OutputField,long l_StepOrder){
	char fname[200];
	int written,ii;
	double x_min,y_min,z_min,dx,dy,dz,t,f_laser_center,v_inchannel,v_invacuum,t1;
	double *x=new double[l_Field_Nx];
	double *y=new double[l_Field_Ny];
	double *z=new double[l_Field_Nz];
	long *l_local_N=new long[i_numprocs];
	long l_1=0;
	MPI_Status M_status;

	hid_t       file, fdataset;                   /* File and dataset            */
	hid_t       fdataspace;                       /* Dataspace handles           */
	hsize_t     dimsf[3], dimsfi[3];              /* Dataset dimensions          */
	herr_t      status;                           /* Error checking              */
	hid_t rank1 = 1;
	hid_t rank2 = 2;
	hid_t rank3 = 3;
	
    t=l_StepOrder*f_dt;
	v_inchannel=1-0.5*(ExternalFields::p_ExternalFields->f_Density_Plasma+1.0/(PI*ExternalFields::p_ExternalFields->f_Laser_FocusW)/(PI*ExternalFields::p_ExternalFields->f_Laser_FocusW));
	v_invacuum=1-0.5*(1.0/(PI*ExternalFields::p_ExternalFields->f_Laser_FocusW)/(PI*ExternalFields::p_ExternalFields->f_Laser_FocusW));
	t1=(ExternalFields::p_ExternalFields->f_ExternalField_End-ExternalFields::p_ExternalFields->f_Laser_LtoFocus)/v_inchannel;
	if(t<t1)
		f_laser_center=v_inchannel*t+ExternalFields::p_ExternalFields->f_Laser_LtoFocus;
	else
        f_laser_center=ExternalFields::p_ExternalFields->f_ExternalField_End+(t-t1)*v_invacuum;

	if(i_OutputField ==1){
		x_min=-f_xLength/2.0;
		y_min=-f_yLength/2.0;
		z_min=-f_zLength/2.0+f_laser_center;
	}
	else if(i_OutputField ==2){
		x_min=-f_xLength/2.0+f_x0Output;
		y_min=-f_yLength/2.0+f_y0Output;
		z_min=-f_zLength/2.0+f_z0Output;
	};
	if(l_Field_Nx>1)  
		dx=f_xLength/(l_Field_Nx-1);
	else
		dx=0.0;
    if(l_Field_Ny>1)  
		dy=f_yLength/(l_Field_Ny-1);
	else
		dy=0.0;
	if(l_Field_Nz>1)  
		dz=f_zLength/(l_Field_Nz-1);
	else
		dz=0.0;
	for(long i=0;i<l_Field_Nx;i++)
		x[i]=x_min+i*dx;
	for(long i=0;i<l_Field_Ny;i++)
		y[i]=y_min+i*dy;
	for(long i=0;i<l_Field_Nz;i++)
		z[i]=z_min+i*dz;

	for(ii=0;ii<i_numprocs;ii++)
		l_local_N[ii]=0;
    ii=0;l_1=0;
	while(l_1<l_Field_Nx*l_Field_Ny*l_Field_Nz){
		l_local_N[ii%i_numprocs]++;
		ii++;
		l_1++;
	}                                      // here records the number processed by each CPU
	for(int ii=1;ii<i_numprocs;ii++)
		l_local_N[ii]+=l_local_N[ii-1];  
	for(int ii=i_numprocs-1;ii>=0;ii--){
		l_local_N[ii]-=l_local_N[0];      // here recordes the start number processed by each CPU
	}
	if(i_my_rank==0){
		cout<<"FieldDistrbution is writting into file: FieldDis_at_"<<l_StepOrder<<".h5"<<endl;
		Domain::p_D->out_Flog<<"FieldDistrbution is writting into file: FieldDis_at_"<<l_StepOrder<<".h5"<<endl;
   
		sprintf(fname,"FieldDis_at_%d.h5",l_StepOrder);
		Domain::p_D->out_Flog << "SAVE FIELD DISTRIBUTION: Opening file " << fname << "\n";
		Domain::p_D->out_Flog.flush();

		double *fdata = new double[l_Field_Nx*l_Field_Ny*l_Field_Nz];

		file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		assert (file >= 0);

		dimsfi[0] = 1;
		written = WriteHDF5RecordDouble(file, "Time", dimsfi[0], &t);
		assert (written >= 0);

		dimsfi[0] = l_Field_Nx;
		written = WriteHDF5RecordDouble(file, "x", dimsfi[0], x);
		assert (written >= 0);
		
		dimsfi[0] = l_Field_Ny;
		written = WriteHDF5RecordDouble(file, "y", dimsfi[0], y);
		assert (written >= 0);

		dimsfi[0] = l_Field_Nz;
		written = WriteHDF5RecordDouble(file, "z", dimsfi[0], z);
		assert (written >= 0);

	    dimsfi[0] = 1;
		written = WriteHDF5RecordLong(file, "Nx", dimsfi[0], &l_Field_Nx);
		assert (written >= 0);	    
		dimsfi[0] = 1;
		written = WriteHDF5RecordLong(file, "Ny", dimsfi[0], &l_Field_Ny);
		assert (written >= 0);
		dimsfi[0] = 1;
		written = WriteHDF5RecordLong(file, "Nz", dimsfi[0], &l_Field_Nz);
		assert (written >= 0);

		dimsf[0] = l_Field_Nx*l_Field_Ny*l_Field_Nz;  


		// Creating the dataspace                                         
		fdataspace = H5Screate_simple(rank1, dimsf, NULL); 
		assert (fdataspace >= 0);

		// Creating the dataset within the dataspace =======================                      
		fdataset = H5Dcreate(file, "Ex", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
		assert (fdataset >= 0);
		status = H5Dclose(fdataset);
		assert (status >= 0);

		// Creating the dataset within the dataspace =======================                      
		fdataset = H5Dcreate(file, "Ez", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
		assert (fdataset >= 0);
		status = H5Dclose(fdataset);
		assert (status >= 0);

		// Creating the dataset within the dataspace =======================                      
		fdataset = H5Dcreate(file, "By", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
		assert (fdataset >= 0);
		status = H5Dclose(fdataset);
		assert (status >= 0);

		// Close the datefile                                            
		status = H5Fclose(file);
		assert (status >= 0);

		Domain::p_D->out_Flog << "SAVE FIELD DISTRIBUTION: Closing file " << fname << "\n";
		Domain::p_D->out_Flog.flush();

	}

	MPI_Barrier(MPI_COMM_WORLD);
    //------------------------------ Here parallelly calculates the Calculated output field ----------------------
	if(i_my_rank<i_numprocs-1 && l_local_N[i_my_rank+1]-l_local_N[i_my_rank]!=0)
		ExternalFields::p_ExternalFields->Cal_OutField(l_StepOrder, t, x_min, y_min, z_min, dx, dy, dz, l_Field_Nx, l_Field_Ny, l_Field_Nz, l_local_N[i_my_rank], l_local_N[i_my_rank+1]-1);
	if(i_my_rank==i_numprocs-1 && l_Field_Nx*l_Field_Ny*l_Field_Nz-l_local_N[i_my_rank]!=0)
		ExternalFields::p_ExternalFields->Cal_OutField(l_StepOrder, t, x_min, y_min, z_min, dx, dy, dz, l_Field_Nx, l_Field_Ny, l_Field_Nz, l_local_N[i_my_rank],l_Field_Nx*l_Field_Ny*l_Field_Nz-1);

#ifdef PARALLELHDF5
	//------------------------------ Here parallelly outputs the Calculated output field ----------------------
	if(i_my_rank<i_numprocs-1 && l_local_N[i_my_rank+1]-l_local_N[i_my_rank]!=0)
		ExternalFields::p_ExternalFields->Write_Field_Parallel(l_StepOrder, t, x_min, y_min, z_min, dx, dy, dz, l_Field_Nx, l_Field_Ny, l_Field_Nz, l_local_N[i_my_rank], l_local_N[i_my_rank+1]-1);
	if(i_my_rank==i_numprocs-1 && l_Field_Nx*l_Field_Ny*l_Field_Nz-l_local_N[i_my_rank]!=0)
		ExternalFields::p_ExternalFields->Write_Field_Parallel(l_StepOrder, t, x_min, y_min, z_min, dx, dy, dz, l_Field_Nx, l_Field_Ny, l_Field_Nz, l_local_N[i_my_rank],l_Field_Nx*l_Field_Ny*l_Field_Nz-1);
#else
	//------------------------------ Here outputs the Calculated output field one by one ----------------------
	long l_Message=1;
	if(i_my_rank>0){
		MPI_Recv(&l_Message,1,MPI_LONG,i_my_rank-1,i_my_rank-1,MPI_COMM_WORLD,&M_status);
	}
	if(i_my_rank<i_numprocs-1 && l_local_N[i_my_rank+1]-l_local_N[i_my_rank]!=0)
		ExternalFields::p_ExternalFields->Write_Field(l_StepOrder, t, x_min, y_min, z_min, dx, dy, dz, l_Field_Nx, l_Field_Ny, l_Field_Nz, l_local_N[i_my_rank], l_local_N[i_my_rank+1]-1);
	if(i_my_rank==i_numprocs-1 && l_Field_Nx*l_Field_Ny*l_Field_Nz-l_local_N[i_my_rank]!=0)
		ExternalFields::p_ExternalFields->Write_Field(l_StepOrder, t, x_min, y_min, z_min, dx, dy, dz, l_Field_Nx, l_Field_Ny, l_Field_Nz, l_local_N[i_my_rank],l_Field_Nx*l_Field_Ny*l_Field_Nz-1);
	if(i_my_rank<i_numprocs-1){
		MPI_Send(&l_Message,1,MPI_LONG,i_my_rank+1,i_my_rank,MPI_COMM_WORLD);                           
	}
#endif

}

void Domain::OutputBeamDiagnosis(long l_StepOrder){
	;
}

void Domain::OutputFieldEvolution(int i_my_rank){
	char FileName[200];
	FILE *FieldEvolution;
	sprintf(FileName,"FieldTempEvolution.dat");
	if(i_my_rank==0){
		double Ex,Ey,Ez,Bx,By,Bz,t;
		FieldEvolution=fopen(FileName,"w");
		for(long i=0;i<l_Step;i++){
			t=i*f_dt;
			ExternalFields::p_ExternalFields->Set_Field(t,f_xPosition,f_yPosition,f_zPosition,Ex,Ey,Ez,Bx,By,Bz);
			fprintf(FieldEvolution,"%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n",t,Ex,Ey,Ez,Bx,By,Bz);
		}
		fclose(FieldEvolution);
	}
}

//---------------------------- Domain::~Domain() -----------------------
Domain::~Domain()
{
   if(out_Flog)
      out_Flog.close();

   if(p_File)
      fclose(p_File);

   if(p_MovieFile)
      fclose(p_MovieFile);
};


