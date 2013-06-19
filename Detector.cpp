#include <iostream>
using namespace std;

#include <stdlib.h>
#include <math.h>
#include "myhdfshell.h"
#include "vdsr.h"

Detector* Detector::p_Detector = NULL;
Detector::Detector(char *infile, int rank): NList("Detector")
{

	Detector::p_Detector = this;
	i_my_rank=rank;
	AddEntry("f_Omegamax", &f_omega_max, 1.0e3);
	AddEntry("f_Omegamin", &f_omega_min, 0.1);
	AddEntry("i_OmegaLog", &i_OmegaLog, 1);
	AddEntry("l_nOmegabins", &l_nEbins, 100);
	AddEntry("f_Thetamax", &f_theta_max, 180.);
	AddEntry("f_Thetamin", &f_theta_min, 0.);
	AddEntry("l_nThetabins", &l_nThetabins, 100);
    AddEntry("f_Phimax", &f_phi_max, 360.);
	AddEntry("f_Phimin", &f_phi_min, 0.);
	AddEntry("l_nPhibins", &l_nPhibins, 100);
	AddEntry("f_Rdistance", &f_R_dp, 100.);
	AddEntry("f_Xlength", &f_Xlength, 1.);
	AddEntry("f_Ylength", &f_Ylength, 1.);
	AddEntry("f_DCenterX", &f_DCenterX, 0.);
	AddEntry("f_DCenterY", &f_DCenterY, 0.);
	AddEntry("i_RealDetector", &i_RealDetector, 0);
	AddEntry("i_InitialEndIntegration", &i_InitialEndIntegration, 0);


	FILE *p_File = NULL;


   if (rank==NIL)
   {
      p_File = fopen(infile,"rt");     
      if (p_File == NULL)
      {
		  cout<<"Detector no input file: "<<infile<<" "<<endl;
		  Domain::p_D->out_Flog<<"Detector no input file: "<<infile<<" "<<endl;
          exit (-1);
      }
   }

   if (p_File)
   {
      rewind(p_File); 
      read(p_File); 
   }
   if(i_my_rank==0){
	    cout<<"Detector: Input file "<<infile<<" is found!"<<endl;
		cout<<"l_nPhibins="<<l_nPhibins<<" l_nEbins="<<l_nEbins<<" iThetabins="<<l_nThetabins<<endl;
        cout<<"Theta_min="<<f_theta_min<<" Theta_max="<<f_theta_max<<" phi_min="<<f_phi_min<<" phi_max="<<f_phi_max<<" unit of Degree"<<endl;
		Domain::p_D->out_Flog<<"Detector: Input file "<<infile<<" is found!"<<endl;
		Domain::p_D->out_Flog<<"l_nPhibins="<<l_nPhibins<<" l_nEbins="<<l_nEbins<<" iThetabins="<<l_nThetabins<<endl;
        Domain::p_D->out_Flog<<"Theta_min="<<f_theta_min<<" Theta_max="<<f_theta_max<<" phi_min="<<f_phi_min<<" phi_max="<<f_phi_max<<" unit of Degree"<<endl;
   }
#ifdef V_MPI
   CBuffer *buf = new CBuffer;
   buf->reset();
   pack_nls(buf);
   Domain::p_D->BroadCast(buf);
   if (rank) unpack_nls(buf);
   delete buf;
#endif
   f_total_emission_energy=f_receive_emission_energy=0.0;
   if(i_RealDetector==0){
	   f_theta_max=f_theta_max/180.0*3.1415926;
	   f_theta_min=f_theta_min/180.0*3.1415926;
       f_phi_max=f_phi_max/180.0*3.1415926;
       f_phi_min=f_phi_min/180.0*3.1415926;
   }
   else{
	   f_theta_max=(f_DCenterX+f_Xlength/2.0)/f_R_dp*3.1415926;
	   f_theta_max=(f_DCenterX-f_Xlength/2.0)/f_R_dp*3.1415926;
	   f_phi_max=(f_DCenterY+f_Ylength/2.0)*3.1415926;
	   f_phi_max=(f_DCenterY-f_Ylength/2.0)*3.1415926;
   }
   if(l_nThetabins ==1) {f_dtheta=f_theta_max-f_theta_min;}
   else {f_dtheta=(f_theta_max-f_theta_min)/(l_nThetabins-1);}
   if(l_nPhibins == 1) {f_dphi=f_phi_max-f_phi_min;}
   else{f_dphi=(f_phi_max-f_phi_min)/(l_nPhibins-1);}
   l_Pixel_number=l_nThetabins*l_nPhibins;
   Pixel_Array= new Pixel*[l_Pixel_number];
   for(long i=0;i<l_nPhibins;i++)                                 // the pixel is first increased in the theta direction
	   for(long j=0;j<l_nThetabins;j++)
		Pixel_Array[i*l_nThetabins+j]= new Pixel(f_theta_min+j*f_dtheta,f_phi_min+i*f_dphi);
   	f_Omega=new double[l_nEbins];
	if (l_nEbins == 1) {
		f_Omega[0]=f_omega_min;
		f_domega=f_omega_max-f_omega_min;
	}
	else{
		if (i_OmegaLog==1){
			f_domega=(log(f_omega_max)-log(f_omega_min))/log(10.0)/(l_nEbins-1);
			for (long i=0;i<l_nEbins;i++){
				f_Omega[i]=f_omega_min*pow(10.0,i*f_domega);
			}
		}
		else{
			f_domega=(f_omega_max-f_omega_min)/(l_nEbins-1);
			for (long i=0;i<l_nEbins;i++){
				f_Omega[i]=f_omega_min+i*f_domega;  
			}
		}
	}

    p_UnitsCGS= new UnitsCGS();

	if(i_my_rank==0){
		cout<<"PhotonPlasma="<<p_UnitsCGS->f_PhotonPlasma_eV<<"eV."<<endl;
		cout<<"***** Detector has been initialized! *****\n"<<endl;
		Domain::p_D->out_Flog<<"PhotonPlasma="<<p_UnitsCGS->f_PhotonPlasma_eV<<"eV."<<endl;
		Domain::p_D->out_Flog<<"***** Detector has been initialized! *****\n"<<endl;
		
	}
	l_Emission_number=l_Recorder_number=0;
};


void Detector::Run(int i_Radiation_Type)  // i_Radiation_Type decides the method for the radiation calculation
{   
	if(i_my_rank==0){
		cout<<"Detector is working ..."<<endl<<"Total Step numbers:"<<Domain::p_D->l_Step<<endl;
		Domain::p_D->out_Flog<<"Detector is working ..."<<endl<<"Total Step numbers:"<<Domain::p_D->l_Step<<endl;
	}

	for(long j=0;j<Domain::p_D->pNumber;j++){            // for every particle do radiation calculation
		if(i_my_rank==0){
			cout<<"--- Particle loop: "<<1.0*j/Domain::p_D->pNumber*100<<"% calculation finished!"<<endl;
			Domain::p_D->out_Flog<<"--- Particle loop: "<<1.0*j/Domain::p_D->pNumber*100<<"% calculation finished!"<<endl;
		}
	
		for(long i=0;i<l_Pixel_number;i++)
			Pixel_Array[i]->CleanImRe_IntegrationXYZ();

        //********************************** Following for the radiation calculation *****************************
		if(i_Radiation_Type==1){                        // Classical Radiation Calculation method
			if(i_InitialEndIntegration==1){
				if(i_my_rank==0) {
					cout<<"InitialEndModification has been used!"<<endl;
					Domain::p_D->out_Flog<<"InitialEndModification has been used!"<<endl;
					;
				}
				for(long i=0;i<l_Pixel_number;i++)                                              
					Pixel_Array[i]->InitialEndModification(Domain::p_D->pa[j]);	
			}

			/*
			Projector_Rotation();	
			for(long i=0;i<l_Pixel_number;i++)
				Detector::p_Detector->Pixel_Array[i]->GetIntensity();
			for(long i=0;i<l_Pixel_number;i++)
				Pixel_Array[i]->CleanImRe_IntegrationXYZ();
			*/
	
			for(long k=0;k<Domain::p_D->l_Step-1;k++){                            // sum over time step
				
				if(i_my_rank==0 && k%(Domain::p_D->l_Step/10)==0){
					cout<<"     Step loop: "<<1.0*k/Domain::p_D->l_Step*100<<"% calculation finished!"<<endl;
					Domain::p_D->out_Flog<<"     Step loop: "<<1.0*k/Domain::p_D->l_Step*100<<"% calculation finished!"<<endl;
				}
				
				Domain::p_D->pa[j]->Interpolation(k); 
				for(long i=0;i<l_Pixel_number;i++)                            // calculate for this particle at this time at different pixels
					Pixel_Array[i]->ProjectSinglePa(Domain::p_D->pa[j],k);	
			} // finished integration on time and particles
    
			Projector_Rotation();
			for(long i=0;i<l_Pixel_number;i++){
				Detector::p_Detector->Pixel_Array[i]->GetIntensity();
			}
		}    // end of i_Radiation_Type==1

		if(i_Radiation_Type==2){          // Monte Carlo Compton Scattering
			for(long k=0;k<Domain::p_D->l_Step-1;k++){                            // sum over time step
				
				if(i_my_rank==0 && k%(Domain::p_D->l_Step/25)==0){
					cout<<"     Step loop: "<<1.0*k/Domain::p_D->l_Step*100<<"% calculation finished!"<<endl;
					Domain::p_D->out_Flog<<"     Step loop: "<<1.0*k/Domain::p_D->l_Step*100<<"% calculation finished!"<<endl;
				}
								
				Domain::p_D->pa[j]->CP_MC_ComptonRadiation(k);
			} 

			Domain::p_D->pa[j]->CP_MC_ProjectSinglePhoton();  // at the end of particle trajectory, project all the photons to the pixels
			cout<<"Particle "<<j<<" has emitted "<<Domain::p_D->pa[j]->l_emission_number<<" photons, and "<<Domain::p_D->pa[j]->l_Recorder<<" has been recorded!"<<endl;
			l_Emission_number+=Domain::p_D->pa[j]->l_emission_number;
			l_Recorder_number+=Domain::p_D->pa[j]->l_Recorder;


			//Domain::p_D->pa[j]->Write_XP(100); // temporal here for debug output particle trajectory
		}  // end of i_Radiation_Type==2

		if(i_Radiation_Type==3){          // Monte Carlo Nonlinear Compton Scattering
			for(long k=0;k<Domain::p_D->l_Step-1;k++){                            // sum over time step
				
				if(i_my_rank==0 && k%(Domain::p_D->l_Step/25)==0){
					cout<<"     Step loop: "<<1.0*k/Domain::p_D->l_Step*100<<"% calculation finished!"<<endl;
					Domain::p_D->out_Flog<<"     Step loop: "<<1.0*k/Domain::p_D->l_Step*100<<"% calculation finished!"<<endl;
				}
								
				Domain::p_D->pa[j]->CP_LBNL_ComptonRadiation(k);
			} 

			Domain::p_D->pa[j]->CP_LBNL_ProjectSinglePhoton();  // at the end of particle trajectory, project all the photons to the pixels
			cout<<"Particle "<<j<<" has emitted "<<Domain::p_D->pa[j]->l_emission_number<<" photons, and "<<Domain::p_D->pa[j]->l_Recorder<<" has been recorded!"<<endl;
			l_Emission_number+=Domain::p_D->pa[j]->l_emission_number;
			l_Recorder_number+=Domain::p_D->pa[j]->l_Recorder;


			//Domain::p_D->pa[j]->Write_XP(100); // temporal here for debug output particle trajectory
		}  // end of i_Radiation_Type==3


		// *********************************** end of Radiation Calculation for particle j *********************************
		cout<<"End of radiation Calculation for Particle "<<j<<endl;
	}

#ifdef V_MPI
	Global_Intensity_reduction();
#endif

	if(i_my_rank==0){
		cout<<"***** Detector has been switched off *****\n"<<endl;
		Domain::p_D->out_Flog<<"***** Detector has been switched off *****\n"<<endl;
		if(i_Radiation_Type==2){
			cout<<"Total emission energy is:"<<f_total_emission_energy<<" eV, received emission energy is:"<<f_receive_emission_energy<<" eV"<<endl;
			cout<<"Totally "<<l_Emission_number_tot<<" photons have been emitted, and "<<l_Recorder_number_tot<<" of them have been recordered by this Detector!"<<endl;
			cout<<"In average each macro electron emits "<<l_Emission_number_tot*1.0/Domain::p_D->l_Pnumber<<" photons"<<endl;
			cout<<"In average each real electron emits "<<l_Emission_number_tot*1.0/Domain::p_D->l_Pnumber/Domain::p_D->pa[0]->Charge2Number_weight<<" photons"<<endl;
			cout<<"Classically in average each real electron emits N=alpha*N_beta*a_beta^2 ~= "<<Domain::p_D->pa[0]->f_averageemissionnumber_classical<<" photons"<<endl;
			cout<<"Here the artificial emission amplify factor is "<<Domain::p_D->f_Amplify_Factor<<endl;
		}
	}
}




void Detector::Projector_Rotation(){
	for(long i=0;i<l_Pixel_number;i++){      // here finish the         n X n X ...      term                                 process for different pixels
		double Ndirection[3],fIntensity[3];

	    Ndirection[0]=sin(Pixel_Array[i]->theta)*cos(Pixel_Array[i]->phi);
		Ndirection[1]=sin(Pixel_Array[i]->theta)*sin(Pixel_Array[i]->phi);
		Ndirection[2]=cos(Pixel_Array[i]->theta);
        
		for(long ii=0;ii<l_nEbins;ii++){
			fIntensity[0]=Ndirection[1]*Pixel_Array[i]->Re_IntegrationZ[ii]-Ndirection[2]*Pixel_Array[i]->Re_IntegrationY[ii];        //Ix=ny*Iz-Nz*Iy
			fIntensity[1]=Ndirection[2]*Pixel_Array[i]->Re_IntegrationX[ii]-Ndirection[0]*Pixel_Array[i]->Re_IntegrationZ[ii];        //Iy=Nz*Ix-Nx*Iz
			fIntensity[2]=Ndirection[0]*Pixel_Array[i]->Re_IntegrationY[ii]-Ndirection[1]*Pixel_Array[i]->Re_IntegrationX[ii];        //Iz=Nx*Iy-Ny*Ix
			Pixel_Array[i]->Re_IntegrationX[ii]=fIntensity[0];
			Pixel_Array[i]->Re_IntegrationY[ii]=fIntensity[1];
			Pixel_Array[i]->Re_IntegrationZ[ii]=fIntensity[2];

			fIntensity[0]=Ndirection[1]*Pixel_Array[i]->Im_IntegrationZ[ii]-Ndirection[2]*Pixel_Array[i]->Im_IntegrationY[ii];
			fIntensity[1]=Ndirection[2]*Pixel_Array[i]->Im_IntegrationX[ii]-Ndirection[0]*Pixel_Array[i]->Im_IntegrationZ[ii];
			fIntensity[2]=Ndirection[0]*Pixel_Array[i]->Im_IntegrationY[ii]-Ndirection[1]*Pixel_Array[i]->Im_IntegrationX[ii];
			Pixel_Array[i]->Im_IntegrationX[ii]=fIntensity[0];
			Pixel_Array[i]->Im_IntegrationY[ii]=fIntensity[1];
			Pixel_Array[i]->Im_IntegrationZ[ii]=fIntensity[2];
		}

		for(long ii=0;ii<l_nEbins;ii++){
			fIntensity[0]=Ndirection[1]*Pixel_Array[i]->Re_IntegrationZ[ii]-Ndirection[2]*Pixel_Array[i]->Re_IntegrationY[ii];        //Ix=Ny*Iz-Nz*Iy
			fIntensity[1]=Ndirection[2]*Pixel_Array[i]->Re_IntegrationX[ii]-Ndirection[0]*Pixel_Array[i]->Re_IntegrationZ[ii];        //Iy=Nz*Ix-Nx*Iz
			fIntensity[2]=Ndirection[0]*Pixel_Array[i]->Re_IntegrationY[ii]-Ndirection[1]*Pixel_Array[i]->Re_IntegrationX[ii];        //Iz=Nx*Iy-Ny*Ix
			Pixel_Array[i]->Re_IntegrationX[ii]=fIntensity[0];
			Pixel_Array[i]->Re_IntegrationY[ii]=fIntensity[1];
			Pixel_Array[i]->Re_IntegrationZ[ii]=fIntensity[2];

			fIntensity[0]=Ndirection[1]*Pixel_Array[i]->Im_IntegrationZ[ii]-Ndirection[2]*Pixel_Array[i]->Im_IntegrationY[ii];
			fIntensity[1]=Ndirection[2]*Pixel_Array[i]->Im_IntegrationX[ii]-Ndirection[0]*Pixel_Array[i]->Im_IntegrationZ[ii];
			fIntensity[2]=Ndirection[0]*Pixel_Array[i]->Im_IntegrationY[ii]-Ndirection[1]*Pixel_Array[i]->Im_IntegrationX[ii];
			Pixel_Array[i]->Im_IntegrationX[ii]=fIntensity[0];
			Pixel_Array[i]->Im_IntegrationY[ii]=fIntensity[1];
			Pixel_Array[i]->Im_IntegrationZ[ii]=fIntensity[2];
		}
	}
}

void Detector::WriteRadiation(){
   char fname[200];
   long ldump=l_nPhibins*l_nThetabins*l_nEbins;
   int written = 0;
   
   f_PhotonPlasma_eV=p_UnitsCGS->f_PhotonPlasma_eV;

   
   if(i_my_rank==0){
		cout<<"Radiation is writting into a file!"<<endl;
		Domain::p_D->out_Flog<<"Radiation is writting into a file!"<<endl;
   }
   hid_t       file, fdataset;         /* File and dataset            */
   hid_t       fdataspace;   /* Dataspace handles           */
   hsize_t     dimsf[3], dimsfi[3];              /* Dataset dimensions          */
   herr_t      status;                /* Error checking              */
   hid_t rank1 = 1;
   hid_t rank2 = 2;
   hid_t rank3 = 3;

   sprintf(fname,"vdsr_synchrotron.h5");
   Domain::p_D->out_Flog << "SAVE SYNCHROTRON: Opening file " << fname << "\n";
   Domain::p_D->out_Flog.flush();



   double *fdata = new double[l_nEbins*l_nThetabins*l_nPhibins];


   file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   assert (file >= 0);

   dimsfi[0] = 1;
   written = WriteHDF5RecordLong(file, "ParticleNumber", dimsfi[0], &(Domain::p_D->l_Pnumber));
   assert (written >= 0);
   written = WriteHDF5RecordDouble(file, "Omegamin", dimsfi[0], &f_omega_min);
   assert (written >= 0);
   written = WriteHDF5RecordDouble(file, "Omegamax", dimsfi[0], &f_omega_max);
   assert (written >= 0);
   written = WriteHDF5RecordDouble(file, "Phimax", dimsfi[0], &f_phi_max);
   assert (written >= 0);
   written = WriteHDF5RecordDouble(file, "Phimin", dimsfi[0], &f_phi_min);
   assert (written >= 0);
   written = WriteHDF5RecordDouble(file, "Thetamax", dimsfi[0], &f_theta_max); 
   assert (written >= 0);
   written = WriteHDF5RecordDouble(file, "Thetamin", dimsfi[0], &f_theta_min);
   assert (written >= 0);
   written = WriteHDF5RecordLong(file, "l_nEbins", dimsfi[0], &l_nEbins);
   assert (written >= 0);
   written = WriteHDF5RecordLong(file, "l_nThetabins", dimsfi[0], &l_nThetabins);
   assert (written >= 0);
   written = WriteHDF5RecordLong(file, "l_nPhibins", dimsfi[0], &l_nPhibins);
   assert (written >= 0);

   //============ E axis ============
   dimsfi[0] = l_nEbins;
   double *Xaxis = new double[dimsfi[0]];
   for (long i=0; i<dimsfi[0]; i++) {
      Xaxis[i] = f_Omega[i]*f_PhotonPlasma_eV;
   }
   written = WriteHDF5RecordDouble(file, "Frequency_eV", dimsfi[0], Xaxis);
   assert (written >= 0);

   delete[] Xaxis;

   //============ Theta axis ============
   dimsfi[0] = l_nThetabins;
   double *Yaxis = new double[dimsfi[0]];
   for (long i=0; i<dimsfi[0]; i++) {
      Yaxis[i] = f_theta_min+i*f_dtheta;
   }
   written = WriteHDF5RecordDouble(file, "Theta", dimsfi[0], Yaxis);
   assert (written >= 0);
   delete[] Yaxis;

   //============ Z axis ============

   dimsfi[0] = l_nPhibins;
   double *Zaxis = new double[dimsfi[0]];
   for (long i=0; i<dimsfi[0]; i++) {
      Zaxis[i] = f_phi_min+f_dphi*i;
   }
   written = WriteHDF5RecordDouble(file, "Phi", dimsfi[0], Zaxis);
   assert (written >= 0);
   delete[] Zaxis;

   //============== Output the radiation data ======================
   dimsf[0] = l_nPhibins;  
   dimsf[1] = l_nThetabins;
   dimsf[2] = l_nEbins;


   // Creating the dataspace                                         
   fdataspace = H5Screate_simple(rank3, dimsf, NULL); 
   assert (fdataspace >= 0);

   // Creating the dataset within the dataspace =======================                      
   fdataset = H5Dcreate(file, "Synchrotron3DX", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
   assert (fdataset >= 0);

   for(long i=0;i<l_nPhibins;i++)
	   for(long j=0;j<l_nThetabins;j++)
		   for(long k=0;k<l_nEbins;k++){
			   fdata[i*l_nThetabins*l_nEbins+j*l_nEbins+k]=Pixel_Array[i*l_nThetabins+j]->ReturnIntensityX(k);
		   }

   // Writing the data to the dataset                                
   status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, fdata);                                                          // 090515 here the last work is originally fdata written by pukhov, it should be fdataAccumulated
   assert (status >= 0);

   status = H5Dclose(fdataset);
   assert (status >= 0);

    // Creating the dataset within the dataspace =======================                     
   fdataset = H5Dcreate(file, "Synchrotron3DY", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
   assert (fdataset >= 0);

   for(long i=0;i<l_nPhibins;i++)
	   for(long j=0;j<l_nThetabins;j++)
		   for(long k=0;k<l_nEbins;k++){
			   fdata[i*l_nThetabins*l_nEbins+j*l_nEbins+k]=Pixel_Array[i*l_nThetabins+j]->ReturnIntensityY(k);
		   }

   // Writing the data to the dataset                                
   status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, fdata);                                                          // 090515 here the last work is originally fdata written by pukhov, it should be fdataAccumulated
   assert (status >= 0);

   status = H5Dclose(fdataset);
   assert (status >= 0);



    // Creating the dataset within the dataspace  ==========================                    
   fdataset = H5Dcreate(file, "Synchrotron3DZ", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
   assert (fdataset >= 0);

   for(long i=0;i<l_nPhibins;i++)
	   for(long j=0;j<l_nThetabins;j++)
		   for(long k=0;k<l_nEbins;k++){
			   fdata[i*l_nThetabins*l_nEbins+j*l_nEbins+k]=Pixel_Array[i*l_nThetabins+j]->ReturnIntensityZ(k);
		   }

   // Writing the data to the dataset                                
   status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, fdata);                                                          // 090515 here the last work is originally fdata written by pukhov, it should be fdataAccumulated
   assert (status >= 0);

   status = H5Dclose(fdataset);
   assert (status >= 0);


     // Creating the dataset within the dataspace =======================                 this output is not related to the polarization     
   fdataset = H5Dcreate(file, "Synchrotron3D", H5T_IEEE_F64LE, fdataspace,H5P_DEFAULT);
   assert (fdataset >= 0);

   for(long i=0;i<l_nPhibins;i++)
	   for(long j=0;j<l_nThetabins;j++)
		   for(long k=0;k<l_nEbins;k++){
			   fdata[i*l_nThetabins*l_nEbins+j*l_nEbins+k]=Pixel_Array[i*l_nThetabins+j]->ReturnIntensityX(k)+Pixel_Array[i*l_nThetabins+j]->ReturnIntensityY(k)+Pixel_Array[i*l_nThetabins+j]->ReturnIntensityZ(k);
		   }

   // Writing the data to the dataset                                
   status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, fdata);                                                          // 090515 here the last work is originally fdata written by pukhov, it should be fdataAccumulated
   assert (status >= 0);

   status = H5Dclose(fdataset);
   assert (status >= 0);


   // Close the datefile                                            
   status = H5Fclose(file);
   assert (status >= 0);
  
   // Close the dataspace                                              
   status = H5Sclose(fdataspace);
   assert (status >= 0);

   delete[] fdata;
}

void Detector::Global_radiation_reduction(){
    double *Send_Data,*Receive_Data;
	long Send_bufsize=l_nEbins*l_nThetabins,Tcount=0,Tnumber=l_nPhibins*l_nThetabins*l_nEbins,Pixelnumber,omeganumber;
	int count=0;
	Send_Data=new double[Send_bufsize];
	Receive_Data=new double[Send_bufsize];
	if(i_my_rank==0){
		cout<<"Global_radiation_reduction begins..."<<endl;
		Domain::p_D->out_Flog<<"Global_radiation_reduction begins..."<<endl;
	}

	//----------------------------------- Send Im_IntegrationX
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->Im_IntegrationX[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){ 
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins; 
							Pixel_Array[Pixelnumber]->Im_IntegrationX[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------
	
		//----------------------------------- Send Im_IntegrationY
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->Im_IntegrationY[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins;
							Pixel_Array[Pixelnumber]->Im_IntegrationY[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------
    //----------------------------------- Send Im_IntegrationZ
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->Im_IntegrationZ[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins;
							Pixel_Array[Pixelnumber]->Im_IntegrationZ[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------


	//----------------------------------- Send Re_IntegrationX
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->Re_IntegrationX[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins;
							Pixel_Array[Pixelnumber]->Re_IntegrationX[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------
		//----------------------------------- Send Re_IntegrationY
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->Re_IntegrationY[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins;
							Pixel_Array[Pixelnumber]->Re_IntegrationY[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------
    //----------------------------------- Send Re_IntegrationZ
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->Re_IntegrationZ[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins;
							Pixel_Array[Pixelnumber]->Re_IntegrationZ[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------




	delete[] Send_Data;
	delete[] Receive_Data;
}



void Detector::Global_Intensity_reduction(){
    double *Send_Data,*Receive_Data;
	long Send_bufsize=l_nEbins*l_nThetabins,Tcount=0,Tnumber=l_nPhibins*l_nThetabins*l_nEbins,Pixelnumber,omeganumber;
	int count=0,ierr_temp;
	Send_Data=new double[Send_bufsize];
	Receive_Data=new double[Send_bufsize];
    if(i_my_rank==0)
		cout<<"Global_Intensity_reduction begins..."<<endl;

    ierr_temp=MPI_Reduce(&l_Emission_number,&l_Emission_number_tot,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    ierr_temp=MPI_Reduce(&l_Recorder_number,&l_Recorder_number_tot,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	//----------------------------------- Send IntensityX
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->IntensityX[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){ 
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins; 
							Pixel_Array[Pixelnumber]->IntensityX[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------
	
//----------------------------------- Send IntensityY
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->IntensityY[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){ 
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins; 
							Pixel_Array[Pixelnumber]->IntensityY[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------

	//----------------------------------- Send IntensityZ
	for(long i=0;i<l_nPhibins;i++)
		for(long j=0;j<l_nThetabins;j++)
			for(long k=0;k<l_nEbins;k++){
				Send_Data[count++]=Pixel_Array[i*l_nThetabins+j]->IntensityZ[k];
				Tcount++;
				if(count==Send_bufsize||Tcount==Tnumber){
					int ierr=MPI_Reduce(Send_Data,Receive_Data,count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(i_my_rank==0){
						Tcount-=count;
						for(long l=0;l<count;l++){ 
							omeganumber=Tcount%l_nEbins;
							Pixelnumber=Tcount/l_nEbins; 
							Pixel_Array[Pixelnumber]->IntensityZ[omeganumber]=Receive_Data[l];
							Tcount++;
						}
					}
					count=0;
				}
			}
	Tcount=0;
	//---------------------------



	delete[] Send_Data;
	delete[] Receive_Data;
}