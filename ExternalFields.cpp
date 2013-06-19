#include <iostream>
using namespace std;

#include <stdlib.h>
#include "vdsr.h"
#include <math.h>

ExternalFields* ExternalFields::p_ExternalFields = NULL;
ExternalFields::ExternalFields(char *infile, int rank): NList("ExternalFields")
{
	ExternalFields::p_ExternalFields = this;
   	i_my_rank=rank;
	
	AddEntry("f_Lambda_L", &f_Lambda_L0, 0.8);   // here in the unit of um
	AddEntry("f_Density_Normalize", &f_Density_Normalize, 1.0);
	AddEntry("f_Density_Plasma", &f_Density_Plasma, 0.01);
	AddEntry("i_Particle_Field", &i_Particle_Field, 1);
	AddEntry("f_Gamma_wake", &f_gamma_wake, 10.0);
	AddEntry("f_r_bubble", &f_r_bubble, 0.5);
	AddEntry("f_ExternalField_Intensity", &f_ExternalField_Intensity, 1.0);
	AddEntry("f_ExternalField_End", &f_ExternalField_End, 1.0e8);
	AddEntry("f_ExternalField_Duration", &f_ExternalField_Duration, 1.0);
	AddEntry("f_Undu_B", &f_Undu_B, 1.0);
	AddEntry("f_Undulator_a", &f_Undulator_a, 1.0);
	AddEntry("f_Undulator_b", &f_Undulator_b, 0.0);
	AddEntry("i_Particle_Trajectory", &i_Particle_Trajectory,0);
	AddEntry("i_Radiation_Reaction", &i_Radiation_Reaction,0);
	AddEntry("i_DampingFlag", &i_DampingFlag,0);
	AddEntry("f_BeginDampingTime", &f_BeginDampingTime,1.0);
	AddEntry("f_EndDampingTime", &f_EndDampingTime,1.0);

	AddEntry("f_Laser_FocusX", &f_Laser_FocusX,0.0);
	AddEntry("f_Laser_FocusY", &f_Laser_FocusY,0.0);
	AddEntry("f_Laser_FocusZ", &f_Laser_FocusZ,0.0);
	AddEntry("f_Laser_FocusW", &f_Laser_FocusW,1.0);
	AddEntry("f_Laser_Mode", &f_Laser_Mode,2.0);
	AddEntry("f_Laser_LtoFocus",&f_Laser_LtoFocus,0.0);
	AddEntry("f_Laser_IniPhase",&f_Laser_IniPhase,0.0);

	AddEntry("f_Laser_Kx", &f_Laser_Kx,0.0);
	AddEntry("f_Laser_Ky", &f_Laser_Ky,0.0);
	AddEntry("f_Laser_Kz", &f_Laser_Kz,0.0);

	AddEntry("i_Laser_Polarization", &i_Laser_Polarization,0);
	AddEntry("i_AzCalculation", &i_AzCalculation,0);

	AddEntry("f_Laser_FocusX_2", &f_Laser_FocusX_2,0.0);
	AddEntry("f_Laser_FocusY_2", &f_Laser_FocusY_2,0.0);
	AddEntry("f_Laser_FocusZ_2", &f_Laser_FocusZ_2,0.0);
	AddEntry("f_Laser_FocusW_2", &f_Laser_FocusW_2,1.0);
	AddEntry("f_Laser_Mode_2", &f_Laser_Mode_2,2.0);
	AddEntry("f_Laser_LtoFocus_2",&f_Laser_LtoFocus_2,0.0);
	AddEntry("f_Laser_IniPhase_2",&f_Laser_IniPhase_2,0.0);

	AddEntry("f_Laser_Kx_2", &f_Laser_Kx_2,0.0);
	AddEntry("f_Laser_Ky_2", &f_Laser_Ky_2,0.0);
	AddEntry("f_Laser_Kz_2", &f_Laser_Kz_2,0.0);
    AddEntry("f_Lambda_L2", &f_Lambda_L2, 0.8);   // here in the unit of um

	AddEntry("i_Laser_Polarization_2", &i_Laser_Polarization_2,0);
	AddEntry("f_ExternalField_Intensity_2", &f_ExternalField_Intensity_2, 1.0);
	AddEntry("f_ExternalField_Duration_2", &f_ExternalField_Duration_2, 1.0);

	FILE *p_File = NULL;
	

   if (rank==NIL)
   {
      p_File = fopen(infile,"rt");     
      if (p_File == NULL)
      {
		  cout<<"ExternalFields no input file: "<<infile<<" "<<endl;
		  // following shoule be noted,since ExternalFields class is defined earlier than domain class
		  //Domain::p_D->out_Flog<<"ExternalFields no input file: "<<infile<<" and exit!\n";
          exit (-1);
      }
   }

   if (p_File)
   {
      rewind(p_File); 
      read(p_File); 
   }
   // now we can only process the laser axis is in the x-z plane, that is Ky=0.0.
   f_Laser_Ky=0.0;
   f_Laser_FocusY=0.0;
   if(f_Laser_Kx!=0)
	    f_Laser_Kx=sqrt(1.0-f_Laser_Kz*f_Laser_Kz)*f_Laser_Kx/fabs(f_Laser_Kx);
   //----------------------------------------------------------------------------
   f_Length_Normalize=f_Lambda_L0*sqrt(1.0/f_Density_Normalize); // unit in um
   f_Laser1_GaussianEnergy=1.375e18*f_ExternalField_Intensity*f_ExternalField_Intensity/f_Lambda_L0/f_Lambda_L0*f_ExternalField_Duration/sqrt(2.0)*f_Laser_FocusW*f_Laser_FocusW/2.0*PI*sqrt(PI)*f_Lambda_L0*f_Lambda_L0*1.0e-8*f_Lambda_L0*1.0e-4/3.0e10; // unit in J
   if(i_my_rank==0){
		cout<<"Lambda_L="<<f_Lambda_L0<<"um n_e="<<f_Density_Normalize<<"n_c"<<" Field is in form of "<<i_Particle_Field<<endl;
		cout<<"Length_Normalize="<<f_Length_Normalize<<"um n_e_Normalize="<<f_Density_Normalize*1.1e21*(1.0/f_Lambda_L0)*(1.0/f_Lambda_L0)<<"/cm^3"<<endl;
		cout<<"Kx="<<f_Laser_Kx<<" Ky="<<f_Laser_Ky<<" Kz="<<f_Laser_Kz<<endl;
		// following shoule be noted,since ExternalFields class is defined earlier than domain class
		//Domain::p_D->out_Flog<<"Lambda_L="<<f_Lambda_L0<<"um n_e="<<f_Density_Normalize<<"n_c"<<" Field is in form of "<<i_Particle_Field<<"\n";
		//Domain::p_D->out_Flog<<"Length_Normalize="<<f_Length_Normalize<<"um n_e_Normalize="<<f_Density_Normalize*1.1e21*(1.0/f_Lambda_L0)*(1.0/f_Lambda_L0)<<"/cm^3"<<"\n";
		//Domain::p_D->out_Flog<<"Kx="<<f_Laser_Kx<<" Ky="<<f_Laser_Ky<<" Kz="<<f_Laser_Kz<<"\n";
		if(i_DampingFlag==1){
			if(i_Particle_Trajectory==0)
				cout<<"Particle velocity Damping is used. BeginDampingTime="<<f_BeginDampingTime<<"  EndDampingTime="<<f_EndDampingTime<<endl;
			else
				cout<<"Warning, Particle Trajectory="<<i_Particle_Trajectory<<", but i_DampingFlag="<<i_DampingFlag<<endl; 
		}
   }
#ifdef V_MPI
   CBuffer *buf = new CBuffer;
   buf->reset();
   pack_nls(buf);
   Domain::p_D->BroadCast(buf);
   if (rank) unpack_nls(buf);
   delete buf;
#endif
};

void ExternalFields::Set_Field(double t,double x,double y,double z,double &Ex,double &Ey,double &Ez,double &Bx, double &By, double &Bz){
	 double Ex1,Ey1,Ez1,Bx1,By1,Bz1;
	 int i_FieldType;
	 Ex=Ey=Ez=Bx=By=Bz=0.0;
	 if(i_Particle_Field<100){
		 i_FieldType=i_Particle_Field;
		 Cal_Field(t,x,y,z,Ex1,Ey1,Ez1,Bx1,By1,Bz1,i_FieldType);
		 Ex+=Ex1;Ey+=Ey1;Ez+=Ez1;Bx+=Bx1;By+=By1;Bz+=Bz1;		 
	 }
	 else if(i_Particle_Field<10000){                                     // two kinds of field
		 i_FieldType=i_Particle_Field/100;
		 Cal_Field(t,x,y,z,Ex1,Ey1,Ez1,Bx1,By1,Bz1,i_FieldType);
		 Ex+=Ex1;Ey+=Ey1;Ez+=Ez1;Bx+=Bx1;By+=By1;Bz+=Bz1;
         i_FieldType=i_Particle_Field%100;
		 Cal_Field(t,x,y,z,Ex1,Ey1,Ez1,Bx1,By1,Bz1,i_FieldType);
		 Ex+=Ex1;Ey+=Ey1;Ez+=Ez1;Bx+=Bx1;By+=By1;Bz+=Bz1;
	 }
	 else if (i_Particle_Field<1000000){                                 // three kinds of field
		 i_FieldType=i_Particle_Field/10000;
		 Cal_Field(t,x,y,z,Ex1,Ey1,Ez1,Bx1,By1,Bz1,i_FieldType);
		 Ex+=Ex1;Ey+=Ey1;Ez+=Ez1;Bx+=Bx1;By+=By1;Bz+=Bz1;
         i_FieldType=(i_Particle_Field-i_FieldType*10000)/100;
		 Cal_Field(t,x,y,z,Ex1,Ey1,Ez1,Bx1,By1,Bz1,i_FieldType);
		 Ex+=Ex1;Ey+=Ey1;Ez+=Ez1;Bx+=Bx1;By+=By1;Bz+=Bz1;
		 i_FieldType=i_Particle_Field-i_FieldType%100;
		 Cal_Field(t,x,y,z,Ex1,Ey1,Ez1,Bx1,By1,Bz1,i_FieldType);
		 Ex+=Ex1;Ey+=Ey1;Ez+=Ez1;Bx+=Bx1;By+=By1;Bz+=Bz1;
	 }
	 else {
	     cout<<"cannot process this kind of laser field, exit -1 ......"<<endl;
		 exit(-1);
	 }

}

void ExternalFields::Cal_Field(double t,double x,double y,double z,double &Ex,double &Ey,double &Ez,double &Bx, double &By, double &Bz, int i_FieldType){
	Ex=Ey=Ez=Bx=By=Bz=0.0;
	switch(i_FieldType){
		case 0:                   // in vacum, no external field
		    Ex=Ey=Ez=Bx=By=Bz=0.0;
	        break;

		case 1:                   // test ion channel
		default:
		{   double ne;
			ne=f_Density_Plasma; 
		    Bx=By=Bz=Ez=0.0;
	        Ex=PI*x*ne;
	        Ey=PI*y*ne;
		}
	        break;
		case 2:                  //here test the blowout bubble	
		{
			double f_beta_wake,ne;
			f_beta_wake=sqrt(1-1/f_gamma_wake/f_gamma_wake);
			ne=f_Density_Plasma;
			if(x*x+y*y+(z-f_beta_wake*t)*(z-f_beta_wake*t)<f_r_bubble*f_r_bubble){
				Ex=PI*x/2.0*ne;
				Ey=PI*y/2.0*ne;
				Ez=PI*(z-f_beta_wake*t)*ne;
				Bx=PI*y/2.0*ne;
				By=-PI*x/2.0*ne;
			}
			else{
                Ex=Ey=Ez=Bx=By=Bz=0.0;
			}
		}
			break;
		case 3:                 // simple circularly polarized laser field, plane wave
		{
			double a0,tao0;
			a0=f_ExternalField_Intensity;
			tao0=f_ExternalField_Duration;
			Ex=a0*sin(2*PI*(z-t))*exp(-(z-t)*(z-t)/(tao0*tao0));
			By=Ex;
			Ey=a0*cos(2*PI*(z-t))*exp(-(z-t)*(z-t)/(tao0*tao0));
			Bx=-Ey;
			Ez=Bz=0.0;
			}
			break;	
		case 4:                 // Nature Focus Undulator
			{
			double Lambda_undulator;
			Lambda_undulator=1.0;
            By= f_Undu_B*cosh(2*PI*y/Lambda_undulator)*cos(2*PI*z/Lambda_undulator);
			Bz=-f_Undu_B*sinh(2*PI*y/Lambda_undulator)*sin(2*PI*z/Lambda_undulator);
			Ex=Ey=Ez=Bx=0.0;
		}
			break;
		case 5:                 // Linear polarized laser field, plane wave
		{
			double a0,tao0;
			a0=f_ExternalField_Intensity;
			tao0=f_ExternalField_Duration;
			switch(i_Laser_Polarization){
				case 0:
				default:
					Ex=a0*sin(2*PI*(z-t))*exp(-(z-t-f_Laser_LtoFocus)*(z-t-f_Laser_LtoFocus)/(tao0*tao0));
					By=Ex;
					Ey=0.0;
					Bx=0.0;
					Ez=Bz=0.0;
					break;
				case 1:
					Ey=a0*sin(2*PI*(z-t))*exp(-(z-t-f_Laser_LtoFocus)*(z-t-f_Laser_LtoFocus)/(tao0*tao0));
					Bx=-Ey;
					Ex=0.0;
					By=0.0;
					Ez=Bz=0.0;
					break;
				case 2:
					Ex=sqrt(0.5)*a0*sin(2*PI*(z-t))*exp(-(z-t-f_Laser_LtoFocus)*(z-t-f_Laser_LtoFocus)/(tao0*tao0));
					By=Ex;
					Ey=sqrt(0.5)*a0*sin(2*PI*(z-t)+PI/2.0)*exp(-(z-t-f_Laser_LtoFocus)*(z-t-f_Laser_LtoFocus)/(tao0*tao0));
					Bx=-Ey;
					Ez=Bz=0.0;
					break;
			}	
		}
			break;	
		case 6:                    // tightly focused long laser pulse
		{
			double a0,tao0;
            a0=f_ExternalField_Intensity;
			tao0=f_ExternalField_Duration;
			double xprime,yprime,zprime,tprime,f_envelopet,f_enveloper,f_rprime2,z_r,w_zprime,ipsilong,ksi,kxi,phi,phiG,phiP,phiR,phi0,Exprime,Eyprime,Ezprime,Bxprime,Byprime,Bzprime;
			xprime=f_Laser_Kz*(x-f_Laser_FocusX)-f_Laser_Kx*(z-f_Laser_FocusZ);  // current x position related to the focus in the laser coordinate system
			yprime=y;
			zprime=f_Laser_Kx*(x-f_Laser_FocusX)+f_Laser_Kz*(z-f_Laser_FocusZ);  // current z position related to the focus in the laser coordinate system
			f_rprime2=yprime*yprime+xprime*xprime;
			z_r=PI*f_Laser_FocusW*f_Laser_FocusW;
			w_zprime=f_Laser_FocusW*sqrt(1+(zprime/z_r)*(zprime/z_r));
            tprime=t+f_Laser_LtoFocus;                                           // current laser peak distance to the focus
			f_envelopet=exp(-(zprime-tprime)*(zprime-tprime)/tao0/tao0);         // temporal shape factor
			phiG=atan(zprime/z_r);
			phiR=PI*f_rprime2*zprime/(z_r*z_r+zprime*zprime);
			phiP=2*PI*(zprime-tprime);
			phi0=f_Laser_IniPhase;
			phi=phi0+phiP-phiR+phiG;                                   // the phase factor
            f_enveloper=a0*f_Laser_FocusW/w_zprime*exp(-f_rprime2/w_zprime/w_zprime); 
			ipsilong=f_Laser_FocusW/z_r;
			ksi=xprime/f_Laser_FocusW;
			kxi=yprime/f_Laser_FocusW;
			switch(i_Laser_Polarization){
				case 0:
				default:
					Exprime=f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi,phiG)/4.0));
					Eyprime=f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG));
					Ezprime=f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					Bxprime=0.0;
					Byprime=f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Bzprime=f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					break;
				case 1:
					Exprime=-f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG));
					Eyprime=f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi,phiG)/4.0));
					Ezprime=f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					Bxprime=-f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Byprime=0.0;
                    Bzprime=f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					break;
				case 2:
                    Exprime=sqrt(0.5)*f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi,phiG)/4.0))-sqrt(0.5)*f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG));
					Eyprime=sqrt(0.5)*f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG)/4.0));
					Ezprime=sqrt(0.5)*f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG));
					Bxprime=-sqrt(0.5)*f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Byprime=sqrt(0.5)*f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Bzprime=sqrt(0.5)*f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					break;
			}
			Ex= Exprime*f_Laser_Kz+Ezprime*f_Laser_Kx; 
			Ey= Eyprime;  
			Ez=-Exprime*f_Laser_Kx+Ezprime*f_Laser_Kz;
			Bx= Bxprime*f_Laser_Kz+Bzprime*f_Laser_Kx;
			By= Byprime;
			Bz=-Bxprime*f_Laser_Kx+Bzprime*f_Laser_Kz;
            double ne;
			ne=f_Density_Plasma; 
	        Ex+=PI*x*ne;
	        Ey+=PI*y*ne;
		}
			break;
		case 7:                                 // laser field calculated from analytical formula
		{
            Ex=Ey=Ez=Bx=By=Bz=0.0;			
            double f_betag_laser,f_Laser_FocusW2,r2,ne,kexi,f_beta_wake,a0,tao0,t_invacuum,t_inplasma,laserbeta_plasma,laserbeta_vacuum,laser_center;
			double f_lambda_p,f_ramp_z=1.0,f_radius,f_dr=1.0,f_ramp_r=1.0,temp1,temp2,temp3,f_dx=0.001,f_dy=0.001,f_dz=0.001,f_dt=0.001;

			f_lambda_p=sqrt(1.0/f_Density_Plasma);
			// ------------------------------------------ Here add the bubble field	-----------------------------------------------------------------		
			if(z<f_ExternalField_End){
				ne=f_Density_Plasma; 
				f_betag_laser=1-0.5*(ne+1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
				f_beta_wake=f_betag_laser;        // wake phase velocity
				if(z>f_ExternalField_End-f_lambda_p)
					f_ramp_z=0.5*(1.0+cos(PI*(z+f_lambda_p-f_ExternalField_End)/f_lambda_p));
				f_radius=sqrt(x*x+y*y+(z-f_beta_wake*t)*(z-f_beta_wake*t));
				if(f_radius>f_r_bubble-f_dr)
					f_ramp_r=0.5*(1.0+cos(PI*(f_radius+f_dr-f_r_bubble)/f_dr));
				f_ramp_z*=f_ramp_r;
				
				if(f_radius*f_radius<f_r_bubble*f_r_bubble && z<f_ExternalField_End){
					Ex+=PI*x/2.0*ne*f_ramp_z;
					Ey+=PI*y/2.0*ne*f_ramp_z;
					Ez+=PI*(z-f_beta_wake*t)*ne*f_ramp_z;
					Bx+=PI*y/2.0*ne*f_ramp_z;
					By+=-PI*x/2.0*ne*f_ramp_z;
					;
				}					
			}

			// ------------------------------------------- Here add the laser field -----------------------------------------
            ne=f_Density_Plasma;
			laserbeta_plasma=1-0.5*(ne+1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
			laserbeta_vacuum=1-0.5*1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW);
			t_inplasma=(f_ExternalField_End-(f_Laser_LtoFocus+f_ExternalField_Duration/2.0))/laserbeta_plasma;
			t_invacuum=(f_ExternalField_End-(f_Laser_LtoFocus-f_ExternalField_Duration/2.0))/laserbeta_plasma;

			if(t<t_inplasma){   // laser is totally in plasma
				laser_center=f_Laser_LtoFocus+laserbeta_plasma*t;
				kexi=z-(laser_center-f_ExternalField_Duration/2.0);
				if(kexi>0.0 && kexi<f_ExternalField_Duration){
					a0=f_ExternalField_Intensity;
					tao0=f_ExternalField_Duration;
					r2=x*x+y*y;
					f_Laser_FocusW2=f_Laser_FocusW*f_Laser_FocusW;
					temp1=-a0*exp(-r2/f_Laser_FocusW2)*((0.5*(1-cos(2*PI*kexi/tao0)))*sin(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z));
					temp2=a0*exp(-r2/f_Laser_FocusW2)*laserbeta_plasma/(2.0*tao0)*sin(2*PI*kexi/tao0)*cos(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z); // this is a higher order modified term
					Ex+=temp1+temp2;     
					temp1=-a0*exp(-r2/f_Laser_FocusW2)*((0.5*(1-cos(2*PI*kexi/tao0)))*sin(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z)*(1-ne/2.0-1.0/(2*PI*PI*f_Laser_FocusW2)));
					temp2=a0/(2.0*tao0)*exp(-r2/f_Laser_FocusW2)*sin(2*PI*kexi/tao0)*cos(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z);            // this is a higher order modified term
					temp3=-1.0/(2*PI)*((az(x+f_dx,y,z,t)-az(x-f_dx,y,z,t))/(2*f_dx));
					By+=temp1+temp2+temp3;

			        Ey+=0.0;
			        Ez+=-1.0/(2*PI)*(az(x,y,z,t+f_dt)-az(x,y,z,t-f_dt))/(2*f_dt);
			        Bx+=1.0/(2*PI)*(az(x,y+f_dy,z,t)-az(x,y-f_dy,z,t))/(2*f_dy);
			        Bz+=1.0/(2*PI)*2*y/f_Laser_FocusW2*a0*exp(-r2/f_Laser_FocusW2)*((0.5*(1-cos(2*PI*kexi/tao0)))*cos(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z));
				}
			}
			else if(t>t_invacuum){                   // laser is totally in vacuum    
				double z_r,f_Laser_W,f_ExLaser,z1;
				tao0=f_ExternalField_Duration/laserbeta_plasma*laserbeta_vacuum;  
				a0=f_ExternalField_Intensity*sqrt(f_ExternalField_Duration/tao0); 
				laser_center=f_ExternalField_End+tao0/2.0+laserbeta_vacuum*(t-t_invacuum);
				kexi=z-(laser_center-tao0/2.0); 
				if(kexi>0.0 && kexi<tao0){
					ne=0.0; 
					f_betag_laser=laserbeta_vacuum;
					r2=x*x+y*y;
					f_Laser_FocusW2=f_Laser_FocusW*f_Laser_FocusW;
					z_r=PI*f_Laser_FocusW2;
					z1=z-f_ExternalField_End;
					f_Laser_W=f_Laser_FocusW*sqrt(1+z1*z1/z_r/z_r);
					f_ExLaser=-a0*f_Laser_FocusW/f_Laser_W*exp(-r2/f_Laser_W/f_Laser_W)*((0.5*(1-cos(2*PI*kexi/tao0)))*sin(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*f_ExternalField_End+PI*r2*z1/(z1*z1+z_r*z_r)-atan(z1/z_r)));
					Ex+=f_ExLaser;
					By+=f_ExLaser; 
				}
			}
			else{   //  laser is part in vacuum, part in channel   t_inplasma<= t <=t_invacuum
				double laserleft,laserright;
				laserleft=f_Laser_LtoFocus-f_ExternalField_Duration/2.0+laserbeta_plasma*t;
				laserright=f_ExternalField_End+laserbeta_vacuum*(t-t_inplasma);
				if(z>laserleft && z<=f_ExternalField_End){
					laser_center=f_Laser_LtoFocus+laserbeta_plasma*t;
					kexi=z-(laser_center-f_ExternalField_Duration/2.0);
					a0=f_ExternalField_Intensity;
					tao0=f_ExternalField_Duration;
					r2=x*x+y*y;
					f_Laser_FocusW2=f_Laser_FocusW*f_Laser_FocusW;
					Ex+=-a0*exp(-r2/f_Laser_FocusW2)*((0.5*(1-cos(2*PI*kexi/tao0)))*sin(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z));
					By+=-a0*exp(-r2/f_Laser_FocusW2)*((0.5*(1-cos(2*PI*kexi/tao0)))*sin(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*z))*(1-ne/2.0-1.0/(2*PI*PI*f_Laser_FocusW2));
				}
				else if(z>f_ExternalField_End && z<=laserright){
						double z_r,f_Laser_W,f_ExLaser,z1;
						tao0=f_ExternalField_Duration/laserbeta_plasma*laserbeta_vacuum;
						a0=f_ExternalField_Intensity*sqrt(f_ExternalField_Duration/tao0);
						laser_center=f_ExternalField_End+tao0/2.0+laserbeta_vacuum*(t-t_invacuum);
						kexi=z-(laser_center-tao0/2.0);
						ne=0.0; 
						f_betag_laser=laserbeta_vacuum;
						r2=x*x+y*y;
						f_Laser_FocusW2=f_Laser_FocusW*f_Laser_FocusW;
						z_r=PI*f_Laser_FocusW2;
						z1=z-f_ExternalField_End;
						f_Laser_W=f_Laser_FocusW*sqrt(1+z1*z1/z_r/z_r);
						f_ExLaser=-a0*f_Laser_FocusW/f_Laser_W*exp(-r2/f_Laser_W/f_Laser_W)*(sin(PI*(kexi)/tao0)*sin(2*PI*(z-t)-(PI*ne+1/(PI*f_Laser_FocusW2))*f_ExternalField_End+PI*r2*z1/(z1*z1+z_r*z_r)-atan(z1/z_r)));
	                    Ex+=f_ExLaser;
						By+=f_ExLaser;						
				}		
			}
		}
			break;
		case 8:                                        // calculate EB from ax,az
		{
			Ex=Ey=Ez=Bx=By=Bz=0.0;			
            double f_betag_laser,f_Laser_FocusW2,ne,kexi,f_beta_wake,a0,tao0,t_invacuum,t_inplasma,laserbeta_plasma,laserbeta_vacuum,laser_center;
			double f_dt=0.001,f_dx=0.001,f_dy=0.001,f_dz=0.001,f_lambda_p,f_ramp_z=1.0,f_radius,f_dr=1.0,f_ramp_r=1.0;
			
			if(f_Density_Plasma!=0.0)
				f_lambda_p=sqrt(1.0/f_Density_Plasma);
			else
                f_lambda_p=10.0;
			
			// ------------------------------------------ Here add the bubble field	-----------------------------------------------------------------		
			if(z<f_ExternalField_End){
				ne=f_Density_Plasma; 
				f_betag_laser=1-0.5*(ne+1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
				f_beta_wake=f_betag_laser;        // wake phase velocity
				if(z>f_ExternalField_End-f_lambda_p)
					f_ramp_z=0.5*(1.0+cos(PI*(z+f_lambda_p-f_ExternalField_End)/f_lambda_p));
				f_radius=sqrt(x*x+y*y+(z-f_beta_wake*t)*(z-f_beta_wake*t));
				if(f_radius>f_r_bubble-f_dr)
					f_ramp_r=0.5*(1.0+cos(PI*(f_radius+f_dr-f_r_bubble)/f_dr));
				f_ramp_z*=f_ramp_r;
				
				if(f_radius*f_radius<f_r_bubble*f_r_bubble && z<f_ExternalField_End){
					Ex+=PI*x/2.0*ne*f_ramp_z;
					Ey+=PI*y/2.0*ne*f_ramp_z;
					Ez+=PI*(z-f_beta_wake*t)*ne*f_ramp_z;
					Bx+=-PI*y/2.0*ne*f_ramp_z;
					By+=-PI*x/2.0*ne*f_ramp_z;
					;
				}	
				
			}	
			
			// ------------------------------------------- Here add the laser field from the vector potential of ax(x,y,z,t), az(x,y,z,t) ------------
			Ex+=-1.0/(2*PI)*(ax(x,y,z,t+f_dt)-ax(x,y,z,t-f_dt))/(2*f_dt);
			Ey+=0.0;
			Ez+=-1.0/(2*PI)*(az(x,y,z,t+f_dt)-az(x,y,z,t-f_dt))/(2*f_dt);
			Bx+=1.0/(2*PI)*(az(x,y+f_dy,z,t)-az(x,y-f_dy,z,t))/(2*f_dy);
			By+=1.0/(2*PI)*((ax(x,y,z+f_dz,t)-ax(x,y,z-f_dz,t))/(2*f_dz))-1.0/(2*PI)*((az(x+f_dx,y,z,t)-az(x-f_dx,y,z,t))/(2*f_dx));
			Bz+=-1.0/(2*PI)*(ax(x,y+f_dy,z,t)-ax(x,y-f_dy,z,t))/(2*f_dy);
		}
			break;
            
        case 9:
        {       // homogenous magnetic field
            Bx+=f_Undu_B;
            Ex=Ey=Ez=By=Bz=0.0;
        }
            break;
            
		case 41:         // simple Undulator
		{
		    double Lambda_undulator;
			Lambda_undulator=1.0;
            By=f_Undu_B*cos(2*PI*z/Lambda_undulator);
			Ex=Ey=Ez=Bx=Bz=0.0;
		}
			break;
		case 42:         // two faces focusing Undulator
		{
			double Lambda_undulator;
			Lambda_undulator=1.0;
			Bx=-f_Undu_B*(f_Undulator_a*y+f_Undulator_b*y*cos(2*2*PI*z/Lambda_undulator));
			By= f_Undu_B*(cos(2*PI*z/Lambda_undulator)-f_Undulator_a*x-f_Undulator_b*x*cos(2*2*PI*z/Lambda_undulator));
			Bz=-f_Undu_B*(2*PI*y/Lambda_undulator*sin(2*PI*z/Lambda_undulator)-2*f_Undulator_b*2*PI*y/Lambda_undulator*x*sin(2*2*PI*z/Lambda_undulator));
			Ex=Ey=Ez=0.0;
		}
			break;
		case 16:        // two tightly focused long laser pulses,   studing for electron dynamics in colliding pulses
		{
			double a0,tao0;
			//case 16, for the first laser pulse
            a0=f_ExternalField_Intensity;
			tao0=f_ExternalField_Duration;
			double xprime,yprime,zprime,tprime,f_envelopet,f_enveloper,f_rprime2,z_r,w_zprime,ipsilong,ksi,kxi,phi,phiG,phiP,phiR,phi0,Exprime,Eyprime,Ezprime,Bxprime,Byprime,Bzprime;
			xprime=f_Laser_Kz*(x-f_Laser_FocusX)-f_Laser_Kx*(z-f_Laser_FocusZ);  // current x position related to the focus in the laser coordinate system
			yprime=y;
			zprime=f_Laser_Kx*(x-f_Laser_FocusX)+f_Laser_Kz*(z-f_Laser_FocusZ);  // current z position related to the focus in the laser coordinate system
			f_rprime2=yprime*yprime+xprime*xprime;
			z_r=PI*f_Laser_FocusW*f_Laser_FocusW;
			w_zprime=f_Laser_FocusW*sqrt(1+(zprime/z_r)*(zprime/z_r));
            tprime=t+f_Laser_LtoFocus;                                           // current laser peak distance to the focus
			f_envelopet=exp(-(zprime-tprime)*(zprime-tprime)/tao0/tao0);         // temporal shape factor
			phiG=atan(zprime/z_r);
			phiR=PI*f_rprime2*zprime/(z_r*z_r+zprime*zprime);
			phiP=2*PI*(zprime-tprime);
			phi0=f_Laser_IniPhase;
			phi=phi0+phiP-phiR+phiG;                                   // the phase factor
            f_enveloper=a0*f_Laser_FocusW/w_zprime*exp(-f_rprime2/w_zprime/w_zprime); 
			ipsilong=f_Laser_FocusW/z_r;
			ksi=xprime/f_Laser_FocusW;
			kxi=yprime/f_Laser_FocusW;
			switch(i_Laser_Polarization){
				case 0:
				default:
					Exprime=f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi,phiG)/4.0));
					Eyprime=f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG));
					Ezprime=f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					Bxprime=0.0;
					Byprime=f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Bzprime=f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					break;
				case 1:
					Exprime=-f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG));
					Eyprime=f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi,phiG)/4.0));
					Ezprime=f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					Bxprime=-f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Byprime=0.0;
                    Bzprime=f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					break;
				case 2:
                    Exprime=sqrt(0.5)*f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi,phiG)/4.0))-sqrt(0.5)*f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG));
					Eyprime=sqrt(0.5)*f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG)/4.0));
					Ezprime=sqrt(0.5)*f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi+TWOPI/4.0,phiG));
					Bxprime=-sqrt(0.5)*f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Byprime=sqrt(0.5)*f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW,w_zprime,phi,phiG);
					Bzprime=sqrt(0.5)*f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW,w_zprime,phi,phiG));
					break;
			}
			Ex+= Exprime*f_Laser_Kz+Ezprime*f_Laser_Kx; 
			Ey+= Eyprime;  
			Ez+=-Exprime*f_Laser_Kx+Ezprime*f_Laser_Kz;
			Bx+= Bxprime*f_Laser_Kz+Bzprime*f_Laser_Kx;
			By+= Byprime;
			Bz+=-Bxprime*f_Laser_Kx+Bzprime*f_Laser_Kz;
		
            //     case 16, for the second laser pulse                the second pulse length, width are both normalized to the fundament unit of the code as the first pulse
			double omega2over1=f_Lambda_L0/f_Lambda_L2;
		    a0=f_ExternalField_Intensity_2*omega2over1;
			tao0=f_ExternalField_Duration_2;
			xprime=f_Laser_Kz_2*(x-f_Laser_FocusX_2)-f_Laser_Kx_2*(z-f_Laser_FocusZ_2);  // current x position related to the focus in the laser coordinate system
			yprime=y;
			zprime=f_Laser_Kx_2*(x-f_Laser_FocusX_2)+f_Laser_Kz_2*(z-f_Laser_FocusZ_2);  // current z position related to the focus in the laser coordinate system
			f_rprime2=yprime*yprime+xprime*xprime;
			z_r=omega2over1*PI*f_Laser_FocusW_2*f_Laser_FocusW_2;
			w_zprime=f_Laser_FocusW_2*sqrt(1+(zprime/z_r)*(zprime/z_r));
            tprime=t+f_Laser_LtoFocus_2;                                           // current laser peak distance to the focus
			f_envelopet=exp(-(zprime-tprime)*(zprime-tprime)/tao0/tao0);         // temporal shape factor
			phiG=atan(zprime/z_r);
			phiR=omega2over1*PI*f_rprime2*zprime/(z_r*z_r+zprime*zprime);
			phiP=omega2over1*2*PI*(zprime-tprime);
			phi0=f_Laser_IniPhase_2;
			phi=phi0+phiP-phiR+phiG;                                   // the phase factor
            f_enveloper=a0*f_Laser_FocusW_2/w_zprime*exp(-f_rprime2/w_zprime/w_zprime); 
			ipsilong=f_Laser_FocusW_2/z_r;
			ksi=xprime/f_Laser_FocusW_2;
			kxi=yprime/f_Laser_FocusW_2;
			switch(i_Laser_Polarization_2){
				case 0:
				default:
					Exprime=f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW_2,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW_2,w_zprime,phi,phiG)/4.0));
					Eyprime=f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW_2,w_zprime,phi,phiG));
					Ezprime=f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG));
					Bxprime=0.0;
					Byprime=f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG);
					Bzprime=f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG));
					break;
				case 1:
					Exprime=-f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW_2,w_zprime,phi,phiG));
					Eyprime=f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW_2,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW_2,w_zprime,phi,phiG)/4.0));
					Ezprime=f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG));
					Bxprime=-f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG);
					Byprime=0.0;
                    Bzprime=f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG));
					break;
				case 2:
                    Exprime=sqrt(0.5)*f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW_2,w_zprime,phi,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW_2,w_zprime,phi,phiG)/4.0))-sqrt(0.5)*f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW_2,w_zprime,phi+TWOPI/4.0,phiG));
					Eyprime=sqrt(0.5)*f_envelopet*f_enveloper*ksi*kxi*(ipsilong*ipsilong*Sn(2,f_Laser_FocusW_2,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*(Sn(0,f_Laser_FocusW_2,w_zprime,phi+TWOPI/4.0,phiG)+ipsilong*ipsilong*(ksi*ksi*Sn(2,f_Laser_FocusW_2,w_zprime,phi+TWOPI/4.0,phiG)-f_rprime2*f_rprime2/f_Laser_FocusW/f_Laser_FocusW*Sn(3,f_Laser_FocusW_2,w_zprime,phi+TWOPI/4.0,phiG)/4.0));
					Ezprime=sqrt(0.5)*f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*ksi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi+TWOPI/4.0,phiG));
					Bxprime=-sqrt(0.5)*f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG);
					Byprime=sqrt(0.5)*f_envelopet*f_enveloper*Sn(0,f_Laser_FocusW_2,w_zprime,phi,phiG);
					Bzprime=sqrt(0.5)*f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG))+sqrt(0.5)*f_envelopet*f_enveloper*kxi*(ipsilong*Cn(1,f_Laser_FocusW_2,w_zprime,phi,phiG));
					break;
			}
			Ex+= Exprime*f_Laser_Kz_2+Ezprime*f_Laser_Kx_2; 
			Ey+= Eyprime;  
			Ez+=-Exprime*f_Laser_Kx_2+Ezprime*f_Laser_Kz_2;
			Bx+= Bxprime*f_Laser_Kz_2+Bzprime*f_Laser_Kx_2;
			By+= Byprime;
			Bz+=-Bxprime*f_Laser_Kx_2+Bzprime*f_Laser_Kz_2;
		}
			break;
	}
}

double ExternalFields::DampingCoefficient(double t){
	double f_coefficient,f_tEnd;
	f_tEnd=Domain::p_D->l_Step*Domain::p_D->GetTs();
	if(t<f_BeginDampingTime && t>0){
		f_coefficient=(cos(PI*(1.0-t/f_BeginDampingTime))+1.0)/2.0;
	}else if(t>f_tEnd-f_EndDampingTime){
		f_coefficient=(cos(PI*(t-(f_tEnd-f_EndDampingTime))/f_EndDampingTime)+1.0)/2.0;
	}
	else{
		f_coefficient=1.0;
	}
	return(f_coefficient);
}


double ExternalFields::ax(double x,double y,double z,double t){                                       //110422 add the ax inside the vacuum
	double beta_g,f_ax,f_laser_center,f_laser_front,f_laser_end,kexi,r2,beta_gv,t_front2boundary,f_newLength,f_newa0,z1,f_distance2front,rs,zr;
	r2=x*x+y*y;
	f_ax=0.0;
	beta_g=1-0.5*(f_Density_Plasma+1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
	beta_gv=1-0.5*(1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
	if(z<=f_ExternalField_End){                                                                     // ax inside plasma
		f_laser_center=beta_g*t+f_Laser_LtoFocus;
		f_laser_front=f_laser_center+f_ExternalField_Duration/2.0;
		f_laser_end=f_laser_center-f_ExternalField_Duration/2.0;
		kexi=z-(f_laser_center-f_ExternalField_Duration/2.0);
		if(z>=f_laser_end && z<=f_laser_front)
			f_ax=f_ExternalField_Intensity*exp(-(r2)/f_Laser_FocusW/f_Laser_FocusW)*(0.5*(1-cos(2*PI*kexi/f_ExternalField_Duration)))*cos(2*PI*(z-t)-(PI*f_Density_Plasma+1.0/PI/f_Laser_FocusW/f_Laser_FocusW)*z);
	}
	else{                                                                                            // ax in vacuum
		t_front2boundary=(f_ExternalField_End-(f_Laser_LtoFocus+f_ExternalField_Duration/2.0))/beta_g;  // Here we have not considered the pulse lengthen
		f_laser_front=f_ExternalField_End+(t-t_front2boundary)*beta_gv;
		f_laser_end=f_laser_front-f_ExternalField_Duration;
		f_newLength=f_ExternalField_Duration*beta_gv/beta_g;
		f_newa0=f_ExternalField_Intensity*sqrt(f_ExternalField_Duration/f_newLength);
		if(z>f_laser_end && z<f_laser_front){
			kexi=z-(f_laser_front-f_newLength);		 // position related to the pulse_end
			z1=z-f_ExternalField_End;	
			zr=PI*f_Laser_FocusW*f_Laser_FocusW;
			rs=f_Laser_FocusW*sqrt(1+z1*z1/zr/zr);
			f_ax=f_newa0*exp(-(r2)/rs/rs)*f_Laser_FocusW/rs*(0.5*(1-cos(2*PI*kexi/f_newLength)))*cos(2*PI*(z-t)-atan(z1/zr)+z1/zr*(r2/rs/rs)-(PI*f_Density_Plasma+1.0/PI/f_Laser_FocusW/f_Laser_FocusW)*f_ExternalField_End);
		}
	}
	return(f_ax);
}
double ExternalFields::az(double x,double y,double z,double t){
	double f_az,dz,dx,beta_g,f_ax,f_laser_center,f_laser_front,f_laser_end,kexi,r2,beta_gv,t_front2boundary,f_newLength,f_newa0,z1,f_distance2front,rs,zr,temp1;
    beta_g=1-0.5*(f_Density_Plasma+1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
	beta_gv=1-0.5*(1.0/(PI*f_Laser_FocusW)/(PI*f_Laser_FocusW));
	t_front2boundary=(f_ExternalField_End-(f_Laser_LtoFocus+f_ExternalField_Duration/2.0))/beta_g;  
	if(t<t_front2boundary)
		f_laser_front=beta_g*t+f_Laser_LtoFocus+f_ExternalField_Duration/2.0;
	else
        f_laser_front=f_ExternalField_End+beta_gv*(t-t_front2boundary);

	f_laser_end=f_laser_front-f_ExternalField_Duration;

	dz=0.05;dx=0.001;
	f_az=0.0;z1=f_laser_end;
	if(z>f_laser_end && z<=f_laser_front){
		switch (i_AzCalculation){
			default:
			case 0:
			{ 
				while(z1<z){
					f_az+=-(ax(x+dx,y,z1,t)-ax(x-dx,y,z1,t))/(2*dx)*dz;             // az calculated from differentiation and integration, time costing
					z1+=dz;
				}
				if(z1-dz<z)
					f_az+=-(ax(x+dx,y,z1,t)-ax(x-dx,y,z1,t))/(2*dx)*(z-(z1-dz));
			}
			break;
			case 1:
			{
				r2=x*x+y*y;
				temp1=(PI*f_Density_Plasma+1.0/PI/f_Laser_FocusW/f_Laser_FocusW);
				if(z<=f_ExternalField_End){
					f_laser_center=beta_g*t+f_Laser_LtoFocus;
					f_laser_front=f_laser_center+f_ExternalField_Duration/2.0;
					f_laser_end=f_laser_center-f_ExternalField_Duration/2.0;
					kexi=z-(f_laser_center-f_ExternalField_Duration/2.0);
					if(z>=f_laser_end && z<=f_laser_front)
						f_az=1.0/(2*PI-temp1)*2*x/f_Laser_FocusW/f_Laser_FocusW*f_ExternalField_Intensity*exp(-(r2)/f_Laser_FocusW/f_Laser_FocusW)*(0.5*(1-cos(2*PI*kexi/f_ExternalField_Duration)))*sin(2*PI*(z-t)-temp1*z);
				}
				else{
					t_front2boundary=(f_ExternalField_End-(f_Laser_LtoFocus+f_ExternalField_Duration/2.0))/beta_g;  // Here we have not considered the pulse lengthen
					f_laser_front=f_ExternalField_End+(t-t_front2boundary)*beta_gv;
					f_laser_end=f_laser_front-f_ExternalField_Duration;
					f_newLength=f_ExternalField_Duration*beta_gv/beta_g;
					f_newa0=f_ExternalField_Intensity*sqrt(f_ExternalField_Duration/f_newLength);
					if(z>f_laser_end && z<f_laser_front){
						kexi=z-(f_laser_front-f_newLength);		 // position related to the pulse_end
						z1=z-f_ExternalField_End;	
						zr=PI*f_Laser_FocusW*f_Laser_FocusW;
						rs=f_Laser_FocusW*sqrt(1+z1*z1/zr/zr);
						f_az=1.0/(2*PI)*2*x/rs/rs*f_newa0*exp(-(r2)/rs/rs)*f_Laser_FocusW/rs*(0.5*(1-cos(2*PI*kexi/f_newLength)))*sin(2*PI*(z-t)-atan(z1/zr)+z1/zr*((r2)/rs/rs+1)-(PI*f_Density_Plasma+1/PI/f_Laser_FocusW/f_Laser_FocusW)*f_ExternalField_End);
					}
				}
			}
			break;
		}
	}
	return(f_az);
}


void ExternalFields::Write_Field(long l_StepOrder, double t, double x_start, double y_start, double z_start, double dx, double dy, double dz, long l_Field_Nx, long l_Field_Ny, long l_Field_Nz, long local_N_Start, long local_N_End){
	double x,y,z,f_dt,f_dx,f_dy,f_dz,Ex,Ey,Ez,Bx,By,Bz;
	char fname[200],fdatasetname[200];
	long i,j,k;

	hid_t       file, fdataset;         /* File and dataset            */
	hid_t       fdataspace,memspace;   /* Dataspace handles           */
	hsize_t     dimsf[1], dimsfmem[1];              /* Dataset dimensions          */
	herr_t      status;                /* Error checking              */
	hid_t rank1 = 1;
	hid_t rank2 = 2;
	hid_t rank3 = 3;
	hsize_t     count[1],offset[1],count_out[1],offset_out[1];             


   dimsf[0] = l_Field_Nx*l_Field_Ny*l_Field_Nz;
   dimsfmem[0] = local_N_End-local_N_Start+1;
   
   /* Creating the dataspace                                         */
   fdataspace = H5Screate_simple(rank1, dimsf, NULL); 
   assert (fdataspace >= 0);

   memspace = H5Screate_simple(rank1, dimsfmem, NULL); 
   assert (memspace >= 0);


	offset[0] = local_N_Start; 
    count[0] = local_N_End-local_N_Start+1;  
	offset_out[0]=0;
	count_out[0]=local_N_End-local_N_Start+1;

	//Open the file and the dataset.

    sprintf(fname,"FieldDis_at_%d.h5",l_StepOrder);
	file = H5Fopen (fname, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file <= 0 && i_my_rank==0) {
         Domain::p_D->out_Flog << "FieldDis re-open: no such file " << fname << "\n";
         Domain::p_D->out_Flog.flush();
         cout << "FieldDis re-open: no such file " << fname << "\n";
    }
	else{
		if(i_my_rank==0){
		//cout<<"my_rank="<<i_my_rank<<"  local_N_Start="<<local_N_Start<<"  local_N_End="<<local_N_End<<" ExternalFieldsOutputFieldDis begin!"<<endl;
			Domain::p_D->out_Flog << "FieldDis re-open: Find file " << fname << "\n";
			Domain::p_D->out_Flog.flush();
		}
	}
	assert (file >= 0);



	sprintf(fdatasetname,"Ex");
	fdataset = H5Dopen (file, fdatasetname);
	//fdataspace = H5Dget_space (fdataset);   

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, Data_ArrayEx);      

    assert (status >= 0);
	status = H5Dclose(fdataset);
	assert (status >= 0);

    sprintf(fdatasetname,"Ez");
	fdataset = H5Dopen (file, fdatasetname);
	fdataspace = H5Dget_space (fdataset);   

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, Data_ArrayEz);                                                         
    assert (status >= 0);
	status = H5Dclose(fdataset);
	assert (status >= 0);


	sprintf(fdatasetname,"By");
	fdataset = H5Dopen (file, fdatasetname);
	fdataspace = H5Dget_space (fdataset);   

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, Data_ArrayBy);                                                         
    assert (status >= 0);
	status = H5Dclose(fdataset);
	assert (status >= 0);


	status = H5Fclose(file);
	assert (status >= 0);

	delete[] Data_ArrayEx;
	delete[] Data_ArrayEz;
	delete[] Data_ArrayBy;
}


void ExternalFields::Write_Field_Parallel(long l_StepOrder, double t, double x_start, double y_start, double z_start, double dx, double dy, double dz, long l_Field_Nx, long l_Field_Ny, long l_Field_Nz, long local_N_Start, long local_N_End){
	double x,y,z,f_dt,f_dx,f_dy,f_dz,Ex,Ey,Ez,Bx,By,Bz;
	char fname[200],fdatasetname[200];
	long i,j,k;

	hid_t       file, fdataset;         /* File and dataset            */
	hid_t       fdataspace,memspace;   /* Dataspace handles           */
	hsize_t     dimsf[1], dimsfmem[1];              /* Dataset dimensions          */
	herr_t      status;                /* Error checking              */
	hid_t rank1 = 1;
	hid_t rank2 = 2;
	hid_t rank3 = 3;
	hsize_t     count[1],offset[1],count_out[1],offset_out[1];         

	hid_t acc_tpl1;  // File access templates
	herr_t ret;   // Generic return value

#ifdef PARALLELHDF5
	acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
	assert(acc_tpl1 >=0);
	ret = H5Pset_fapl_mpio(acc_tpl1,MPI_COMM_WORLD,MPI_INFO_NULL);
    assert(ret >=0);
#endif


   dimsf[0] = l_Field_Nx*l_Field_Ny*l_Field_Nz;
   dimsfmem[0] = local_N_End-local_N_Start+1;
   
   /* Creating the dataspace                                         */
   fdataspace = H5Screate_simple(rank1, dimsf, NULL); 
   assert (fdataspace >= 0);

   memspace = H5Screate_simple(rank1, dimsfmem, NULL); 
   assert (memspace >= 0);


	offset[0] = local_N_Start; 
    count[0] = local_N_End-local_N_Start+1;  
	offset_out[0]=0;
	count_out[0]=local_N_End-local_N_Start+1;

	//Open the file and the dataset.

    sprintf(fname,"FieldDis_at_%d.h5",l_StepOrder);

#ifdef PARALLELHDF5
	file = H5Fopen (fname, H5F_ACC_RDWR, acc_tpl1);
#else 
    file = H5Fopen (fname, H5F_ACC_RDWR, H5P_DEFAULT);
#endif

    if (file <= 0 && i_my_rank==0) {
         Domain::p_D->out_Flog << "FieldDis re-open: no such file " << fname << "\n";
         Domain::p_D->out_Flog.flush();
         cout << "FieldDis re-open: no such file " << fname << "\n";
    }
	else{
		if(i_my_rank==0){
		//cout<<"my_rank="<<i_my_rank<<"  local_N_Start="<<local_N_Start<<"  local_N_End="<<local_N_End<<" ExternalFieldsOutputFieldDis begin!"<<endl;
			Domain::p_D->out_Flog << "FieldDis re-open: Find file " << fname << "\n";
			Domain::p_D->out_Flog.flush();
		}
	}
	assert (file >= 0);

#ifdef PARALLELHDF5
	ret = H5Pclose(acc_tpl1);
	assert (ret >= 0);
#endif



	sprintf(fdatasetname,"Ex");
	fdataset = H5Dopen (file, fdatasetname);
	//fdataspace = H5Dget_space (fdataset);   

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, Data_ArrayEx);      

    assert (status >= 0);
	status = H5Dclose(fdataset);
	assert (status >= 0);

    sprintf(fdatasetname,"Ez");
	fdataset = H5Dopen (file, fdatasetname);
	fdataspace = H5Dget_space (fdataset);   

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, Data_ArrayEz);                                                         
    assert (status >= 0);
	status = H5Dclose(fdataset);
	assert (status >= 0);


	sprintf(fdatasetname,"By");
	fdataset = H5Dopen (file, fdatasetname);
	fdataspace = H5Dget_space (fdataset);   

    //Define hyperslab in the dataset.   
    status = H5Sselect_hyperslab(fdataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      assert (status >= 0);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
      assert (status >= 0);
    status = H5Dwrite(fdataset, H5T_NATIVE_DOUBLE, memspace, fdataspace,
        H5P_DEFAULT, Data_ArrayBy);                                                         
    assert (status >= 0);
	status = H5Dclose(fdataset);
	assert (status >= 0);


	status = H5Fclose(file);
	assert (status >= 0);

	delete[] Data_ArrayEx;
	delete[] Data_ArrayEz;
	delete[] Data_ArrayBy;
}


void ExternalFields::Cal_OutField(long l_StepOrder, double t, double x_start, double y_start, double z_start, double dx, double dy, double dz, long l_Field_Nx, long l_Field_Ny, long l_Field_Nz, long local_N_Start, long local_N_End){
	double x,y,z,f_dt,f_dx,f_dy,f_dz,Ex,Ey,Ez,Bx,By,Bz;

	long i,j,k;
	Data_ArrayEx= new double[local_N_End-local_N_Start+1];
	Data_ArrayEz= new double[local_N_End-local_N_Start+1];
	Data_ArrayBy= new double[local_N_End-local_N_Start+1];

	f_dx=f_dy=f_dz=f_dt=0.001;
	for(long l=local_N_Start;l<=local_N_End;l++){
		i=l%l_Field_Nz;
		k=l/(l_Field_Nz*l_Field_Nx);
		j=(l-l_Field_Nz*l_Field_Nx*k)/l_Field_Nz; 
		x=x_start+j*dx; 
		y=y_start+k*dy;
		z=z_start+i*dz; 
		Set_Field(t,x,y,z,Ex,Ey,Ez,Bx,By,Bz);
		Data_ArrayEx[l-local_N_Start]=Ex;
		Data_ArrayEz[l-local_N_Start]=Ez;
		Data_ArrayBy[l-local_N_Start]=By;
	}
}


double ExternalFields::CP_MC_LaserField_GetNp(){	// only laser 1 is used, Gaussian spatial and temporal shape
   return(pow(PI/2.0,2.5)*f_Lambda_L0*f_Lambda_L0*1.0e-8/UnitsCGS::p_UnitsCGS->f_r_e/UnitsCGS::p_UnitsCGS->f_Electron_ComptonWavelength*f_ExternalField_Intensity*f_ExternalField_Intensity*f_ExternalField_Duration*f_Laser_FocusW*f_Laser_FocusW);
}

double ExternalFields::CP_MC_LaserField_GetPhi(double x,double y,double z,double t){
	double a0,tao0;
    //case 16, for the first laser pulse
            
	tao0=f_ExternalField_Duration;
	double xprime,yprime,zprime,tprime,f_rprime2,z_r,w_zprime,phi,phiG,phiP,phiR,phi0;
	xprime=f_Laser_Kz*(x-f_Laser_FocusX)-f_Laser_Kx*(z-f_Laser_FocusZ);  // current x position related to the focus in the laser coordinate system
	yprime=y;
	zprime=f_Laser_Kx*(x-f_Laser_FocusX)+f_Laser_Kz*(z-f_Laser_FocusZ);  // current z position related to the focus in the laser coordinate system
	f_rprime2=yprime*yprime+xprime*xprime;
	z_r=PI*f_Laser_FocusW*f_Laser_FocusW;
	w_zprime=f_Laser_FocusW*sqrt(1+(zprime/z_r)*(zprime/z_r));
    tprime=t+f_Laser_LtoFocus;                                           // current laser peak distance to the focus
	phiG=atan(zprime/z_r);
	phiR=PI*f_rprime2*zprime/(z_r*z_r+zprime*zprime);
	phiP=2*PI*(zprime-tprime);
	phi0=f_Laser_IniPhase;
	phi=phi0+phiP-phiR+phiG;                                             // the phase factor
	
    return(phi);
}
double ExternalFields::CP_MC_LaserField_GetlocalEnvelope_nplocal(double x,double y,double z,double t){
	double a0,tao0,tp1;
	a0=f_ExternalField_Intensity;
	tao0=f_ExternalField_Duration;
	double xprime,yprime,zprime,tprime,f_envelopet,f_enveloper,f_rprime2,z_r,w_zprime,ipsilong,ksi,kxi,phi,phiG,phiP,phiR,phi0,Exprime,Eyprime,Ezprime;
	xprime=f_Laser_Kz*(x-f_Laser_FocusX)-f_Laser_Kx*(z-f_Laser_FocusZ); // current x position related to the focus in the laser coordinate system
	yprime=y;
	zprime=f_Laser_Kx*(x-f_Laser_FocusX)+f_Laser_Kz*(z-f_Laser_FocusZ);  // current z position related to the focus in the laser coordinate system
	f_rprime2=yprime*yprime+xprime*xprime;
	z_r=PI*f_Laser_FocusW*f_Laser_FocusW;
	w_zprime=f_Laser_FocusW*sqrt(1+(zprime/z_r)*(zprime/z_r));
    tprime=t+f_Laser_LtoFocus;                                           // current laser peak distance to the focus
	f_envelopet=exp(-(zprime-tprime)*(zprime-tprime)/tao0/tao0);         // temporal shape factor
    f_enveloper=f_Laser_FocusW/w_zprime*exp(-f_rprime2/w_zprime/w_zprime); 
	tp1=f_envelopet*f_enveloper; 
	tp1=tp1/(z_r*tao0*sqrt(PI/2.0)/2.0*f_Lambda_L0*1e-4*f_Lambda_L0*1e-4*f_Lambda_L0*1e-4); 
	if(i_Laser_Polarization==0 || i_Laser_Polarization==1)
		return(tp1);
	else if(i_Laser_Polarization==2)
		return(tp1*2.0);
	else{
		cout<<"linear polarized laser pulse has been used!!!"<<endl;
		return(tp1);
	}
}

double ExternalFields::CP_MC_PhotonK_Selection(){   // Mento Carlo samping of y=x-1 satisfy exp(-y^2/sigema^2) sigema=lambda0/piL,here L is the laser pulse length defined in electric field
	double sigema,y,u,v,w;
	sigema=1.0/((f_ExternalField_Duration/f_Density_Normalize)*PI);
Sampling: 

	loopg:
    u=((double)rand()/(double)RAND_MAX);
    v=((double)rand()/(double)RAND_MAX);
    w=(2.0*u-1)*(2.0*u-1)+(2.0*v-1)*(2.0*v-1);
    if(w>1.0) goto loopg;
    y=u*sqrt(-2.0*log(w)/w)*sigema;
	u=((double)rand()/(double)RAND_MAX);
	if(u<0.5)
		y=-y;
	if(y>-1){
		return(y+1);
	}
	else goto Sampling;
}
