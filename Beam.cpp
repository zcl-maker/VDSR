#include <iostream>
using namespace std;

#include <stdlib.h>
#include "vdsr.h"
#include <math.h>

Beam* Beam::p_Beam = NULL;
Beam::Beam(char *infile, int rank, long Pnumber, int numprocs,int flag): NList("Beam")
{
	long i,j;
    
    FILE *Finitial_position;
	Beam::p_Beam = this;
	i_my_rank=rank;
	l_Pnumber=Pnumber;
	f_sampling_data=new double[l_Pnumber];
	i_numprocs=numprocs;
    l_myPnumber=l_Pnumber/numprocs;
	l_start_position=l_myPnumber*i_my_rank;                               // first particle order

	if(numprocs-i_my_rank<l_Pnumber%numprocs+1){                          // Here means the last few processes will possibly take one particle more than the first few processes.
		l_myPnumber++; 
		l_start_position+=l_Pnumber%numprocs-(numprocs-i_my_rank);
	}
    
	// Generally goto statements are not very safe, I will redo with if ... {} (rsg) 
	// if(flag==0) goto endBeam;                                          // Here means the particle initial positions are read from an existed file, do not use beam initialization

    if(flag!=0)
    {
        f_using_data=new double[l_myPnumber*6];                              // used for save x0,y0,z0,px0,py0,pz0
        AddEntry("f_Qcharge", &f_Qcharge, 1.0e2);
        AddEntry("f_Energy", &f_Energy, 1.0e3);
        AddEntry("f_EnergySpread", &f_EnergySpread, 0.1);
        AddEntry("f_RadiusX", &f_RadiusX, 20.);
        AddEntry("f_RadiusY", &f_RadiusY, 10.);
        AddEntry("f_RadiusZ", &f_RadiusZ, 10.);
        AddEntry("f_RadiusR", &f_RadiusR, 0.0);      //Here the default value should set to be 0, otherwise ShapeX,Y,Z cannot be decided from the input file
        AddEntry("f_CenterX", &f_CenterX, 0.);
        AddEntry("f_CenterY", &f_CenterY, 0.);
        AddEntry("f_CenterZ", &f_CenterZ, 0.);
        AddEntry("f_CenterR", &f_CenterR, 0.);
        AddEntry("i_ShapeX", &i_ShapeX, 0);
        AddEntry("i_ShapeY", &i_ShapeY, 0);
        AddEntry("i_ShapeZ", &i_ShapeZ, 0);
        AddEntry("i_ShapeR", &i_ShapeR, 0);
        AddEntry("f_RadiusPX", &f_RadiusPX, 20.);
        AddEntry("f_RadiusPY", &f_RadiusPY, 10.);
        AddEntry("f_RadiusPZ", &f_RadiusPZ, 10.);
        AddEntry("f_RadiusPR", &f_RadiusPR, 0.);
        AddEntry("f_CenterPX", &f_CenterPX, 0.);
        AddEntry("f_CenterPY", &f_CenterPY, 0.);
        AddEntry("f_CenterPZ", &f_CenterPZ, 0.);
        AddEntry("f_CenterPR", &f_CenterPR, 0.);
        AddEntry("i_ShapePX", &i_ShapePX, 0);
        AddEntry("i_ShapePY", &i_ShapePY, 0);
        AddEntry("i_ShapePZ", &i_ShapePZ, 0);
        AddEntry("i_ShapePR",&i_ShapePR,0);
        AddEntry("f_ExpTemp",&f_ExpTemp,16);
        AddEntry("f_Exp_min",&f_Exp_min,24.463);
        AddEntry("f_Exp_max",&f_Exp_max,208.434);


        FILE *p_File = NULL;


        if (rank==NIL)
        {
            p_File = fopen(infile,"rt");     
            if (p_File == NULL)
            {
                cout<<"Beam no input file: "<<infile<<" "<<endl;
                exit (-1);
            }
        }

        if (p_File)
        {
            rewind(p_File); 
            read(p_File); 
        }  

        if(i_my_rank==0){

            cout<<"     Beam: Beam Charge="<<f_Qcharge<<"pC Energy="<<f_Energy<<"MeV EnergySpread="<<f_EnergySpread*100<<"%"<<endl;
            cout<<"     RadiusX="<<f_RadiusX<<"  RadiusY="<<f_RadiusY<<"  RadiusZ="<<f_RadiusZ<<"\n     CenterX="<<f_CenterX<<"  CenterY="<<f_CenterY<<"  CenterZ="<<f_CenterZ<<" unit of lambdap"<<endl;
            cout<<"     RadiusPX="<<f_RadiusPX<<"  RadiusPY="<<f_RadiusPY<<"  RadiusPZ="<<f_RadiusPZ<<"\n     CenterPX="<<f_CenterPX<<"  CenterPY="<<f_CenterPY<<"  CenterPZ="<<f_CenterPZ<<" unit of mec"<<endl;
            cout<<"     ShapePX="<<i_ShapePX<<"  ShapePY="<<i_ShapePY<<"  ShapePZ="<<i_ShapePZ<<" ShapeR="<<i_ShapeR<<endl;
            cout<<"\n****** Beam is creating ******\n\n"<<endl;
            Domain::p_D->out_Flog<<"     Beam: Beam Charge="<<f_Qcharge<<"pC Energy="<<f_Energy<<"MeV EnergySpread="<<f_EnergySpread*100<<"%"<<endl;
            Domain::p_D->out_Flog<<"     RadiusX="<<f_RadiusX<<"  RadiusY="<<f_RadiusY<<"  RadiusZ="<<f_RadiusZ<<"\n     CenterX="<<f_CenterX<<"  CenterY="<<f_CenterY<<"  CenterZ="<<f_CenterZ<<" unit of lambdap"<<endl;
            Domain::p_D->out_Flog<<"     RadiusPX="<<f_RadiusPX<<"  RadiusPY="<<f_RadiusPY<<"  RadiusPZ="<<f_RadiusPZ<<"  RadiusPR"<<f_RadiusPR<<"\n     CenterPX="<<f_CenterPX<<"  CenterPY="<<f_CenterPY<<"  CenterPZ="<<f_CenterPZ<<"  CenterPR="<<f_CenterPR<<" unit of mec"<<endl;
            Domain::p_D->out_Flog<<"     ShapePX="<<i_ShapePX<<"  ShapePY="<<i_ShapePY<<"  ShapePZ="<<i_ShapePZ<<" ShapeR="<<i_ShapeR<<"  ShapePR="<<i_ShapePR<<endl;
            Domain::p_D->out_Flog<<"\n****** Beam is creating ******\n\n";
        }   
#ifdef V_MPI
        CBuffer *buf = new CBuffer;
        buf->reset();
        pack_nls(buf);
        Domain::p_D->BroadCast(buf);
        if (rank) unpack_nls(buf);
        delete buf;

#endif
        if(i_my_rank==0){
            Finitial_position=fopen("Initial_Particle_XP.dat","w");
        }

        double *f_sampling_Theta; // define here, because it will be used for f_RadiusPR>0
        f_sampling_Theta=new double[l_Pnumber];
        /********************* Begin of Position Sampling ***************************/
        if (f_RadiusR>0){          // X,Y distributions are decided from r,\theta distribution
            Sampling(f_CenterR,f_RadiusR,l_Pnumber,i_ShapeR);  
            double *f_sampling_R;
            f_sampling_R=new double[l_Pnumber];
            for(j=0;j<l_Pnumber;j++){
                f_sampling_R[j]=f_sampling_data[j];
            }
            Sampling(PI,PI,l_Pnumber,3);
            for(j=0;j<l_Pnumber;j++)
                f_sampling_Theta[j]=f_sampling_data[j];

            for(j=0;j<l_Pnumber;j++){
                f_sampling_data[j]=fabs(f_sampling_R[j])*cos(f_sampling_Theta[j]);
            }
            if(i_my_rank==0){
                for(j=0;j<l_Pnumber;j++)
                    fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
                fprintf(Finitial_position,"\n");
            }
            for(i=0;i<l_myPnumber;i++)
                f_using_data[0+i*6]=f_sampling_data[l_start_position+i];

            for(j=0;j<l_Pnumber;j++)
                f_sampling_data[j]=fabs(f_sampling_R[j])*sin(f_sampling_Theta[j]);
            if(i_my_rank==0){
                for(j=0;j<l_Pnumber;j++)
                    fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
                fprintf(Finitial_position,"\n");
            }
            for(i=0;i<l_myPnumber;i++)
                f_using_data[1+i*6]=f_sampling_data[l_start_position+i];
            delete f_sampling_R;
        }
        else{                      // X,Y distributions are decided by themselves	   
            Sampling(f_CenterX,f_RadiusX,l_Pnumber,i_ShapeX);
            if(i_my_rank==0){
                for(j=0;j<l_Pnumber;j++)
				fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
			fprintf(Finitial_position,"\n");
            }

		for(i=0;i<l_myPnumber;i++)
			f_using_data[0+i*6]=f_sampling_data[l_start_position+i];   
        p_Beam->Sampling(f_CenterY,f_RadiusY,l_Pnumber,i_ShapeY);
		if(i_my_rank==0){
			for(j=0;j<l_Pnumber;j++)
				fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
			fprintf(Finitial_position,"\n");
		}
		for(i=0;i<l_myPnumber;i++)
			f_using_data[1+i*6]=f_sampling_data[l_start_position+i];
	}                          // end of else

	p_Beam->Sampling(f_CenterZ,f_RadiusZ,l_Pnumber,i_ShapeZ);
	if(i_my_rank==0){
		for(j=0;j<l_Pnumber;j++)
			fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
		fprintf(Finitial_position,"\n");
	}
	for(i=0;i<l_myPnumber;i++)
		f_using_data[2+i*6]=f_sampling_data[l_start_position+i];

    /********************* End of Position Sampling ***************************/
    /********************* Begin Momentum Sampling ***************************/
        
    // Reseeding the random number generator (rsg), not physical in the bubble regime due to symmetry
    // The particles should not spiral
    //    srand(13124131);    // let's reseed the random number generator // rsg
    //    Sampling(PI,PI,l_Pnumber,3);
    //    for(j=0;j<l_Pnumber;j++)
    //        f_sampling_Theta[j]=f_sampling_data[j];
    // end reseeding
    
        
    

    
    if (f_RadiusPR>0){          // PX,PY distributions are decided from PR, in radius direction
		p_Beam->Sampling(f_CenterPR,f_RadiusPR,l_Pnumber,i_ShapePR);

		double *f_sampling_PR;
		f_sampling_PR=new double[l_Pnumber];
        // p_perp
		for(j=0;j<l_Pnumber;j++){
		    f_sampling_PR[j]=f_sampling_data[j];
		}

        double *f_sampling_thetasign;                         // use this to make sure the electrons initially |px/py|=|x/y|
		f_sampling_thetasign=new double[l_Pnumber];
	    Sampling(0.5,0.5,l_Pnumber,3);
	    for(j=0;j<l_Pnumber;j++)
            f_sampling_thetasign[j]=(f_sampling_data[j]>0.5)?1.0:0.0;
        
		// px
		for(j=0;j<l_Pnumber;j++){			
			f_sampling_data[j]=fabs(f_sampling_PR[j])*cos(f_sampling_Theta[j]+f_sampling_thetasign[j]*PI);
	    }

		if(i_my_rank==0){
			for(j=0;j<l_Pnumber;j++)
				fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
			fprintf(Finitial_position,"\n");
		}
		for(i=0;i<l_myPnumber;i++)
			f_using_data[3+i*6]=f_sampling_data[l_start_position+i];

        // py
		for(j=0;j<l_Pnumber;j++){
			f_sampling_data[j]=fabs(f_sampling_PR[j])*sin(f_sampling_Theta[j]+f_sampling_thetasign[j]*PI);
	    }
		
		if(i_my_rank==0){
			for(j=0;j<l_Pnumber;j++)
				fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
			fprintf(Finitial_position,"\n");
		}
		for(i=0;i<l_myPnumber;i++)
			f_using_data[4+i*6]=f_sampling_data[l_start_position+i];
		delete f_sampling_PR;
	}
	else{
		p_Beam->Sampling(f_CenterPX,f_RadiusPX,l_Pnumber,i_ShapePX);
		if(i_my_rank==0){
			for(j=0;j<l_Pnumber;j++)
				fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
			fprintf(Finitial_position,"\n");
		}
		for(i=0;i<l_myPnumber;i++)
			f_using_data[3+i*6]=f_sampling_data[l_start_position+i];

		p_Beam->Sampling(f_CenterPY,f_RadiusPY,l_Pnumber,i_ShapePY);
		if(i_my_rank==0){
			for(j=0;j<l_Pnumber;j++)
				fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
			fprintf(Finitial_position,"\n");
		}
		for(i=0;i<l_myPnumber;i++)
			f_using_data[4+i*6]=f_sampling_data[l_start_position+i];
	}

	p_Beam->Sampling(f_CenterPZ,f_RadiusPZ,l_Pnumber,i_ShapePZ);
	if(i_my_rank==0){
		for(j=0;j<l_Pnumber;j++)
			fprintf(Finitial_position,"%15.5e",f_sampling_data[j]);
		fprintf(Finitial_position,"\n");
	}
	for(i=0;i<l_myPnumber;i++)
		f_using_data[5+i*6]=f_sampling_data[l_start_position+i];
    delete f_sampling_Theta;

	/********************** End of Momentum Sampling ************************/

	if(i_my_rank==0)
		fclose(Finitial_position);           // end of initial position initialization
//endBeam: ;   // don't need it anymore (rsg)
    }
	if(i_my_rank==0){
		cout<<"****** Beam initial position and momentum distribution is initialized! *****"<<"\n End of Beam initialization"<<endl;
		Domain::p_D->out_Flog<<"****** Beam initial position and momentum distribution is initialized! *****"<<"\n End of Beam initialization"<<endl;
	}
};


void Beam::Sampling(double f_Center,double f_Radius,long l_number, int i_shape){
    long i,j;
	switch(i_shape){
		 case 0:{        // Gaussian distribution:    exp[-(x-f_CenterX)^2/(2*f_RadiusX^2)]
             double x1, x2, w;
             j=0;
			 do{
				do {
					x1 = 2.0 * rand()/RAND_MAX - 1.0;
					x2 = 2.0 * rand()/RAND_MAX - 1.0;
					//x1 = 2.0 * F_rand() - 1.0;
					//x2 = 2.0 * F_rand() - 1.0;
					w = x1 * x1 + x2 * x2;
				} while ( w >= 1.0 );
				w = sqrt( (-2.0 * log( w ) ) / w );
				f_sampling_data[j++] = f_Center+x1*w*f_Radius;
			 } while (j<l_number);
			 break;
			 }
		 case 1:{         // Square distribution
			 double df=2.0*f_Radius/l_number;
			 f_Center=f_Center+df/2.0;
			 for(i=0;i<l_number;i++)
				 f_sampling_data[i]=f_Center-f_Radius+df*i;
			 break;
			 }
		 case 2:{         // Constant
			 for(i=0;i<l_number;i++)
				 f_sampling_data[i]=f_Center;
			 break;
			 }
		case 3:{         // Uniform but random selection
			 for(i=0;i<l_number;i++)
				 f_sampling_data[i]=f_Center-f_Radius+2.0*rand()/RAND_MAX*f_Radius;
			 break;
			 }
	   	case 4:{         // exponential distribution for E_k
			 double *f_Exp_bin,*f_Number_bin,f_Total_bin=0,f_dExp;
             int Binnumber=l_number/32,l_Bin_order=0;
			 double f_Pbin_order=0,f_Exp_binCenter;
			 f_Number_bin=new double[Binnumber];
			 f_Exp_bin=new double[Binnumber];
			 f_dExp=(f_Exp_max-f_Exp_min)/(Binnumber-1.0);  
			 for(int j=0;j<Binnumber;j++){
				f_Exp_bin[j]=f_Exp_min+j*f_dExp;
				f_Exp_binCenter=f_Exp_bin[j]+0.5*f_dExp;
				f_Number_bin[j]=(f_Exp_binCenter/sqrt(1+f_Exp_binCenter*f_Exp_binCenter))*exp(-1/f_ExpTemp*(sqrt(1+f_Exp_binCenter*f_Exp_binCenter)-1));
				f_Total_bin+=f_Number_bin[j];
			 }
			 for(int j=0;j<Binnumber;j++)
				f_Number_bin[j]=f_Number_bin[j]/f_Total_bin*l_number;
             
			 for(int i=0;i<l_number;i++){
				f_sampling_data[i]=f_Exp_bin[l_Bin_order]+1.0*(i-f_Pbin_order)/(f_Number_bin[l_Bin_order])*f_dExp;
				if(i>f_Pbin_order+f_Number_bin[l_Bin_order]){
					f_Pbin_order+=f_Number_bin[l_Bin_order];
					l_Bin_order++;
				}
			 }
			 delete f_Number_bin;
			 delete f_Exp_bin;
			 break;
			 }
			   			   
		case 5:{        // radial Gaussian distribution: r0 satisfies distribution as:   r*exp[-r^2/(f_Radius^2)], f_RadiusX=f_RadiusY=f_Radius/sqrt(2)
             double x1, x2, w, f_RadiusX,r0;
			 f_RadiusX=f_Radius/sqrt(2.0);  // to sample r0, we sample x and y first then use r0=sqrt(x^2+y^2) to get r0, x,y should satisfy distribution like exp(-x^2/2f_RadiusX^2)
             j=0;
			 do{
				do {
					x1 = 2.0 * rand()/RAND_MAX - 1.0;
					x2 = 2.0 * rand()/RAND_MAX - 1.0;
					//x1 = 2.0 * F_rand() - 1.0;
					//x2 = 2.0 * F_rand() - 1.0;
					w = x1 * x1 + x2 * x2;
				} while ( w >= 1.0 );
				w = sqrt( (-2.0 * log( w ) ) / w );
                r0=sqrt(x1*x1+x2*x2)*w*f_RadiusX;
				f_sampling_data[j++] = f_Center+r0;
			 } while (j<l_number); 
			 break;
			 }
			 
	}
};
	 
double Beam::F_rand()
{
  long a=16807, m=2147483647, q=127773, r=2836;
  long low, high;
  static long seed=31207321;
  double fnumb;
  
  high=seed/q;
  low=seed-q*high;
  seed=a*low-r*high;
  if(seed<=0) seed=seed+m;
  fnumb=(double)seed/2147483646.0;
  return(fnumb);
}

