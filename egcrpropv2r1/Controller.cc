//Controller.cc/

//local headers/
#include "Controller.h"
#include "InputParameters.h"
#include "Particle.h"
#include "Magnetic.h"
#include "Cosmology.h"
#include "Dynamics.h"
#include "Outputfiles.h"
#include "Mathematics.h"
#include "Geometry.h"

#include "IntegratorRK4th.h"
#include "IntegratorBoris.h"
#include "IntegratorVay.h"
#include "IntegratorNewEuler.h"

//constructor of the class/
Controller::Controller ( ){

}

//destructor of the class/
Controller::~Controller ( ){

}

//methods of the Controller class
double Controller::LoadingInputParameters(){//1
	
	time_start = clock();//initiating the CPU time	
	
	//@@@@@ATENTION@@@@ATENTION@@@@@ATENTION@@@@ATENTION@@@@@ATENTION@@@
	mpreal::set_default_prec(mpfr::digits2bits(digits)); //setup default precision in bits for all subsequent computations
	
	InputParameters* Input=new InputParameters();Particle* Part=new Particle();
	Magnetic*Magn=new Magnetic();Cosmology* Cosm=new Cosmology();Outputfiles* Outp=new Outputfiles();
		
	Input->ReadingInputParameters();
	Part->SettingParticleParameters();/*setting the particle's parameters*/
	Magn->SettingMagneticField();/*if MF_model is 1 or 2 loading the tensors Theta and Phi from the file theta_phi_tensors.dat*/		
	Cosm->SettingCosmologicalParameters();
	if(loop==0)Outp->OpeningSummary();
	
	//the standard Unix random number generator 
	srand( (unsigned)time(NULL) ); //setting the 'random seed' based on time	
	
	//the GSL random number generator
	strcpy (RANDGENERATOR, gsl_rng_name (rng));	
    if(r_seed==0){//2
		gsl_rng_set(rng, time(NULL)); //setting the 'random seed' based on time
	}//2
	else gsl_rng_set(rng, mySeedDIRETIONS); //setting the 'random seed' based on time
	
	delete Input;delete Part;
	delete Magn;delete Cosm;delete Outp;
	return 1.0;
}//1
double Controller::SettingInitialConditions(){//1 
   
   	Geometry* Geom=new Geometry();Dynamics* Dyns=new Dynamics();   	 		 
	
	E0 = (JbyEeV)*(E_EeV)*Z_aux; //transforming the energy from [EeV] to [J] 
	E=E0;//energy initialization/

	if(strcmp(INTMETHOD, "Runge-Kutta4th") == 0){//2 
		p0_mod=E0/(c); //particle's momentum: p = E/c				
		GammaRK4th(&Gamma,p0_mod, m0);//particle's Lorentz Factor
		Dyns->CyclotronFrequency(); //particle's cyclotron frequency
		Dyns->StepTime(); //step time: Dt = xi*(2*pi)/omega_B
		r_L=p0_mod/((m0*(Gamma))*(omega_B)); //particle's Larmor radius [Mpc]
	}//2		
	else if(strcmp(INTMETHOD,"Boris") == 0){//2 
		u0_mod=E0/(m0*c); //u = Gamma*velocity = p/m0		
		GammaBoris(&Gamma,u0_mod);				
		Dyns->CyclotronFrequency(); //particle's cyclotron frequency
		Dyns->StepTime(); //step time: Dt = xi*(2*pi)/omega_B 
		r_L=u0_mod/((omega_B)*(Gamma));		
	}//2
	else if(strcmp(INTMETHOD,"Vay") == 0){//2
		BetaVay(&beta,&E0); //particle's velocity parameter: beta = sqrt(1-(m0*c²/E)²);	
		LorentzFactorVay(&Gamma,&beta);//(void)printf("Gamma(beta) = %.16e\n", Gamma);
		v0_mod = beta*c; 		
		u0_mod=((v0_mod)*(Gamma)).toDouble(); //u = Gamma*velocity = p/m0			
		/*Dyns->Gamma_BorisnewLeapFrog(&Gamma,u0_mod);*/ //(void)printf("Gamma(u0) = %.16e\n", Gamma);
		Dyns->CyclotronFrequency(); //particle's cyclotron frequency
		Dyns->StepTime(); //step time: Dt = xi*(2*pi)/omega_B 
		r_L=(v0_mod/((omega_B))).toDouble();
	}//2
	else if(strcmp(INTMETHOD,"NewEuler") == 0){//2 //New Euler method / 		
		BetaNewEuler(&beta,&E0);
		v0_mod = beta*c; //initial velocity modulus: v = beta*c 
		LorentzFactorNewEuler(&Gamma,&beta);		
		Dyns->CyclotronFrequency(); //particle's cyclotron frequency
		Dyns->StepTime(); //step time: Dt = xi*(2*pi)/omega_B  		
		r_L=((v0_mod)/(omega_B)).toDouble();	
	}//2
								
	delete Geom;delete Dyns;							
	return 1.0;
}//1

double Controller::SettingRange(){//1 //setting the range of the initial conditions
   
	Cosmology* Cosm=new Cosmology();Outputfiles* Outp=new Outputfiles();
	
	Cosm->SettingRedShift();
	
	//checking if the initial energy is in the energy range for EGCRProp run
	if((E0*(1+z0)) > E_max ){//1															
		(void)printf("** Wrong choice for requested energy  E*(1+z) [EeV] = %lf (z= %lf) > 3.0 ZeV  **\n",(E0*(1+z0)/(JbyEeV)),z0);exit(EXIT_FAILURE);					
	}//1
	if( (E0*(1+z0)) < E_min){//1															
		(void)printf("** Wrong choice for requested energy E*(1+z) [EeV] = %lf (z= %lf) < 0.1 EeV **\n", ((E0*(1+z0))/(JbyEeV)),z0);exit(EXIT_FAILURE);			
	}//1
	
	D_max = (mbyMpc)*D_Mpc;//the maximum distance/							
	Outp->PrintingSummaryBegin(); //printing the summary at the begin of the simulations				s				
	Outp->OpeningDynamicsParameters(); //opening the output files	
	SettingAverage(); //setting the average values
				
	delete Cosm;delete Outp;	 							
	return 1.0;
}//1
double Controller::SettingParticlePosition(){//1
   
	Mathematics * Math = new Mathematics();Outputfiles* Outp=new Outputfiles();   	
	
	if(N_iter<=N_traj){//1
		Outp->OpeningTrajectory(); //opening the trajectory file		
	}//1   
				    				
	for(i=0;i<=2;i++){//2 //initializing the position at the origin and the maximum and minimum values for the plotting											
		x0[i]    = Origin;
		O[i]     = x0[i];		
		x_max[i] = -100*x0[i];
		x_min[i] = 100*x0[i];			
	}//2	
	
	//random directions for the initial p,v or u vectors/	
	//GSL's generator
	theta =(Math->RandNum_gsl(0,pi,rng));
	phi = (Math->RandNum_gsl(0,pi*2,rng));	
	
	//Standard generator
	//theta =(Math->RandNum_std(0,pi));
	//phi = (Math->RandNum_std(0,pi*2));
		
	//////DANGER///DANGER///DANGER//DANGER///DANGER///DANGER///DANGER///
	////testing
	//theta =pi;
	//phi = 0;	
	//////DANGER///DANGER///DANGER//DANGER///DANGER///DANGER///DANGER///
	
	if(strcmp(INTMETHOD, "Runge-Kutta4th") == 0){//2 //initial momentum/			
		p0[0] = (p0_mod)*sin(theta)*cos(phi);		
		p0[1] = (p0_mod)*sin(theta)*sin(phi);								
		p0[2] = (p0_mod)*cos(theta);								
		for(i=0;i<=2;i++){//3 //reference momentum/								
			p0_ref[i]= p0[i];				
		}//3		
	}//2					
	if(strcmp(INTMETHOD,"NewEuler") == 0 || strcmp(INTMETHOD,"Vay") == 0){//2 //initial velocity/					
		v0[0] = (v0_mod)*sin(theta)*cos(phi);		
		v0[1] = (v0_mod)*sin(theta)*sin(phi);								
		v0[2] = (v0_mod)*cos(theta);								
		for(i=0;i<=2;i++){//3 //reference velocity								
			v0_ref[i]= v0[i];				
		}//3		
	}//2			 								
	if(strcmp(INTMETHOD,"Boris") == 0 || strcmp(INTMETHOD,"Vay") == 0){//2 //initial u vector: u0/				
		u0[0] = (u0_mod)*sin(theta)*cos(phi);		
		u0[1] = (u0_mod)*sin(theta)*sin(phi);								
		u0[2] = (u0_mod)*cos(theta);							
		for(i=0;i<=2;i++){//3 //reference u vector/								
			u0_ref[i]= u0[i];				
		}//3
		
	}//2			
						
	delete Math;delete Outp;								
	return 1.0;
}//1 
double Controller::SettingDynamicsParameters(){//1 //initializing the particle's energy, redshift, traveled distance and time
   	 
	E=E0;
	z=z0;			
	for(i=0;i<=2;i++){//2 //initializing the position at the origin 											
		x[i] = x0[i];									
	}//2
	D=0.0;
	t= 0.0;	
	step=0;
	
	return 1.0;
}//1
double Controller::SettingAverage(){//1//initializing  the variables  for the average calculation		
	
	E_avg=0.0;
	D_avg=0.0;
	D0_avg=0.0;
	DC_avg=0.0;
	DbyD0_avg=0.0;
	defl_avg=0.0;
	
	//initializing the counts of the cutoffs in redshift and energy
	zcut_counts=0;
	E1cut_counts=0;			
	E2cut_counts=0;
		
	return 1.0;
}//1

double Controller::ParticleTracking(){//1
	
	Outputfiles* Outp=new Outputfiles();Magnetic* Magn=new Magnetic();Cosmology* Cosm=new Cosmology();Dynamics* Dyns=new Dynamics();   
	
	SettingParticlePosition();//initializing: initial position with random  orientation//																				
	SettingDynamicsParameters();// E,z,x,D,t,step;	
				
	do{//Loop: trajectory//
		Outp->PrintingTrajectory();//printing the trajectories/
		Magn->MagneticField(); //(void)printf("%.10e %.10e %.10e\n",B[0],B[1], B[2]);
	
		if(strcmp(INTMETHOD, "Runge-Kutta4th") == 0){//2			
			Dyns->PusherRK4th (p, p0,x,x0,&Gamma,&E, B, q, m0, &D, &delta_x,delta_t);									
			
			for (j=0;j<=2;j++){//3	
					p0[j] = p[j];
					x0[j]=x[j];	
			}//3	
		
		}//2
		else if(strcmp(INTMETHOD,"Boris") == 0){//2
			
			Dyns->PusherBoris(u,u0,x,x0,&Gamma,&E,B,q,m0,&D,&delta_x,delta_t);

			for (j=0;j<=2;j++){//3	
				u0[j] = u[j];
				x0[j]=x[j];							
			}//3
		}//2									
		else if(strcmp(INTMETHOD,"Vay") == 0){//2
				
		Dyns->PusherVay(u,u0,v,v0,x,x0,&Gamma,&E,B,q,m0,&D,&delta_x,delta_t);

		//cout << "BETA FACTOR = " << beta.toString()<< "\n";

			for (j=0;j<=2;j++){//3	
				v[j]  = u[j]/Gamma;
				v0[j] = v[j];
				u0[j] = u[j];
				x0[j] = x[j];
			}//3
		}//2
		else{//2			
		
			Dyns->PusherNewEuler(v,v0,x,x0,&Gamma,&beta,a_L,&E,B,q,m0,&D,&delta_x,delta_t);
			//cout << "BETA FACTOR = " << beta.toString()<< "\n";

			for (j=0;j<=2;j++){//3	
				v0[j] = v[j];	
				x0[j]=x[j];		
			}//3	
		}//2
		Dyns->CyclotronFrequency(); //particle's cyclotron frequency
		Dyns->StepTime(); //step time: Dt = xi*(2*pi)/omega_B
		Dyns->Energy();
		
		Cosm->CalculatingCosmologicalParameters();//Calculating z, DC
	
		//REVISARRRRRRRRRRRRRR
		if(strcmp(ENERGYLOSS, "ON") == 0){CuttingTracking();}
		//(void)printf("travelled_cells %d, step %d\n",travelled_cells, step);
		StoppingTracking();
	
	if(cut!=0) D_aux = (mbyMpc)*D_Mpc;//cut in redshift or energy// 
	}//Loop: trajectory//
	while(D_aux<=D_max);////deflectionXenergy//	
	//while(step<=100);////step of the simulation//
	//while(t<=t_H);//cut in Hubble time// 
	CalculatingDynamicsParameters();// E, D, D0, DC, D/D0, defl//
	Outp->ClosingTrajectory();
						
	delete Outp;delete Magn;delete Cosm;delete Dyns;				
	return 1.0;
}//1

double Controller::CalculatingDynamicsParameters(){//1
   	 
   	 Geometry * Geom = new Geometry();   	 		 
   	 Outputfiles* Outp=new Outputfiles();

   	 D0    = Geom->DistPP(x,O);
   	 DbyD0 = D/Geom->DistPP(x,O);	
		if(strcmp(INTMETHOD, "Runge-Kutta4th") == 0){//2			
		defl  = (Geom->Angle2Vec(p,p0_ref));
	}//2
	else if(strcmp(INTMETHOD,"NewEuler") == 0){//2			
		defl  = (Geom->Angle2Vec(v,v0_ref)).toDouble();
	}//2
	else if(strcmp(INTMETHOD,"Boris") == 0 ||  strcmp(INTMETHOD,"Vay") == 0){//2			
		defl  = (Geom->Angle2Vec(u,u0_ref));
	}//2	
	Outp->PrintingDynamicsParameters();

	E_avg+=E;
	D_avg+=D;
	D0_avg+=D0;
	DC_avg+=D_C;
	DbyD0_avg+=DbyD0;
	defl_avg+=defl;
	//SettingValues();				
	delete Geom;delete Outp;										
	return 1.0;
}//1
double Controller::CuttingTracking(){//1 
	cut=0;
	if(z<0.0){//1 //cut in redshift/				
		z=0.0;
		//(void)printf("the redishift < 0.0, z = %.5e (DC = %.3lf)\n", z, (D_C/Mpc));
		zcut_counts++;									
		cut=1;
	}//1
	if((E*(1+z))/(JbyEeV) > E_max ){//1															
		//(void)printf("** Wrong choice for requested energy  E*(1+z) [EeV] > 3.0 ZeV  **\n");
		//(void)printf("** E*(1+z) [EeV] = %lf (z= %lf) **\n", (E*(1+z)/(JbyEeV)),z);		
		E1cut_counts++;
		cut=1;	
	}//1
	if( (E*(1+z))/(JbyEeV) < E_min){//1															
		//(void)printf("** Wrong choice for requested energy E*(1+z) [EeV] < 0.1 EeV **\n");
		//(void)printf("** E*(1+z) [EeV] = %lf (z= %lf) **\n", ((E*(1+z))/(JbyEeV)),z);		
		E2cut_counts++;
		cut=1;	
	}//1 
						
	return 1.0;
}//1
double Controller::StoppingTracking(){//1 //stop
	if(strcmp(DSTOP, "traveled distance") == 0){//2
		D_aux=D;
	}//2			
	else if(strcmp(DSTOP,"trajectory radius") == 0){//2
		D_aux=(mbyMpc)*D_C;
	}//2
	
	t += delta_t;
	step++;	 					
	return 1.0;
}//1
double Controller::CalculatingAverage(){//1

	Dynamics* Dyns=new Dynamics();Outputfiles* Outp=new Outputfiles();
	
	E_avg = E_avg/N_evts;		
	D_avg = D_avg/N_evts;		
	D0_avg = D0_avg/N_evts;
	DC_avg =DC_avg/N_evts;		
	DbyD0_avg = DbyD0_avg/N_evts;
	Dampl=(DbyD0_avg - 1)*100;
	
	defl_avg = defl_avg/N_evts;
	Dyns->SpreadingAngle();
	Outp->PrintingSummaryBody();
	Outp->ClosingDynamicsParameters();
	if(loop==0)Outp->PrintingSummaryEnd();
			
	delete Dyns;delete Outp;
	return 1.0;
}//1






