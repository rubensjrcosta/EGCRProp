//Dynamics.cc/

//local headers/
#include "Dynamics.h"
#include "Geometry.h"
#include "Mathematics.h"

//constructor of the class/
Dynamics::Dynamics ( )       
{
  
}

//destructor of the class/
Dynamics::~Dynamics ( )
{
        
}

//methods of the Dynamics class//
double Dynamics::PusherRK4th(double p[3],double p0[3],double x[3], double x0[3],double* Gamma,double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t){Geometry * Geom = new Geometry();
							
	double *u,*u0;
	
	u0 = new double[ODE_num];
		
	//initializing the position and momentum
	for(i=0;i<=2;i++){//2
		u0[i]=x0[i];	
		u0[i+3]=p0[i];
	}//2
		
	// Otherwise, advance to time t0 + dt, and have RK4 estimate the solution u there.
	u = IntegrationRK4th(Gamma,u0, Delta_t, uprimeRK4th);
	
	//shifting the data to prepare for another step.
	for (i=0;i<ODE_num; i++){u0[i] = u[i];}
	
	for (i=0;i<=2;i++){//2
		x[i] = u0[i];
		p[i] = u0[i+3]; 
	}//2
	
	*Delta_x = Geom->DistPP(x,x0);
	*D +=*Delta_x;
	
	GammaRK4th(Gamma,Geom->Mod(p),m0);
	
	delete Geom;delete [] u0;delete [] u;return 1.0;
}//1

double Dynamics::PusherBoris(double u[3],double u0[3],double x[3], double x0[3],double* Gamma,double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t){Geometry * Geom = new Geometry();

	UpdateVelocityBoris(u,u0,Gamma,Delta_t);
	PushParticleBoris(x,x0,u,Gamma,Delta_t);
	
	*D += Geom->DistPP(x,x0);
	*Delta_x = Geom->DistPP(x,x0);
	
	GammaBoris(Gamma,Geom->Mod(u0));	
		
 delete Geom;return 1.0;
}//1

double Dynamics::PusherVay(double u[3],double u0[3],mpfr::mpreal v[3],mpfr::mpreal v0[3],double x[3], double x0[3],double* Gamma,double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t){Geometry * Geom = new Geometry();
	
	UpdateVelocityVay(u,u0,v0,Gamma,B,Delta_t);
	PushParticle_Vay(x,x0,v0,Delta_t);
	
	*D += Geom->DistPP(x,x0);
	*Delta_x = Geom->DistPP(x,x0);
	
	//GammaVay(Gamma,Geom->Mod(u0));	
	
	delete Geom;return 1.0;
}//1  

double Dynamics::PusherNewEuler(mpfr::mpreal v[3],mpfr::mpreal v0[3],double x[3], double x0[3],double* Gamma,mpfr::mpreal beta[3],mpfr::mpreal a_L[3],double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t){Geometry * Geom = new Geometry();
	
	
	LorentzAcceleration(a_L,v,B);
	UpdateVelocityNewEuler(v,v0,Gamma,a_L,Delta_t);
	PushParticleNewEuler(x,x0,v0,Gamma,a_L,Delta_t);
	BetaNewEuler(beta,E);
	LorentzFactorNewEuler(Gamma,beta); 
				
	
	*D += Geom->DistPP(x,x0);
	*Delta_x = Geom->DistPP(x,x0);
	
	//GammaVay(Gamma,Geom->Mod(u0));	
	
	delete Geom;return 1.0;
}//1  

double Dynamics::Energy(){//1//particle's energy calculator/

	Mathematics * Math =  new Mathematics();Geometry * Geom = new Geometry();
	
	double XI_aux;	
	double sigl=-1;
		
	if(strcmp(ENERGYLOSS, "OFF") == 0){E = (Gamma)*m0*pow(c,2);}
		else{//2
		if(strcmp(BTRACK,"ON") == 0){sigl = -sigl;}	
				
		if((E) <= 0.016375556){//3
			E = E*(1.0 + sigl*(delta_x/(mbyMpc))*(pow(XI_extr[0],-1)/pow((1+(z)),3)+ pow(XI_ad,-1)/pow((1+(z)),-1.5)));	
			Gamma=E/(m0*pow(c,2));		
		}//3	
		else if((E)*(1+(z)) >= 487.0040579066){//3
			(void)printf("energy cut %lf\n", E/(JbyEeV));
			E = E*(1.0 + sigl*(delta_x/(mbyMpc))*(pow(XI_extr[1],-1)/pow((1+(z)),3)+ pow(XI_ad,-1)/pow((1+(z)),-1.5)));
		}//3		
		else{//3	
			for (i=0;i<=line_loss-2;i++){//4						
				if((E)*(1+(z))>=E_z0[i] && (E)*(1+(z)) <=E_z0[i+1]){//5
					XI_aux= Math->Interpolate(E,E_z0[i] , XI_z0[i], E_z0[i+1],  XI_z0[i+1]);
						//(void)printf("qq %lf\n", (1.0/(E_z0[i]-E_z0[i+1])));
						//(void)printf("E %lf\nEloss+1 %lf Xloss+1 %lf\n  Eloss %lf Xloss%lf\n%lf",*E, E_z0[i+1], XI_z0[i+1],E_z0[i], XI_z0[i],XI_aux);
				}//5
			}//4
			E = E*(1.0 + sigl*(delta_x/(mbyMpc))*(pow(XI_aux,-1)/pow((1+(z)),3)+ pow(XI_ad,-1)/pow((1+(z)),-1.5)));
			Gamma=E/(m0*pow(c,2));		

		}//3
	}//2
	
	//(void)printf("%.12f\n", (*E)*(JbyEeV));		

	delete Math;delete Geom;
	return 1.0;
}//1
double Dynamics::SpreadingAngle(){//1 //calculating the sum of values of the particle's energy, redshift, traveled distance and time
   	
   	//the deflection angle after crossing a single cell/ 
	defl_1 = (l_coh/r_L)*(degbyrad)*2/3;
	if(defl_1>=180.0)defl_1 =180.0-fmod(defl_1,180.0);
		
//@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE//
	//N_aux = (DC_avg*(mbyMpc))/l_coh;		
	N_aux = (D_avg)/l_coh;		
	//(void)printf("N_aux %.3lf\n",N_aux);
//@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE@@@@@@@@IMPORTANTE//

	//the total deflection angle after crossing  N cell/
	defl_N = defl_1*(sqrt(N_aux));					
	if(defl_N>=180.0)defl_N =180.0-fmod(defl_N,180.0);
	
	return 1.0;
}//1
double Dynamics::CyclotronFrequency(){//1 //particle's cyclotron frequency
	
	Geometry * Geom = new Geometry();
	
	omega_B = (Geom->ScalMod(q)/(m0*(Gamma)))*(Geom->Mod(B)); 
					
	delete Geom;					
	return 1.0;
}//1
double Dynamics::StepTime(){//1 //stop
	
	delta_t=t_scal*((pi*2)/(omega_B)); //step time: Dt = xi*(2*pi)/omega_B
						
	return 1.0;
}//1


