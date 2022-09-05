
//local headers/
#include "Cosmology.h"
#include  "Geometry.h"

//constructor of the class/
Cosmology::Cosmology ( ){

}

//destructor of the class/
Cosmology::~Cosmology ( ){

}

//cosmological parameters/
struct cosmology_parameters { double H0; double OmegaM; double OmegaL; double OmegaK; };

//cosmological function/
double cosmology_function(double z, void * p) {//1
	
	struct cosmology_parameters * parameters = (struct cosmology_parameters *)p;
     
	double H0 = (parameters->H0);
	double OmegaM = (parameters->OmegaM);	
	double OmegaL = (parameters->OmegaL);
	double OmegaK = (parameters->OmegaK);

   return  1.0/(pow(OmegaL + pow(1.0+z,3.0)*OmegaM  + pow(1.0+z,2.0)*OmegaK, 0.5));
}//1

//methods of the Cosmology class/
double Cosmology::SettingCosmologicalParameters(){//1

	H0 = 70000.0;//Hubble constante [m/s/Mpc]
	OmegaM = 0.3;//total matter density/
	OmegaL = 0.7;//dark energy density/	
	OmegaK = 1.0 - (OmegaM+OmegaL);//density parameter (the Universe curvature)/
	cosmo_parameters[0]=H0,cosmo_parameters[1]=OmegaM,cosmo_parameters[2]=OmegaL,cosmo_parameters[3]=OmegaK;

	HubbleTime(&t_H, H0);
	HubbleDistance(&D_H, H0);
		
return 1.0;
}//1

//Hubble time/
double Cosmology::HubbleTime(double* t_H, double H0){//1

*t_H = (1.0/H0)*(mbyMpc);

return 1.0;
}//1

//Hubble distance/
double Cosmology::HubbleDistance(double* D_H, double H0){//1

*D_H = c/H0;

return 1.0;
}//1

double Cosmology::SettingRedShift(){//1
	
	if(b_track == 0){//2 //setting the redshift z at the shot point				
		RedShift(&z0, &D_H, D_Mpc);					
	}//2			
	else{//2					
		z0=0.0;	
	}//2
	
return 1.0;
}//1

//redshift/
double Cosmology::RedShift(double* z, double* D_H, double radial_dist){//1

*z = radial_dist/(*D_H);

return 1.0;
}//1


double Cosmology::CalculatingCosmologicalParameters(){//1
	
	Geometry 	* Geom = new Geometry();

	RedShift(&z, &D_H, Geom->DistPP(x,O)/(mbyMpc));			
	ComovDist(&D_C,  &D_H,&z,Geom->DistPP(x,O)/(mbyMpc), cosmo_parameters);
	if(b_track==0){/*2*/z=fabs(z0-z);}/*2*/

	delete Geom;
return 1.0;
}//1


//comoving distance/
double Cosmology::ComovDist(double* D_C,double* D_H, double* z,double radial_dist, double cosmo_parameters[4]){//1
	
	double H0 =cosmo_parameters[0];
	double OmegaM = cosmo_parameters[1] ;
	double OmegaL = cosmo_parameters[2] ;
	double OmegaK = cosmo_parameters[3];
	
	RedShift(z,D_H,radial_dist);
	
	gsl_integration_workspace * work_ptr = gsl_integration_workspace_alloc (1000);
  
	double lower_limit = 0.0, abs_error = 1.0e-8, rel_error = 1.0e-8, Iz, error;
    
	
	gsl_function FUNCTION;
	
	struct cosmology_parameters parameters = { H0, OmegaM, OmegaL, OmegaK};  
 
  
	FUNCTION.function = &cosmology_function;
	FUNCTION.params = &parameters;

    gsl_integration_qags(&FUNCTION, lower_limit, *z, abs_error, rel_error, 1000, work_ptr, &Iz, &error);

	*D_C = *D_H*Iz;
	
	gsl_integration_workspace_free(work_ptr);

  return 1.0;
}//1

//transverse comoving distance/
double Cosmology::TransverseDist(double* D_M, double* D_C, double* D_H,  double* z, double radial_dist, double cosmo_parameters[4]){//1

	double H0 =cosmo_parameters[0];
	double OmegaM = cosmo_parameters[1] ;
	double OmegaL = cosmo_parameters[2] ;
	double OmegaK = cosmo_parameters[3];
	
	RedShift(z,D_H,radial_dist);
	ComovDist(D_C,D_H, z,radial_dist, cosmo_parameters);

	if(OmegaK>0){//2
		*D_M = (*D_H)*pow(sqrt(OmegaK),-1)*sinh(sqrt(OmegaK)*(*D_C)/(*D_H));
	}//2
	else if (OmegaK==0){//2
		*D_M=*D_C;
	}//2	
	else{//2
		*D_M = (*D_H)*pow(sqrt(OmegaK),-1)*sin(sqrt(OmegaK)*(*D_C)/(*D_H));
	}//2
	return 1.0;
}//1

//angular diameter distance/
double Cosmology::AngularDist(double* D_A, double* D_M, double* D_C, double* D_H,double* z, double radial_dist, double cosmo_parameters[4]) {//1

	double H0 =cosmo_parameters[0];
	double OmegaM = cosmo_parameters[1] ;
	double OmegaL = cosmo_parameters[2] ;
	double OmegaK = cosmo_parameters[3];
    
    RedShift(z,D_H,radial_dist);
    ComovDist(D_C,D_H, z,radial_dist, cosmo_parameters);
    TransverseDist(D_M, D_C,D_H,z,radial_dist, cosmo_parameters);
	
    *D_A = (*D_M)*pow(1.0+(*z),-1);
    
    return 1.0;
}//1

//luminosity distance/
double Cosmology::LuminosityDist( double* D_L,  double* D_M, double* D_C, double* D_H,double* z, double radial_dist, double cosmo_parameters[4]) {//1
    
    RedShift(z,D_H,radial_dist);
    ComovDist(D_C,D_H, z,radial_dist, cosmo_parameters);
    TransverseDist(D_M, D_C,D_H,z,radial_dist, cosmo_parameters);
    
     *D_L = (*D_M)*(1.0+(*z));
     
    return 1.0;
}//1
