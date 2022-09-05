double UpdateVelocityBoris(double u[3], double u0[3], double* Gamma, double Delta_t);
double PushParticleBoris(double x[3], double x0[3], double u[3], double* Gamma,double Delta_t);
double GammaBoris(double* Gamma, double u_mod);

double UpdateVelocityBoris(double u[3], double u0[3], double* Gamma, double Delta_t){Geometry * Geom = new Geometry();

	double  u_prime[3], vec_s[3], u0xtau[3],  u_pxtau[3];
	double tau[3],tau2;
		
	for (i=0;i<=2;i++){//2
		tau[i] = (Delta_t/2)*(q*(B[i])/(m0*(*Gamma))); 
	}//2	
	
	tau2=pow(Geom->Mod(tau),2);
		    		
	Geom->CrossProd(u0xtau,u0,tau);
	
	for (i=0;i<=2;i++){//2
		u_prime[i] = (u0[i] + u0xtau[i]);
	}//2	
	
	Geom->CrossProd(u_pxtau,u_prime,tau);
		
	for (i=0;i<=2;i++){//2
		u[i] = u0[i] + (2/(1 + tau2))*(u_pxtau[i]);			    
	}//2 
	
	delete Geom;return 1.0;	
}//1   

double PushParticleBoris(double x[3], double x0[3], double u[3], double* Gamma,double Delta_t){Geometry * Geom = new Geometry();
	
	for (i=0;i<=2;i++){//2
		x[i] = x0[i] + u[i]*(Delta_t/(*Gamma));		
	}//2
		
	delete Geom;return 1.0;	
}//1

double GammaBoris(double* Gamma, double u_mod){  
	
	*Gamma = sqrt(1+pow(u_mod/c,2));
	
	return 1.0;
}//1

