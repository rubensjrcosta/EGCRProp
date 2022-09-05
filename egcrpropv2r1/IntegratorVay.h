double UpdateVelocityVay(double u[3], double u0[3],  mpfr::mpreal v0[3], double* Gamma,double B[3], double Delta_t);
double PusherVay(double u[3],double u0[3], mpfr::mpreal v[3], mpfr::mpreal v0[3], double x0[3], double x[3], double* Gamma,double B[3], double q, double m0,double* D, double* Delta_x,double Delta_t);
double GammaVay(double* Gamma, double u_mod);
double BetaVay(mpfr::mpreal* beta, double* E);
double LorentzFactorVay(double* Gamma,mpfr::mpreal* beta); 

double UpdateVelocityVay(double u[3], double u0[3],  mpfr::mpreal v0[3], double* Gamma,double B[3], double Delta_t){Geometry * Geom = new Geometry();Dynamics* Dyns=new Dynamics();

	double  u_prime[3], vec_s[3], s,  u_pxt[3];
	double  Gamma_prime, tau_star,sigma;
	double tau_prime[3], tau[3],taup2,u_pdott;	
	double u_minus[3]; 
	mpfr::mpreal taup_aux[3],v0xtaup[3];

	mpfr::mpreal one = 1.0;	
		
	//calculating the vector tau' and tau^2 
	for (i=0;i<=2;i++){//2
		tau_prime[i] = (q*(B[i])/(m0))*(Delta_t/2); 
	    taup_aux[i]=tau_prime[i]*one;
	}//2
	taup2=pow(Geom->Mod(tau_prime),2);

	//calculating the vector u'	
	Geom->CrossProd(v0xtaup,v0,taup_aux);	 
	for (i=0;i<=2;i++){//2
	u_prime[i] = u0[i] + (v0xtaup[i]).toDouble();
	}//2	

	GammaVay(&Gamma_prime,Geom->Mod(u_prime)); 
	tau_star=Geom->ScalProd(u_prime,tau_prime)/(c);
	sigma = pow(Gamma_prime,2) - taup2;
	
	*Gamma = sqrt((sigma+sqrt(pow(sigma,2)+4*(taup2+pow(tau_star,2))))/2);	
	
	for (i=0;i<=2;i++){//2
	tau[i] = tau_prime[i]/(*Gamma);
	}//2			   
	
	u_pdott=Geom->ScalProd(u_prime,tau);	
	
	s=1/(1+ pow(Geom->Mod(tau),2));
	Geom->CrossProd(u_pxt,u_prime,tau);	
	
	for (i=0;i<=2;i++){//2
		u[i] = s*(u_prime[i] + (u_pdott)*tau[i] + u_pxt[i]);
	}//2
	
	delete Geom;delete Dyns;return 1.0;
}//1   
double PushParticle_Vay(double x[3], double x0[3], mpfr::mpreal v0[3],double Delta_t){//1

	Geometry * Geom = new Geometry();
	
	mpfr::mpreal one = 1.0;

	for (i=0;i<=2;i++){//2
		x[i] = (x0[i]*one + v0[i]*Delta_t).toDouble();		
	}//2
	
	delete Geom;return 1.0;	
}//1

double GammaVay(double* Gamma, double u_mod){//1  
	
	*Gamma = sqrt(1+pow(u_mod/c,2));
	
	return 1.0;
}//1

double BetaVay(mpfr::mpreal* beta, double* E){//1 

    mpfr::mpreal one = 1.0;    
    	
	*beta = sqrt(one - pow((m0*pow(c,2))/(*E),2));

	return 1.0;
}//1

double LorentzFactorVay(double* Gamma, mpfr::mpreal* beta){//1

	 mpfr::mpreal one = 1.0;
	
	*Gamma = (one/sqrt(1-pow(*beta,2))).toDouble();

	return 1.0;
}//1
