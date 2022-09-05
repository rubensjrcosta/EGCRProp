//methods for the integration of the equation of motion  with the Runge-Kutta 4th method
double LorentzAcceleration(mpfr::mpreal a_L[3],  mpfr::mpreal v[3],double B[3]);
double UpdateVelocityNewEuler(mpfr::mpreal v[3], mpfr::mpreal v0[3], double* Gamma,  mpfr::mpreal  a_L[3],  double Delta_t);
double PushParticleNewEuler(double x[3],double x0[3], mpfr::mpreal v0[3],double* Gamma, mpfr::mpreal  a_L[3], double Delta_t);
double BetaNewEuler(mpfr::mpreal* beta, double* E);
double LorentzFactorNewEuler(double* Gamma,mpfr::mpreal* beta); 

double LorentzAcceleration(mpfr::mpreal a_L[3], mpfr::mpreal v[3], double B[3]){Geometry * Geom = new Geometry();
		
	mpfr::mpreal B_aux[3];
	mpfr::mpreal vxB[3];
		
	for (i=0;i<=2;i++){B_aux[i] =B[i];}	
	
	Geom -> CrossProd(vxB,v,B_aux);
	
	for (i=0;i<=2;i++){a_L[i]=q*(vxB[i])/((Gamma)*m0);}
	
	delete Geom;return 1.0;
}//1 


double UpdateVelocityNewEuler(mpfr::mpreal v[3], mpfr::mpreal v0[3], double* Gamma,  mpfr::mpreal  a_L[3],  double Delta_t){Geometry * Geom = new Geometry();
	
	mpfr::mpreal v_prime[3],v_norm[3];	
		
	for (i=0;i<=2;i++){//2
	
		v_prime[i]= v0[i] + (a_L[i]*Delta_t);
	
	}//2	
	Geom->Norm(v_norm,v_prime);
	for (i=0;i<=2;i++){//2
	
		v[i]=v_norm[i]*(Geom->Mod(v0));
	
	}//2
	
	delete Geom;return 1.0;	
}//1
double PushParticleNewEuler(double x[3],double x0[3], mpfr::mpreal v0[3],double* Gamma, mpfr::mpreal  a_L[3],  double Delta_t){Geometry * Geom = new Geometry();
	
	for (i=0;i<=2;i++){//2		
		
		x[i] = x0[i] + (v0[i]*Delta_t + (a_L[i]*(pow(Delta_t,2)))/2).toDouble();
	
	}//2	

	delete Geom;return 1.0;
}//1

double BetaNewEuler(mpfr::mpreal* beta, double* E){//1//beta = sqrt(1-(m0*c²/E)²); 

    mpfr::mpreal one = 1.0;    
    	
	*beta = sqrt(one - pow((m0*pow(c,2))/(*E),2));

	return 1.0;
}//1

double LorentzFactorNewEuler(double* Gamma, mpfr::mpreal* beta){//1//Gamma = 1/sqrt(1-beta²)     

	 mpfr::mpreal one = 1.0;
	
	*Gamma = (one/sqrt(1-pow(*beta,2))).toDouble();

	return 1.0;
}//1
