double *uprimeRK4th(double t, double u[]);
double *IntegrationRK4th(double* Gamma,double u0[], double dt, double *f(double t,  double u[]));
double GammaRK4th(double* Gamma, double p_mod, double m0);

int ODE_num = 6;
double t0 =0.0;
	
double *uprimeRK4th(double t,double u[]){Geometry* Geom = new Geometry();
	double *uprime,p[3],pxB[3];
	double Energy;
	
	//set of six ordinary differential equations	
	uprime = new double[ODE_num];
	
	for(i=3;i<=5;i++){//2
		p[i-3]=u[i];	
	}//2	
		
	for(i=0;i<=2;i++){//2
		uprime[i] =  (u[i+3]*c)/Geom->Mod(p);
	}//2
	
	Geom->CrossProd(pxB,p,B);
	for(i=3;i<=5;i++){//2
		uprime[i] = (q*c)*(pxB[i-3])/Geom->Mod(p);
	}//2

	delete Geom;			
	return uprime;
}

double *IntegrationRK4th(double* Gamma,double u0[], double dt, double *f (double t, double u[] ) ){//1
	double *f0, *f1,*f2,*f3;
	double t1,t2,t3;
	double *u,*u1,*u2,*u3;
	
	//Geting four sample values of the derivative.
	f0 = f (t0, u0 );
	
	t1 = t0 + dt / 2.0;
	u1 = new double[ODE_num];
	for ( i = 0; i < ODE_num; i++ ){//2
		u1[i] = u0[i] + dt * f0[i] / 2.0;
	}//2
	f1 = f (t1,u1 );
	
	t2 = t0 + dt / 2.0;
	u2 = new double[ODE_num];
	for ( i = 0; i < ODE_num; i++ ){//2
		u2[i] = u0[i] + dt * f1[i] / 2.0;
	}//2
	f2 = f (t2,u2);
	
	t3 = t0 + dt;
	u3 = new double[ODE_num];
	for ( i = 0; i < ODE_num; i++ ){//2
		u3[i] = u0[i] + dt * f2[i];
	}//2
	f3 = f (t3,u3 );
	
	//Combining  to estimate the solution.
	u = new double[ODE_num];
	for ( i = 0; i < ODE_num; i++ ){//2
		u[i] = u0[i] + dt * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]) / 6.0;
	}//2	
	//  Free memory.
	delete [] f0;
	delete [] f1;
	delete [] f2;
	delete [] f3;
	delete [] u1;
	delete [] u2;
	delete [] u3;
	
	return u;
}//1
  
double GammaRK4th(double* Gamma, double p_mod, double m0){    	 	    

	*Gamma = p_mod/(m0*(c));
	
	return 1.0;
}//1
