//Mathematics.cc/

//local headers/
#include "Mathematics.h"
#include  "Geometry.h"

//constructor of the class Mathematics/
Mathematics::Mathematics ( ){

}

//destructor of the class Mathematics/
Mathematics::~Mathematics ( ){

}

//methods of the Mathematics class/
double Mathematics::RandNum_std(const double rand_min, const double rand_max){//1 //random number generator/

	Geometry * Geom = new Geometry();
	double rand_num; 
	double delta;//delta values to rand number
	
	delta = rand_max - rand_min;
	rand_num = (delta*(Geom->ScalMod(sin(rand()))) + rand_min);
	
	delete Geom;
	return (rand_num);
}//1

double Mathematics::RandNum_gsl(const double rand_min, const double rand_max,gsl_rng* r){//1 //random number generator/

	Geometry * Geom = new Geometry();
	double rand_num; 
	double delta;//delta values to rand number
  	double sigma=1.0;
	delta = rand_max - rand_min;
	double u;
	
	if(r_dist==1)  {u= gsl_rng_uniform (r);}	
	else{u= gsl_ran_gaussian(r,sigma);}	

	rand_num = (delta*(Geom->ScalMod(sin(u))) + rand_min);
	
	delete Geom;
	return (rand_num);
}//1


//numerical interpolation/ 
double Mathematics::Interpolate(double x_int, double x_inf, double y_inf, double x_sup,  double y_sup){//1

	double wi= sqrt(pow(1.0/(x_int-x_inf),2));
	double ws=sqrt(pow(1.0/(x_sup - x_int),2));	
	double y_int=(y_inf*wi+y_sup*ws)/(wi+ws);
	
	return y_int;
}//1




////****************************************************************************80
//double Mathematics::rk4 ( double t0, double u0, double dt, double f ( double t, double u ) ){//1
	//double f0;
	//double f1;
	//double f2;
	//double f3;
	//double t1;
	//double t2;
	//double t3;
	//double u;
	//double u1;
	//double u2;
	//double u3;
	////
	////  Get four sample values of the derivative.
	////
	//f0 = f ( t0, u0 );
	
	//t1 = t0 + dt / 2.0;
	//u1 = u0 + dt * f0 / 2.0;
	//f1 = f ( t1, u1 );
	
	//t2 = t0 + dt / 2.0;
	//u2 = u0 + dt * f1 / 2.0;
	//f2 = f ( t2, u2 );
	
	//t3 = t0 + dt;
	//u3 = u0 + dt * f2;
	//f3 = f ( t3, u3 );
	////
	////  Combine to estimate the solution at time T0 + DT.
	////
	//u = u0 + dt * ( f0 + 2.0 * f1 + 2.0 * f2 + f3 ) / 6.0;

  //return u;
//}//1

//double* Mathematics::rk4vec (double B[3],double q,double t0, int ODE_num, double u0[], double dt, 
//double *f (double B[3],double q, double t, int ODE_num, double u[] ) ){//1
  
	//double *f0;
	//double *f1;
	//double *f2;
	//double *f3;
	//int i;
	//double t1;
	//double t2;
	//double t3;
	//double *u;
	//double *u1;
	//double *u2;
	//double *u3;
	////
	////  Get four sample values of the derivative.
	////
	//f0 = f (B,q,t0, ODE_num, u0 );
	
	//t1 = t0 + dt / 2.0;
	//u1 = new double[ODE_num];
	//for ( i = 0; i < ODE_num; i++ ){//2
		//u1[i] = u0[i] + dt * f0[i] / 2.0;
	//}//2
	//f1 = f (B,q, t1, ODE_num, u1 );
	
	//t2 = t0 + dt / 2.0;
	//u2 = new double[ODE_num];
	//for ( i = 0; i < ODE_num; i++ ){//2
		//u2[i] = u0[i] + dt * f1[i] / 2.0;
	//}//2
	//f2 = f (B,q, t2, ODE_num, u2 );
	
	//t3 = t0 + dt;
	//u3 = new double[ODE_num];
	//for ( i = 0; i < ODE_num; i++ ){//2
		//u3[i] = u0[i] + dt * f2[i];
	//}//2
	//f3 = f (B,q, t3, ODE_num, u3 );
	////
	////  Combine them to estimate the solution.
	////
	//u = new double[ODE_num];
	//for ( i = 0; i < ODE_num; i++ ){//2
		//u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
	//}//2
	////
	
	////Free memory.
	////
	//delete [] f0;
	//delete [] f1;
	//delete [] f2;
	//delete [] f3;
	//delete [] u1;
	//delete [] u2;
	//delete [] u3;
	
	//return u;
//}//1
////****************************************************************************80
