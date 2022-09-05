#ifndef __MATHEMATICS_H
#define __MATHEMATICS_H

class Mathematics
{
private:
//double ;

public:

//! Constructor of the class.
Mathematics ( );

//! Destructor of the class.
~Mathematics ( );

//methods of the Mathematics class/
double RandNum_std(const double rand_min, const double rand_max); 
double RandNum_gsl(const double rand_min, const double rand_max,gsl_rng* r);
double Interpolate(double x_int, double x_inf, double y_inf, double x_sup,  double y_sup);  

//double rk4 ( double t0, double u0, double dt, double f ( double t, double u ) );
//double *rk4vec (double B[3],double q,double t0, int ODE_num, double u0[], double dt, 
//double *f (double B[3],double q,double t, int ODE_num, double u[] ) );
  
};

#endif // ____MATHEMATICS_H
