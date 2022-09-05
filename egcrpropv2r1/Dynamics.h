#ifndef __DYNAMICS_H
#define __DYNAMICS_H
  
class Dynamics
{
private:
//double ;
    
public:
  
//! Constructor of the class.
Dynamics( );

//! Destructor of the class.
~Dynamics ( );

//methods of the Dynamics class/
double PusherRK4th(double p[3],double p0[3],double x[3], double x0[3],double* Gamma,double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t);
double PusherBoris(double u[3],double u0[3],double x[3], double x0[3],double* Gamma,double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t);
double PusherVay(double u[3],double u0[3],mpfr::mpreal v[3],mpfr::mpreal v0[3],double x[3], double x0[3],double* Gamma,double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t);
double PusherNewEuler(mpfr::mpreal v[3],mpfr::mpreal v0[3],double x[3], double x0[3],double* Gamma,mpfr::mpreal beta[3],mpfr::mpreal a_L[3],double* E, double B[3], double q,double m0,double* D, double* Delta_x,double Delta_t);

double Energy();
double SpreadingAngle();
double CyclotronFrequency();
double StepTime();
};    
#endif // __DYNAMICS
