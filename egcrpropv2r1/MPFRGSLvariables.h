#ifndef __MPFRGSLVARIABLES_H
#define __MPFRGSLVARIABLES_H

using mpfr::mpreal;  

const int digits = 39; //required precision of computations in MPFR C++ library

//mpfr variables ---------------------------------------------------------------/
mpreal beta;//particle's velocity parameter/
mpreal v0[3],v0_mod,v[3], v0_ref[3];//particle's velocity [m/s]/
mpreal F[3]; //Lorentz force [Kg*(m/s^2)]/
mpreal a_L[3];//Lorentz acceleration [m/s^2]

//gsl variables ----------------------------------------------------------------/
gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);// GSL's Mersenne Twister generator
//gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);// GSL's Taus generator
unsigned long mySeedTENSORS = 123;
unsigned long mySeedDIRETIONS = 123;
			
#endif // __MPFRGSLVARIABLES_H
