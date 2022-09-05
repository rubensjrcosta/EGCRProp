#ifndef __COSMOLOGY_H
#define __COSMOLOGY_H

class Cosmology
{
private:
//double ;

public:

//! Constructor of the class.
Cosmology ( );

//! Destructor of the class.
~Cosmology ( );

//methods of the Cosmology class/
double SettingCosmologicalParameters();
double HubbleTime(double* t_H, double H0);
double HubbleDistance(double* D_H, double H0);
double SettingRedShift();
double RedShift(double* z, double* D_H, double radial_dist);
double CalculatingCosmologicalParameters();
double ComovDist(double* D_C,double* D_H, double* z,double radial_dist, double cosmo_parameteres[4]);
double TransverseDist(double* D_M, double* D_C, double* D_H,  double* z, double radial_dist, double cosmo_parameters[4]);
double AngularDist(double* D_A, double* D_M,  double* D_C,double* D_H,double* z, double radial_dist, double cosmo_parameters[4]);
double LuminosityDist( double* D_L,  double* D_M, double* D_C, double* D_H,double* z, double radial_dist, double cosmo_parameters[4]);

};
#endif // ____COSMOLOGY_H
