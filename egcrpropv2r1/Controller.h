#ifndef __CONTROL_H
#define __CONTROL_H

class Controller
{
private:

public:

//! Constructor of the class.
Controller ( );

//! Destructor of the class.
~Controller ( );

//methods of the Controller class/
double LoadingInputParameters(); 
double SettingInitialConditions(); 
double SettingParticlePosition();

double SettingRange();
double SettingAverage();
double CalculatingAverage(); 

double SettingDynamicsParameters();
double CalculatingDynamicsParameters();
double ParticleTracking();

double CuttingTracking();
double StoppingTracking();

};

#endif // __CONTROL

