#ifndef __PARTICLE_H
#define __PARTICLE_H

class Particle
{
private:
//double ;

public:

//! Constructor of the class.
Particle ( );

//! Destructor of the class.
~Particle ( );

//methods of the Particle class/
float SettingParticleParameters();
float SettingEnergyLossesParameters();
float Freeing();
};

#endif // ____PARTICLE_H
