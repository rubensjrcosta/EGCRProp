#ifndef __MAGNETIC_H
#define __MAGNETIC_H

class Magnetic
{
private:

public:

//! Constructor of the class.
Magnetic ( );

//! Destructor of the class.
~Magnetic ( );

//methods of the Magnetic class/
float SettingMagneticField();
float LoadingMagneticTensors();
float GridPositioning();
float MagneticField();
float RegularComponent();
float SphericalStructure();
float CubicStructure();

};

#endif // __MAGNETIC
