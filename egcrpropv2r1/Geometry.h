#ifndef __GEOMETRY_H
#define __GEOMETRY_H

using mpfr::mpreal;  
using std::cout;
using std::endl;

class Geometry
{
private:
//mpreal ;

public:

//! Constructor of the class.
Geometry ( );

//! Destructor of the class.
~Geometry ( );

//methods of the Geometry class
mpreal DistPP(mpreal A[3], mpreal B[3]);//the distance between two points/
double DistPP(double A[3], double B[3]);
mpreal CrossProd(mpreal vec3[3], mpreal vec1[3],mpreal vec2[3]);//cross product/
double CrossProd(double vec3[3], double vec1[3],double vec2[3]);
mpreal Norm(mpreal versor[3], mpreal vector[3]);//normalization of a vector/
double Norm(double versor[3], double vector[3]);
mpreal AddVec(mpreal sum[3],mpreal vec1[3], mpreal vec2[3]);//vector addition/
//mpreal MultScal(mpreal prod[3], mpreal vec1[3], mpreal factor);//multiplication of vectors by a scalar/
mpreal ScalMod(mpreal scalar);//modulus of a scalar/
double ScalMod(double scalar);
mpreal Mod(mpreal vector[3]);//modulus of a vector/
double Mod(double vector[3]);
mpreal ScalProd(mpreal vec1[3], mpreal vec2[3]);//scalar product/
double ScalProd(double vec1[3], double vec2[3]);
mpreal Angle2Vec(mpreal vec1[3], mpreal vec2[3]);//the angle between vectors/
double Angle2Vec(double vec1[3], double vec2[3]);


};

#endif // __GEOMETRY
