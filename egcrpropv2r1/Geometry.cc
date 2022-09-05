//Geometry.cc/
#include "Geometry.h"

//constructor of the class/
Geometry::Geometry ( ){

}

//destructor of the class/
Geometry::~Geometry ( ){

}

//the distance between two points/
mpreal Geometry::DistPP(mpreal A[3], mpreal B[3]){//1
	mpreal delta;
	mpreal dist;
    int i;
    
    delta=0.0;
    for(i=0;i<=2;i++){//2
		delta=delta+((A[i]-B[i])*(A[i]-B[i]));
	}//2
     
     dist=sqrt(delta);
     return(dist);
     
}//1

//the distance between two points double/
double Geometry::DistPP(double A[3], double B[3]){//1
	double delta;
	double dist;
    int i;
    
    delta=0.0;
    for(i=0;i<=2;i++){//2
		delta=delta+((A[i]-B[i])*(A[i]-B[i]));
	}//2
     
     dist=sqrt(delta);
     return(dist);
     
}//1
                 
//cross product/
mpreal Geometry::CrossProd(mpreal vec3[3], mpreal vec1[3],mpreal vec2[3]){//1
	  
	vec3[0]=(vec1[1]*vec2[2]-vec1[2]*vec2[1]);
	vec3[1]=(vec1[2]*vec2[0]-vec1[0]*vec2[2]);
	vec3[2]=(vec1[0]*vec2[1]-vec1[1]*vec2[0]);

return 1.0;
}//1

//cross product double/
double Geometry::CrossProd(double vec3[3], double vec1[3], double  vec2[3]){//1
	  
	vec3[0]=(vec1[1]*vec2[2]-vec1[2]*vec2[1]);
	vec3[1]=(vec1[2]*vec2[0]-vec1[0]*vec2[2]);
	vec3[2]=(vec1[0]*vec2[1]-vec1[1]*vec2[0]);

return 1.0;
}//1

//normalization of a vector/
mpreal Geometry::Norm(mpreal versor[3], mpreal vector[3]){//1
	
	mpreal modulus;
	int i;
	
	modulus=Mod(vector);
	
	if (modulus!=0.){//2
		
		for(i=0;i<=2;i++){//3			
			versor[i]=vector[i]/modulus;
		}//3
     }//2
     return(1);
}//1

//normalization of a vector double/
double Geometry::Norm(double versor[3], double vector[3]){//1
	
	double modulus;
	int i;
	
	modulus=Mod(vector);
	
	if (modulus!=0.){//2
		
		for(i=0;i<=2;i++){//3			
			versor[i]=vector[i]/modulus;
		}//3
     }//2
     return(1);
 }//1
      
//vector addition/
mpreal Geometry::AddVec(mpreal sum[3],mpreal vec1[3], mpreal vec2[3]){//1
	int i;
    for(i=0;i<=2;i++){//2
		sum[i]=vec1[i]+vec2[i];
	}//2
	  return(1);
}//1
    
////multiplication of vectors by a scalar/
//mpreal Geometry::MultScal(mpreal prod[3], mpreal vec1[3], mpreal factor){//1
	//int i;
    //for(i=0;i<=2;i++){//2
		//prod[i]=vec1[i]*factor;
	//}//2
//}//1
    
//modulus of a scalar/
mpreal Geometry::ScalMod(mpreal scalar){//1		
	mpreal scal_mod=sqrt(scalar*scalar);
	return(scal_mod);
}//1

//modulus of a scalar/
double Geometry::ScalMod(double scalar){//1		
	double scal_mod=sqrt(scalar*scalar);
	return(scal_mod);
}//1

//modulus of a vector/
mpreal Geometry::Mod(mpreal vector[3]){//1		
	 mpreal modulus=sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);    
     return(modulus);
}//1

//modulus of a vector double/
double Geometry::Mod(double vector[3]){//1		
	 double modulus=sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);    
     return(modulus);
}//1

//scalar product/
mpreal Geometry::ScalProd(mpreal vec1[3], mpreal vec2[3]){//1
	int i;
	mpreal scal_prod;		
	scal_prod=0;
	
	for(i=0;i<=2;i++){//2
		scal_prod+=vec1[i]*vec2[i];
	}//2
    return(scal_prod);
}//1

//scalar product double/
double Geometry::ScalProd(double vec1[3], double vec2[3]){//1
	int i;
	double scal_prod;		
	scal_prod=0;
	
	for(i=0;i<=2;i++){//2
		scal_prod+=vec1[i]*vec2[i];
	}//2
    return(scal_prod);
}//1

//the angle between two vectors/
mpreal Geometry::Angle2Vec(mpreal vec1[3], mpreal vec2[3]){//1	

	
	mpreal scalar,mod1,mod2,theta;
    scalar=ScalProd(vec1,vec2);
    mod1=Mod(vec1);
    mod2=Mod(vec2);
     
    theta=(acos(scalar/(mod1*mod2)))*(degbyrad);
    
    return(theta);
}//1

//the angle between two vectors double/
double Geometry::Angle2Vec(double vec1[3], double vec2[3]){//1	
	
	double scalar,mod1,mod2,theta;
    scalar=ScalProd(vec1,vec2);
    mod1=Mod(vec1);
    mod2=Mod(vec2);
     
    theta=(acos(scalar/(mod1*mod2)))*(degbyrad);
 
    return(theta);
}//1
