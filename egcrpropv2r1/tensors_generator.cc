/*
    File: tensors_generator.cc 
	Purpose: Program for the generate the tensors of the Extragalactic Magnetic Field  					 
	Author: R. P. Costa Jr.
	Started: 19/03/2015		
	Atualized:26/07/2017

*/

//local headers -------------------------------------------------------/
#include "includes.h"
#include "defines.h"
#include "globalvariables.h"
#include "MPFRGSLvariables.h"

//classes -------------------------------------------------------------/
#include "Mathematics.cc"
#include "Geometry.cc"

FILE *magnetic_tensors_out;	

int main (void){ 

	//the standard Unix random number generator 
	srand( (unsigned)time(NULL) ); //setting the 'random seed' based on time	
	
	//the GSL random number generator
    //(r_seed==0)gsl_rng_set(rng, time(NULL)); //setting the 'random seed' based on time
	gsl_rng_set(rng, mySeedTENSORS); //setting the random seed: 123 as default (mySeedTENSORS in MPFRGSLvariables.h line 17)
	
	Mathematics * Math = new Mathematics();
	
	sprintf(string,"magnetic_tensors_SEED%ld.dat",mySeedTENSORS);
	magnetic_tensors_out=fopen(string,"w");
	if (magnetic_tensors_out==NULL){(void)printf("Can't open %s\n",string);return(0);}
	
	for (i=0;i<=Unv_scl;i++){//1
		for (j=0;j<=Unv_scl;j++){//2
			for (k=0;k<=Unv_scl;k++){//3
				theta=Math->RandNum_gsl(0.0,(pi),rng);
				phi=Math->RandNum_gsl(0.0,2*(pi),rng);  
				(void)fprintf(magnetic_tensors_out,"%lf %lf\n",theta,phi); 
			}//3
		}//2			
	}//1	
	(void)fclose(magnetic_tensors_out);

	sprintf(archive,"TENSORSGENERATOR-SUCCESSFULLYCOMPLETED_SEED%ld.dat",mySeedTENSORS);
	successfully_completed=fopen(archive,"a");
	if (successfully_completed==NULL){(void)printf("\nCan't open %s\n",archive);return(0);}
	
	(void)fprintf(successfully_completed,"******************************************************************************************************************************\n");


return (0);

}


 
 
		


