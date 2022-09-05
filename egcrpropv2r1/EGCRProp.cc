//****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCR//
//EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp//

/*
    File: Prog.cc 
	Purpose: Program for ultra-high energy charged particles propagation 
	in the intergalactic medium assumed to have a cellular structure of 
	coherent magnetic fields.			 
	Authors: R. P. Costa Jr and M. A. Leigui de Oliveira
	Started: 19/03/2015		
	Atualized:04/01/2019

*/
//****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCR//
//EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp//

#include "EGCRProp.h"
		
int main (void){//Main//
	    
	Controller* Control=new Controller();		
	
	Control->LoadingInputParameters();
	Control->SettingInitialConditions();  	
	Control->SettingRange();	
	 
	for(N_iter=1;N_iter<=N_evts;N_iter++)
	{//2
		Control->ParticleTracking();		
	}//2	    								
	Control->CalculatingAverage();	
		
	(void)printf("mySeedDIRETIONS = %.ld, mySeedTENSORS = %.ld\n\n", mySeedDIRETIONS, mySeedTENSORS);
			

	delete Control;		
	return 0;
}//Main//
