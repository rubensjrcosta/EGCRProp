//Outputfiles.cc/

//local headers/
#include "Outputfiles.h"

//constructor of the class/
Outputfiles::Outputfiles ( ){

}

//destructor of the class/
Outputfiles::~Outputfiles ( ){

}

////methods of the Outputfiles class/
float Outputfiles::OpeningSummary(){//1
	
	sprintf(archive,"EGCRProp-summary-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_Nevts%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,l_coh/(mbyMpc),N_evts);
	summary=fopen(archive,"w");
	if(summary==NULL){(void)printf("\nCan't open %s\n",archive);return(0);}
		
	return 1.0;
}//1

float Outputfiles::PrintingSummaryBegin(){//1 
	
	(void)fprintf(summary,"**************************************************************************************************************************************\n");
	(void)fprintf(summary,"************************************************* DATE : %s  TIME: %s *************************************************\n",__DATE__,__TIME__);
	(void)fprintf(summary,"**************************************************************************************************************************************\n");
	(void)fprintf(summary,"***************** INPUT PARAMETERS - EGCRPROP RUN ********************************** INPUT PARAMETERS - EGCRPROP RUN *****************\n\n");
	(void)fprintf(summary,"%s (Z=%d,A=%d, rest mass [MeV/c^2] = %.3lf, total charge [C] = %.3e) \n",PARTICLENAME,Z,A,m0/(MeVc2bykg),q);		
	(void)fprintf(summary,"%s [EeV] = %.2e, maximum distance [Mpc] = %.2e\n",EMRINDEX,(E_EeV),D_Mpc);	
	(void)fprintf(summary,"energy loss is %s\n",ENERGYLOSS);		
	(void)fprintf(summary,"Magnetic Field Model: %s\n",MFMODEL);
	(void)fprintf(summary,"magnetic field strength [nG] = %.1lf, coherence length [Mpc] = %.1lf\n",B_r/(TbynG),l_coh/(mbyMpc));	
	if(strcmp(MFMODEL,"Spherical Structure") == 0) (void)fprintf(summary,"skin depth [%%] = %.1lf, number of nearby cells to be considered = %d\n",s_depth*100,n_cells);
	(void)fprintf(summary,"backtracking is %s, stopping distance is the %s, magnetic rigidity is %s\n",BTRACK,DSTOP,MRIGIDITY);											
	(void)fprintf(summary,"number of launched particles %d, number of plotted trajectories %d \n", N_evts,N_traj);		
	if(strcmp(INTMETHOD,"NewEuler") == 0 || strcmp(INTMETHOD,"NewLeapFrog") == 0)(void)fprintf(summary,"equation of motion integration: %s method (MPFR digits %d), time scale = %.2e\n",INTMETHOD,digits,t_scal);
	else(void)fprintf(summary,"equation of motion integration: %s method, time scale = %.2e\n",INTMETHOD,t_scal);
	(void)fprintf(summary,"lauching directions: %s distribution, ",RANDDISTRIBUTION);	
	if(r_seed==0)(void)fprintf(summary,"seed = %s\n\n", SEED);
    else if (r_seed==1)(void)fprintf(summary,"seed = %ld\n", mySeedDIRETIONS);
    (void)fprintf(summary,"magnetic tensors generator: uniform distribution, seed = %ld\n", mySeedTENSORS);
	(void)fprintf(summary,"random number generator:  GSL's %s\n\n",RANDGENERATOR);	
    
	(void)fprintf(summary,"**************************************************************************************************************************************\n");
	(void)fprintf(summary,"**************************************************************************************************************************************\n");
				
	PrintingDynamicsQuantities();	
	
	return 1.0;
}//1

float Outputfiles::PrintingDynamicsQuantities(){//1
	
	(void)fprintf(summary,"DYNAMICS PARAMETERS:\n\n");
	//(void)fprintf("Gamma = %.3e, omega_B [s⁻¹]: %.3e, delta_t [s]: %.3e, \nomega_B*delta_t: %.3e, r_L [Mpc]: %.3e\n\n", Gamma,omega_B,delta_t,omega_B*delta_t,r_L/(mbyMpc));
	(void)fprintf(summary,"D_max (%s) [Mpc] = %.2e\n\n", DSTOP,D_Mpc);		
	(void)fprintf(summary,"E%s%s [EeV] = %.2e (z = %.3e)\n",EINDEX,MRINDEX, E0/((JbyEeV)*Z_aux),z0);	

	return 1.0;
}//1

float Outputfiles::OpeningDynamicsParameters(){//1
	
	if(strcmp(DSTOP,"traveled distance")== 0) sprintf(archive,"EGCRProp-D-defl-Enrg-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,(l_coh/(mbyMpc)),E0/((JbyEeV)*Z_aux),D_Mpc,N_evts);
	else  sprintf(archive,"EGCRProp-D0-defl-Enrg-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,(l_coh/(mbyMpc)),E0/((JbyEeV)*Z_aux),D_Mpc,N_evts);
	EGCRProp_out=fopen(archive,"a");
	if (EGCRProp_out==NULL){(void)printf("\nCan't open %s\n",archive);return(0);}	

	return 1.0;
}//1

float Outputfiles::PrintingDynamicsParameters(){//1
	
	
	if(strcmp(DSTOP,"traveled distance")== 0)(void)fprintf(EGCRProp_out,"%.15lf %.15lf %.15lf\n", D/(mbyMpc),defl,(E)/((JbyEeV)*Z_aux));
	else (void)fprintf(EGCRProp_out,"%.15lf %.15lf %.15lf\n", D_C,defl,(E)/((JbyEeV)*Z_aux));

	return 1.0;
}//1

float Outputfiles::OpeningTrajectory(){//1
	
	sprintf(archive,"EGCRProp-traj-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d_%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,(l_coh/(mbyMpc)),E0/((JbyEeV)*Z_aux),D_Mpc,N_evts,N_iter);
	trajectory_out=fopen(archive,"w");
	if (trajectory_out==NULL){(void)printf("\nCan't open %s\n",archive);return(0);}
							
	return 1.0;
}//1

float Outputfiles::PrintingTrajectory(){//1
	
	if(N_iter<=N_traj){//1								
		(void)fprintf(trajectory_out,"%.10lf %.10lf %.10lf\n",(x0[0]-O[0])/(mbyMpc),(x0[1]-O[1])/(mbyMpc),(x0[2]-O[2])/(mbyMpc));					
	}//1	
					
	return 1.0;
}//1

float Outputfiles::PrintingSummaryBody(){//1
	
	if(strcmp(BTRACK,"OFF") == 0){strcpy (EINDEX, E_index[0]);}	
	else {strcpy (EINDEX, E_index[1]);}
	
	if(strcmp(ENERGYLOSS, "OFF") == 0){strcpy (EINDEX, E_index[2]);}

	(void)fprintf(summary,"E%s%s [EeV] = %.2e (z = %.3e)\n\n", EINDEX,MRINDEX,(E_avg/((JbyEeV)*Z_aux)),z);
	(void)fprintf(summary,"Dampl_avg (D0 = %.2e Mpc) [%%] = %.3e\n\n",DC_avg,Dampl);	
	(void)fprintf(summary,"defl_avg [degree] = %.3e\n\n", (defl_avg));

	
	//if(zcut_counts != 0){//1
		//(void)fprintf(summary,"cut in the redshift: %d events of %d\n",zcut_counts,N_evts);
		//(void)printf("cut in the redshift: %d events of %d\n",zcut_counts,N_evts);
	//}//1
	
	//if(Ecut_counts1 != 0){//1
		//(void)fprintf(summary,"cut in the minimum energy: %d events of %d\n",Ecut_counts1,N_evts);
		//(void)printf("cut in the minimum energy: %d events of %d\n",Ecut_counts1,N_evts);
	//}//1

	//if(Ecut_counts2 != 0){//1
		//(void)fprintf(summary,"cut in the maximum energy: %d events of %d\n",Ecut_counts2,N_evts);
		//(void)printf("cut in the maximum energy: %d events of %d\n",Ecut_counts2,N_evts);
	//}//1
	
	if(strcmp(BTRACK,"OFF") == 0){strcpy (EINDEX,  E_index[1]);}
	else {strcpy (EINDEX,  E_index[0]);}
	
	if(strcmp(ENERGYLOSS, "OFF") == 0){strcpy (EINDEX,  E_index[2]);}
		
	return 1.0;
}//1

float Outputfiles::PrintingSummaryEnd(){//1
	
	Particle* Part=new Particle();
	PrintingCosmologicalParameters();
	
	(void)fprintf(summary,"\ntotal CPU time elapsed [s]: %.3e\n",((double) (clock()- time_start)) / CLOCKS_PER_SEC);	
	(void)fprintf(summary,"**************************************************************************************************************************************\n");
	(void)fprintf(summary,"************************************************* DATE : %s  TIME: %s *************************************************\n",__DATE__,__TIME__);
	(void)fprintf(summary,"**************************************************************************************************************************************\n\n");	
	
	SuccessfullyCompleted();													
	ClosingSummary();
	//Part->Freeing();
	
	delete Part;	
	return 1.0;
}//1

float Outputfiles::PrintingCosmologicalParameters(){//1
	
	(void)fprintf(summary,"\n**************************************************************************************************************************************\n");
	(void)fprintf(summary,"**************************************************************************************************************************************\n");	
	(void)fprintf(summary,"COSMOLOGICAL PARAMETERS:\n\n");	
	(void)fprintf(summary,"H0 (Hubble constant) [(Km/s)/Mpc]: %.1lf, Hubble time: [s] = %.2e, Hubble distance: [Mpc] = %.2e\n", (cosmo_parameters[0])/1000, t_H,D_H);	
	(void)fprintf(summary,"total matter density: %.1lf, dark energy density: %.1lf, space curvature: %.1lf\n",cosmo_parameters[1],cosmo_parameters[2],cosmo_parameters[3]); 	
    
	return 1.0;
}//1

float Outputfiles::SuccessfullyCompleted(){//1
	
	sprintf(archive,"EGCRProp-SUCCESSFULLYCOMPLETED-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_Nevts%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,(l_coh/(mbyMpc)),N_evts);
	successfully_completed=fopen(archive,"a");
	if (successfully_completed==NULL){(void)printf("\nCan't open %s\n",archive);return(0);}
	
	(void)fprintf(successfully_completed,"******************************************************************************************************************************\n");
	(void)fprintf(successfully_completed,"********************************************* DATE : %s  TIME: %s *********************************************\n",__DATE__,__TIME__);
	(void)fprintf(successfully_completed,"******************************************************************************************************************************\n");
	(void)fprintf(successfully_completed,"*********** SUCCESSFULLY COMPLETED ************** SUCCESSFULLY COMPLETED ************** SUCCESSFULLY COMPLETED ***************\n");
	(void)fprintf(successfully_completed,"******************************************************************************************************************************\n\n");

	return 1.0;
}//1

float Outputfiles::ClosingDynamicsParameters(){//1
	
	(void)fclose(EGCRProp_out);
		
	return 1.0;
}//1

float Outputfiles::ClosingSummary(){//1
	
	(void)fclose(summary);
	(void)fclose(successfully_completed);
		
	return 1.0;
}//1

float Outputfiles::ClosingTrajectory(){//1
	
	if(N_traj!=0 && N_iter<=N_traj)(void)fclose(trajectory_out);	
		
	return 1.0;
}//1
