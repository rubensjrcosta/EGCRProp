//InputParameters.cc/

//local headers/
#include "InputParameters.h"

FILE  *parameters_in;

//constructor of the class/
InputParameters::InputParameters ( ){

}

//destructor of the class/
InputParameters::~InputParameters ( ){

}

//methods of the Input Parameters class
float InputParameters::ReadingInputParameters(){//1 
 	
	double parameters;	
	parameters_in=fopen("input_parameters.dat","r");
	if(parameters_in==NULL){(void)printf("\n Can't open input_parameters.dat\n\n");exit(EXIT_FAILURE);return(0);}
	parameters=fgetc(parameters_in);
	(void)ungetc(parameters,parameters_in);
	
	for (i=1;i<=10;i++){//2
		(void)fgets(line,sizeof(line),parameters_in);			
		if(i==1){sscanf(line,"%d",&P_spec);}
		if(i==2){sscanf(line,"%lf %lf",&E_EeV, &D_Mpc);}
		if(i==3){sscanf(line,"%d",&E_loss);}
		if(i==4){sscanf(line,"%d",&MF_model);}			
		if(i==5){//3	
			if(MF_model==1 || MF_model==2){//4	
				sscanf(line,"%lf %lf",&B_nG,&l_coh );
			}//4			
			else{//4	
				sscanf(line,"%lf",&B_nG);				
			}//4
		}//3
		if(i==6){sscanf(line,"%lf %d",&s_depth,&n_cells);}
		if(i==7){sscanf(line,"%d %d %d", &b_track,&d_stop,&m_rigidity);}
		if(i==8){sscanf(line,"%d %d",&N_evts, &N_traj);}		
		if(i==9){sscanf(line,"%d %lf",&i_method,&t_scal);}	
		if(i==10){sscanf(line,"%d %d",&r_dist, &r_seed);}	
	}//2
	parameters=fgetc(parameters_in);		
	(void)fclose(parameters_in);
	CheckingInputParameters();
	
	return 1.0;
}//1

float InputParameters::CheckingInputParameters(){//1  
 	
 	if(P_spec != 1 && P_spec != 2 && P_spec != 3 && P_spec != 4 && P_spec != 5){(void)printf("** Wrong choice for requested particle specie: options available (1 proton;2 helium;3 oxygen;4 silicon;5 iron) **\n");exit(EXIT_FAILURE);}
	
	if(E_EeV < E_min || E_EeV > E_max){(void)printf("** Wrong choice for requested particle's initial energy: E_EeV needs to be 0.1 EeV <= E_EeV <= 30.0 ZeE **\n");exit(EXIT_FAILURE);}
	
	if(D_Mpc<=0.0){(void)printf("** Wrong choice for requested maximum distance: D_max needs to be > 0.0 **\n");exit(EXIT_FAILURE);}
	
		if(E_loss==0){//2
		strcpy (ENERGYLOSS, onoff_Eloss[0]);
		strcpy (EINDEX,  E_index[2]);				 		 			 				
	}//2		
	else if(E_loss==1){//2		
		strcpy (ENERGYLOSS, onoff_Eloss[1]);			 		 		
	}//2		
	else{(void)printf("** Wrong choice for requested energy losses: options available (0 off; 1 on) **\n");exit(EXIT_FAILURE);}
	
	if(MF_model== 1){//2
		strcpy (MFMODEL, mf_model[0]);			
	}//2
	else if(MF_model == 2){//2
		s_depth=0,n_cells=0;
		strcpy (MFMODEL, mf_model[1]);	
	}//2
	else if(MF_model == 3){//2
		s_depth=0,n_cells=0;
		strcpy (MFMODEL, mf_model[2]);	
	}//2
	//else if(MF_model==XXXXXX){//2
			//strcpy (MFMODEL, mf_model[XXXXXX]);	
	//}//2
	else{(void)printf("** Wrong choice for requested magnetic field model: options available (1 Cellular Structure,2 Cubic Domain,3 Uniform Magnetic Field) **\n");exit(EXIT_FAILURE);}

	if(B_nG<=0.0){(void)printf("** Wrong choice for requested magnetic field strength: B_nG needs to be > 0.0 **\n");exit(EXIT_FAILURE);}
	
	if(MF_model==1 || MF_model==2){//2
		if(l_coh<=0){(void)printf("** Wrong choice for requested coherence length: l_coh needs to be > 0.0 **\n");exit(EXIT_FAILURE);}	
	}//2
	if(MF_model==1){//2
		
		if(s_depth<=0){(void)printf("** Wrong choice for requested skin depth: s_depth needs to be > 0%% **\n");exit(EXIT_FAILURE);}		
		
		if(n_cells<2 || n_cells>27){(void)printf("** Wrong choice for requested nearby cells to be considered: n_cells needs to be 2 <= n_cells <= 27 **\n");exit(EXIT_FAILURE);}	
	}//2	
	
	if(b_track == 0){//2
		if(strcmp(ENERGYLOSS, "ON") == 0)strcpy (EINDEX,  E_index[1]);			 	
		strcpy (BTRACK, onoff_btrack[0]);			 		 			
		//the redshift at the initial position is known only with backtracking ON (where z0=0.0) or with the stopping distance  2 (the trajectory radius)
		if(d_stop == 1){(void)printf("** Wrong choice: the  backtracking OFF works only with d_stop = 2 **\n");exit(EXIT_FAILURE);}			
	}//2
	else if(b_track == 1){//2
		if(strcmp(ENERGYLOSS, "ON") == 0)strcpy (EINDEX,  E_index[0]);	
		strcpy (BTRACK, onoff_btrack[1]);			 	
	}//2
	else{(void)printf("** Wrong choice for requested backtracking: options available (0 off;1 on) **\n");exit(EXIT_FAILURE);}	
	
	if(d_stop == 1){//2
			strcpy (DSTOP, D_stop[0]);			 					
		}//2        
	else if(d_stop == 2){//2	
		strcpy (DSTOP, D_stop[1]);			 					 					
	}//2
	else{(void)printf("** Wrong choice for requested stopping distance: options available (1 traveled distance;2 trajectory radius) **\n");exit(EXIT_FAILURE);}	
	
	if(m_rigidity == 0){//2
		strcpy (MRINDEX, MR_index[0]);			 	
		strcpy (EMRINDEX, EMR_index[0]);			 	
		strcpy (MRIGIDITY, onoff_Mr[0]);			 			
		}//2
	else if(m_rigidity == 1){//2
		strcpy (MRINDEX, MR_index[1]);			 	
		strcpy (EMRINDEX, EMR_index[1]);			 	
		strcpy (MRIGIDITY, onoff_Mr[1]);			 	
	}//2
	else{(void)printf("** Wrong choice for requested magnetic rigidity: options available (0 off;1 on) **\n");exit(EXIT_FAILURE);}

	if(N_evts <= 0){(void)printf("** Wrong choice for requested number of particles to be launched: N_evts needs to be > 0 **\n");exit(EXIT_FAILURE);}
	
	if(N_traj < 0 || N_traj > N_evts ){(void)printf("** Wrong choice for requested number of plotted trajectories: N_traj needs to be >= 0 and <= N_evts **\n");exit(EXIT_FAILURE);}
						
		 if(i_method== 1){//2
			strcpy (INTMETHOD, int_method[0]);			 	
		}//2
	else if(i_method == 2){//2
			strcpy (INTMETHOD, int_method[1]);	
	}//2
	else if(i_method == 3){//2
			strcpy (INTMETHOD, int_method[2]);	
	}//2
	else if(i_method == 4){//2
			strcpy (INTMETHOD, int_method[3]);	
	}//2
	else{(void)printf("** Wrong choice for requested equation of motion integration: options available (1 Runge-Kutta 4th,2 New Euler,3 Boris,4 New LeapFrog) **\n");exit(EXIT_FAILURE);}	
	
	if(t_scal<=0.0){(void)printf("** Wrong choice for requested time scale: t_scal needs to be > 0.0 **\n");exit(EXIT_FAILURE);}
		
	if(r_dist==1)strcpy(RANDDISTRIBUTION, rand_dist[0]);		
	else if(r_dist==2)(RANDDISTRIBUTION, rand_dist[1]);	
	else {(void)printf("** Wrong choice for requested rand_dist: options available (1 uniform;2 gaussian, with the parameter sigma=1) **\n");exit(EXIT_FAILURE);}	
	
	if(r_seed==0){strcpy(SEED, "BASEDONTIME");}
	else if(r_seed!=1){(void)printf("** Wrong choice for requested rand_seed: options available (0 based on time, 1 fixed (123 as default - MPFRGSLvariables.h line 16)) **\n");exit(EXIT_FAILURE);}

	return 1.0;
}//1








