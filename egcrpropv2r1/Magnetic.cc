//Magnetic.cc/

//local headers/
#include "Magnetic.h"
#include  "Geometry.h"

FILE  *magnetic_tensors;

double Theta, Phi;//auxiliary tensors theta and phi 
double loc[3],loc_aux[3];
int locint_aux[3],locint_swap[3];		
int ijk_aux[3];
int loc_int[3];

const int nearbycell_num = 27;//number of nearby cells
int nearbycell;
double dist_cell[nearbycell_num-1];//distance between particle's position and one of the cell of the 27 nearby cells
double theta_cell[nearbycell_num-1];//theta tensor of one cell of the 27 nearby cells
double phi_cell[nearbycell_num-1];//phi tensor of one cell of the 27 nearby cells											
double dist_inverse_sum;//the sum  of the distance inverse 
double Unv_sclplusOne;//double number = 171
int Unv_sclplusOne_int;//integer number = 171
int jj, min;
double swap;
		
//the constructor of the class/
Magnetic::Magnetic ( ){

}

//the destructor of the class/
Magnetic::~Magnetic ( ){

}

//the methods of the class/
float Magnetic::SettingMagneticField(){//1

	if(	strcmp(MFMODEL,"Spherical Structure") == 0 || 	strcmp(MFMODEL,"Cubic Structure") == 0){//2
		/*number of cells = (Universe Sacale +1)^3*/;
		Unv_sclplusOne=Unv_scl+1;Unv_sclplusOne_int=int(Unv_sclplusOne);
				
		B_r = (TbynG)*B_nG;//transforming the magnectic field strength to Tesla [T]
		l_coh = l_coh*(mbyMpc);//transforming the coherence length to [m]			
		s_depth = s_depth/100;//transforming the skin depth to [%]			
			
		Origin=l_coh*((Unv_sclplusOne_int)*(10*(int(D_Mpc))+1));
	
		for(i=0;i<=2;i++){//2 //initializing the position at the Universe's	center										
			x0[i] = Origin;
			locint_swap[i]=Origin/(l_coh);			
		}//2		
		LoadingMagneticTensors();
		GridPositioning();
		if(	strcmp(MFMODEL,"Spherical Structure") == 0){SphericalStructure();}	
		if(	strcmp(MFMODEL,"Cubic Structure") == 0){CubicStructure();}	
	
	}//2		
//else if(MF_model==XXXXXX){//11 //dipole magnetic field/
	//for (i=0;i<=2;i++){//2 
		//m_mag=????///		
	//}//2
//}//11		
	else{//2		
		B_reg = (TbynG)*B_nG;//transforming the magnectic field strength to Tesla [T]
		Origin=0.0;		
		RegularComponent();
	}//2
		
	return 1.0;
}//1

float Magnetic::LoadingMagneticTensors(){//1 

	double tensors;
	sprintf(string,"magnetic_tensors_SEED%ld.dat",mySeedTENSORS);
	magnetic_tensors=fopen(string,"r");
	if (magnetic_tensors==NULL){(void)printf("\n Can't open %s\n\n", string);exit(EXIT_FAILURE);return(0);}
	tensors=fgetc(magnetic_tensors);
	(void)ungetc(tensors,magnetic_tensors);
	
	//loading the tensors Theta and Phi/
	for (i=0;i<=Unv_scl;i++){//2
		for (j=0;j<=Unv_scl;j++){//3
			for (k=0;k<=Unv_scl;k++){//4
				(void)fgets(line,sizeof(line),magnetic_tensors);
				sscanf(line,"%lf %lf",&Theta_Tensor[i][j][k],&Phi_Tensor[i][j][k]); 	
				//(void)printf("%lf %lf\n",Theta_Tensor[i][j][k],Phi_Tensor[i][j][k]); 
			}//4
		}//3
	}//2
	tensors=fgetc(magnetic_tensors);	
	(void)fclose(magnetic_tensors);
	
	return 1.0;
}//1

float Magnetic::GridPositioning(){//1

	Geometry * Geom = new Geometry();

	//calculation the cell position 	
	for(i=0;i<=2;i++){//2		
		loc[i] = x0[i]/(l_coh);	
	//(void)printf("Mpc %lf, ", x0[i]/(l_coh));																

	}//2
	
	loc_int[0] =  int(loc[0]);
	if(Geom->ScalMod(loc[0] - loc_int[0]) >= 0.5){ loc_int[0] = loc_int[0] + 1;}		
	locint_aux[0] = -loc_int[0] + int(fmod(loc_int[0],Unv_sclplusOne));
			
	loc_int[1] =  int(loc[1]);
	if(Geom->ScalMod(loc[1] - loc_int[1]) >=  0.5){ loc_int[1] = loc_int[1] + 1;}
	locint_aux[1] = -loc_int[1] + int(fmod(loc_int[1],Unv_sclplusOne));

	loc_int[2] =  int(loc[2]);
	if(Geom->ScalMod(loc[2] - loc_int[2]) >= 0.5){ loc_int[2]= loc_int[2] + 1;}	
	locint_aux[2] = -loc_int[2] + int(fmod(loc_int[2],Unv_sclplusOne));

	//(void)printf("\nloc_swap = %d %d %d\n", locint_swap[0],locint_swap[1],locint_swap[2]);		
	//(void)printf("loc(real) = %d %d %d\n", loc_int[0],loc_int[1],loc_int[2]);	
	//(void)printf("loc_aux   = %d %d %d\n", locint_aux[0],locint_aux[1],locint_aux[2]);	
	//(void)printf("loc**       = %d %d %d\n\n", loc_int[0]+locint_aux[0],loc_int[1]+locint_aux[1],loc_int[2]+locint_aux[2]);	
	
	if(locint_swap[0]!=loc_int[0]||locint_swap[1]!=loc_int[1]||locint_swap[2]!=loc_int[2]){//2
		for (i=0;i<=2;i++){//3		
			locint_swap[i]=loc_int[i];
		}//3
		travelled_cells++;
		//(void)printf("travelled_cells %d, step %d\n",travelled_cells, step);
		//(void)printf("loc_swap = %d %d %d\n", locint_swap[0],locint_swap[1],locint_swap[2]);		
		//(void)printf("loc(real) = %d %d %d\n", loc_int[0],loc_int[1],loc_int[2]);	
	}//2
															

	loc_aux[0]=(l_coh)*loc_int[0];
	loc_aux[1]=(l_coh)*loc_int[1];
	loc_aux[2]=(l_coh)*loc_int[2];	
		
	delete Geom;
	
	return 1.0;
}//1

float Magnetic::MagneticField(){//1

	Geometry * Geom = new Geometry();

	if(	strcmp(MFMODEL,"Spherical Structure") == 0 || 	strcmp(MFMODEL,"Cubic Structure") == 0){//2

		GridPositioning();
	
		if(	strcmp(MFMODEL,"Spherical Structure") == 0){//3 //cellular structure
			SphericalStructure();
		}//3
		else if(strcmp(MFMODEL,"Cubic Structure") == 0){//3 //cubic dimain
			CubicStructure();	
		}//3

	}//2	
//else if(MF_model==XXXXXX){//11 //dipole magnetic field/
	//for (i=0;i<=2;i++){//2 
		//B[i] = - 2*( (Geom->ScalProd(m_mag,x))*(x[i]))/(pow(Geom->Mod(x),5)) - ( m_mag[i] )/( pow(Geom->Mod(x),3));		
	//}//2
//}//11	
	else{//2 //uniform magnetic field/
		RegularComponent();	
	}//2					
	delete Geom;
	return 1.0;
}//1

float Magnetic::RegularComponent(){//1
	
	B[0]=0;
	B[1]=0;
	B[2]=B_reg;
	
	return 1.0;
}//1

float Magnetic::SphericalStructure(){//1
	
	Geometry * Geom = new Geometry();
		
	if(Geom->DistPP(x0,loc_aux)<= (1-s_depth)*(l_coh/2)){//2
				
		//(void)printf("-----------------inside----------------------\n");
		//(void)printf("loc       = %d %d %d\n", loc_int[0]+locint_aux[0],loc_int[1]+locint_aux[1],loc_int[2]+locint_aux[2]);																
		//(void)printf("dist/lcoh      = %.3lf\n", (Geom->DistPP(x0,loc_aux)/l_coh));									
		////(void)printf("dist/Mpc      = %.3lf\n", (Geom->DistPP(x0,loc_aux)/(mbyMpc)));									

		Theta = Theta_Tensor[(loc_int[0]+locint_aux[0])][(loc_int[1]+locint_aux[1])][(loc_int[2]+locint_aux[2])];
		Phi = Phi_Tensor[(loc_int[0]+locint_aux[0])][(loc_int[1]+locint_aux[1])][(loc_int[2]+locint_aux[2])];
		//(void)printf("%.3lf %.3lf\n", Theta, Phi);

		B[0]=B_r*sin(Theta)*cos(Phi);
		B[1]=B_r*sin(Theta)*sin(Phi);
		B[2]=B_r*cos(Theta);
	
	}//2 
	else{//2
	
		//(void)printf("-------------------outside-------------------\n");
		//(void)printf("loc       = %d %d %d\n", loc_int[0]+locint_aux[0],loc_int[1]+locint_aux[1],loc_int[2]+locint_aux[2]);																
		//(void)printf("dist/lcoh      = %.5lf\n", (Geom->DistPP(x0,loc_aux)/l_coh));									
		////(void)printf("dist /Mpc     = %.5lf\n", (Geom->DistPP(x0,loc_aux)/(mbyMpc)));									

		nearbycell=0;
		
		for (i=loc_int[0]-1;i<=loc_int[0]+1;i++){//2
				
			ijk_aux[0]=Geom->ScalMod(i+locint_aux[0]);
			
			if(i+locint_aux[0]==-1)
				ijk_aux[0]=Unv_scl;
			
			if(i%Unv_sclplusOne_int==0)
				ijk_aux[0]=0;
							
			for (j=loc_int[1]-1;j<=loc_int[1]+1;j++){//3
				
				ijk_aux[1]=Geom->ScalMod(j+locint_aux[1]);
				
				if(j+locint_aux[1]==-1)
					ijk_aux[1]=Unv_scl;
					
				if(j%Unv_sclplusOne_int==0)
					ijk_aux[1]=0;
							
				for (k=loc_int[2]-1;k<=loc_int[2]+1;k++){//4
					
					ijk_aux[2]=Geom->ScalMod(k+locint_aux[2]);
			
					if(k+locint_aux[2]==-1)
					ijk_aux[2]=Unv_scl;
					
					if(k%Unv_sclplusOne_int==0)
						ijk_aux[2]=0;
																		
					//(void)printf(" loc [%d] = %d %d %d\n", nearbycell,i+locint_aux[0],j+locint_aux[1],k+locint_aux[2]);									

					loc_aux[0]= l_coh*i;
					loc_aux[1]=l_coh*j;
					loc_aux[2]=l_coh*k;
                    dist_cell[nearbycell]=Geom->DistPP(x0,loc_aux)/l_coh;
					
					//(void)printf("dist [%d] = %.4lf \n", nearbycell,dist_cell[nearbycell]);
					//(void)printf("loc = %d %d %d\n\n", ijk_aux[0],ijk_aux[1],ijk_aux[2]);									
				
					theta_cell[nearbycell] = Theta_Tensor[ijk_aux[0]][ijk_aux[1]][ijk_aux[2]];
					phi_cell[nearbycell] = Phi_Tensor[ijk_aux[0]][ijk_aux[1]][ijk_aux[2]];
					
					//(void)printf("%lf %lf\n", Theta_Tensor[ijk_aux[0]][ijk_aux[1]][ijk_aux[2]],Phi_Tensor[ijk_aux[0]][ijk_aux[1]][ijk_aux[2]]); 

					nearbycell++;					
				}//4
			}//3		
	}//2
									
	for (nearbycell = 0; nearbycell <=(nearbycell_num-1) ; nearbycell++){//1 
	    
	    min = nearbycell; 
	    for (jj = (nearbycell+1); jj <= (nearbycell_num-1); jj++) {//2 
	      if(dist_cell[jj] < dist_cell[min]){//3 
	        min = jj; 
	      }//3 
	    }//2
	
		if (nearbycell != min){//2 
	      swap = dist_cell[nearbycell]; 
	      dist_cell[nearbycell] = dist_cell[min]; 
	      dist_cell[min] = swap;
	      
	      swap = theta_cell[nearbycell]; 
	      theta_cell[nearbycell] = theta_cell[min]; 
	      theta_cell[min] = swap;
	      
	      swap = phi_cell[nearbycell]; 
	      phi_cell[nearbycell] = phi_cell[min]; 
	      phi_cell[min] = swap; 
	    }//2 
	  }//1 
	
	dist_inverse_sum =0.0;
	Theta=0.0;
	Phi=0.0;
	for (nearbycell=0;nearbycell<=(n_cells-1);nearbycell++){//6
										
		//(void)printf("dist min[%d] = %.4lf \n", nearbycell,dist_cell[nearbycell]);
	
		dist_inverse_sum +=pow(dist_cell[nearbycell],-1);
												
	    Theta+=(theta_cell[nearbycell])/dist_cell[nearbycell];
		Phi+=(phi_cell[nearbycell])/dist_cell[nearbycell];
						
		}//6					
		
		Theta=Theta/dist_inverse_sum;
		Phi= Phi/dist_inverse_sum;
						
		B[0]=B_r*sin(Theta)*cos(Phi);
		B[1]=B_r*sin(Theta)*sin(Phi);
		B[2]=B_r*cos(Theta);
				
	}//2
	delete Geom;
	return 1.0;
}//1

float Magnetic::CubicStructure(){//1
	
	Theta = Theta_Tensor[(loc_int[0]+locint_aux[0])][(loc_int[1]+locint_aux[1])][(loc_int[2]+locint_aux[2])];
	Phi = Phi_Tensor[(loc_int[0]+locint_aux[0])][(loc_int[1]+locint_aux[1])][(loc_int[2]+locint_aux[2])];
	//(void)printf("%.3lf %.3lf\n", Theta, Phi);
	B[0]=B_r*sin(Theta)*cos(Phi);
	B[1]=B_r*sin(Theta)*sin(Phi);
	B[2]=B_r*cos(Theta);
	
	return 1.0;
}//1
