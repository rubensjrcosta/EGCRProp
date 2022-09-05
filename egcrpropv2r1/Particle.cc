//Particle.cc/

//local headers/
#include "Particle.h"

FILE *energy_loss_in;

//the constructor of the class/
Particle::Particle ( ){

}

//the destructor of the class/
Particle::~Particle ( ){

}

//the methods of the class/
float Particle::SettingParticleParameters(){//1
		
	if(b_track==0){//2
		for(i=0;i<=4;i++){//3
			strcpy (PARTICLE_NAME[i], part_name[i]);
		}//3
	}//2
	else{//2		
		for(i=0;i<=4;i++){//3
			strcpy (PARTICLE_NAME[i], antipart_name[i]);
		}//3
	}//2
	
	if(P_spec == 1){//2 //PROTON/	
		Z=1; //atomic number/
		A=1;//mass number/
		strcpy (PARTICLENAME, PARTICLE_NAME[0]);//particle name/				
		strcpy (NUCLEINAME, nuclei_name[0]);//particle loss file name/				
		strcpy (NUCLEUSSYMBOL, nuclei_symbol[0]);							
		line_loss=133;//size of the vectors E_z0 and XI_z0/		
		color =1;//Black (color for trajectory ploting: ROOT MACRO)
		kcolor=920;//kGray					
	}//2
	else if(P_spec == 2) {//2 //HELIUM/
		Z=2; 
		A=4;
		strcpy (PARTICLENAME, PARTICLE_NAME[1]);
		strcpy (NUCLEINAME, nuclei_name[1]);							
		strcpy (NUCLEUSSYMBOL, nuclei_symbol[1]);							
		line_loss=129;
		color =4;//Blue			
		kcolor =600;//kBlue
	}//2
	else if(P_spec == 3) {//2 //OXYGEN/
		Z=8; 
		A=16;		
		strcpy (PARTICLENAME, PARTICLE_NAME[2]);
		strcpy (NUCLEINAME, nuclei_name[2]);		
		strcpy (NUCLEUSSYMBOL, nuclei_symbol[2]);											 	
		line_loss=99;
		color =6;//Violet
		kcolor=616;//kMagenta
	}//2
	else if(P_spec == 4) {//2 //SILICON/
		Z=14; 		
		A=28;		
		strcpy (PARTICLENAME, PARTICLE_NAME[3]);
		strcpy (NUCLEINAME, nuclei_name[3]);	
		strcpy (NUCLEUSSYMBOL, nuclei_symbol[3]);											
		line_loss=93;
		color =3;//Green	
		kcolor=416;//kGreen 
	}//2
	else {//2 //IRON/	
		Z=26; 
		A=56;		
		strcpy (PARTICLENAME, PARTICLE_NAME[4]);
		strcpy (NUCLEINAME, nuclei_name[4]);
		strcpy (NUCLEUSSYMBOL, nuclei_symbol[4]);										 	
		line_loss=106;
		color =2;//Red
		kcolor=632;//kRed			
	}//2
	
	q=e*(Z); //total charge [C]/
	if(b_track==1) q =-q;//propagating antiparticle/
	
	m0 = (m_p)*(A);//particle mass ~ nuclear mass ~ A*(proton mass)/
	
	if(strcmp(MRIGIDITY,"OFF") == 0){Z_aux=1;}
	else{Z_aux=Z;}//the energy is multiplied by the nuclear charge number Z/
	
	if(strcmp(ENERGYLOSS, "ON") == 0){SettingEnergyLossesParameters();}//loading the particle's energy loss parameters/		
				
	return 1.0;
}//1

float Particle::SettingEnergyLossesParameters(){//1

		double particle_loss;		
		XI_ad = 4282.749;
		sprintf(archive,"%s_loss.dat",NUCLEINAME);
		energy_loss_in=fopen(archive,"r");
		if(energy_loss_in==NULL){(void)printf("\n Can't open %s_loss.dat\n\n",NUCLEINAME);exit(EXIT_FAILURE);return(0);}
		particle_loss=fgetc(energy_loss_in);
		(void)ungetc(particle_loss,energy_loss_in);

		E_z0=(double*) malloc(line_loss*sizeof(double));if(!E_z0){printf("**Alocation Error E_z0\n**");exit(-1);}
		XI_z0=(double*) malloc(line_loss*sizeof(double));if(!XI_z0){printf("**Alocation Error XI_z0\n**");exit(-1);}
					
		for (i=0;i<=line_loss-1;i++){//2			
			(void)fgets(line,sizeof(line),energy_loss_in);
			sscanf(line,"%lf %lf",&E_z0[i],&XI_z0[i]);
		}//2
		particle_loss=fgetc(energy_loss_in);		
		(void)fclose(energy_loss_in);
		
		//extremes X_loss and E_loss
		XI_extr[0] = XI_z0[0];XI_extr[1] = XI_z0[line_loss-1];
		 E_extr[0] =  E_z0[0]; E_extr[1] =  E_z0[line_loss-1];			
		 	
return 1.0;
}//1

float Particle::Freeing(){//1
if(strcmp(ENERGYLOSS, "ON") == 0){free(E_z0);free(XI_z0);};
return 1.0;
}//1
