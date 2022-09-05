{//1
gROOT->Reset();

//local headers -------------------------------------------------------/
#include "Riostream.h"
#include "defines.h"
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include "globalvariables.h"

//classes -------------------------------------------------------------/
#include "InputParameters.cc"
#include "Particle.cc"
 
FILE *energy_defl_out,*distance_defl_out;

InputParameters* Input=new InputParameters();Particle *Part = new Particle();

Input->ReadingInputParameters();
Part->SettingParticleParameters();	

sprintf(string,"EGCRProp-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d_%d.root",b_track,d_stop, m_rigidity,E_loss,A,B_nG,l_coh,E_EeV,D_Mpc,N_evts,N_iter);
TFile *EGCRProp_ROOT = new TFile(string,"RECREATE");
	
for(i=0;i<=2;i++){//2 											
	x_max[i] = -(D_Mpc)*pow(10,10);
	x_min[i] =  (D_Mpc)*pow(10,10);	
}//2

//plotting the trajectories
if(N_traj != 0){//2	
	
	TCanvas *cavanstraj = new TCanvas("cavanstraj", "Graph Draw Options",65,24,1040,776);gBenchmark->Start("canvas1");
	cavanstraj->Range(-0.131737,-1.719298,1.132958,3.473684);
	cavanstraj->ToggleEventStatus();cavanstraj->ToggleToolTips();
	cavanstraj->SetGridx();cavanstraj->SetGridy();
	cavanstraj->SetTickx(1);cavanstraj->SetTicky(1);
	cavanstraj->SetTopMargin(0.1875);
			
	for(N_iter = 1;N_iter <= N_traj;N_iter++){//Loop iteration

		TPolyLine3D *trajectory = new TPolyLine3D(1000000);
		sprintf(archive,"EGCRProp-traj-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d_%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,l_coh,E_EeV,D_Mpc,N_evts,N_iter);
		trajectory_out=fopen(archive,"r");
		if(trajectory_out==NULL){(void)printf("Can't open EGCRProp-traj-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d_%d.dat",b_track,d_stop, m_rigidity,E_loss,A,B_nG,l_coh,E_EeV,D_Mpc,N_evts,N_iter);return(0);}

		line_count=0;getting=0;
		do{//3
			fgets(line,100,trajectory_out);
			getting = sscanf(line,"%lf %lf %lf",&x[0],&x[1],&x[2]);
			if(getting=1){//4
				//calculating the values of the maximum and the minimum
				for (i=0;i<=2;i++){//5  
					if(x[i]>=x_max[i])x_max[i]=x[i];				
					if(x[i]<=x_min[i])x_min[i]=x[i];							
				}//5  
				trajectory->SetPoint(line_count,x[0],x[1],x[2]);
				line_count++;	
			}//4
		}//3
		while(!feof(trajectory_out)); 
		fclose(trajectory_out);
		trajectory->SetLineColor(color);
		trajectory->Draw();

	}//Loop iteration
	
	TView *view = TView::CreateView(1);
	view->SetRange(x_min[0],x_min[1],x_min[2],x_max[0],x_max[1],x_max[2]);
	TPaveText *pt = new TPaveText(0.1245174,0.775,0.4604247,0.9694444,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextAlign(12);
	pt->SetTextFont(102);
	pt->SetTextSize(0.00996016);
	
	if(strcmp(BTRACK,"OFF") == 0){sprintf(string,"#scale[3]{#color[1]{%s}}",NUCLEINAME);}
	else{sprintf(string,"#scale[3]{#color[1]{anti%s nuclei}}",NUCLEINAME);}text = pt->AddText(string);
	sprintf(string,"#scale[3]{#color[1]{b_{track} = %s, M_{r} = %s, E_{loss} = %s}}",BTRACK,MRIGIDITY,ENERGYLOSS);text = pt->AddText(string);
	sprintf(string,"#scale[3]{#color[1]{E_{%s}%s = %.1lf EeV, D = %.1lf Mpc (%s)}}",EINDEX,MRINDEX,E_EeV,D_Mpc,DSTOP);text = pt->AddText(string);
	sprintf(string,"#scale[3]{#color[1]{B_{r} = %.1lf nG, l_{coh} = %.1lf Mpc}}",B_nG,l_coh);text = pt->AddText(string);
	sprintf(string,"#scale[3]{#color[1]{N_{evts} = %d, N_{traj} = %d}}",N_evts,N_traj);text = pt->AddText(string);
	text->SetTextColor(1);
	pt->Draw();
	
	TAxis3D *rulers = new TAxis3D();
	rulers->SetXTitle("X [Mpc]");rulers->SetYTitle("Y [Mpc]");rulers->SetZTitle("Z [Mpc]");  
	rulers->SetAxisColor(1);rulers->SetLabelColor(1);rulers->SetLabelSize(0.03);
	rulers->SetLabelOffset(-0.005);rulers->SetTitleOffset(1.5);
	rulers->GetXaxis()->SetTitleSize(0.035);
	rulers->GetYaxis()->SetTitleSize(0.035);
	rulers->GetZaxis()->SetTitleSize(0.035);
	rulers->Draw();  

	sprintf(string,"EGCRProp-traj-%d%d%d_Eloss%d_A%d_B%.1f_lcoh%.1f_E%.1f_D%.1f_Nevts%d_%d.png",b_track,d_stop, m_rigidity,E_loss,A,B_nG,l_coh,E_EeV,D_Mpc,N_evts,N_traj);
	cavanstraj->Print(string);
	cavanstraj->cd();
	cavanstraj->SetSelected(cavanstraj);
	cavanstraj->ToggleToolBar();cavanstraj->ToggleEventStatus();cavanstraj->ToggleEditor();
	cavanstraj->Write();gBenchmark->Show("canvas1");

}//2
	
}//1





