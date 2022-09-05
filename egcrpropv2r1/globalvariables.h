#ifndef __GLOBALVARIABLES_H
#define __GLOBALVARIABLES_H

//files/
FILE *summary,*EGCRProp_out,*trajectory_out,*successfully_completed; 

//global variables ----------------------------------------------------/
clock_t time_start;//time CPU control/

const double E_min =  0.1*(JbyEeV);//minimum energy 0.1 EeV/
const double E_max = 30.0*(JbyZeV);//maximum energy 30 ZeV/
	
static char archive [170],line[100];//files variables/
int N_evts;//number of particles launched/ 
int N_traj;//number of trajectories ploted/
int b_track;//backtracking mode/
int d_stop;//stopping distance mode/
int m_rigidity; //magnetic rigidity mode/
int E_loss;//energy loss mode/  
int P_spec;//particle specie/
int MF_model;//magnetic field model/
double B_nG;//the magnetic field strength [nG]/
double l_coh;//the coherence length [Mpc]/
double s_depth;//the skin depth in % of coherence length/
int n_cells;//the number of nearby cells to be considered/
double t_scal;//time scale/
double E_EeV;//Energy in EeV/
double D_Mpc ;//Distance in Mpc/
int i_method;//equation of motion integration method/
int r_dist;//random distribution/
int r_seed;char SEED[50];//random seed/
int A; //particle's mass number/ 
int Z,Z_aux; //particle's atomic number/
double q; //particle's charge [C]/
double m0; //particle's rest mass [Kg]/
double cosmo_parameters[4];//cosmological parameters/
double H0 ;//Hubble constante [m/s/Mpc]/
double OmegaM ;//total matter density/
double OmegaL ;//dark energy density/
double OmegaK ;//density parameter (the Universe curvature)/
double t_H ;//Hubble time [s]/
double D_H ;//Hubble distance [Mpc]/
double D_C;//comoving distance [Mpc]/
double DC_avg;//average(_avg) comoving distance [Mpc]/
double z0,z,z_aux;//redshift/
double delta_z;
int i,j,k;//counters/
int step;//equation of motion step/
int cicles;//number of cicles
int N_iter;//number of iterations/ 
double D,D_max,D_avg,D_aux; //traveled distance [m]/
double D0,D0_avg;//rectilinear distance [Mpc]/

double t; //time [s]/
double delta_t; //time step [s]/
double B[3];//Magnetic Field [T]/
double B_r;//random component magnitude of the magnetic field [nG]/
double B_reg;//regular component magnitude of the magnetic field [nG]/
double *E_z0,*XI_z0;//Energy at z=0, Energy loss length at z=0/
double XI_ad;//Adiabatic energy loss length/
double XI_extr[2],E_extr[2];//E_z0 and XI_z0 extremes/
int line_loss;
int line_count;
static double Theta_Tensor[171][171][171],Phi_Tensor[171][171][171]; //magnetic field tensors/
double theta,phi;//theta and phi random numbers (spherical coordinates)/
const int Unv_scl=170;//Universe scale (maximum 170): Universe volume Unv_scl^{3} in Mpc/
int travelled_cells;
double u0[3],u0_mod,u[3],u0_ref[3];//particle's vector u [(kg*m/)s]/
double p0[3],p0_mod,p[3],p0_ref[3];//particle's momentum [(kg*m/)s]/
double E0,E,E_avg,E_dev;//particle's energy [J]/
double Gamma;//Lorentz Factor/
double omega_B;//cyclotron frequency [s^-1]/
double r_L,rL_error; //Larmor radius [Mpc]/
double x[3],x0[3];//particle's position [m]/
double delta_x;//displacement in each step [m]/
double O[3],Origin;//trajectory origin [m]/
 
double DbyD0,DbyD0_avg;//distances ratio (D)
double Dampl;//distance amplification/
double defl,defl_avg,defl_dev,defl_err;//simulated  deflection angle/ 
double defl_1,defl_N;//analytical deflection angle for one (1) and (N) N cells [degrees]/ 
double  N_aux;	

int zcut_counts;
int E1cut_counts, E2cut_counts;
int cut;

char D_stop[2][50]={"traveled distance","trajectory radius"}, DSTOP[50];
char E_index[3][50]={"A","I",""}, EINDEX[50];
char MR_index[2][50]={"","/Z"}, MRINDEX[50];
char EMR_index[2][50]={"particle's energy","particle's magnetic rigidity"}, EMRINDEX[50];
char onoff_btrack[2][50]={"OFF","ON"},BTRACK[50];
char onoff_Mr[2][50]={"OFF","ON"}, MRIGIDITY[50];
char onoff_Eloss[2][50]={"OFF","ON"}, ENERGYLOSS[50];
char mf_model[3][50]={"Spherical Structure","Cubic Structure","Uniform Magnetic Field"}, MFMODEL[50];
char int_method[4][50]={"Runge-Kutta4th","Boris","Vay","NewEuler"}, INTMETHOD[50];
char part_name[5][100]={"proton ","helium nuclei, ⁴He:","oxygen nuclei, ¹⁶O:","silicon nuclei, ²⁸Si:","iron nuclei, ⁵⁶Fe:"}; 
char antipart_name[5][100]={"antiproton","antihelium nuclei, ⁴He⁻²:","antioxygen nuclei, ¹⁶O⁻⁸:","antisilicon nuclei, ²⁸Si⁻¹⁴:","antiiron nuclei, ⁵⁶Fe⁻²⁶:"}; 
char PARTICLE_NAME[5][100],PARTICLENAME[100];
char nuclei_symbol[5][10]={"proton","He","O","Si","Fe"},NUCLEUSSYMBOL[10]; 
char nuclei_name[5][100]={"proton","helium","oxygen","silicon","iron"},NUCLEINAME[100]; 
char rand_dist[2][100]={"uniform","gaussian"}, RANDDISTRIBUTION[100];
char RANDGENERATOR[50];
//printing parameters
int color,kcolor;
int getting;
char string[100];
double x_min[3],x_max[3];

//histogram parameters
double  defl_max,defl_min; 
double  enrg_max,enrg_min;
int bins = 300;

//loops 
int loop=0;
int energy_count, distance_count;

int distance_count_MAX=2;
int energy_count_MAX=2;
int particle_count_MAX=2;

double Distance[7] = {3.8,17.1,35.0,50.0,60.0,75.0,100.0};//3.8 (Centaurus A distance cite[]), 17.1 (Virgo A distance cite[])

double Energy[6] = {1000.0,300.0,100.0,30.0,10.0,3.0};
////double Energy[10] = {58.0, 63.0, 64.0, 66.0, 69.0, 70.0, 79.0, 80.0, 84.0, 148.0};/*ten nearby events of Centaurus A from Auger Observatory cite[]*/
////double Energy[18] = {84.7, 74.5, 74.0, 73.1, 72.5, 69.5, 66.7, 64.8,62.7,61.5,60.7,60.0,59.5,58.8,58.6,53.3,52.4,52.1};//ten nearby events of Centaurus A from Auger Observatory [AUGER2015]

			
#endif // __GLOBALVARIABLES_H
