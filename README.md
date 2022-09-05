//****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp//
//****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp//

###################In order to use this software, you will need to have installed:##########################################

GMP library (version >= 4.2.1) - https://gmplib.org/ 
MPFR library(version >=2.3.1) - http://mpfr.org 
GNU Scientific Library (GSL) - https://www.gnu.org/software/gsl/

##############################################Installing:####################################################################

sudo apt-get install libgmp3-dev

sudo apt-get install libmpfr-dev libmpfr-doc libmpfr4 libmpfr4-dbg

sudo apt-get install libgsl0ldbl

##############################################Compiling:#####################################################################

g++ -std=c++11 EGCRProp.cc -o EGCRProp -lmpfr -lgsl -lgslcblas -lgmp

###########################The parameters_input.dat file accepts the following options:######################################

P_spec 
E_EeV 				D_Mpc
E_loss
MF_model
B_nG 				l_coh
s_depth				n_cells
b_track				d_stop			m_rigidity
N_evts				N_traj 
i_method 			t_scal
r_dist				r_seed

###########################################The parameters input description:################################################

P_spec 				= particle species (1 proton; 2 helium; 3 oxygen; 4 silicon; 5 iron)
E_EeV 				= the particle's energy in EeV
D_Mpc 				= the distance in Mpc
E_loss				= the energy losses (0 off; 1 on)
MF_model 			= magnetic field model (1 Cellular Structure,2 Cubic Domain, 3 Uniform Magnetic Field)
B_nG				= magnetic field strength in nG
l_coh 				= coherence length in Mpc
s_depth 			= the skin depth in % of coherence length
n_cells  			= the number of nearby cells to be considered
b_track 			= backtracking (0 off; 1 on) !!!backtracking off works only for the stopping distance mode 1
d_stop 				= stopping distance (1 cut in distance traveled; 2 cut in trajectory radius) 
m_rigidity 			= magnetic rigidity (0 off; 1 on)
N_evts 				= number of launched particles
N_traj 				= number of plotted trajectories 
i_method 			= equation of motion integration (1 Runge-Kutta 4th, 2 Boris, 3 Vay*, 4 New Euler*)
t_scal 				= time scale
r_dist				= random distribution (1 uniform, 2 gaussian, with the parameter sigma=1)
r_seed 				= random seed (0 based on time, 1 fixed (123 as default, mySeed (line 17, file MPFRGSLvariables.h))

# * the number of decimal places is defined by variable digits, line 6 in the file MPFRGSLvariables.h
####################The three-dimensional viewer ande histogram generation compilation:#####################################

root MACRO_ROOT.cpp
		  
####################The tests were made from Ubuntu 14.04.5 LTS with the compiler g++ 4.8.4 ################################

//****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp//
//****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp****EGCRProp//

