**EGCRProp**<a name="cite_ref-1"></a>[<sup>[1]</sup>](#cite_note-1) a is a Monte Carlo code written in C ++ for the propagation of Ultra-High Energy
Cosmic Rays (UHECRs) through Extragalactic Magnetic Fields (EMFs). In the simulation, it is
calculated the three-dimensional trajectory of nuclei (from proton to iron) through unitary cells of
coherent magnetic fields, as well as the particle’s dynamic parameters: the trajectory radius or the
traveled distance, the energy, and the deflection angle, are recorded after each event to be further
analyzed.

#### In order to use this software, you will need to have installed:

GMP library (version >= 4.2.1) - https://gmplib.org/ <br>
MPFR library(version >=2.3.1) - http://mpfr.org <br>
GNU Scientific Library (GSL) - https://www.gnu.org/software/gsl/ <br>

#### Installing:

sudo apt-get install libgmp3-dev <br>

sudo apt-get install libmpfr-dev libmpfr-doc libmpfr4 libmpfr4-dbg <br>

sudo apt-get install libgsl0ldbl <br>

#### Compiling:

g++ -std=c++11 EGCRProp.cc -o EGCRProp -lmpfr -lgsl -lgslcblas -lgmp <br>

#### The parameters_input.dat file accepts the following options:

P_spec <br>
E_EeV 				D_Mpc <br>
E_loss <br>
MF_model <br>
B_nG 				l_coh <br>
s_depth				n_cells <br>
b_track				d_stop			m_rigidity <br>
N_evts				N_traj <br>
i_method 			t_scal <br>
r_dist				r_seed <br>

#### The parameters input description:

P_spec 				= particle species (1 proton; 2 helium; 3 oxygen; 4 silicon; 5 iron) <br>
E_EeV 				= the particle's energy in EeV <br>
D_Mpc 				= the distance in Mpc <br>
E_loss				= the energy losses (0 off; 1 on) <br>
MF_model 			= magnetic field model (1 Cellular Structure,2 Cubic Domain, 3 Uniform Magnetic Field) <br>
B_nG				= magnetic field strength in nG <br>
l_coh 				= coherence length in Mpc <br>
s_depth 			= the skin depth in % of coherence length <br>
n_cells  			= the number of nearby cells to be considered <br>
b_track 			= backtracking (0 off; 1 on) !!!backtracking off works only for the stopping distance mode 1 <br>
d_stop 				= stopping distance (1 cut in distance traveled; 2 cut in trajectory radius) <br>
m_rigidity 			= magnetic rigidity (0 off; 1 on) <br>
N_evts 				= number of launched particles <br>
N_traj 				= number of plotted trajectories <br> 
i_method 			= equation of motion integration (1 Runge-Kutta 4th, 2 Boris, 3 Vay*, 4 New Euler<a name="cite_ref-2"></a>[<sup>[2]</sup>](#cite_note-2)) <br>
t_scal 				= time scale <br>
r_dist				= random distribution (1 uniform, 2 gaussian, with the parameter sigma=1) <br>
r_seed 				= random seed (0 based on time, 1 fixed (123 as default, mySeed (line 17, file MPFRGSLvariables.h)) <br>

#### The three-dimensional viewer ande histogram generation compilation:

root MACRO_ROOT.cpp <br>
		  
#### The tests were made from Ubuntu 14.04.5 LTS with the compiler g++ 4.8.4 

<a name="cite_note-1"></a>1. [](#cite_ref-1)Please, acknowledge the use of the EGCRProp code by citing R.P. Costa Junior, and M.A. Leigui de Oliveira. "A Numerical Model for the Propagation of Ultra-High Energy Cosmic Rays through Extragalactic Magnetic Fields", Proc. of the 35th ICRC (ICRC2017), PoS(ICRC2017)477.

<a name="cite_note-2"></a>2. [](#cite_ref-2)The number of decimal places is defined by variable digits, line 6 in the file MPFRGSLvariables.h <br>
