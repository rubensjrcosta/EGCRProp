#ifndef __DEFINES_H
#define __DEFINES_H

//physical constants from Particle Physics Booklet 2012 ---------------------------------------------/
const double c   =	299792458.0;  			 //speed of light [m/s]
const double m_p =	pow(10,-27)*1.672621777; //proton mass 	   [Kg]      
const double m_e =	pow(10,-31)*9.10938291;  //electron mass   [Kg]      
const double e   =	pow(10,-19)*1.602176565; //charge unit      [C]  
const double pi  =	3.141592653589793238; 		 //pi number
const double pc  =	pow(10,16)*3.0856776; 	 //definition of parsec in [m]


//units conversions ---------------------------------------------------------------------------------/

const double mbypc  =	pow(10,16)*3.0856776; 	 //parsec 	   [pc] to meter [m]
const double mbykpc =	pow(10,3)*(mbypc); 		 //kiloparsec [kpc] to meter [m]       
const double mbyMpc =	pow(10,6)*(mbypc); 		 //megaparsec [Mpc] to meter [m]
const double mbyGpc =	pow(10,9)*(mbypc); 		 //gigaparsec [Gpc] to meter [m]

const double JbyeV  =	pow(10,-19)*1.602176565; //electronvolt       [eV] to joule [J] 
const double JbykeV =	pow(10, 3)*(JbyeV);	     //kiloelectronvolt  [keV] to joule [J] 
const double JbyMeV =	pow(10, 6)*(JbyeV); 	 //megaelectronvolt  [MeV] to joule [J] 
const double JbyGeV =	pow(10, 9)*(JbyeV); 	 //gigaelectronvolt  [GeV] to joule [J] 
const double JbyTeV =	pow(10,12)*(JbyeV); 	 //teraelectronvolt  [TeV] to joule [J] 
const double JbyPeV =	pow(10,15)*(JbyeV); 	 //petaelectronvolt  [PeV] to joule [J] 
const double JbyEeV =	pow(10,18)*(JbyeV); 	 //exaelectronvolt   [EeV] to joule [J] 
const double JbyZeV =	pow(10,21)*(JbyeV); 	 //zettaelectronvolt [ZeV] to joule [J] 

const double TbyG 	=	pow(10, -4); 	//Gauss 	    [G]  to Tesla [T] 
const double TbymuG =	pow(10,-10); 	//microGauss  [muG]  to Tesla [T] 
const double TbynG 	=	pow(10,-13); 	//nanoGauss    [nG]  to Tesla [T] 
const double TbypG 	=	pow(10,-16); 	//picoGauss    [pG]  to Tesla [T] 

const double MeVc2bykg =	pow(10,-30)*1.78266181; //Mega-electron-Volts/c2 [MeV/c2] to Kilogram [kg]
const double degbyrad  =	pow(pi,-1)*180; 		//radians to degrees 

//Use always parentheses: (J_EeV), (MeVc2bykg), ..., instead of J_EeV, MeVc2bykg ... !!!/
			
#endif // __DEFINES_H


