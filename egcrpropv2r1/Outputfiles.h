#ifndef __OUTPUTFILES_H
#define __OUTPUTFILES_H

class Outputfiles
{
private:

public:

//! Constructor of the class.
Outputfiles ( );

//! Destructor of the class.
~Outputfiles ( );

//methods of the Outputfiles class/
float OpeningSummary();
float PrintingSummaryBegin(); 
float PrintingSummaryBody();
float PrintingSummaryEnd();
float PrintingDynamicsQuantities();
float PrintingCosmologicalParameters();
float ClosingSummary();

float OpeningDynamicsParameters();
float PrintingDynamicsParameters();
float ClosingDynamicsParameters();

float OpeningTrajectory();
float PrintingTrajectory();
float ClosingTrajectory();

float SuccessfullyCompleted();
				
};

#endif // __OUTPUTFILES_

