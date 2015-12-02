/***************************************************************************
                          AdaptTimeStep.cpp  -  description
                             -------------------
    begin                : jeu mar 25 2004
    copyright            : (C) 2004 by Olivier Roussel and Alexei Tsigulin
    email                : roussel@ict.uni-karlsruhe.de; lpsoft@mail.ru
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// It is included a check of remaining time when the cicle finish for the LTS-case
// 18/07/2012-CIRM


#include "Carmen.h"

void AdaptTimeStep()
{
	int 	RemainingIterations;
	real 	RemainingTime;
	real	OldTimeStep = TimeStep;
	real 	S = RKFSafetyFactor*(9.*exp(-(real)IterationNo)+1.);

	// Security : do nothing if ConstantTimeStep is true

	if (ConstantTimeStep)
		return;
	

	// Compute remaining time

	RemainingTime = PhysicalTime-ElapsedTime; //MD
	
	// For multiresolution computations with local time stepping, only every 2^a(L-2) iterations

	if (Multiresolution && TimeAdaptivity && (IterationNo%(1<<(TimeAdaptivityFactor*(ScaleNb-4))) != 0)){
		if (RemainingTime <= 0.) //MD
			{
			IterationNb = IterationNo;
			}
			else if (RemainingTime < TimeStep)
			{
			TimeStep = RemainingTime;
			IterationNb = IterationNo + 1;
	                TimeAdaptivity=false; //MD
			}
		return;
}

	// Compute remaining time

	// RemainingTime = PhysicalTime-ElapsedTime;

	
	if (StepNb != 3)
	{
		// Only for Navier-Stokes
		if (EquationType != 6) return;
		
		// In this case, use time adaptivity based on CFL
		TimeStep = CFL*SpaceStep/EigenvalueMax;

		if (Viscosity != 0.)
        		TimeStep = Min(TimeStep, CFL*0.5*Re*SpaceStep*SpaceStep);	
	}
	else
	{
		// For the Runge-Kutta-Fehlberg 2(3) method only
		TimeStep *= exp ( log(RKFAccuracyFactor/RKFError)/3.);
	}
			
	// Limiter
	
	if (TimeStep > (1.+0.5*S)*OldTimeStep)
		TimeStep = (1.+0.5*S)*OldTimeStep;

	if (TimeStep < (1.-0.5*S)*OldTimeStep)
		TimeStep = (1.-0.5*S)*OldTimeStep;
		
	// Recompute IterationNb
	
	if (RemainingTime <= 0.)
	{
		IterationNb = IterationNo;
	}
	else if (RemainingTime < TimeStep)
	{
		TimeStep = RemainingTime;
		IterationNb = IterationNo + 1;
                TimeAdaptivity=false; //MD
	}
	else
	{
		RemainingIterations =  (int)(RemainingTime/TimeStep);
		IterationNb = IterationNo + RemainingIterations;
	}

	return;

}
