/***************************************************************************
                          InitTimeStep.cpp  -  description
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

#include "Carmen.h"

void InitTimeStep()
{

	// --- Init TimeStep --------------------------------------------------------------

	if (TimeStep == 0)
	{
		switch(EquationType)
		{
			// NAVIER-STOKES
  			case 6:
				TimeStep = CFL*SpaceStep/EigenvalueMax;

				if (Viscosity != 0.)
        				TimeStep = Min(TimeStep, CFL*0.5*Re*SpaceStep*SpaceStep);

				break;

		};
	}
		
	// --- Compute number of iterations -----------------------------------------------

	if (PhysicalTime != 0. && IterationNb == 0)
		IterationNb =  (int)(ceil(PhysicalTime/TimeStep));

	// --- Compute Refresh ------------------------------------------------------------

	if (Refresh == 0)
		Refresh = (int)(ceil(IterationNb/(RefreshNb*1.)));

	// --- Compute PrintEvery ---------------------------------------------------------

	if ((PrintEvery == 0)&&(ImageNb != 0))
		PrintEvery = (int)(ceil(IterationNb/(ImageNb*1.)));

	// --- Compute iterations for print  ----------------------------------------------

	if (PrintTime1 != 0.)
		PrintIt1 = (int)(ceil(PrintTime1/TimeStep));

	if (PrintTime2 != 0.)
		PrintIt2 = (int)(ceil(PrintTime2/TimeStep));

	if (PrintTime3 != 0.)
		PrintIt3 = (int)(ceil(PrintTime3/TimeStep));

	if (PrintTime4 != 0.)
		PrintIt4 = (int)(ceil(PrintTime4/TimeStep));

	if (PrintTime5 != 0.)
		PrintIt5 = (int)(ceil(PrintTime5/TimeStep));

	if (PrintTime6 != 0.)
		PrintIt6 = (int)(ceil(PrintTime6/TimeStep));

	// --- Compute FV reference time -----------------------------------------------------------

	if (Multiresolution)
		FVTimeRef = CPUTimeRef(IterationNbRef,ScaleNbRef);

}
