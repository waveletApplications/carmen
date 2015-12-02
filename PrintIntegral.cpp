/***************************************************************************
                          PrintIntegral.cpp  -  description
                             -------------------
    begin                : Wed Jun 27 2001
    copyright            : (C) 2001 by Olivier Roussel & Alexei Tsigulin
    email                : roussel@ict.uni-karlsruhe.de, lpsoft@mail.ru
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
/*
______________________________________________________________________________________________

Print integral values into file "FileName"
______________________________________________________________________________________________

*/

void PrintIntegral(const char* FileName)
{
	// --- Local variables ---

	real 	t;				// time
	FILE	*output;  // output file
	int i;					// counter
	real Volume=1.;	// Total volume	

	// --- Open file ---

	if ( (IterationNo == 0) ? (output = fopen(FileName,"w")) : (output = fopen(FileName,"a")) )
	{
	// HEADER

		if (IterationNo == 0)
		{

			fprintf(output, "#");
			fprintf(output, TXTFORMAT, " Time");

 			switch(EquationType)
			{
				case 1:
				case 2:
					fprintf(output, TXTFORMAT, "Error Max");
					fprintf(output, TXTFORMAT, "Error L2");
					fprintf(output, TXTFORMAT, "Error L1");
					if (Dimension == 1)
					{
             					fprintf(output, TXTFORMAT, "Momentum: exact");
						fprintf(output, TXTFORMAT, "computed");
						fprintf(output, TXTFORMAT, "Energy: exact");
						fprintf(output, TXTFORMAT, "computed");
					}
					break;

				case 3:
				case 5:
					fprintf(output, TXTFORMAT, "Flame velocity ");
					break;

				case 4:
					fprintf(output, TXTFORMAT, "Reaction rate");
					fprintf(output, TXTFORMAT, "Radius");
					fprintf(output, TXTFORMAT, "Flame velocity");
					break;

				case 6:
					fprintf(output, TXTFORMAT, "CFL");
					fprintf(output, TXTFORMAT, "Momentum");
					fprintf(output, TXTFORMAT, "Energy");
					fprintf(output, TXTFORMAT, "Enstrophy");
					//fprintf(output, TXTFORMAT, "Int. Vorticity");
					//fprintf(output, TXTFORMAT, "Baroclinic Ef.");
					//if (CMin[1]==3)
					//	fprintf(output, TXTFORMAT, "Force");
					break;
			};

			if (Multiresolution)
			{
				fprintf(output, TXTFORMAT, "Memory comp.");
				fprintf(output, TXTFORMAT, "CPU comp.");
				if (ExpectedCompression != 0. || CVS)
					fprintf(output, TXTFORMAT, "Tolerance");
				if (CVS)
					fprintf(output, TXTFORMAT, "Av. Pressure");
					
			}

			if (!ConstantTimeStep)
			{
				if (StepNb == 3) fprintf(output, TXTFORMAT, "RKF Error");
				fprintf(output, TXTFORMAT,"Next time step");
				fprintf(output, "%13s  ", "IterationNo");
				fprintf(output, "%13s  ", "IterationNb");
			}

			fprintf(output,"\n");

		}

		if (ConstantTimeStep)
			t=IterationNo*TimeStep;
		else
			t = ElapsedTime;		

		fprintf(output, FORMAT, t);

		switch (EquationType)
		{
			// ADVECTION, BURGERS
			// Integral values are errors
			case 1:
			case 2:
				fprintf(output, FORMAT, ErrorMax);				
				fprintf(output, FORMAT, ErrorL2);			
				fprintf(output, FORMAT, ErrorMid);
				if (Dimension == 1)
				{
					fprintf(output, FORMAT, ExactMomentum);
					fprintf(output, FORMAT, GlobalMomentum);
					fprintf(output, FORMAT, ExactEnergy);
					fprintf(output, FORMAT, GlobalEnergy);
				}	
				break;

			// FLAME FRONT AND INTERACTION FLAME-CURL
			// Integral value is flame velocity
			case 3:
			case 5:
				fprintf(output, FORMAT, FlameVelocity);
				break;

			// FLAME BALL
			// Integral value are gloabal reaction rate, average radius and flame velocity
			case 4:
				if (IterationNo > 3*Refresh)
					AverageRadius = 0.5*(AverageRadius+2*PreviousAverageRadius-PreviousAverageRadius2);

				if (IterationNo > 2*Refresh)
					FlameVelocity = (AverageRadius-PreviousAverageRadius2)/(2*TimeStep*Refresh);
				else
					FlameVelocity = 0.;

        			PreviousAverageRadius2 = PreviousAverageRadius;
				PreviousAverageRadius	 = AverageRadius;

				fprintf(output, FORMAT, GlobalReactionRate);
				fprintf(output, FORMAT, AverageRadius);
				if (IterationNo > 2*Refresh)
					fprintf(output, FORMAT, FlameVelocity);
				else
					fprintf(output, FORMAT, 0.);

				break;

			// NAVIER-STOKES

			case 6:

				// --- Compute total volume ---

				for (i=1; i<= Dimension; i++)
					Volume *= fabs(XMax[i]-XMin[i]);

				// Print CFL
				fprintf(output, FORMAT,EigenvalueMax*TimeStep/SpaceStep);

				// Print momentum and energy
				fprintf(output, FORMAT, GlobalMomentum);
				fprintf(output, FORMAT, GlobalEnergy);
				fprintf(output, FORMAT, GlobalEnstrophy);
				//fprintf(output, FORMAT, IntVorticity);
				//fprintf(output, FORMAT, BaroclinicEffect);
				
				//if (CMin[1]==3)
				//	fprintf(output, FORMAT, ForceX);
				break;
    };

		if (Multiresolution)
		{
			fprintf(output, FORMAT, (1.*CellNb)/(1<<(ScaleNb*Dimension)));
			fprintf(output, FORMAT, CPUTime.CPUTime()/(IterationNo*FVTimeRef));

			if (ExpectedCompression != 0.)
				fprintf(output, FORMAT, Tolerance);
		
    			if (CVS)
			{
				fprintf(output, FORMAT, ComputedTolerance(ScaleNb));
				fprintf(output, FORMAT, QuantityAverage.value(1));
			}
    	}

		if (!ConstantTimeStep)
		{
			if (StepNb == 3) fprintf(output, FORMAT, RKFError);
			fprintf(output, FORMAT, TimeStep);
			fprintf(output, "%13i  ", IterationNo);
			fprintf(output, "%13i  ", IterationNb);
		}

		fprintf(output, "\n");			
		fclose(output);
	}
	else
	{
		cout << "PrintIntegral.cpp: In method `void PrintIntegral(Node*, char*)´:\n";
		cout << "PrintIntegral.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [PrintIntegral.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
