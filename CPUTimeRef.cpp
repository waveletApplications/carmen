/***************************************************************************
                          CPUTimeRef.cpp  -  description
                             -------------------
    begin                : Sun Jul 7 2002
    copyright            : (C) 2002 by Olivier Roussel and Alexei Tsigulin
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

/*
____________________________________________________________________________

Compute the reference CPU time with the finite volume solver.
The output is the CPU time for 1 iteration.
____________________________________________________________________________

*/
double CPUTimeRef(int iterations, int scales)
{
	// --- Local variables ---------------------------------------------------

	int OldIterationNb=0;
	int OldScaleNb=0;
	real OldTimeStep=0.;
	bool ConstantTimeStepOld=ConstantTimeStep;
		
	double result=0.;

  	Timer CPURef;
	FineMesh* MeshRef;

	// --- Execution ---------------------------------------------------------

	// Toggle on : Compute reference CPU time

	ComputeCPUTimeRef = true;
	
	// Toggle off : Constant time step
	
	ConstantTimeStep = true;

	// backup values of IterationNb and ScaleNb

  	OldIterationNb 	= IterationNb;
	OldScaleNb 	= ScaleNb;
	OldTimeStep 	= TimeStep;

	// use reference values
	IterationNb 	= iterations;
	ScaleNb     	= scales;
	TimeStep	= 0.;

	one_D=1; two_D=1;
	if (Dimension >= 2) one_D=1<<ScaleNb;
	if (Dimension == 3) two_D=1<<ScaleNb;
 
	// init mesh
  	MeshRef = new FineMesh;

  	// Iterate on time

	for (IterationNo = 1; IterationNo <= IterationNb; IterationNo++)
	{
		// start timer
		CPURef.start();

		// Compute time evolution
		TimeEvolution(MeshRef);

		// check CPU Time
	  	CPURef.check();

		// stop timer
		CPURef.stop();
  	}
		
	// Compute CPUTimeRef
	result = CPURef.CPUTime();
	result *= 1./IterationNb;
	result *= 1<<(Dimension*(OldScaleNb-ScaleNb));

	// delete MeshRef
	delete MeshRef;
	
	// restore values of IterationNb and ScaleNb
  	IterationNb = OldIterationNb;
	ScaleNb = OldScaleNb;
	TimeStep = OldTimeStep;
	IterationNo = 0;

	one_D=1; two_D=1;
	if (Dimension >= 2) one_D=1<<ScaleNb;
	if (Dimension == 3) two_D=1<<ScaleNb;
 
	// Toggle off : Compute reference CPU time

	ComputeCPUTimeRef = false;

	// Restore the value of ConstantTimeStep
	
	ConstantTimeStep = ConstantTimeStepOld;

	return result;
}
