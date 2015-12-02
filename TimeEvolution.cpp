/***************************************************************************
                          TimeEvolution.cpp  -  description
                             -------------------
    begin                : Thu Jun 7 2001
    copyright            : (C) 2001 by Olivier Roussel & Alexei Tsigulin
    email                : roussel@ict.uni-karlsruhe.de, lpsoft@mail.ru
 ***************************************************************************/

#include "Carmen.h"
/*
______________________________________________________________________________________________

Time evolution for finite volume with multiresolution
______________________________________________________________________________________________

*/
void TimeEvolution(Node* Root)
{

		// --- Smooth data --- 

		if (EquationType==6 && SmoothCoeff != 0.) 
			Root->smooth();
		
		// --- Store cell-average values of leaves ---
		Root->store();
		
		// --- Refresh tree structure ---
		RefreshTree(Root);
	
		// --- Only for Navier-Stokes ---

		if (EquationType==6) 
		{
			// For the first iteration, compute initial gradient values
			
			if (IterationNo == 1 && Viscosity !=0.) Root->computeGradient();
			
			// For McMP3 and mixed McMP3 / McCormack, store gradient values into temporary ones
			
			if (SchemeNb > 5 && Viscosity != 0.)
				Root->storeGrad();
		}

		for (StepNo = 1; StepNo <= StepNb; StepNo++)
		{
			// --- Compute divergence ---
			Root->computeDivergence();

			// --- Runge-Kutta step ---
			Root->RungeKutta();
			
			// --- Compute gradient of conservative quantities (Navier-Stokes only) ---
			if (EquationType == 6 && Viscosity !=0.) Root->computeGradient();
				
			// --- Refresh tree structure ---
			RefreshTree(Root);
		}

		// --- Check stability ---
		Root->checkStability();

		// --- Compute integral values ---
		Root->computeIntegral();
		
		// --- Compute total number of cells and leaves ---

		TotalCellNb += CellNb;
		TotalLeafNb += LeafNb;

		// --- Compute elapsed time and adapt time step ---

		ElapsedTime += TimeStep;
		if (!ConstantTimeStep) AdaptTimeStep();

}
/*
______________________________________________________________________________________________

Time evolution for finite volume on fine grid (no multiresolution)
______________________________________________________________________________________________

*/
void TimeEvolution(FineMesh* Root)
{

		// --- Store cell-average values into temporary ---
		Root->store();
		
		// --- Compute velocity gradient (Navier-Stokes only) ---
		if (EquationType == 6) 
		{
			if (IterationNo == 1 && Viscosity !=0.)	{
          		Root->computeGradient(1); //first compute gradient for neighbours cells (Input Parameter = 1)
  	  		CPUExchange(Root, SendGrad); //then send it to another processors
          		Root->computeGradient(0); //now, compute gradien for the internal cells while inter-CPU exchanges
                        //are doing in the background (Input Parameter = 0)
          		CommTimer.start();        //Start Communication Timer, only for perfomance analysize
#if defined PARMPI
          		if (MPIRecvType == 1) MPI_Waitall(4*Dimension,req,st);     //Waiting while inter-CPU exchanges are finished
#endif
          		CommTimer.stop();         //Stop Communication Timer
        	}
		
		// For McMP3 and mixed McMP3 / McCormack, store gradient values into temporary ones
			
		if (SchemeNb > 5 && Viscosity !=0.) Root->storeGrad();
		}

		for (StepNo = 1; StepNo <= StepNb; StepNo++)
		{
			// --- Compute divergence for neighbour cells ---
			//The same conception with background computations, see upper...
			Root->computeDivergence(1);
			// --- Runge-Kutta step for neighbour cells ---
			Root->RungeKutta(1);
			// --- Start inter-CPU exchanges ---
			CPUExchange(Root, SendQ);
			// --- Compute divergence for internal cells ---
			Root->computeDivergence(0);
			// --- Runge-Kutta step for internal cells ---
			Root->RungeKutta(0);
#if defined PARMPI
			CommTimer.start();    //Communication Timer Start
			//Waiting while inter-CPU exchanges are finished
			if (MPIRecvType == 1)                     //for nonblocking recive...
				MPI_Waitall(4*Dimension,req,st);     
			CommTimer.stop();
#endif
  
			
	// --- Compute velocity gradient (Navier-Stokes only) ---
	if (EquationType==6 && Viscosity !=0.) 
	{
          	Root->computeGradient(1); //For neighbour cells
  	  	CPUExchange(Root, SendGrad);
          	Root->computeGradient(0); //for internal...
          	CommTimer.start();

#if defined PARMPI
          	if (MPIRecvType == 1) MPI_Waitall(4*Dimension,req,st);
#endif
          CommTimer.stop();
         }
  		}

		// --- Check stability ---
		Root->checkStability();

		// --- Compute integral values ---
		Root->computeIntegral();

		// --- Compute elapsed time and adapt time step ---

		if (!ComputeCPUTimeRef)
		{
			ElapsedTime += TimeStep;
			if (!ConstantTimeStep) AdaptTimeStep();
		}

		// --- Compute time-average values ---

		if (TimeAveraging)
			Root->computeTimeAverage();

}
