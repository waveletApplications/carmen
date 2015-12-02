/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Thu Jun 7 2001, May 18 , 2011 Margarete Domingues InitTimeStep before PrintIntegral in the initialization 
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Carmen.h"

/*
______________________________________________________________________________________________

Main
______________________________________________________________________________________________

*/
int main(int argc, char *argv[])
{
	// --- Init Cluster variable ---------------------------------------------------------

	// carmen 0 => local execution i.e. show time on screen
	// carmen 1 => cluster execution i.e. refresh Performance.dat (default)

#if defined PARMPI
	// --- MPI Runtime system initialization
	// size - total number of processors
	// rnak - current CPU
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  	size=1;
	rank=0;
#endif
  
	if (argc == 2)
		Cluster = atoi(argv[1]);

	// --- Print messages on screen- -----------------------------------------------------

	cout << "carmen: begin execution.\n\n";
	printf("Carmen %4.2f \n",CarmenVersion);
	cout << "Copyright (C) 2000-2005 by Olivier Roussel.\n";
	cout << "All rights reserved.\n\n";
	
#if defined PARMPI
	//Synchronize all parallel branches
  	MPI_Barrier(MPI_COMM_WORLD);
#endif

  	CPUTime.start();
  
  // --- Create first node of the tree structure ---------------------------------------

	Node 	*Mesh=0;
	FineMesh *FMesh=0;

	
	// --- Init global values (See Parameters.h and Parameters.cpp) ----------------------------

	cout << "carmen: init computation ...\n";
	InitParameters();

  // --- Debug output information for parallel execution -------------------------------------
  
#if defined PARMPI                                                                            
  	if (Multiresolution) 
  	{
		printf("\nParallel Multiresolution solver not implemented yet!\n");
		exit(0);
	}
	
  	printf("My Rank=%d\n",rank);
  
  	// --- Each CPU print his coordinates in the virtual processor cart ------------------------
	printf("Cart_i = %d;     Cart_j = %d;     Cart_k = %d;\n",coords[0],coords[1],coords[2]);
	
  	// --- Each CPU print his computation domain
  	printf("Xmin = %lf;     XMax = %lf;\n",XMin[1],XMax[1]);
	printf("Ymin = %lf;     YMax = %lf;\n",XMin[2],XMax[2]);
	printf("Zmin = %lf;     ZMax = %lf;\n",XMin[3],XMax[3]);
	
  	// --- And the local scale number ----------------------------------------------------------
  	printf("ScaleNb = %d\n",ScaleNb);
#endif
  
	// --- Allocate ----------------------------------------------------------------------------

	if (Multiresolution)
		Mesh = new Node;
	else
		FMesh = new FineMesh;

	// --- Init tree structure -----------------------------------------------------------------

	if (Multiresolution)
	{
		InitTree(Mesh);
		RefreshTree(Mesh);
	}

	// --- Compute initial integral values and init time step ----------------------------------

	if (Multiresolution)
		Mesh->computeIntegral();
	else
		FMesh->computeIntegral();
	
	// -- Write integral values --
	
	
			
	// -- Compute initial time step --
	InitTimeStep();

	if (rank==0) PrintIntegral("Integral.dat");

	// --- Save initial values into files ------------------------------------------------------

	if (PrintEvery == 0)
	{
		if (Multiresolution)
			View(Mesh, "Tree_0.dat", "Mesh_0.dat", "Average_0.dat");
		else
			View(FMesh,"Average_0.dat");
	}
	
	// --- When PrintEvery != 0, save initial values into specific name format ---

	if (PrintEvery != 0)
	{
		if (Multiresolution)
			ViewEvery(Mesh, 0);
		else
			ViewEvery(FMesh, 0);
	}

  	// --- Parallel execution only --------------------------------------------
  	// --- Save to disk DX header for ouput files -----------------------------
  	// --- This file is needed for the external postprocessing (merging files from the different processors)

#if defined PARMPI

	real 	tempXMin[4];
	real	tempXMax[4];
	
  	// --- Save original task parameters for the parallel execution
  	int tempScaleNb=ScaleNb;

  	// --- Simulate sequantial running
  	ScaleNb=AllTaskScaleNb;

  	for (int i=0;i<4;i++) 
  	{
  		tempXMin[i]=XMin[i];
		tempXMax[i]=XMax[i];
		// --- Simulate sequantial running
		XMin[i]=AllXMin[i];
		XMax[i]=AllXMax[i];
  	}

  	// --- Write header with parameters, as we have run sequantial code
  	if (rank==0) FMesh->writeHeader("header.txt");

  	// Restore variables
  	for (int i=0;i<4;i++) 
  	{
		XMin[i]=tempXMin[i];
   		XMax[i]=tempXMax[i];
	}
	
  	ScaleNb=tempScaleNb;

#endif  
  
  	// --- Done ---

	cout << "carmen: done.\n";

	// --- Write solver type ---

	if (Multiresolution)
		cout << "carmen: multiresolution (MR) solver.\n";
	else
		cout << "carmen: finite volume (FV) solver.\n";

	// --- Write number of iterations ---

	if (IterationNb == 1)
		cout << "carmen: compute 1 iteration ...\n";
  	else
		cout << "carmen: compute " << IterationNb << " iterations ...\n";
	
	printf("\n\n\n");

	// --- Begin time iteration ----------------------------------------------------------------

	for (IterationNo = 1; IterationNo <= IterationNb; IterationNo++)
	{

		// --- Time evolution procedure ---
		if (Multiresolution)
			TimeEvolution(Mesh);
		else
			TimeEvolution(FMesh);

		// --- Remesh ---
		 if (Multiresolution) Remesh(Mesh);

		// --- Check CPU Time ---
	  	CPUTime.check();

		// --- Write information every (Refresh) iteration ---
		if ((IterationNo-1)%Refresh == 0)
		{                                          
			// - Write integral values -
			if (rank==0) PrintIntegral("Integral.dat");

			if (Cluster == 0)
				ShowTime(CPUTime);  // Show time on screen
		  //else
			if (rank==0) Performance("carmen.prf"); // Refresh file "carmen.prf"
    		}

		// --- Backup data every (10*Refresh) iteration ---
		if ((IterationNo-1)%(10*Refresh) == 0 && UseBackup)
		{
			if (Multiresolution)
				Backup(Mesh);
			else
				Backup(FMesh);
    		}

		// --- Print solution if IterationNo = PrintIt1 to PrintIt6 ---
		if (Multiresolution)
			ViewIteration(Mesh);
		else
			ViewIteration(FMesh);

		// --- Print solution if IterationNo is a multiple of PrintEvery ---

		if (PrintEvery != 0)
		{
     			if (IterationNo%PrintEvery == 0)
			{
				if (Multiresolution)
					ViewEvery(Mesh, IterationNo);
				else
					ViewEvery(FMesh, IterationNo);
			}
		}

	// --- End time iteration ------------------------------------------------------------------
	}	

	// --- Backup final data  ------------------------------------------------------------------

	IterationNo--;

	if (UseBackup)
  {
		if (Multiresolution)
			Backup(Mesh);
		else
			Backup(FMesh);
	}

	// --- Write integral values ---------------------------------------------------------------
	
	if (rank==0) PrintIntegral("Integral.dat");
	
	IterationNo++;

	// --- Save values into file ---------------------------------------------------------------

	if (Multiresolution)
		View(Mesh, "Tree.dat", "Mesh.dat", "Average.dat");
	else
		View(FMesh, "Average.dat");

	cout << "\ncarmen: done.\n";
	
	// --- Analyse performance and save it into file -------------------------------------------
				
	if (rank==0) Performance("carmen.prf");	

	// --- End ---------------------------------------------------------------------------------

	if (Multiresolution)
		delete Mesh;
	else
		delete FMesh;


#if defined PARMPI

  	//free memory for the MPI runtime variables
	delete[] disp;
	delete[] blocklen;
	int sz;
	MPI_Buffer_detach(&MPIbuffer,&sz);
	// for (int i = 0; i < 4*Dimension; i++)	MPI_Request_free(&req[i]);
	MPI_Finalize();
	
#endif
    
  	cout <<"carmen: end execution.\n";
	return EXIT_SUCCESS;
}
