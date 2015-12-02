/***************************************************************************
                          View.cpp  -  description
                             -------------------
    begin                : Thu Jun 7 2001
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
#include <cstdlib>
#include <cstring>

using namespace std;



/*
______________________________________________________________________________________________

Vizualization for multiresolution
______________________________________________________________________________________________

*/
void View(Node* Root, const char* TreeFileName, const char* MeshFileName, const char* AverageFileName)
{
	char buf[256];
	int iaux;

	// write tree (debugging only)
	if (debug) Root->writeTree(TreeFileName);

	// Compute velocity gradient (only for Navier-Stokes equation)

	if (EquationType==6) Root->computeGradient();

	// write mesh for graphic visualisation
	if (Dimension != 1)
	{
		Root->writeHeader(MeshFileName);
		Root->writeAverage(MeshFileName);

		// Compress data (if parameter ZipData is true)
		if (ZipData)
		{
			sprintf(buf,"gzip %s",MeshFileName);
			iaux=system(buf);
  	}
  }
	else
		Root->writeMesh(MeshFileName);


	// write cell-averages in multiresolution representation (1D) or on fine grid (2-3D)
	if (Dimension != 1)
	{
		Root->writeFineGrid(AverageFileName,ScaleNb+PrintMoreScales);

		// Compress data
		if (ZipData)
		{
			sprintf(buf,"gzip %s",AverageFileName);
			iaux=system(buf);
		}
	}
	else
	{
		Root->writeHeader(AverageFileName);
		Root->writeAverage(AverageFileName);
  }
}
/*
______________________________________________________________________________________________

Vizualization for finite volume on fine grid (no multiresolution)
______________________________________________________________________________________________

*/
void View(FineMesh* Root, const char* AverageFileName)
{
	char buf[256];
	int iaux;


  char CPUFileName[255];
#if defined PARMPI
	sprintf(CPUFileName,"%d_%d_%d_%s",coords[0],coords[1],coords[2],AverageFileName);
#else
	strcpy(CPUFileName, AverageFileName);
#endif

  // Compute velocity gradient (only for Navier-Stokes equation)
  
	if (EquationType==6) 
	{
    		Root->computeGradient(1);
    		Root->computeGradient(0);
   	}

	// write header for graphic visualization
	Root->writeHeader(CPUFileName);

	// write cell-average values for graphic visualization
	Root->writeAverage(CPUFileName);

  // Compress data
	if (Dimension != 1)
	{
		if (ZipData)
		{
      			sprintf(buf,"gzip %s",CPUFileName);
			iaux=system(buf);
		}
	}

		// --- Write time-average values into file ---

		if (TimeAveraging)
			Root->writeTimeAverage("TimeAverage.dat");

}
