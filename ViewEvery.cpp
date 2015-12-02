/***************************************************************************
                          ViewEvery.cpp  -  description
                             -------------------
    begin                : Fri Jun 14 2002
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
__________________________________________________________________________

Print solution every PrintEvery iteration
__________________________________________________________________________

Multiresolution solver
--------------------------------------------------------------------------
*/
void ViewEvery(Node* Root, int arg)
{
	char AverageName[256];		// File name for AverageNNN.dat
	char MeshName[256];       // File name for MeshNNN.dat
	char AverageFormat[256];	// File format for AverageNNN.dat
	char MeshFormat[256];			// File format for MeshNNN.dat

	sprintf(AverageFormat, "Average%s0%ii.dat", "%", DigitNumber(IterationNb));
	sprintf(AverageName, AverageFormat, arg);				
	sprintf(MeshFormat, "Mesh%s0%ii.dat", "%", DigitNumber(IterationNb));
	sprintf(MeshName, MeshFormat,arg);

	View(Root, "Tree.dat", MeshName, AverageName);
}

/*
--------------------------------------------------------------------------
Finite volume solver
--------------------------------------------------------------------------
*/

void ViewEvery(FineMesh* Root, int arg)
{
	char AverageName[256];		// File name for AverageNNN.dat
	char AverageFormat[256];	// File format for AverageNNN.dat

	sprintf(AverageFormat, "Average%s0%ii.dat", "%", DigitNumber(IterationNb));
	sprintf(AverageName, AverageFormat, arg);				

	View(Root, AverageName);
}
