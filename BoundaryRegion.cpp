/***************************************************************************
                          BoundaryRegion.cpp  -  description
                             -------------------
    begin                : lun jan 24 2005
    copyright            : (C) 2005 by Olivier Roussel
    email                : roussel@ict.uni-karlsruhe.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/*
____________________________________________________________________________

Inits the function BoundaryRegion

User file "carmen.bc" must contain the following lines, e.g.

if (x <= XMin)
	Region = Inflow;

if (x >= XMax)
	Region = Outflow;

if (y <= YMin)
	Region = AdiabaticSolid;

if (y >= YMax)
	Region = IsothermalSolid;

When no file "carmen.bc" exists, it is created as empty.

Accepted values for Region are:

Fluid;
Inflow;
Outflow;
FreeSlipSolid;
AdiabaticSolid;
IsothermalSolid;
____________________________________________________________________________

NOTE: Only for Navier-Stokes (EquationType = 6)
____________________________________________________________________________

IMPORTANT: The boundaries are inside the computational domain.
Ex: If the computational domain in x is [1:10], then inflow sould be
[1:1.1] and outflow [9.9:10].
____________________________________________________________________________
*/

#include "Carmen.h"

int BoundaryRegion(const Vector& X)
{
	real x=0., y=0., z=0.;

  	int Fluid;
	int Inflow;
	int Outflow;
	int FreeSlipSolid;
	int AdiabaticSolid;
	int IsothermalSolid;
	
	int Region;
  	
	// Only in UseBoundaryRegions = true
	
	if (!UseBoundaryRegions) return 0;
	
	// --- Init values ---
	
	Fluid   	= 0;
	Inflow  	= 1;
	Outflow 	= 2;
	FreeSlipSolid 	= 3;
	AdiabaticSolid 	= 4;
	IsothermalSolid = 5;
	
	Region = Fluid;
	
	x = X.value(1);
  	y = (Dimension > 1)? X.value(2):0.;
	z = (Dimension > 2)? X.value(3):0.;

	#include "carmen.bc"

	if (EquationType != 6)
		Region = Fluid;

	return Region;
}
