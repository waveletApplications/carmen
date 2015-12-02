/***************************************************************************
                          AnalyticAverage.cpp  -  description
                             -------------------
    begin                : Tue Nov 13 2001
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
/*
______________________________________________________________________________________________

DESCRIPTION
______________________________________________________________________________________________

This function return the cell-average value of the analytic solution of the advection-diffusion
or viscous Burgers equation.
______________________________________________________________________________________________

1D version
______________________________________________________________________________________________
*/

#include "Carmen.h"

real AnalyticAverage(real x, real dx, real t)
{
	// --- Local variables


	if (EquationType == 1 && Dimension == 1)
//			return InitAverage(x, 0., 0.).value(1);
			return .5*(Advection(x-0.5*dx,t) + Advection(x+0.5*dx,t));
                        //return Advection(x,t);

	else if (EquationType == 2 && Dimension == 1)
			return .5*( Burgers(x-0.5*dx,t) + Burgers(x+0.5*dx,t) );
//			return Burgers(x,t);

	else
	{
		cout << "AnalyticAverage.cpp: In method `void AnalyticAverage(real, real, real)´:\n";
		cout << "AnalyticAverage.cpp: equation type and/or dimension mismatch \n";
		cout << "carmen: *** [AnalyticAverage.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
/*
______________________________________________________________________________________________

2D version
______________________________________________________________________________________________

*/

real AnalyticAverage(real x, real dx, real y, real dy, real t)
{
	// --- Local variables

	real result = 0.;
	int i,j;

	if (EquationType == 1 && Dimension == 2)
	{
		for (i=-1;i<=1;i+=2)
		for (j=-1;j<=1;j+=2)
			result += InitAverage(x+0.5*i*dx, y+0.5*j*dy, 0.).value(1);

		result *= 0.25;
		return result;
	}
	else if (EquationType == 6 && Dimension == 2)
  {
		for (i=-1;i<=1;i+=2)
		for (j=-1;j<=1;j+=2)
		{
			real y2 = y+0.5*j*dy;
			result += (y2 <= -1. || y2 >= 1.) ? 0.:1.-y2*y2;
    }
		result *=0.25;
		return result;
	}
	else
	{
		cout << "AnalyticAverage.cpp: In method `void AnalyticAverage(real, real, real, real, real)´:\n";
		cout << "AnalyticAverage.cpp: equation type and/or dimension mismatch \n";
		cout << "carmen: *** [AnalyticAverage.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}

/*
______________________________________________________________________________________________

3D version
______________________________________________________________________________________________

*/

real AnalyticAverage(real x, real dx, real y, real dy, real z, real dz, real t)
{
	// --- Local variables

	real result = 0.;
	int i,j,k;

	if (EquationType == 1 && Dimension == 3)
	{
		for (i=-1;i<=1;i+=2)
		for (j=-1;j<=1;j+=2)
		for (k=-1;k<=1;k+=2)
			result += MoveGauss(x+0.5*i*dx, y+0.5*j*dy, z+0.5*k*dz, t);

		result *= 0.125;
		return result;
	}
	else if (EquationType == 6 && Dimension == 3)
	{
		for (i=-1;i<=1;i+=2)
		for (j=-1;j<=1;j+=2)
		for (k=-1;k<=1;k+=2)
    {
			real z2 = z+0.5*k*dz;
			result +=(z2 <= -1. || z2 >= 1.) ? 0.:1.-z2*z2;
    }

		result += 0.125;
		return result;
  }
	else
	{
		cout << "AnalyticAverage.cpp: In method `void AnalyticAverage(real, real, real, real, real)´:\n";
		cout << "AnalyticAverage.cpp: equation type and/or dimension mismatch \n";
		cout << "carmen: *** [AnalyticAverage.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
