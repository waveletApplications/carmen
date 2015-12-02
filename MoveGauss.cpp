/***************************************************************************
                          MoveGauss.cpp  -  description
                             -------------------
    begin                : Wed Dec 12 2001
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

This function return the analytic solution of 2D advection-diffusion for
the initial condition :
	
			u(x) = infinity if x = 0, u(x) = 0 elsewhere.

The solution is a moving gaussian.
______________________________________________________________________________________________

*/

#include "Carmen.h"
/*
______________________________________________________________________________________________

1D case
______________________________________________________________________________________________

*/

real MoveGauss(const real x, const real t)
{
		real t0;
		real r2;

		t0 = Re/(4.*pi);
		r2 = (x-Celerity*t)*(x-Celerity*t);

 		return Re/(4.*pi*(t+t0))*exp(-0.25*Re*r2/(t+t0));
}
/*
______________________________________________________________________________________________

2D case
______________________________________________________________________________________________

*/

real MoveGauss(const real x, const real y, const real t)
{
		real t0;
		real r2;

		t0 = Re/(4.*pi);
		r2 = (x-Celerity*t)*(x-Celerity*t)+y*y;

 		return Re/(4.*pi*(t+t0))*exp(-0.25*Re*r2/(t+t0));
}
/*
______________________________________________________________________________________________

3D case
______________________________________________________________________________________________

*/

real MoveGauss(const real x, const real y, const real z, const real t)
{
		real t0;
		real r2;

		t0 = Re/(4.*pi);
		r2 = (x-Celerity*t)*(x-Celerity*t)+y*y+z*z;

 		return Re/(4.*pi*(t+t0))*exp(-0.25*Re*r2/(t+t0));
}
