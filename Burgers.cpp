/***************************************************************************
                          Burgers.cpp  -  description
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

/*
______________________________________________________________________________________________

DESCRIPTION
______________________________________________________________________________________________

This function return the analytic solution of viscous Burgers equation for
the initial condition :
	
			u(x) = 1 if x < 0, u(x) = 0 elsewhere.
______________________________________________________________________________________________

*/

#include "Carmen.h"

real Burgers(real x, real t)
{
	real t2;

	t2 = t;

	if (t2 == 0.)
		return Step(x);
	else
		return .5*( 1 - tanh(0.25*Re*(x-0.5*t2)) );
}

