/***************************************************************************
                          Advection.cpp  -  description
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

This function return the analytic solution of the advection-diffusion equation for
the initial condition :

			u(x) = 1 if x < 0, u(x) = 0 elsewhere

______________________________________________________________________________________________

*/
#include "Carmen.h"

real Advection(real x, real t)
{
	real t0 = 0.1;

	return 0.5*erfc(0.5*(x-(t+t0))*sqrt(Re/(t+t0)));
}
