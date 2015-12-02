/***************************************************************************
                          CurlVelocity.cpp  -  description
                             -------------------
    begin                : Fri May 31 2002
    copyright            : (C) 2002 by Olivier Roussel and Alexei Tsigulin
    email                : roussel@ict.uni-karlsruhe.de ; lpsoft@mail.ru
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
______________________________________________________________________________________________

DESCRIPTION

 This function returns the function velocity (EquationType = 5 only)
______________________________________________________________________________________________

*/
real CurlVelocity(real x, real y, real t, int AxisNo)
{
	real result=0;

	real r=0;           		// radius
	real v=0;								// velocity
	real x1=0, y1=0, t1=0; 	// to use in case of change of variables

	real x0 = 0.;
	real y0 = 0.;
	real t0 = 0.1;

	x1 = x - x0;
	y1 = y - y0;
	t1 = t + t0;

	r = sqrt(x1*x1+y1*y1);

	v = (r==0.)?0.:Circulation/(2*pi*r)*(1-exp(-Re*r*r/(4*t1)));
	
  switch(AxisNo)
	{
		case 1:
			result = (r==0)?0:(-y1*v/r);
			break;

    case 2:
			result = (r==0)?0:(x1*v/r);
			break;

		default:
			result = 0;
	};

	return result;	
}

