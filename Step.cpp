/***************************************************************************
                          Step.cpp  -  description
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

This function returns
	
			u(x) = 	1   if x < 0
							0   if x > 0
							1/2 if x = 0
______________________________________________________________________________________________

*/


#include "Carmen.h"

real Step(real x)
{
	if (x < 0.)
		return 1.;
	else if (x > 0.)
		return 0.;
	else
		return .5;
}
