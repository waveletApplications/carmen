/***************************************************************************
                          InitAverage.cpp  -  description
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

Vector InitAverage(real x, real y, real z)
{
	// --- Local variables ---

	Vector Result(QuantityNb);
	real	 *Q;
  Q = new real [QuantityNb+1];
	int		 n;

	// --- Init Q ---

	for (n = 1; n <= QuantityNb; n++)
		Q[n]=0.;

	// --- Use definition of initial Q contained in file 'initial' ---

	#include "carmen.ini"

	// --- Fill vector Result and return it ---

	for (n = 1; n <= QuantityNb; n++)
		Result.setValue(n, Q[n]);

  delete[] Q;
    
	return Result;
}
