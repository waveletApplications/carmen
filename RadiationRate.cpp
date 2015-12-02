/***************************************************************************
                          RadiationRate.cpp  -  description
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

real RadiationRate(const real T)
{
	real T1, T2;  // Temporary temperatures

   T1 = 1./Alpha - 1;
   T2 = T + T1;
   return Sigma * ( T2*T2*T2*T2 - T1*T1*T1*T1 );
}
