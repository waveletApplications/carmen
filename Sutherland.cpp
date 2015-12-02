/***************************************************************************
                          Sutherland.cpp  -  description
                             -------------------
    begin                : Wed Apr 30 2003
    copyright            : (C) 2003 by Olivier Roussel and Alexei Tsigulin
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

/*
_____________________________________________________________________________

This function computes the dimensionless viscosity following Sutherland's law
_____________________________________________________________________________

*/
#include "Carmen.h"

real Sutherland(real T)
{
	return exp(1.5*log(T))*(1+110/TRef)/(T+110/TRef);
}

