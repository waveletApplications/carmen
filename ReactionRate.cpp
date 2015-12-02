/***************************************************************************
                          ReactionRate.cpp  -  description
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

real ReactionRate(const real T, const real Y)
{
	real result;

	result = Ze*Ze/(2*Le)*Y*exp( Ze*(T-1) / (1+Alpha*(T-1)) );
	return (result < 0.)? 0.:result;
}
