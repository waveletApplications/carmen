/***************************************************************************
                          DigitNumber.cpp  -  description
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

 This function returns the number of digits of an integer
______________________________________________________________________________________________

*/
int DigitNumber(int arg)
{
	int result;
	int i;

	result = 0;
	i = arg;

	while(i != 0)
	{
		i/=10;
		result++;
	}

	return result;
}
