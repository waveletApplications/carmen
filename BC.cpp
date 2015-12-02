/***************************************************************************
                          BC.cpp  -  description
                             -------------------
    begin                : ven nov 5 2004
    copyright            : (C) 2004 by Olivier Roussel and Alexei Tsigulin
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


#include "Carmen.h"
/*
______________________________________________________________________________________________

DESCRIPTION
______________________________________________________________________________________________

This function returns the position of i, taking into account boundary conditions.
For example: for periodic boundary conditions, N=256 elements, if i=-1, return 255.
For Neuman boundary conditions, N=256 elements, if i=-1, return 0.
______________________________________________________________________________________________
                     
*/



int BC(int i, int AxisNo, int N)
{
  int result=-999999;

	if (AxisNo > Dimension)
		return 0;

#if defined PARMPI
	if (CMin[AxisNo] == 3) result=i;    //Periodic
	else
	{
		if (i<0) if ((coords[0]==0 && AxisNo==1) || (coords[1]==0 && AxisNo==2) || (coords[2]==0 && AxisNo==3)) result=-i-1;  //Neumann
		if (i>=N)		if ((coords[0]==CartDims[0]-1 && AxisNo==1) ||	(coords[1]==CartDims[1]-1 && AxisNo==2)
										 || (coords[2]==CartDims[2]-1 && AxisNo==3)) result=2*N-i-1;
	  if (result==-999999) result=i;    //Not the boundary, simple cell from anouther CPU
	}


#else
	if (CMin[AxisNo] == 3)
		result = (i+N)%N;                        // Periodic
	else
	  result = ((i+N)/N==1)? i : (2*N-i-1)%N;  // Neumann
#endif

	return result;
}
