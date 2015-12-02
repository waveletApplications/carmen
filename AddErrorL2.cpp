/***************************************************************************
                          AddErrorL2.cpp  -  description
                             -------------------
    begin                : lun fév 9 2004
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

# include "Carmen.h"

/***************************************************************************
 * Adds a L2-error <i>NewError</i> to <i>N-1</i> L2-errors stored in <i>TotalError</i>.
 *
 */

real AddErrorL2(int N, real TotalError, real NewError)
{
	return sqrt( ( (N-1)*TotalError*TotalError + NewError*NewError ) / N );
}
