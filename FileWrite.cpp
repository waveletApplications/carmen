/***************************************************************************
                          FileWrite.cpp  -  description
                             -------------------
    begin                : mer oct 20 2004
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

/***************************************************************************
 * Writes in binary or ASCII mode the real number <i>arg</i> into the file <i>f</i>.
 * The global parameter <i>DataIsBinary</i> determines this choice.
 *
 */
int FileWrite(FILE *f, const char* format, real arg)
{
	int result;
	real x;
	
	x = arg;

	if (DataIsBinary)
		result = fwrite(&x,sizeof(real),1,f);
	else
		result = fprintf(f, format, x);

  return result;
}
