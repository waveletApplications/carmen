/***************************************************************************
                          PreProcessor.h  -  description
                             -------------------
    begin                : Fri Oct 26 2001
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

#ifndef XYZ987 
#include "carmen.pre"
#endif
   
/*
______________________________________________________________________________________________

 DO NOT CHANGE THESE LINES
______________________________________________________________________________________________

*/

#ifdef DOUBLE_PRECISION
	#define real double
	#define FORMAT "%23.16e  "
	#define TXTFORMAT "%-23s  "
	#define REAL "double"
  #define BACKUP_FILE_FORMAT "%lf"
#else
	#define real float
	#define FORMAT "%15.8e  "
	#define TXTFORMAT "%-15s  "
	#define REAL "float"
  #define BACKUP_FILE_FORMAT "%f"
#endif

// ALWAYS USED
#define byte unsigned char

#ifdef DOUBLE_PRECISION
#define MPI_Type MPI_DOUBLE
#else
#define MPI_Type MPI_FLOAT
#endif




