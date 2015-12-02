/***************************************************************************
                          ShowTime.cpp  -  description
                             -------------------
    begin                : Fri Jun 14 2002
    copyright            : (C) 2002 by Olivier Roussel and Alexei Tsigulin
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

void ShowTime(Timer arg)
{
//	double ftime;									// real time
	double ctime;									// CPU time
//	unsigned int 	ttime, rtime; 	// total and remaining  real time (in seconds)
	unsigned int  tctime, rctime; // total and remaining CPU time (in seconds)

	int day, hour, min, sec;
	unsigned int rest;

	// --- Write total and remaining estimated time -----------------------

//	ftime  = arg.GetRealTime();
	ctime  = arg.CPUTime();
//	ttime  = (unsigned int)((ftime*IterationNb)/IterationNo);
//	rtime  = (unsigned int)((ftime*(IterationNb-IterationNo))/IterationNo);
	tctime = (unsigned int)((ctime*IterationNb)/IterationNo);
	rctime = (unsigned int)((ctime*(IterationNb-IterationNo))/IterationNo);

	// --- Show total time ------------------------------------------------

	rest =  tctime;
  day  =  rest/86400;
	rest %= 86400;
	hour =  rest/3600;
	rest %= 3600;
	min	 =  rest/60;
	rest %= 60;
	sec  =  rest;

	printf("\033[1A\033[1A");

	if (tctime >= 86400)
		printf("Total CPU time     (estimation) : %5d day %2d h %2d min %2d s\n", day, hour, min, sec);

	if ((tctime < 86400)&&(tctime >= 3600))
		printf("Total CPU time     (estimation) : %2d h %2d min %2d s        \n", hour, min, sec);

	if ((tctime < 3600)&&(tctime >= 60))
		printf("Total CPU time     (estimation) : %2d min %2d s              \n", min, sec);

	if (tctime < 60)
		printf("Total CPU time     (estimation) : %2d s                      \n", sec);

	// --- Show remaining time ------------------------------------------------

	rest =  rctime;
  day  =  rest/86400;
	rest %= 86400;
	hour =  rest/3600;
	rest %= 3600;
	min	 =  rest/60;
	rest %= 60;
	sec  =  rest;

	if (rctime >= 86400)
		printf("Remaining CPU time (estimation) : %5d day %2d h %2d min %2d s\n", day, hour, min, sec);

	if ((rctime < 86400)&&(rctime >= 3600))
		printf("Remaining CPU time (estimation) : %2d h %2d min %2d s        \n", hour, min, sec);

	if ((rctime < 3600)&&(rctime >= 60))
		printf("Remaining CPU time (estimation) : %2d min %2d s              \n", min, sec);

	if (rctime < 60)
		printf("Remaining CPU time (estimation) : %2d s                      \n", sec);

}
