/***************************************************************************
                          Timer.h  -  description
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

#ifndef TIMER_H
#define TIMER_H
#include <time.h>
#include <sys/times.h>
#include "PreProcessor.h"

/***************************************************************************
 * An object Timer gives information on the CPU time of long-time computations.
 *
 */
class Timer
{
/*
______________________________________________________________________________________________

PUBLIC METHODS
______________________________________________________________________________________________

*/
/***************************************************************************
 * Constructor. Initialize timer.
 *
 */
public:
	Timer();
								
/***************************************************************************
 * Resets time and start.
 *
 */
	void   resetStart();

/***************************************************************************
 * Adds CPU time and real time to their buffers and resets. For long
 * computations, it is recommended to do this operation at least once per
 * iteration.
 *
 */
	void	 check();

/***************************************************************************
 * Starts timer.
 *			
 */
	void start();

/***************************************************************************
 * Stop timer and, if asked, returns CPU time from previous start in seconds.
 *
 */
	double stop();

/***************************************************************************
 * Returns CPU time from previous start in seconds.
 *
 */
	double CPUTime();

/***************************************************************************
 * Returns real time from previous start in seconds.
 *
 */
	double realTime();

/***************************************************************************
 * Adds time to buffer (only when a computation is restarted).
 *
 */
	void add(double cpuTime, double realTime);


/*
______________________________________________________________________________________________

PRIVATE VARIABLES
______________________________________________________________________________________________

*/
private:

  clock_t 	StartCPUTime;  			// CPU time
  time_t 		StartRealTime;      // real time
  double 		sumCPUtime;         // elapsed CPU time in seconds
  double 		sumRealtime;    		// elapsed real time in seconds

  bool 			TimerOn;            // true if timer is running, false elswhere
};

#endif
