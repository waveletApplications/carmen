/***************************************************************************
                          Timer.cpp  -  description
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

#include "Timer.h"

/*
______________________________________________________________________________________________

	Constructor
______________________________________________________________________________________________

*/
Timer::Timer()
{
  StartCPUTime = clock();
  time(&StartRealTime);

  sumCPUtime = 0.;
  sumRealtime = 0.;
  TimerOn = true;
}
/*
______________________________________________________________________________________________

	Reset time and start
______________________________________________________________________________________________

*/
void Timer::resetStart()
{
  StartCPUTime = clock();
  time(&StartRealTime);

  sumCPUtime = 0.;
  sumRealtime = 0.;
  TimerOn = true;
}
/*
______________________________________________________________________________________________

	Start and, eventually, return CPU time from previous start in seconds
______________________________________________________________________________________________

*/
void Timer::start()
{
  StartCPUTime = clock();
  time(&StartRealTime);

  TimerOn = true;
  return;
}
/*
______________________________________________________________________________________________

	Stop and, eventually, return CPU time from previous start in seconds
______________________________________________________________________________________________

*/
double Timer::stop()
{
	// --- Local variables ---

  clock_t EndCPUTime;        // end CPU time
  time_t EndRealTime;        // end real time
	

	// --- Execution ---

  EndCPUTime = clock();
  time(&EndRealTime);

  sumCPUtime  += double((unsigned long int)EndCPUTime-StartCPUTime)/ (unsigned long int)CLOCKS_PER_SEC;
  sumRealtime += difftime(EndRealTime,StartRealTime);

  TimerOn = false;
  return sumCPUtime;
}
/*
______________________________________________________________________________________________

	Check
______________________________________________________________________________________________

*/
void Timer::check()
{
	// --- Local variables ---

  clock_t EndCPUTime;        // end CPU time	

	// --- Execution ---

  EndCPUTime   =  clock();
  sumCPUtime   += double((unsigned long int)EndCPUTime-StartCPUTime)/ (unsigned long int)CLOCKS_PER_SEC;
  StartCPUTime =  EndCPUTime;

}

/*
______________________________________________________________________________________________

	Return CPU time from previous start in seconds
______________________________________________________________________________________________

*/
double Timer::CPUTime()
{
	// --- Local variables ---
  clock_t EndCPUTime;

 	// --- Execution ---
  if (TimerOn)
  {
		// If timer is running, compute and return time in seconds
    EndCPUTime = clock();
    return sumCPUtime+double((unsigned long int)EndCPUTime-StartCPUTime)/ (unsigned long int)CLOCKS_PER_SEC;
  }
  else
		// else return time previously computed between last start and stop procedures
    return sumCPUtime;
}

/*
______________________________________________________________________________________________

	Return real time from previous start in seconds	
______________________________________________________________________________________________

*/
double Timer::realTime()
{
	// --- Local variables ---

  time_t EndRealTime;

 	// --- Execution ---

  if (TimerOn)
  {
		// If timer is running, compute and return time in seconds
    time(&EndRealTime);
    return sumRealtime+difftime(EndRealTime,StartRealTime);
  }
  else
		// else return time previously computed between last start and stop procedures
    return sumRealtime;
}
/*
______________________________________________________________________________________________

	Add time (Only when computation is restarted)
______________________________________________________________________________________________

*/
void Timer::add(double cpuTime, double realTime)
{
		sumCPUtime=cpuTime;
		sumRealtime=realTime;
}
