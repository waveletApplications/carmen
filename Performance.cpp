/***************************************************************************
                          Performance.cpp  -  description
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

void Performance(const char* FileName)
{
	// --- Local variables ---

	bool EndComputation;			// True if end computation
	int FineCellNb;				// Number of cells on fine grid
	FILE 	*output;			// Pointer to output file

	double realtimefull;			//full real time
	double ftime;				// real time
	double ctime;				// CPU time
	unsigned int 	ttime, rtime; 		// total and remaining  real time (in seconds)
	unsigned int  tctime=0, rctime=0; 	// total and remaining CPU time (in seconds)
	unsigned int	rest;
	int	day, hour, min, sec;

	// --- Init EndComputation

	EndComputation = (IterationNo > IterationNb);

	// --- Compute FineCellNb ---

	FineCellNb = 1<<(ScaleNb*Dimension);	

	// --- Write in file ---

/*  
	char CPUFileName[255];
#if defined PARMPI
	sprintf(CPUFileName,"%d_%d_%d_%s",coords[0],coords[1],coords[2],FileName);
//	strcpy(CPUFileName, FileName);
#else
	strcpy(CPUFileName, FileName);
#endif
*/
    
	if ((output = fopen(FileName,"w")))
	{

		realtimefull=ftime = CPUTime.realTime();
		ctime = CPUTime.CPUTime();

		if (!EndComputation)
		{
			ttime  = (unsigned int)((ftime*IterationNb)/IterationNo);
			rtime  = (unsigned int)((ftime*(IterationNb-IterationNo))/IterationNo);
			tctime = (unsigned int)((ctime*IterationNb)/IterationNo);
			rctime = (unsigned int)((ctime*(IterationNb-IterationNo))/IterationNo);
    }

		fprintf(output, "Dimension             : %12i\n", Dimension);

		if (EndComputation)
			fprintf(output, "Iterations            : %12i\n", IterationNb);
		else
		{
			fprintf(output, "Iterations (total)    : %12i\n", IterationNb);
			fprintf(output, "Iterations (elapsed)  : %12i\n", IterationNo);
			fprintf(output, "In progress           :%13.6f %%\n", 100.*IterationNo/(1.*IterationNb));
    }

		fprintf(output, "Scales (max)          : %12i\n", ScaleNb);
		fprintf(output, "Cells (max)           : %12i\n", (1<<(ScaleNb*Dimension)));

		if (Multiresolution)
		  fprintf(output, "Solver                :           MR\n");
		else
		  fprintf(output, "Solver                :           FV\n");

		//fprintf(output, "Time integration      :     explicit\n");
		fprintf(output, "Time accuracy order   : %12i\n", StepNb);
		fprintf(output, "Time step             :%13.6e s\n", TimeStep);

		if (EndComputation)
		{
			fprintf(output, "Physical time         :%13.6e s\n", TimeStep * IterationNb);
			fprintf(output, "CPU time (s)          :%13.6e s\n", ctime);

			if (Multiresolution)
				fprintf(output, "CPU time / it. x pt   :%13.6e s\n", ctime/TotalLeafNb);
			else
      				fprintf(output, "CPU time / it. x pt   :%13.6e s\n", ctime/((1<<(ScaleNb*Dimension))*IterationNb));

			if (Multiresolution)
			{
				
				fprintf(output, "Leaf compression      :%13.6f %% \n", (100.*TotalLeafNb)/(1.0*IterationNb*FineCellNb));
				fprintf(output, "Memory compression    :%13.6f %% \n", (100.*TotalCellNb)/(1.0*IterationNb*FineCellNb));				
				fprintf(output, "CPU compression       :%13.6f %% \n", (100.*ctime)/(IterationNb*FVTimeRef));
			}
			else
			{
				fprintf(output, "Leaf compression      :%13.6f %% \n", 100.);
				fprintf(output, "Memory compression    :%13.6f %% \n", 100.);
				fprintf(output, "CPU compression       :%13.6f %% \n", 100.);
			}
  	}
		else
		{
			fprintf(output, "Total physical time   :%13.6e s\n", TimeStep * IterationNb);
			fprintf(output, "Elapsed physical time :%13.6e s\n", TimeStep * IterationNo);
			if (Multiresolution)
			{
				fprintf(output, "Leaf compression      :%13.6f %% \n", (100.*TotalLeafNb)/(1.0*IterationNb*FineCellNb));
				fprintf(output, "Memory compression    :%13.6f %% \n", (100.*TotalCellNb)/(1.0*IterationNo*FineCellNb));
				fprintf(output, "CPU compression       :%13.6f %% \n", (100.*ctime)/(IterationNo*FVTimeRef));
			}
			else
			{
				fprintf(output, "Leaf compression      :%13.6f %% \n", 100.);
				fprintf(output, "Memory compression    :%13.6f %% \n", 100.);
				fprintf(output, "CPU compression       :%13.6f %% \n", 100.);
			}
		}

		fprintf(output, "\n");

		if (EndComputation)
		{
			// --- Print final time ------------------------------------------------

			rest =  (unsigned int)(ctime);
  		day  =  rest/86400;
			rest %= 86400;
			hour =  rest/3600;
			rest %= 3600;
			min	 =  rest/60;
			rest %= 60;
			sec  =  rest;
			rest =  (unsigned int)(ctime);

			if (rest >= 86400)
				fprintf(output, "CPU time : %5d day %2d h %2d min %2d s\n", day, hour, min, sec);

			if ((rest < 86400)&&(rest >= 3600))
				fprintf(output, "CPU time : %2d h %2d min %2d s\n", hour, min, sec);

			if ((rest < 3600)&&(rest >= 60))
				fprintf(output, "CPU time : %2d min %2d s\n", min, sec);

			if (rest < 60)
				fprintf(output, "CPU time : %2d s\n", sec);
		}
		else
		{
			// --- Print total time ------------------------------------------------

			rest =  tctime;
  		day  =  rest/86400;
			rest %= 86400;
			hour =  rest/3600;
			rest %= 3600;
			min	 =  rest/60;
			rest %= 60;
			sec  =  rest;

			if (tctime >= 86400)
				fprintf(output, "Total CPU time     (estimation) : %5d day %2d h %2d min %2d s\n", day, hour, min, sec);

			if ((tctime < 86400)&&(tctime >= 3600))
				fprintf(output, "Total CPU time     (estimation) : %2d h %2d min %2d s\n", hour, min, sec);

			if ((tctime < 3600)&&(tctime >= 60))
				fprintf(output, "Total CPU time     (estimation) : %2d min %2d s\n", min, sec);

			if (tctime < 60)
				fprintf(output,"Total CPU time     (estimation) : %2d s\n", sec);

			// --- Print remaining time ------------------------------------------------

			rest =  rctime;
  		day  =  rest/86400;
			rest %= 86400;
			hour =  rest/3600;
			rest %= 3600;
			min	 =  rest/60;
			rest %= 60;
			sec  =  rest;

			if (rctime >= 86400)
				fprintf(output, "Remaining CPU time (estimation) : %5d day %2d h %2d min %2d s\n", day, hour, min, sec);

			if ((rctime < 86400)&&(rctime >= 3600))
				fprintf(output, "Remaining CPU time (estimation) : %2d h %2d min %2d s\n", hour, min, sec);

			if ((rctime < 3600)&&(rctime >= 60))
				fprintf(output, "Remaining CPU time (estimation) : %2d min %2d s\n", min, sec);

			if (rctime < 60)
				fprintf(output, "Remaining CPU time (estimation) : %2d s\n", sec);

		}


#if defined PARMPI
    fprintf(output,  "\n");
    fprintf(output,  "Real time (time() function) :%lf\n",  realtimefull);
		fprintf(output,  "clock() function :%lf\n",  ctime);
		fprintf(output,  "\nCommunications real timer: %lf\n",  CommTimer.realTime());
		fprintf(output,  "Communications clock():%lf\n",  CommTimer.CPUTime());
#endif
    
		fprintf(output, "\n");
		fclose (output);	
	}
	else
	{
		cout << "Performance.cpp: In method `void Performance(Node*, char*)':\n";
		cout << "Performance.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [Performance.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
