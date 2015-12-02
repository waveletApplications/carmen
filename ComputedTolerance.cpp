/***************************************************************************
                          ComputedTolerance.cpp  -  description
                             -------------------
    begin                : May 13 2011 changed by Kai, Sonia, Margarete
    copyright            : (C) 2005 by Olivier Roussel
    email                : roussel@ict.uni-karlsruhe.de
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

real ComputedTolerance(const int ScaleNo)
{
//if ThresholdNorm==0 const Tolerance, else L1 Harten norm

if(ThresholdNorm) return((Tolerance/GlobalVolume)*exp(Dimension*(ScaleNo-ScaleNb+1)*log(2.)));
  else return(Tolerance);

//This on is in the old version of 1.54
real Result=0.;

	// Coefficients for the normalization of the wavelet basis 
	real a = 1.;
	real b = 0.;
	real LocalVolume = GlobalVolume*exp(-log(2.)*Dimension*ScaleNo);
	
	// Choose coefficients a and b in function of ThresholdNorm
	
	switch(ThresholdNorm)
	{
		case 0:
		default: 
			if (EquationType == 6)
			{
				a = 0.5;
				b = 1.;
			}
			else
			{
				a = 1.;
				b = 0.;
			}
			break;
			
		case 1:
			a = 1.;
			b = 0.;
			break;
		
		case 2:
			a = 0.5;
			b = 0.;
			break;
		
		case 3:
			a = 0.5;
			b = 1.;
			break;
	};

	
	if (EquationType == 6) Result = Tolerance*exp(-log(LocalVolume)*(a-b/Dimension)); 
		//Result = Tolerance*exp(-log(LocalVolume)*(a-b/Dimension))*exp(log(2)*(-Dimension*(ScaleNb))); //Meg L1 mod ThresholdNorm ==1 
		//Result = Tolerance*exp(Dimension*(ScaleNo-ScaleNb+1)*log(2.)); //Meg Harten
		//Result = Tolerance*exp((ScaleNo-ScaleNb+1)*log(2.)); //Meg Harten
	else
//		Result = Tolerance*exp(Dimension*(ScaleNo-ScaleNb)*log(2.));
		Result = Tolerance*exp(-log(LocalVolume)*(a-b/Dimension))*exp(-Dimension*(ScaleNb)*log(2.));

//cout<<"ScaleNo, Tolerance0, Tolerance :"<<ScaleNo <<" , "<< Tolerance <<",  "<< Result <<endl;

	return Result;

}
