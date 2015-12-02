/***************************************************************************
                          SchemeAUSM  -  description
                             -------------------
    begin                : ven oct 21 2005
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
/*
______________________________________________________________________________________________

Returns the Euler flux using the AUSM+ method and the Van Albada limiter (family of MUSCL schemes).
______________________________________________________________________________________________

*/

Vector SchemeAUSM(const Cell& Cell1,const Cell& Cell2,const Cell& Cell3,const Cell& Cell4, int AxisNo)
{

	// --- Local variables ---------------------------------------------------------------

	// General variables
	
	Vector	LeftAverage(QuantityNb);	//
	Vector	RightAverage(QuantityNb);	// Conservative quantities
	Vector	Result(QuantityNb);		// Euler flux
		
	// Variables for the AUSM+ scheme
	
	Vector	VL(Dimension), VR(Dimension);	// Left and right velocities
	real	rhoL=0., rhoR=0.;		// Left and right densities
	real	eL=0., eR=0.;			// Left and right energies per unit of mass
	real	HL=0., HR=0.;			// Left and right enthalpies
		
	real	pL =0., pR =0., pC =0.;		// Left, right and centered pressures
	real	aL =0., aR =0., aC =0.;		// Left, right and centered speeds of sound
	real	MaL=0., MaR=0., MaC=0.;		// Left, right and centered Mach numbers
		
	real	M4L=0., M4R=0., P5L=0., P5R=0.; // Weight functions
	real	MassFlux=0.;			// Mass flux
								
	// Variables for the limiter
	
	real   r = 0.;				// Slope quotient
	real   Limiter = 0.;			// Slope limiter
	real   LeftSlope = 0., RightSlope = 0.; // Left and right slopes
	int i;

	bool UseLimiter = true;
			
	// --- Van Albada limiter --------------------------------------------
		
	if (!UseLimiter)
	{
		LeftAverage = Cell2.average();
		RightAverage = Cell3.average();
	}
	else
	{
	for (i=1; i<=QuantityNb; i++)
	{
		// --- Compute left cell-average value ---
			
		if (Cell2.average(i) != Cell1.average(i))
		{
			RightSlope 	= Cell3.average(i)-Cell2.average(i);
			LeftSlope	= Cell2.average(i)-Cell1.average(i);
			r		= RightSlope/LeftSlope;
			Limiter 	= (r > 0) ? (r*r+r)/(1+r*r) : 0.;
				
			LeftAverage.setValue(i, Cell2.average(i) + 0.5*Limiter*LeftSlope);
		}
		else
			LeftAverage.setValue(i, Cell2.average(i));
				
		// --- Compute right cell-average value ---
			
		if (Cell3.average(i) != Cell2.average(i))
		{
			RightSlope 	= Cell4.average(i)-Cell3.average(i);
			LeftSlope	= Cell3.average(i)-Cell2.average(i);
			r		= RightSlope/LeftSlope;
			Limiter 	= (r > 0) ? (r*r+r)/(1+r*r) : 0.;
				
			RightAverage.setValue(i, Cell3.average(i) - 0.5*Limiter*LeftSlope);
		}
		else
			RightAverage.setValue(i, Cell3.average(i));
	}
	}
		
	// --- AUSM scheme -------------------------------------------------------------
		
	// --- Extract left and right natural variables ---
		
	// Left and right densities
	rhoL = LeftAverage.value(1);
	rhoR = RightAverage.value(1);
		
	// Left and right velocities
	for (i=1;i<=Dimension;i++)
	{
		VL.setValue( i, LeftAverage.value(i+1)/rhoL );
		VR.setValue( i, RightAverage.value(i+1)/rhoR );
	} 
		
	// Left and right energies per unit of mass
	eL = LeftAverage.value(Dimension+2)/rhoL;
	eR = RightAverage.value(Dimension+2)/rhoR;
		
	// --- Compute derived variables ---
		
	// Left and right pressures
	pL = (Gamma -1.)*rhoL*( eL - 0.5*(VL*VL) );
	pR = (Gamma -1.)*rhoR*( eR - 0.5*(VR*VR) );
		
	// Left and right enthalpies per unit of mass
	HL = eL + pL/rhoL;
	HR = eR + pR/rhoR;
		
	// Left and right speeds of sound
	aL = sqrt(Gamma*pL/rhoL);
	aR = sqrt(Gamma*pR/rhoR);
		
	// Compute centered speed of sound
	//aC = 0.5*(aL+aR);
	aC = sqrt(aL*aR);
			
	// Left and right Mach numbers
	MaL = VL.value(AxisNo)/aC;
	MaR = VR.value(AxisNo)/aC;
			
	// Compute weight functions
		
	if (Abs(MaL) >= 1.)
	{
		M4L = (MaL > 0.) ? MaL:0.;
		P5L = (MaL > 0.) ? 1. :0.;
	}
	else
	{
		M4L = 0.25 * power2(MaL+1.)		+ 0.125 * power2(MaL*MaL-1.);
		P5L = 0.25 * power2(MaL+1.) * (2.-MaL) 	+ 3./16.* power2(MaL*MaL-1.) * MaL;
	}
	
	if (Abs(MaR) >= 1.)
	{
		M4R = (MaR < 0.) ? MaR:0.;
		P5R = (MaR < 0.) ? 1. :0.;
	}
	else
	{
		M4R = -0.25*power2(MaR-1.)		- 0.125 * power2(MaR*MaR-1.);  
		P5R =  0.25*power2(MaR-1.) * (2.+MaR)	- 3./16.* power2(MaR*MaR-1.) * MaR;
	}
			
	// Compute centered pressure, Mach number and mass flux
		
	pC  = pL*P5L + pR*P5R;
	MaC = M4L + M4R;
		
	MassFlux = (MaC > 0.) ? aC*rhoL*MaC : aC*rhoR*MaC;
		
	// Compute Euler flux
		
	Result.setValue(1, MassFlux);
		
	for (i=1; i<=Dimension; i++)
	{
		Result.setValue(i+1,	0.5*MassFlux*     (VR.value(i)+VL.value(i)) 
				- 	0.5*Abs(MassFlux)*(VR.value(i)-VL.value(i))
				+ 	((AxisNo == i)? pC : 0.));
	}
		
	Result.setValue(Dimension+2, 0.5*MassFlux*(HR+HL) - 0.5*Abs(MassFlux)*(HR-HL));
		
	// --- Return Euler flux ---------------------------------------------------------------
	
	return Result;
}
