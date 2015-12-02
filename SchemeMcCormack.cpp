/***************************************************************************
                          SchemeMcCormack.cpp  -  description
                             -------------------
    begin                : mer fév 25 2004
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

Returns the Euler flux using the 2-4 McCormack numerical scheme.
______________________________________________________________________________________________

*/

Vector SchemeMcCormack(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, int AxisNo)
{
	// --- Local variables ---

	Vector Result(QuantityNb);

	real   	rho1, rho2;  			// Density
	Vector 	V1(Dimension), V2(Dimension);	// Velocity
	real	p1, p2;				// Pressure
	real   	rhoE1, rhoE2; 			// Energy
	real	rhoY1=0.;			// Partial mass
	real	rhoY2=0.;			// Partial mass
	Vector 	F1(QuantityNb), F2(QuantityNb);	// Fluxes
	int  	i;
	int    	it, Direction;
	real	A2 = 7./6.;
	real	A1 = -1./6.;
	
	// --- Chose direction ---
	
      	it = IterationNo%(1<<Dimension);

	if (AxisNo==1)
		Direction = it%2;
	else if (AxisNo==2)
		Direction = (it/2)%2;
	else
		Direction = (it/2)/2;

	Direction = (Direction + (StepNo-1))%2;
	
	// --- Get physical quantities ---
	
	if (Direction == 1)
	{
      		rho1 = Cell4.density();
      		rho2 = Cell3.density();

		V1 = Cell4.velocity();
		V2 = Cell3.velocity();

		rhoE1 = Cell4.energy();
		rhoE2 = Cell3.energy();

		p1 = Cell4.pressure();
		p2 = Cell3.pressure();
		
		if (ScalarEqNb == 1)
		{
			rhoY1=Cell4.average(Dimension+3);
			rhoY2=Cell3.average(Dimension+3);
		}
	}
	else
	{
      		rho1 = Cell1.density();
      		rho2 = Cell2.density();

		V1 = Cell1.velocity();
		V2 = Cell2.velocity();

		rhoE1 = Cell1.energy();
		rhoE2 = Cell2.energy();

		p1 = Cell1.pressure();
		p2 = Cell2.pressure();
	
		if (ScalarEqNb == 1)
		{
			rhoY1=Cell1.average(Dimension+3);
			rhoY2=Cell2.average(Dimension+3);
		}
	}
	
	// --- Correct partial masses ---
	
	if (ScalarEqNb == 1)
	{
		if (rhoY1 < 0.) rhoY1 = 0.;
		if (rhoY2 < 0.) rhoY2 = 0.;
		if (rhoY1 > rho1) rhoY1 = rho1;
		if (rhoY2 > rho2) rhoY2 = rho2;
	}
	
	// --- Compute fluxes ---

	F1.setValue(1, rho1*V1.value(AxisNo));
	F2.setValue(1, rho2*V2.value(AxisNo));

	for (i=1;i<=Dimension;i++)
	{
		F1.setValue(i+1, (AxisNo==i)? F1.value(1)*V1.value(i)+ p1 : F1.value(1)*V1.value(i));
		F2.setValue(i+1, (AxisNo==i)? F2.value(1)*V2.value(i)+ p2 : F2.value(1)*V2.value(i));
      	}

 	F1.setValue(Dimension+2, (rhoE1+p1)*V1.value(AxisNo));
 	F2.setValue(Dimension+2, (rhoE2+p2)*V2.value(AxisNo));

	if (ScalarEqNb == 1)
	{
		F1.setValue(Dimension+3,rhoY1*V1.value(AxisNo));
		F2.setValue(Dimension+3,rhoY2*V2.value(AxisNo)); 
	}
		
	// --- Compute values at the interface (McCormack 2,4) ---

	Result  = A2*F2 + A1*F1;

	return Result;

}

/*
______________________________________________________________________________________________

Old version
______________________________________________________________________________________________


Vector SchemeMcCormack(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, int AxisNo)
{
	// --- Local variables ---

	Vector Result(QuantityNb);

	real   rho1, rho2, rho3, rho4;  																				// Density
	Vector V1(Dimension), V2(Dimension), V3(Dimension), V4(Dimension);			// Velocity
	real	 p1, p2, p3, p4;																									// Pressure
	real   rhoE1, rhoE2, rhoE3, rhoE4; 																			// Energy
	Vector F1(QuantityNb), F2(QuantityNb), F3(QuantityNb), F4(QuantityNb);	// Fluxes
	int  	 i;
	int    it, Direction;
	
	// --- Get physical quantities ---

	// density

      rho1 = Cell1.density();
      rho2 = Cell2.density();
      rho3 = Cell3.density();
      rho4 = Cell4.density();

	// velocity

	V1 = Cell1.velocity();
	V2 = Cell2.velocity();
	V3 = Cell3.velocity();
	V4 = Cell4.velocity();

	// energy

	rhoE1 = Cell1.energy();
	rhoE2 = Cell2.energy();
	rhoE3 = Cell3.energy();
	rhoE4 = Cell4.energy();

	// pressure

	p1 = Cell1.pressure();
	p2 = Cell2.pressure();
	p3 = Cell3.pressure();
	p4 = Cell4.pressure();

	// --- Compute fluxes ---

	F1.setValue(1, rho1*V1.value(AxisNo));
	F2.setValue(1, rho2*V2.value(AxisNo));
	F3.setValue(1, rho3*V3.value(AxisNo));
	F4.setValue(1, rho4*V4.value(AxisNo));

			for (i=1;i<=Dimension;i++)
			{
				F1.setValue(i+1, (AxisNo==i)? rho1*V1.value(AxisNo)*V1.value(i)+p1 : rho1*V1.value(AxisNo)*V1.value(i));
				F2.setValue(i+1, (AxisNo==i)? rho2*V2.value(AxisNo)*V2.value(i)+p2 : rho2*V2.value(AxisNo)*V2.value(i));
				F3.setValue(i+1, (AxisNo==i)? rho3*V3.value(AxisNo)*V3.value(i)+p3 : rho3*V3.value(AxisNo)*V3.value(i));
				F4.setValue(i+1, (AxisNo==i)? rho4*V4.value(AxisNo)*V4.value(i)+p4 : rho4*V4.value(AxisNo)*V4.value(i));
      }

 			F1.setValue(Dimension+2, (rhoE1+p1)*V1.value(AxisNo));
 			F2.setValue(Dimension+2, (rhoE2+p2)*V2.value(AxisNo));
 			F3.setValue(Dimension+2, (rhoE3+p3)*V3.value(AxisNo));
 			F4.setValue(Dimension+2, (rhoE4+p4)*V4.value(AxisNo));

			// --- Compute values at the interface (McCormack 2,4) ---

      it = IterationNo%(1<<Dimension);

			if (AxisNo==1)
				Direction = it%2;
			else if (AxisNo==2)
				Direction = (it/2)%2;
			else
				Direction = (it/2)/2;

			Direction = (Direction + (StepNo-1))%2;

			Result  = (Direction==1)? (7./6.*F3 - 1./6.*F4) : (7./6.*F2 - 1./6.*F1);

			return Result;

}
*/

