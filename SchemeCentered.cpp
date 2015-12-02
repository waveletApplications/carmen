/***************************************************************************
                          SchemeCentered  -  description
                             -------------------
    begin                : Qui Jul 28 2005
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

Returns the Euler flux using the 4th order centered scheme.
______________________________________________________________________________________________

*/

Vector SchemeCentered(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, int AxisNo)
{
	// --- Local variables ---

	Vector Result(QuantityNb);

	real   rho1, rho2, rho3, rho4, rho;  																				// Density
	Vector V1, V2, V3, V4, V;																									// Velocity
	real	 p1, p2, p3, p4, p;	
	real	rhoE;																								// Pressure
	int  	 i;
	
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

	// pressure

	p1 = Cell1.pressure();
	p2 = Cell2.pressure();
	p3 = Cell3.pressure();
	p4 = Cell4.pressure();

	// Scheme centered on the natural variables
	
	rho = 1./12.*(-rho1 + 7.*rho2 + 7.*rho3 -rho4);
	V   = 1./12.*(-V1   + 7.*V2   + 7.*V3   -V4  );
	p   = 1./12.*(-p1   + 7.*p2   + 7.*p3   -p4  );

	// Compute energy per unit of volume

	rhoE = p/(Gamma-1) + 0.5*rho*(V*V);	
	
	// Compute scheme
     
	Result.setValue(1, rho*V.value(AxisNo));

	for (i=1;i<=Dimension;i++)
		Result.setValue(i+1, (AxisNo==i)? rho*V.value(AxisNo)*V.value(i)+p : rho*V.value(AxisNo)*V.value(i));

 	Result.setValue(Dimension+2, (rhoE+p)*V.value(AxisNo));
 

	return Result;

}
