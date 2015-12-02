/***************************************************************************
                          Smagorinsky.cpp  -  description
                             -------------------
    begin                : mar avr 19 2005
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

real Smagorinsky(const Cell& UserCell)
{
	// The Smagorinsky model corresponds to the description given in Vreman et al, 1995.
	
	real 		Result=0.;

	Matrix 	S(Dimension,Dimension);   	// Stress tensor
	real		rho=0.;			// Dimensionless density
	real		dx=1.;			// Cell size
	real		s=0.;			// Characteristic stress
	real		divV=0.;		// Divergence of the velocity

	real 		TwoThird=2./3.;
	
	int i, j;

	// Set to zero if it is not in the smallest scale

	if (UserCell.size(1) > (XMax[1]-XMin[1])/(1<<ScaleNb) ) return 0.;

	// --- Get density and dx ---

	rho = UserCell.density();

	// dx = (Dx*Dy*Dz)^(1/3)
	
	for (i=1; i<=Dimension;i++)
		dx *= UserCell.size(i);
	
	dx = exp(1./Dimension*log(dx));

	// --- Compute divV ---

	for (i=1; i<=Dimension;i++)
		divV += UserCell.gradient(i,i+1);

	// --- Compute S ---

	for (i=1; i<=Dimension;i++)
	for (j=1; j<=Dimension;j++)
		S.setValue(i, j, UserCell.gradient(i,j+1)+UserCell.gradient(j,i+1)-((i==j)? TwoThird*divV: 0.));

	// --- Compute s = sqrt( 2*S(i,j)*S(i,j) ) ---

	for (i=1; i<=Dimension;i++)
	for (j=1; j<=Dimension;j++)
		s+= 0.5*S.value(i,j)*S.value(i,j);

	s = sqrt(s);

	// --- Compute eddy-viscosity ---
 	
	Result = Re*rho*ModelConstant*ModelConstant*dx*dx*s;

	return Result;	
}
