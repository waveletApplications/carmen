/***************************************************************************
                          NormMaxQuantities.cpp  -  description
                             -------------------
    begin                : ven jan 23 2004
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
____________________________________________________________________________

Compute the Linf norm of a vector containing the physical quantities
divided by their characteristic value
____________________________________________________________________________

*/

real NormMaxQuantities(const Vector& V)
{
  int AxisNo=1;
	Vector W(4);
	real MomentumMax=0.;
	real MomentumL2=0.;

	// --- For Navier-Stokes only, divide components by their maximal values ----

  if (EquationType == 6)
	{
		if (CVS)
    		{
			for (AxisNo = 1; AxisNo <= Dimension; AxisNo++)
				MomentumL2 += power2(V.value(AxisNo+1));

			MomentumL2 = sqrt(MomentumL2);
			return MomentumL2;  
		}
		else
		{
			// Density
    			W.setValue(1, V.value(1)/QuantityMax.value(1));
		
    			// Momentum : take the maximum of the momentum in the directions
			// Take the norm of the vector of details
			
			for (AxisNo = 1; AxisNo <= Dimension; AxisNo++)
			{
				MomentumMax = Max( MomentumMax, QuantityMax.value(AxisNo+1) );
				W.setValue(2, W.value(2) + V.value(AxisNo+1)*V.value(AxisNo+1));
			}
			
			W.setValue(2, sqrt(W.value(2))/MomentumMax );	
				
			// Energy
			W.setValue(3 , V.value(Dimension+2)/QuantityMax.value(Dimension+2) );
			
			// Scalar (maximum of Y is 1 => maximum of rhoY is maximum of rho)
			if (ScalarEqNb == 1)
			{
				W.setValue(4, V.value(Dimension+3)/QuantityMax.value(1));
			}
			
		}
	}
	else
		W = V;

  // --- Compute Linf norm --

	return NMax(W);
}
