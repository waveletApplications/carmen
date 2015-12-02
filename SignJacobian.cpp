/***************************************************************************
                          SignJacobian.cpp  -  description
                             -------------------
    begin                : mer mar 10 2004
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

int SignJacobian(const Vector Q1, const Vector Q2, const int i, const int j, const int AxisNo)
{
	  real 		rho=0.;            // density
		real 		E=0.;              // energy per unit of mass
		Vector 	v(Dimension);      // velocity
		real 		result=0.;         // result
		int 		n=0;							 // counter
		int 		sgn=0;							// sign of result

		Vector Q(QuantityNb);

		// --- Rho average method ----------------------------------------------------

		

		// ---- Init rho, v,e --------------------------------------------------------

    rho = Q.value(1);

		for (n=1; n<=Dimension; n++)
			v.setValue(n,Q.value(n+1)/rho);

    E = Q.value(Dimension+2)/rho;

 		// --- Compute term of the jacobian matrix -----------------------------------

		if (j==1)
		{
				// --- Compute df/d(rho) ---
				if (i==1)
					result = 0.;
				else if (i==(Dimension+2))
					result = -v.value(AxisNo)*(Gamma*E-(Gamma-1.)*(v*v));
				else
					result = -v.value(AxisNo)*v.value(i-1) + ((AxisNo==(i-1)) ? 0.5*(Gamma-1.)*(v*v) : 0.);

		}
		else if (j==(Dimension+2))
		{
				// --- Compute df/d(rhoE) ---
				if (i==1)
					result = 0.;
				else if (i==(Dimension+2))
					result = Gamma*v.value(AxisNo);
				else
					result = (AxisNo==(i-1)) ? Gamma-1.:0.;
		}
		else
		{
				// --- Compute dF/d(rhoV) ---
				if (i==1)
					result = (AxisNo==j-1) ? 1.:0.;
				else if (i==(Dimension+2))
					result = -(Gamma-1.)*v.value(j-1)*v.value(AxisNo) + ((AxisNo==(j-1)) ? Gamma*E - 0.5*(Gamma-1.)*(v*v): 0.);
				else
				{
						if (AxisNo==(i-1))
							result = (1.-Gamma)*v.value(j-1) + ((i==j) ? 2.*v.value(j-1) : 0.);
						else
							if (i==j)
								result = v.value(AxisNo);
							else if (AxisNo==(j-1))
								result = v.value(i-1);
							else
								result = 0.;
				}	
		}

		// --- Compute sign ------------------------------------------------------

		if (result > 0.)
			sgn = 1;
		else if (result < 0.)
			sgn = -1;
		else
			sgn = 0;
	
		return sgn;
}
