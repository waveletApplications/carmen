/***************************************************************************
                          Flux.cpp  -  description
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

Vector Flux(Cell& Cell1, Cell& Cell2, Cell& Cell3, Cell& Cell4, int AxisNo)
{
	// --- Local variables ---

	Vector Result(QuantityNb);

	real 	u1,u2, u3, u4, u; 			// Velocity in cells
	real	T1,T2;									// Temperatures in cells
	real	Y1,Y2;									// Concentrations in cells
	real	dx;											// Cell size in direction AxisNo		
	real  x,y,t;	

	int it, Direction;
	
	// --- Get cell size ---

	dx = .5*(Cell2.size(AxisNo)+Cell3.size(AxisNo));

	// --- Compute flux depending on EquationType ---

	switch(EquationType)
	{
		// --- TYPE 1 : LINEAR ADVECTION ---

    case 1:
			u1 = Cell1.average(1);
			u2 = Cell2.average(1);
			u3 = Cell3.average(1);
			u4 = Cell4.average(1);

			// 2nd order ENO scheme for advective terms

			if (AxisNo == 1 && Celerity != 0.)
			{
				if (SchemeNb == 3) // ENO 2
				{
					if (Celerity == 1.)
						u = u2 + 0.5*MinAbs(u3-u2, u2-u1);
					else if (Celerity == -1.)
						u = u3 + 0.5*MinAbs(u4-u3, u3-u2);
					else
						u = 0.; // Error -> Celerity can be 1 or -1 only.
				}
					
				else if (SchemeNb == 2) // Centered 4
					u = (- u1 + 7.*u2 + 7.*u3 - u4)/12.;

				else if (SchemeNb == 20) // Centered 2
					u = 0.5*(u2 + u3);

				else 
				{
      					it = IterationNo%(1<<Dimension);

					if (AxisNo==1)
						Direction = it%2;
					else if (AxisNo==2)
						Direction = (it/2)%2;
					else
						Direction = (it/2)/2;

					Direction = (Direction + (StepNo-1))%2;
					u = (Direction == 1)? (7./6.*u3    - 1./6.*u4   ) :(7./6.*u2    - 1./6.*u1   );  // Mc Cormack (2-4)
					//u = (Direction == 1)? u3:u2;  // Mc Cormack (2-2)
				}
			}
			else
					u = 0.;
		 	
			Result.setValue(1, u - Viscosity*(u3 - u2)/(Re*dx) );

			break;

		// --- TYPE 2 : BURGERS ---

    case 2:
			u1 = Cell1.average(1);
			u2 = Cell2.average(1);
			u3 = Cell3.average(1);
			u4 = Cell4.average(1);

			u = .5*(u2+u3);

			// 2nd order ENO scheme for advective terms

			if (u > 0)
				u = u2 + 0.5*MinAbs(u3-u2, u2-u1);
			else if (u < 0)
				u = u3 + 0.5*MinAbs(u4-u3, u3-u2);

			Result.setValue(1, .5*u*u - Viscosity*(u3 - u2)/(Re*dx) );
			break;

		// --- TYPE 3 : FLAME FRONT ---

		case 3:

			// --- Get temperatures and concentrations ---

  		T1 = Cell2.temperature();
  		Y1 = Cell2.concentration();
  		T2 = Cell3.temperature();
  		Y2 = Cell3.concentration();

			// --- Compute advective and diffusive flux ---

			if (AxisNo == 1)
			{
				Result.setValue(1, .5*FlameVelocity*(T1+T2) -         ( (T2 - T1)/dx ));
				Result.setValue(2, .5*FlameVelocity*(Y1+Y2) - 1./Le * ( (Y2 - Y1)/dx ));
			}
			else
			{
				Result.setValue(1, -         ( (T2 - T1)/dx ));
				Result.setValue(2, - 1./Le * ( (Y2 - Y1)/dx ));
			}
			break;

		// --- TYPE 4 : FLAME BALL ---

		case 4:

			// --- Get temperatures and concentrations ---

  		T1 = Cell2.temperature();
  		Y1 = Cell2.concentration();
  		T2 = Cell3.temperature();
  		Y2 = Cell3.concentration();

			// --- Compute diffusive flux ---

			Result.setValue(1, -         ( (T2 - T1)/dx ));
			Result.setValue(2, - 1./Le * ( (Y2 - Y1)/dx ));

			break;

		// --- TYPE 5 : INTERACTION FLAME - CURL ---

		case 5:

			// --- Get temperatures and concentrations ---

  		T1 = Cell2.temperature();
  		Y1 = Cell2.concentration();
  		T2 = Cell3.temperature();
  		Y2 = Cell3.concentration();
			 x = .5*(Cell2.center(1)+Cell3.center(1));
			 y = .5*(Cell2.center(2)+Cell3.center(2));
			 t = IterationNo*TimeStep;

			// --- Compute advective and diffusive flux ---

				Result.setValue(1, .5*CurlVelocity(x,y,t,AxisNo)*(T1+T2) -         ( (T2 - T1)/dx ));
				Result.setValue(2, .5*CurlVelocity(x,y,t,AxisNo)*(Y1+Y2) - 1./Le * ( (Y2 - Y1)/dx ));
			break;

		// --- TYPE 6 : NAVIER-STOKES ---

		case 6:

			Cell C1, C2, C3, C4;

	    		int BoundaryCell1 = BoundaryRegion(Cell1.center());
  			int BoundaryCell2 = BoundaryRegion(Cell2.center());
  			int BoundaryCell3 = BoundaryRegion(Cell3.center());
  			int BoundaryCell4 = BoundaryRegion(Cell4.center());

			bool UseBoundaryCells = (UseBoundaryRegions && (BoundaryCell1!=0 || BoundaryCell2!=0 || BoundaryCell3!=0 || BoundaryCell4!=0));

			// --- Take into account boundary conditions ---

			if (UseBoundaryCells)
				GetBoundaryCells(Cell1, Cell2, Cell3, Cell4, C1, C2, C3, C4, AxisNo);

			// --- Compute inviscid fluxes ---

			switch(SchemeNb)
			{
				case 1:
				default:
	        			if (UseBoundaryCells)
						Result = SchemeMcCormack(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeMcCormack(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;

				case 2:
	        			if (UseBoundaryCells)
						Result = SchemeCentered(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeCentered(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;

				case 3:
	        			if (UseBoundaryCells)
						Result = SchemeAUSM(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeAUSM(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;
				case 31:
	        			if (UseBoundaryCells)
						Result = SchemeAUSMDV(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeAUSMDV(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;
				case 4:
	        			if (UseBoundaryCells)
						Result = SchemeENO2(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeENO2(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;
					
				case 5:
	        			if (UseBoundaryCells)
						Result = SchemeENO3(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeENO3(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;
					
				case 6:
	        			if (UseBoundaryCells)
						Result = SchemeOsmp3Mc(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeOsmp3Mc(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;
					
				case 7:
	        			if (UseBoundaryCells)
						Result = SchemeOsmp5Mc(C1, C2, C3, C4, AxisNo);
					else
						Result = SchemeOsmp5Mc(Cell1, Cell2, Cell3, Cell4, AxisNo);
					break;
					
				case 8:
	        			if (Cell1.size(AxisNo) > (XMax[AxisNo]-XMin[AxisNo])/(1<<ScaleNb))
					{
						if (UseBoundaryCells)
							Result = SchemeMcCormack(C1, C2, C3, C4, AxisNo);
						else
							Result = SchemeMcCormack(Cell1, Cell2, Cell3, Cell4, AxisNo);
					}
					else
					{
						if (UseBoundaryCells)
							Result = SchemeOsmp3Mc(C1, C2, C3, C4, AxisNo);
						else
							Result = SchemeOsmp3Mc(Cell1, Cell2, Cell3, Cell4, AxisNo);
					}	
					break;				
				};

			// --- Add viscous fluxes ---

			if (Viscosity != 0.)
			{
				if (UseBoundaryCells)
					Result += ViscousFlux(C2,C3,AxisNo);
				else
					Result += ViscousFlux(Cell2,Cell3,AxisNo);
      }
			break;
	};

	return Result;	
}
