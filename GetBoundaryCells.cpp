/***************************************************************************
                          GetBoundaryCells.cpp  -  description
                             -------------------
    begin                : lun jan 24 2005
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
____________________________________________________________________________

Computes the cells C1, C2, C3, C4 in function of the cells Cell1, Cell2,
Cell3, Cell4 taking into account boundary conditions.
____________________________________________________________________________

*/

void GetBoundaryCells(Cell& Cell1, Cell& Cell2, Cell& Cell3, Cell& Cell4,
		Cell& C1, Cell& C2, Cell& C3, Cell& C4, const int AxisNo)
{
  // --- Local variables ---

	int InCell1, InCell2, InCell3, InCell4;     				// Boundary conditions in cells 1, 2, 3, 4
	real P1, P2, P3, P4;							// Pressures in cells 1, 2, 3, 4
	real T1, T2, T3, T4;							// Temperatures in cells 1, 2, 3, 4
	real rho1, rho2, rho3, rho4;						// Densities in cells 1, 2, 3, 4
	Vector V1(Dimension), V2(Dimension), V3(Dimension), V4(Dimension);	// Velocities in cells 1, 2, 3, 4
	real e1, e2, e3, e4;							// Energies in cell 1, 2, 3, 4
	real Y1=0., Y2=0., Y3=0., Y4=0.;					// Partial mass in cell 1, 2, 3, 4
	
	int i;			// Counter
	
 	// --- Init C1, C2, C3, C4 ---

	C1 = Cell1;
	C2 = Cell2;
	C3 = Cell3;
	C4 = Cell4;

	// --- Depending on the boundary region type, transform C1, C2, C3, C4 ---

  	InCell1 = BoundaryRegion(Cell1.center());
  	InCell2 = BoundaryRegion(Cell2.center());
  	InCell3 = BoundaryRegion(Cell3.center());
  	InCell4 = BoundaryRegion(Cell4.center());

	// --- Cell2 IN THE BOUNDARY, Cell3 IN THE FLUID -----------------------------------------

  	if (InCell2 !=0 && InCell3 == 0)
	{
		switch(InCell2)
		{
			// INFLOW
			case 1:
				// Dirichlet on temperature
				T2 = Cell2.temperature();
				T1 = Cell1.temperature();

				// Extrapolate pressure
				P2 = Cell3.oldPressure();
				P1 = P2;
				// P1 = 2*P2 - Cell3.pressure();
				
				// Compute density
				rho2 = Gamma*Ma*Ma*P2/T2;
				rho1 = Gamma*Ma*Ma*P1/T1;
					
				// Dirichlet on momentum
				V2 = (Cell2.density()/rho2)*Cell2.velocity();
				V1 = (Cell1.density()/rho1)*Cell1.velocity();
				
				// Dirichlet on partial mass
				if (ScalarEqNb == 1)
				{
					Y2 = Cell2.average(Dimension+3)/Cell2.density();
					Y1 = Cell1.average(Dimension+3)/Cell1.density();
				}	
				
				// Compute energies
				e2 = P2/((Gamma-1.)*rho2) + 0.5*N2(V2);
				e1 = P1/((Gamma-1.)*rho1) + 0.5*N2(V1);

				// Correct C1 and C2
				C2.setAverage(1, rho2);
				C1.setAverage(1, rho1);
				for (i=1; i<=Dimension; i++)
				{
					C2.setAverage(i+1, rho2*V2.value(i));
					C1.setAverage(i+1, rho1*V1.value(i));
				}
				C2.setAverage(Dimension+2,rho2*e2);
				C1.setAverage(Dimension+2,rho1*e1);
				
				if (ScalarEqNb == 1)
				{
					C2.setAverage(Dimension+3, rho2*Y2);
					C1.setAverage(Dimension+3, rho1*Y1);
				}
				break;
					
			// OUTFLOW : use the old value of the neighbour
			case 2:
				C2.setAverage(Cell3.oldAverage());
				C1.setAverage(Cell3.oldAverage());
				
				// Also change the values in the boundary
				Cell2.setAverage(C2.average());
				Cell1.setAverage(C1.average());
				break;

			// FREE-SLIP WALL : Neuman on all quantities
			case 3:

				C2 = Cell3;
				C1 = Cell4;
				break;

			// ADIABATIC WALL 
			case 4:

				// Dirichlet on velocity
				V2 = Cell2.velocity();
				V1 = Cell1.velocity();

				// Neuman on temperature
				T2 = Cell3.temperature();
				T1 = Cell4.temperature();

				// Neuman on pressure
				P2 = Cell3.pressure();
				P1 = Cell4.pressure();
								
				// Extrapolate pressure
				//P2 = 2*Cell3.pressure()-Cell4.pressure();
				//P1 = P2;

				// Compute densities
				rho2 = Gamma*Ma*Ma*P2/T2;
				rho1 = Gamma*Ma*Ma*P1/T1;
				
				// Compute energies
				e2 = P2/((Gamma-1.)*rho2) + 0.5*N2(V2);
				e1 = P1/((Gamma-1.)*rho1) + 0.5*N2(V1);

				// Correct C1 and C2
				C2.setAverage(1, rho2);
				C1.setAverage(1, rho1);
				for (i=1; i<=Dimension; i++)
				{
					C2.setAverage(i+1, rho2*V2.value(i));
					C1.setAverage(i+1, rho1*V1.value(i));
				}
				C2.setAverage(Dimension+2,rho2*e2);
				C1.setAverage(Dimension+2,rho1*e1);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
				{
					C2.setAverage(Dimension+3, Cell3.average(Dimension+3));
					C1.setAverage(Dimension+3, Cell4.average(Dimension+3));
				}
				
				break;

			// ISOTHERMAL WALL
			case 5:

				// Dirichlet on velocity
				V2 = Cell2.velocity();
				V1 = Cell1.velocity();

				// Dirichlet on temperature
				T2 = Cell2.temperature();
				T1 = Cell1.temperature();

				// Neuman on pressure
				P2 = Cell3.pressure();
				P1 = Cell4.pressure();
				
				// Extrapolate pressure
				//P2 = 2*Cell3.pressure()-Cell4.pressure(); 
				//P1 = P2;

				// Compute densities
				rho2 = Gamma*Ma*Ma*P2/T2;
				rho1 = Gamma*Ma*Ma*P1/T1;
				
				// Compute energies
				e2 = P2/((Gamma-1.)*rho2) + 0.5*N2(V2);
				e1 = P1/((Gamma-1.)*rho1) + 0.5*N2(V1);

				// Correct C1 and C2
				C2.setAverage(1, rho2);
				C1.setAverage(1, rho1);
				for (i=1; i<=Dimension; i++)
				{
					C2.setAverage(i+1, rho2*V2.value(i));
					C1.setAverage(i+1, rho1*V1.value(i));
				}
				C2.setAverage(Dimension+2,rho2*e2);
				C1.setAverage(Dimension+2,rho1*e1);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
				{
					C2.setAverage(Dimension+3, Cell3.average(Dimension+3));
					C1.setAverage(Dimension+3, Cell4.average(Dimension+3));
				}
				
				break;
		};
	return;
	}

  // --- Cell1 IN THE BOUNDARY, Cell2 IN THE FLUID ------------------------------------------------

	if (InCell1 != 0 && InCell2 == 0)
	{
		switch(InCell1)
		{
			// INFLOW 
			case 1:
				// Dirichlet on temperature
				T1 = Cell1.temperature();

				// Extrapolate pressure from old value
				P1 = Cell2.oldPressure();
				
				// Compute density
				rho1 = Gamma*Ma*Ma*P1/T1;
				
				// Dirichlet on momentum
				V1 = (Cell1.density()/rho1)*Cell1.velocity();

				// Dirichlet on partial mass
				if (ScalarEqNb == 1)
					Y1 = Cell1.average(Dimension+3)/Cell1.density();	
				
				// Compute energies
				e1 = P1/((Gamma-1.)*rho1) + 0.5*N2(V1);

				// Correct C1
				C1.setAverage(1, rho1);
				for (i=1; i<=Dimension; i++)
					C1.setAverage(i+1, rho1*V1.value(i));
				C1.setAverage(Dimension+2,rho1*e1);
			
				if (ScalarEqNb == 1)
					C1.setAverage(Dimension+3, rho1*Y1);
				break;
		
			// OUTFLOW : Get old value of the neighbour
			case 2:

				C1.setAverage(Cell2.oldAverage());
				break;
				
			// FREE-SLIP WALL : Neuman on all quantities
			case 3:

				C1 = Cell2;
				break;

			// ADIABATIC WALL
			case 4:

				// Dirichlet on velocity
				V1 = Cell1.velocity();

				// Neuman on temperature
				T1 = Cell2.temperature();
				
				// Neuman on pressure
				P1 = Cell2.pressure();
				
				// Extrapolate pressure
				//P1 = 2*Cell2.pressure()-Cell3.pressure();

				// Compute density
				rho1 = Gamma*Ma*Ma*P1/T1;
				
				// Compute energy
				e1 = P1/((Gamma-1.)*rho1) + 0.5*N2(V1);

				// Correct C1
				C1.setAverage(1, rho1);
				for (i=1; i<=Dimension; i++)
					C1.setAverage(i+1, rho1*V1.value(i));
				C1.setAverage(Dimension+2,rho1*e1);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
					C1.setAverage(Dimension+3, Cell2.average(Dimension+3));
				
				break;

			// ISOTHERMAL WALL 
			case 5:

				// Dirichlet on velocity
				V1 = Cell1.velocity();

				// Dirichlet on temperature
				T1 = Cell1.temperature();

				// Neuman on pressure
				P1 = Cell2.pressure();
				
				// Extrapolate pressure
				//P1 = 2*Cell2.pressure()-Cell3.pressure();

				// Compute density
				rho1 = Gamma*Ma*Ma*P1/T1;
				
				// Compute energies
				e1 = P1/((Gamma-1.)*rho1) + 0.5*N2(V1);

				// Correct C1
				C1.setAverage(1, rho1);
				for (i=1; i<=Dimension; i++)
					C1.setAverage(i+1, rho1*V1.value(i));
				C1.setAverage(Dimension+2,rho1*e1);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
					C1.setAverage(Dimension+3, Cell2.average(Dimension+3));
				
				break;
		};
	return;
	}
	// --- Cell3 IN THE BOUNDARY, Cell2 IN THE FLUID -----------------------------------------

  if (InCell3 !=0 && InCell2 == 0)
	{
		switch(InCell3)
		{
			// INFLOW 
			case 1:
				// Dirichlet on temperature
				T3 = Cell3.temperature();
				T4 = Cell4.temperature();

				// Extrapolate pressure from old value
				P3 = Cell2.oldPressure();
				P4 = P3;
				//P4 = 2*P3 - Cell2.pressure();
				
				// Compute densities
				rho3 = Gamma*Ma*Ma*P3/T3;
				rho4 = Gamma*Ma*Ma*P4/T4;
				
				// Dirichlet on momentum
				V3 = (Cell3.density()/rho3)*Cell3.velocity();
				V4 = (Cell4.density()/rho4)*Cell4.velocity();

				// Dirichlet on partial mass
				if (ScalarEqNb == 1)
				{
					Y3 = Cell3.average(Dimension+3)/Cell3.density();
					Y4 = Cell4.average(Dimension+3)/Cell4.density();
				}
				
				// Compute energies
				e3 = P3/((Gamma-1.)*rho3) + 0.5*N2(V3);
				e4 = P4/((Gamma-1.)*rho4) + 0.5*N2(V4);

				// Correct C1 and C2
				C3.setAverage(1, rho3);
				C4.setAverage(1, rho4);
				for (i=1; i<=Dimension; i++)
				{
					C3.setAverage(i+1, rho3*V3.value(i));
					C4.setAverage(i+1, rho4*V4.value(i));
				}
				C3.setAverage(Dimension+2,rho3*e3);
				C4.setAverage(Dimension+2,rho4*e4);
				
				if (ScalarEqNb == 1)
				{
					C3.setAverage(Dimension+3,rho3*Y3);
					C4.setAverage(Dimension+3,rho4*Y4);
				}
				break;
				
			// OUTFLOW 
			case 2:

				C3.setAverage(Cell2.oldAverage());
				C4.setAverage(Cell2.oldAverage());	
				//C4.setAverage(2*C3.average()-Cell2.average());
				
				// Also change the values in the boundary
				Cell3.setAverage(C3.average());
				Cell4.setAverage(C4.average());
				break;

			// FREE-SLIP WALL : Neuman on all quantities
			case 3:

				C3 = Cell2;
				C4 = Cell1;
				break;

			// ADIABATIC WALL 
			case 4:

				// Dirichlet on velocity
				V3 = Cell3.velocity();
				V4 = Cell4.velocity();

				// Neuman on temperature
				T3 = Cell2.temperature();
				T4 = Cell1.temperature();

				// Neuman on pressure
				P3 = Cell2.pressure();
				P4 = Cell1.pressure();
				
				// Extrapolate pressure
				//P3 = 2*Cell2.pressure()-Cell1.pressure();
				//P4 = P3;

				// Compute densities
				rho3 = Gamma*Ma*Ma*P3/T3;
				rho4 = Gamma*Ma*Ma*P4/T4;
				
				// Compute energies
				e3 = P3/((Gamma-1.)*rho3) + 0.5*N2(V3);
				e4 = P4/((Gamma-1.)*rho4) + 0.5*N2(V4);

				// Correct C3 and C4
				C3.setAverage(1, rho3);
				C4.setAverage(1, rho4);
				for (i=1; i<=Dimension; i++)
				{
					C3.setAverage(i+1, rho3*V3.value(i));
					C4.setAverage(i+1, rho4*V4.value(i));
				}
				C3.setAverage(Dimension+2,rho3*e3);
				C4.setAverage(Dimension+2,rho4*e4);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
				{
					C3.setAverage(Dimension+3, Cell2.average(Dimension+3));
					C4.setAverage(Dimension+3, Cell1.average(Dimension+3));
				}
				
				break;

			// ISOTHERMAL WALL 
			case 5:

				// Dirichlet on velocity
				V3 = Cell3.velocity();
				V4 = Cell4.velocity();

				// Dirichlet on temperature
				T3 = Cell3.temperature();
				T4 = Cell4.temperature();

				// Neuman on pressure
				P3 = Cell2.pressure();
				P4 = Cell1.pressure();
				
				// Extrapolate pressure
				//P3 = 2*Cell2.pressure()-Cell1.pressure();	
				//P4 = P3;

				// Compute densities
				rho3 = Gamma*Ma*Ma*P3/T3;
				rho4 = Gamma*Ma*Ma*P4/T4;
				
				// Compute energies
				e3 = P3/((Gamma-1.)*rho3) + 0.5*N2(V3);
				e4 = P4/((Gamma-1.)*rho4) + 0.5*N2(V4);

				// Correct C3 and C4
				C3.setAverage(1, rho3);
				C4.setAverage(1, rho4);
				for (i=1; i<=Dimension; i++)
				{
					C3.setAverage(i+1, rho3*V3.value(i));
					C4.setAverage(i+1, rho4*V4.value(i));
				}
				C3.setAverage(Dimension+2,rho3*e3);
				C4.setAverage(Dimension+2,rho4*e4);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
				{
					C3.setAverage(Dimension+3, Cell2.average(Dimension+3));
					C4.setAverage(Dimension+3, Cell1.average(Dimension+3));
				}
				
				break;
		};
	return;
	}

  // --- Cell4 IN THE BOUNDARY, Cell3 IN THE FLUID ------------------------------------------------

	if (InCell4 != 0 && InCell3 == 0)
	{
		switch(InCell4)
		{
			// INFLOW 
			case 1:
				// Dirichlet on temperature
				T4 = Cell4.temperature();

				// Extrapolate pressure from old value
				P4 = Cell3.oldPressure();
				
				// Compute density
				rho4 = Gamma*Ma*Ma*P4/T4;
				
				// Dirichlet on momentum
				V4 = (Cell4.density()/rho4)*Cell4.velocity();

				// Dirichlet on partial mass
				if (ScalarEqNb == 1)
					Y4 = Cell4.average(Dimension+3)/Cell4.density();
				
				// Compute energies
				e4 = P4/((Gamma-1.)*rho4) + 0.5*N2(V4);

				// Correct C4
				C4.setAverage(1, rho4);
				for (i=1; i<=Dimension; i++)
					C4.setAverage(i+1, rho4*V4.value(i));
				C4.setAverage(Dimension+2,rho4*e4);
				
				if (ScalarEqNb == 1)
					C4.setAverage(Dimension+3,rho4*Y4);
				break;
		
			// OUTFLOW : Use old cell-average values of the neighbour
			case 2:

				C4.setAverage(Cell3.oldAverage());
				break;

			// FREE-SLIP WALL : Neuman on all quantities
			case 3:

				C4 = Cell3;
				break;

			// ADIABATIC WALL 
			case 4:

				// Dirichlet on velocity
				V4 = Cell4.velocity();

				// Neuman on temperature
				T4 = Cell3.temperature();

				// Neuman on pressure
				P4 = Cell3.pressure();
				
				// Extrapolate pressure
				//P4 = 2*Cell3.pressure()-Cell2.pressure();

				// Compute density
				rho4 = Gamma*Ma*Ma*P4/T4;
				
				// Compute energy
				e4 = P4/((Gamma-1.)*rho4) + 0.5*N2(V4);

				// Correct C4
				C4.setAverage(1, rho4);
				for (i=1; i<=Dimension; i++)
					C4.setAverage(i+1, rho4*V4.value(i));
				C4.setAverage(Dimension+2,rho4*e4);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
					C4.setAverage(Dimension+3, Cell3.average(Dimension+3));
				
				break;

			// ISOTHERMAL WALL 
			case 5:

				// Dirichlet on velocity
				V4 = Cell4.velocity();

				// Dirichlet on temperature
				T4 = Cell4.temperature();

				// Neuman on pressure
				P4 = Cell3.pressure();
				
				// Extrapolate pressure
				//P4 = 2*Cell3.pressure()-Cell2.pressure();

				// Compute density
				rho4 = Gamma*Ma*Ma*P4/T4;
				
				// Compute energies
				e4 = P4/((Gamma-1.)*rho4) + 0.5*N2(V4);

				// Correct C4
				C4.setAverage(1, rho4);
				for (i=1; i<=Dimension; i++)
					C4.setAverage(i+1, rho4*V4.value(i));
				C4.setAverage(Dimension+2,rho4*e4);
				
				// Neuman on partial mass
				if (ScalarEqNb == 1)
					C4.setAverage(Dimension+3, Cell3.average(Dimension+3));
				
				break;
		};
	return;
	}

}
