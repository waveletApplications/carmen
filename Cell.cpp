/***************************************************************************
                          Cell.cpp  -  description
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

/*
______________________________________________________________________________________________
	
Constructor
______________________________________________________________________________________________

*/

Cell::Cell () : 
	X(Dimension), 
	dX(Dimension), 
	Q(QuantityNb), 
	Qs(QuantityNb),
	Qlow(((!ConstantTimeStep && StepNb > 2) || TimeAdaptivity)? QuantityNb:0), 
	Qold((UseBoundaryRegions)? QuantityNb:0),
	D(QuantityNb), 
	Grad((EquationType==6)? Dimension:0, (EquationType==6)? QuantityNb:0),
	Grads((EquationType==6 && SchemeNb > 5)? Dimension:0, (EquationType==6 && SchemeNb > 5)? QuantityNb:0)
{
	// Empty constructor
	;
}
/*
______________________________________________________________________________________________
	
Distructor
______________________________________________________________________________________________

*/
Cell::~Cell()  	
{
	// Empty distructor
}	
/*
______________________________________________________________________________________________

	Set procedures
______________________________________________________________________________________________

*/	
/*
______________________________________________________________________________________________

*/
real Cell::vorticity(const int AxisNo) const
{
	real result;

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `real Cell::vorticity(int)':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

  	switch(Dimension)
	{
		// 1D case : no vorticiy
		case 1:
			result = 0.;
			break;

		// 2D case : the vorticity is scalar
		case 2:
			result = (AxisNo==3)?(gradient(1,3)-gradient(2,2)):0.;  // dv/dx - du/dy
			break;

		// 3D case: the vorticity is a 3D vector
		case 3:
		default:
			switch(AxisNo)
			{
        			case 1:
					result = gradient(2,4)-gradient(3,3); // dw/dy - dv/dz
					break;

        			case 2:
					result = gradient(3,2)-gradient(1,4); // du/dz - dw/dx
					break;

        			case 3:
				default:
					result = gradient(1,3)-gradient(2,2); // dv/dx - du/dy
					break;
      			}
			break;
	};

	return result;
}
/*
______________________________________________________________________________________________

*/
Vector Cell::vorticity() const
{
	// --- Local variables ---

	Vector W(3);

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `Vector Cell::vorticity()':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	switch(Dimension)
	{
		case 1:
			break;

		case 2:
			W.setValue(3,gradient(1,3)-gradient(2,2)); // dv/dx - du/dy
			break;

		case 3:
			W.setValue(1,gradient(2,4)-gradient(3,3)); // dw/dy - dv/dz
			W.setValue(2,gradient(3,2)-gradient(1,4)); // du/dz - dw/dx
			W.setValue(3,gradient(1,3)-gradient(2,2)); // dv/dx - du/dy
      break;
	};

	return W;
}

/*
______________________________________________________________________________________________

*/
real	Cell::pressure() const
{
	// Conservative quantities;

	real rho, rhoE;
	Vector rhoV(Dimension);

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `real Cell::pressure()':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	// Get conservative quantities

      	rho = Q.value(1);

	for (int i=1 ; i<=Dimension ; i++)
		rhoV.setValue(i,Q.value(i+1));

	rhoE = Q.value(Dimension+2);

	// Return pressure

	return (Gamma-1.)*(rhoE - .5*(rhoV*rhoV)/rho);


}

/*
______________________________________________________________________________________________

*/
real	Cell::tempPressure() const
{
	// Conservative quantities;

	real rho, rhoE;
	Vector rhoV(Dimension);

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `real Cell::pressure()':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	// Get conservative quantities

      rho = Qs.value(1);

			for (int i=1 ; i<=Dimension ; i++)
				rhoV.setValue(i,Qs.value(i+1));

			rhoE = Qs.value(Dimension+2);

	// Return pressure

	return (Gamma-1.)*(rhoE - .5*(rhoV*rhoV)/rho);

}

/*
______________________________________________________________________________________________

*/
real	Cell::oldPressure() const
{
	// Conservative quantities;

	real rho, rhoE;
	Vector rhoV(Dimension);

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `real Cell::pressure()':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	// Get conservative quantities

      	rho = Qold.value(1);

	for (int i=1 ; i<=Dimension ; i++)
		rhoV.setValue(i,Qold.value(i+1));

	rhoE = Qold.value(Dimension+2);

	// Return pressure

	return (Gamma-1.)*(rhoE - .5*(rhoV*rhoV)/rho);

}

/*
______________________________________________________________________________________________

*/
real  Cell::temperature() const
{
	// Conservative quantities;

	real rho, p, T;

	// Warning: not for advection-diffusion and Burgers

	if (EquationType < 3)
	{
		cout << "Cell.cpp: In method `real Cell::temperature()':\n";
		cout << "Cell.cpp: EquationType mismatch\n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	if (EquationType < 6)
	{
  	T = Q.value(1);
	}
	else
	{
		// Get conservative quantities

    rho = density();
    p = pressure();
    T = Gamma*Ma*Ma*p/rho;
	}

	// Return temperature

	return T;
}

/*
______________________________________________________________________________________________

*/
real  Cell::tempTemperature() const
{
	// Conservative quantities;

	real rho, p, T;

	// Warning: not for advection-diffusion and Burgers

	if (EquationType < 3)
	{
		cout << "Cell.cpp: In method `real Cell::temperature()':\n";
		cout << "Cell.cpp: EquationType mismatch\n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	if (EquationType < 6)
	{
  	T = Qs.value(1);
	}
	else
	{
		// Get conservative quantities

    		rho = tempDensity();
		p = tempPressure();
    		T = Gamma*Ma*Ma*p/rho;
	}

	// Return temperature

	return T;
}
/*
______________________________________________________________________________________________

*/
Vector Cell::velocity() const
{

	// Local variables

	Vector V(Dimension);
	int i;

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `Vector Cell::velocity()':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	for (i=1; i<=Dimension; i++)
		V.setValue(i,Q.value(i+1)/Q.value(1));

	return	V;
}

/*
______________________________________________________________________________________________

*/
Vector Cell::tempVelocity() const
{

	// Local variables

	Vector V(Dimension);
	int i;

	// Warning: only for Navier-Stokes equations

	if (EquationType != 6)
	{
		cout << "Cell.cpp: In method `Vector Cell::velocity()':\n";
		cout << "Cell.cpp: EquationType mismatch \n";
		cout << "carmen: *** [Cell.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

	for (i=1; i<=Dimension; i++)
		V.setValue(i,Qs.value(i+1)/Qs.value(1));

	return	V;
}

/*
______________________________________________________________________________________________

*/
real	Cell::volume() const
{
	int AxisNo = 1;
	real result = 1.;
	
	for (AxisNo = 1; AxisNo <= Dimension; AxisNo++)
		result *= size(AxisNo);
		
	return result;
	
}

/*
______________________________________________________________________________________________

	Test procedures
______________________________________________________________________________________________

*/
bool Cell::isOverflow() const
{
	// --- Local variables ---

	int n;					// Counter on the quantities

	// --- If one of the values is overflow, return true ---

	for (n = 1; n <= QuantityNb; n++)
		if (average().isNaN()) return true;

	return false;
}
/*
______________________________________________________________________________________________

	Operator
______________________________________________________________________________________________

*/
void Cell::operator=(const Cell& C)
{
	// local cell becomes equal to cell C

	setCenter(C.center());
	setSize(C.size());
	setAverage(C.average());
	setTempAverage(C.tempAverage());
	setDivergence(C.divergence());
	
	if (EquationType==6)
	{
		setGradient(C.gradient());
		setTempGradient(C.tempGradient());
	}
	
	return;	
}

/*
______________________________________________________________________________________________

	Returns true if the cell is inside the boundary
______________________________________________________________________________________________

*/

bool Cell::isInsideBoundary() const
{
	int ei = 1;
	int ej = (Dimension > 1) ? 1:0;
	int ek = (Dimension > 2) ? 1:0;

	int i, j, k;

	Vector Edge(Dimension);

	bool result=false;

	// Loop: if one edge of the cell is in the boundary, return true
	for (i=-ei; i <= ei; i+=2)
	for (j=-ej; j <= ej; j+=2)
	for (k=-ek; k <= ek; k+=2)
	{
		Edge.setValue(1,center(1) + i*0.5*size(1));

		if (Dimension > 1)
			Edge.setValue(2,center(2) + j*0.5*size(2));
		if (Dimension > 2)
			Edge.setValue(3,center(3) + k*0.5*size(3));

		if (BoundaryRegion(Edge) != 0) result = true;
	}

	return result;

}
/*
______________________________________________________________________________________________

	Returns true if the cell is inside the boundary
______________________________________________________________________________________________

*/

bool Cell::isInFluid() const
{
	return (!isInsideBoundary());
}
