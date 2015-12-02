/***************************************************************************
                          PrintGrid.cpp  -  description
                             -------------------
    begin                : Tue Mar 19 2002
    copyright            : (C) 2002 by Olivier Roussel & Alexei Tsigulin
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
PrintGrid::PrintGrid(int L)
{

	int n=0; // Array size

	localScaleNb=L;
	elementNb = n =	(1<<localScaleNb)+1;
	if (Dimension > 1) elementNb *= n;
	if (Dimension > 2) elementNb *= n;

	Q  = new Vector[elementNb];
	Qt = new Vector[elementNb];
  
  int i;
  for (i=0;i<elementNb;i++) 
  {
    Q[i].setDimension(QuantityNb);
    Qt[i].setDimension(QuantityNb);
  }
  
  
}
/*
______________________________________________________________________________________________

Distructor
______________________________________________________________________________________________

*/
PrintGrid::~PrintGrid()
{
	delete[] Q;
	delete[] Qt;
}
/*
______________________________________________________________________________________________

Set procedure
______________________________________________________________________________________________

*/
void PrintGrid::setValue(const int i, const int j, const int k, const Vector& UserAverage)
{
	// --- Local variables ---

	int n=(1<<localScaleNb)+1;	// n = 2^localScaleNb+1

  *(Q + i + n*(j + n*k)) = UserAverage;
}
/*
______________________________________________________________________________________________

Get procedure
______________________________________________________________________________________________

*/
Vector PrintGrid::value(const int i, const int j, const int k) const
{
	// --- Local variables ---

	int n = (1<<localScaleNb)+1; // n = 2^localScaleNb

	return *(Q + i + n*(j + n*k));
}
/*
______________________________________________________________________________________________

*/

real PrintGrid::value(const int i, const int j, const int k, const int QuantityNo) const
{
	// --- Local variables ---

	int n = (1<<localScaleNb)+1; 						// n = 2^localScaleNb

	return (Q + i + n*(j + n*k))->value(QuantityNo);
}

/*
______________________________________________________________________________________________

*/

real PrintGrid::cellValue(const int i, const int j, const int k, const int QuantityNo) const
{
	// --- Local variables ---

	int n = (1<<localScaleNb)+1; 					// n = 2^localScaleNb+1
	int li=0, lj=0, lk=0;									// local i,j,k


	if (CMin[1] == 2)
	  li = ((i+n)/n==1)? i : (2*n-i-1)%n;  // Neumann
	else
		li = (i+n)%n;                        // Periodic

	// -- in y --

	if (Dimension > 1)
	{
		if (CMin[2] == 2)
	 		lj = ((j+n)/n==1)? j : (2*n-j-1)%n; // Neumann
		else
			lj = (j+n)%n;                       // Periodic
	}

	// -- in z --

	if (Dimension > 2)
	{
		if (CMin[3] == 2)
	  	lk = ((k+n)/n==1)? k : (2*n-k-1)%n; // Neumann
		else
			lk = (k+n)%n;												// Periodic
	}

	return (Q + li + n*(lj + n*lk))->value(QuantityNo);
}

/*
______________________________________________________________________________________________

*/

Vector PrintGrid::velocity(const int i, const int j, const int k) const
{
	Vector V(Dimension);

 for (int AxisNo=1; AxisNo <= Dimension; AxisNo++)
		V.setValue( AxisNo, cellValue(i,j,k,AxisNo+1)/cellValue(i,j,k,1) );

	return V;
}

/*
______________________________________________________________________________________________

*/

real PrintGrid::pressure(const int i, const int j, const int k) const
{
	Vector V(Dimension);
	real rho, rhoE;

	rho  = density(i,j,k);
	V    = velocity(i,j,k);
	rhoE = energy(i,j,k);

	return (Gamma-1)*(rhoE - .5*rho*(V*V));
}

/*
______________________________________________________________________________________________

*/

real PrintGrid::temperature(const int i, const int j, const int k) const
{
	real rho, p;

	if (EquationType >=3 && EquationType <=5)
		return value(i,j,k,1);

	rho  = density(i,j,k);
	p    = pressure(i,j,k);

	return Gamma*Ma*Ma*p/rho;
}
/*
______________________________________________________________________________________________

*/

real PrintGrid::concentration(const int i, const int j, const int k) const
{

	if (EquationType >=3 && EquationType <=5)
		return value(i,j,k,2);

	return 0.;
}
/*
______________________________________________________________________________________________

*/

// Returns 0 in 1D, the scalar vorticity in 2D, the vorticity norm in 3D

real PrintGrid::vorticity(const int i, const int j, const int k) const
{
	int L=localScaleNb;
	real dx=0., dy=0., dz=0.;
	real U=0.,V=0.,W=0.;
	real Uy1=0., Uy2=0., Uz1=0., Uz2=0.;
	real Vx1=0., Vx2=0., Vz1=0., Vz2=0.;
	real Wx1=0., Wx2=0., Wy1=0., Wy2=0.;

	int n = (1<<L); 					// n = 2^localScaleNb

	real result=0.;

	if (Dimension == 1)
		return 0.;

	// Compute vorticity components

	dx = (XMax[1]-XMin[1])/n;
	dy = (XMax[2]-XMin[2])/n;

	Vx1 = velocity(BC(i-1,1,n), BC(j  ,2,n),BC(k,3,n),2);
	Vx2 = velocity(BC(i+1,1,n), BC(j  ,2,n),BC(k,3,n),2);
	Uy1 = velocity(BC(i  ,1,n), BC(j-1,2,n),BC(k,3,n),1);
	Uy2 = velocity(BC(i  ,1,n), BC(j+1,2,n),BC(k,3,n),1);

	if (Dimension == 2)
		W = (Vx2-Vx1)/(2.*dx) - (Uy2-Uy1)/(2.*dy);
	else
	{
		dz = (XMax[3]-XMin[3])/(1<<L);

		Uz1 = velocity(BC(i  ,1,n), BC(j  ,2,n),BC(k-1,3,n),1);
		Uz1 = velocity(BC(i  ,1,n), BC(j  ,2,n),BC(k+1,3,n),1);

		Vz1 = velocity(BC(i  ,1,n), BC(j  ,2,n),BC(k-1,3,n),2);
		Vz1 = velocity(BC(i  ,1,n), BC(j  ,2,n),BC(k+1,3,n),2);

		Wx1 = velocity(BC(i-1,1,n), BC(j  ,2,n),BC(k  ,3,n),3);
		Wx2 = velocity(BC(i+1,1,n), BC(j  ,2,n),BC(k  ,3,n),3);

		Wy1 = velocity(BC(i  ,1,n), BC(j-1,2,n),BC(k  ,3,n),3);
		Wy2 = velocity(BC(i  ,1,n), BC(j+1,2,n),BC(k  ,3,n),3);

		U = (Wy2-Wy1)/(2.*dy) - (Vz2-Vz1)/(2.*dz);
		V = (Uz2-Uz1)/(2.*dz) - (Wx2-Wx1)/(2.*dx);
		W = (Vx2-Vx1)/(2.*dx) - (Uy2-Uy1)/(2.*dy);

	}

	switch(Dimension)
	{
		case 2:
			result = W;
			break;

		case 3:
			result = sqrt(U*U+V*V+W*W);
	};

	return result;
}
/*
______________________________________________________________________________________________

Get temporary value (private)
______________________________________________________________________________________________

*/
Vector PrintGrid::tempValue(const int l, const int i, const int j, const int k) const
{
	// --- Local variables ---

	int li=0, lj=0, lk=0; 		// Local i,j,k
	int n=(1<<l);	// n = 2^l
	int N=(1<<localScaleNb); // n=2^localScaleNb

	// --- Periodicity or Neumann ? (Dirichlet not accepted here) ---

	// -- in x --

	if (CMin[1] == 2)
	  li = ((i+n)/n==1)? i : (2*n-i-1)%n;  // Neumann
	else
		li = (i+n)%n;                        // Periodic

	// -- in y --

	if (Dimension > 1)
	{
		if (CMin[2] == 2)
	 		lj = ((j+n)/n==1)? j : (2*n-j-1)%n; // Neumann
		else
			lj = (j+n)%n;                       // Periodic
	}

	// -- in z --

	if (Dimension > 2)
	{
		if (CMin[3] == 2)
	  	lk = ((k+n)/n==1)? k : (2*n-k-1)%n; // Neumann
		else
			lk = (k+n)%n;												// Periodic
	}

	// --- return pointer to cell ---

	return *( Qt + li + (N+1)*(lj + (N+1)*lk) );
}
/*
______________________________________________________________________________________________

	Refresh procedure
______________________________________________________________________________________________

*/
void PrintGrid::refresh()
{
	for (int n=0; n<elementNb; n++)
		*(Qt+n) = *(Q+n);
}
/*
______________________________________________________________________________________________

	Predict procedure
______________________________________________________________________________________________

*/
void PrintGrid::predict(const int l, const int i, const int j, const int k)
{
	// --- Local variables ---

	int pi=1, pj=1, pk=1;	// Parity of i,j,k

	Vector	Result(QuantityNb);

	// --- Init result with the cell-average value of Qt ---

	Result = tempValue(l-1,(i+4)/2-2,(j+4)/2-2,(k+4)/2-2);

	// --- 1D case ---

	pi = (i%2 == 0)?1:-1;
	Result += (pi*-.125) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2,(k+4)/2-2);
	Result -= (pi*-.125) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2,(k+4)/2-2);

	// --- 2D case ---

 	if (Dimension > 1)
	{
		pj = (j%2 == 0)?1:-1;
		Result += (pj*-.125) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2+1,(k+4)/2-2);
		Result -= (pj*-.125) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2-1,(k+4)/2-2);

		Result += (pi*pj*.015625) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2+1,(k+4)/2-2);
		Result -= (pi*pj*.015625) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2-1,(k+4)/2-2);
		Result -= (pi*pj*.015625) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2+1,(k+4)/2-2);
		Result += (pi*pj*.015625) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2-1,(k+4)/2-2);
  }

	// --- 3D case ---

 	if (Dimension > 2)
	{
		pk = (k%2 == 0)?1:-1;
		Result += (pk*-.125) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2,(k+4)/2-2+1);
		Result -= (pk*-.125) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2,(k+4)/2-2-1);

		Result += (pi*pk*.015625) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2,(k+4)/2-2+1);
		Result -= (pi*pk*.015625) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2,(k+4)/2-2-1);
		Result -= (pi*pk*.015625) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2,(k+4)/2-2+1);
		Result += (pi*pk*.015625) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2,(k+4)/2-2-1);

		Result += (pj*pk*.015625) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2+1,(k+4)/2-2+1);
		Result -= (pj*pk*.015625) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2+1,(k+4)/2-2-1);
		Result -= (pj*pk*.015625) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2-1,(k+4)/2-2+1);
		Result += (pj*pk*.015625) * tempValue(l-1,(i+4)/2-2,(j+4)/2-2-1,(k+4)/2-2-1);

		Result += (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2+1,(k+4)/2-2+1);
		Result -= (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2+1,(k+4)/2-2-1);
		Result -= (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2-1,(k+4)/2-2+1);
		Result += (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2+1,(j+4)/2-2-1,(k+4)/2-2-1);
		Result -= (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2+1,(k+4)/2-2+1);
		Result += (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2+1,(k+4)/2-2-1);
		Result += (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2-1,(k+4)/2-2+1);
		Result -= (pi*pj*pk*-.001953125) * tempValue(l-1,(i+4)/2-2-1,(j+4)/2-2-1,(k+4)/2-2-1);
  }

	setValue(i,j,k,Result);
}
/*
______________________________________________________________________________________________

	Compute point-values
______________________________________________________________________________________________

*/
void PrintGrid::computePointValue()
{
	int i=0, j=0, k=0;
	int l=localScaleNb;
	int n=(1<<l);

	switch(Dimension)
	{
		case 1:
			for (i=0;i<=n;i++)
				setValue(i,j,k,.5*(tempValue(l,i-1,j,k)+tempValue(l,i,j,k)));
			break;

    case 2:
			for (i=0;i<=n;i++)
			for (j=0;j<=n;j++)
				setValue(i,j,k, .25*(tempValue(l,i-1,j-1,k)+tempValue(l,i-1,j,k)+tempValue(l,i,j-1,k)+tempValue(l,i,j,k)) );
			break;

		default:
			for (i=0;i<=n;i++)
			for (j=0;j<=n;j++)
			for (k=0;k<=n;k++)
				setValue(i,j,k,.125*(tempValue(l,i-1,j-1,k-1)+tempValue(l,i-1,j,k-1)+tempValue(l,i,j-1,k-1)+tempValue(l,i,j,k-1)
												+tempValue(l,i-1,j-1,k)+tempValue(l,i-1,j,k)+tempValue(l,i,j-1,k)+tempValue(l,i,j,k)) );
			break;
	};
}

