/***************************************************************************
                          TimeAverageGrid.cpp  -  description
                             -------------------
    begin                : ven déc 3 2004
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

Constructors
______________________________________________________________________________________________

*/
TimeAverageGrid::TimeAverageGrid(const int UserScaleNb, const int UserQuantityNb)
{

	// Init SampleNb, LocalScaleNb, LocalQuantityNb, and ElementNb

	SampleNb 				= 0;
	LocalScaleNb    = UserScaleNb;
	LocalQuantityNb = UserQuantityNb;
	ElementNb       = 1 << ( Dimension*LocalScaleNb );

  // Allocate array of time-averages

//	Q  = new Vector[ElementNb](LocalQuantityNb); //!!!No ISO C++ compartable!!!
  	Q  = new Vector[ElementNb];
    int i;
    for (i=0;i<ElementNb;i++) Q[i].setDimension(LocalQuantityNb); 
}

/*
______________________________________________________________________________________________

*/

TimeAverageGrid::TimeAverageGrid(const int UserScaleNb)
{

  // Default quantity number (3 in 1D, 6 in 2D, 10 in 3D)

	switch(Dimension)
	{
		case 1:
			LocalQuantityNb = 3;
			break;

		case 2:
			LocalQuantityNb = 6;
			break;

		case 3:
			LocalQuantityNb = 10;
			break;

	};

	// Init SampleNb, LocalScaleNb, LocalQuantityNb, and ElementNb

	SampleNb 				= 0;
	LocalScaleNb    = UserScaleNb;
	ElementNb       = 1 << ( Dimension*LocalScaleNb );

  // Allocate array of time-averages

//	Q  = new Vector[ElementNb](LocalQuantityNb); !!!No ISO C++ compartable!!!
  	Q  = new Vector[ElementNb];
    int i;
    for (i=0;i<ElementNb;i++) Q[i].setDimension(LocalQuantityNb);
}
/*
______________________________________________________________________________________________

Distructor
______________________________________________________________________________________________

*/
TimeAverageGrid::~TimeAverageGrid()
{
  // Deallocate array of time-averages

	delete[] Q;
}
/*
______________________________________________________________________________________________

Update value at the position ElementNo
______________________________________________________________________________________________

*/
void TimeAverageGrid::updateValue(const int ElementNo, const int QuantityNo, const real UserValue)
{
	real eps = 1./(SampleNb+1.);

	(Q+ElementNo)->setValue(QuantityNo, eps*UserValue + (1-eps)*(Q+ElementNo)->value(QuantityNo));
}

/*
______________________________________________________________________________________________

Update value at the position i,j,k
______________________________________________________________________________________________

*/
void TimeAverageGrid::updateValue(const int i, const int j, const int k, const int QuantityNo, const real UserValue)
{
	// --- Local variables ---

	int n = (1<<LocalScaleNb);
  int ElementNo = i + n*(j + n*k);

	updateValue(ElementNo, QuantityNo, UserValue);
}
/*
______________________________________________________________________________________________

Update value with the vector of conservative quantities at the position i,j,k
______________________________________________________________________________________________

*/
void TimeAverageGrid::updateValue(const int i, const int j, const int k, const Vector arg)
{
	real rho=0.;
	real U=0., V=0., W=0.;

	// Compute density and velocity

	rho = arg.value(1);
	U = arg.value(2)/rho;
	if (Dimension > 1) V = arg.value(3)/rho;
	if (Dimension > 2) W = arg.value(4)/rho;		


  // update values
  
	updateValue(i,j,k,1,rho);

	switch(Dimension)
	{
		case 1:
			updateValue(i,j,k,2,rho*U);
			updateValue(i,j,k,3,rho*U*U);
			break;

		case 2:
			updateValue(i,j,k,2,rho*U);
			updateValue(i,j,k,3,rho*V);
			updateValue(i,j,k,4,rho*U*U);
			updateValue(i,j,k,5,rho*U*V);
			updateValue(i,j,k,6,rho*V*V);
			break;

		case 3:
			updateValue(i,j,k,2,rho*U);
			updateValue(i,j,k,3,rho*V);
			updateValue(i,j,k,4,rho*W);
			updateValue(i,j,k,5,rho*U*U);
			updateValue(i,j,k,6,rho*U*V);
			updateValue(i,j,k,7,rho*V*V);
			updateValue(i,j,k,8,rho*U*W);
			updateValue(i,j,k,9,rho*V*W);
			updateValue(i,j,k,10,rho*W*W);
			break;
	};

}
/*
______________________________________________________________________________________________

Get element No of the Reynolds-stress tensor at the position i,j,k
______________________________________________________________________________________________

*/
real TimeAverageGrid::stress(const int i, const int j, const int k, const int No) const
{
	real rho = density(i,j,k);
	real rhoV1=0., rhoV2=0., rhoV1V2=0.;

	switch(No)
	{
		case 1:
			rhoV1   = value(i,j,k,2);  // rhoU
			rhoV2   = value(i,j,k,2);  // rhoU
			rhoV1V2 = value(i,j,k,Dimension+2);  // rhoUU
			break;

		case 2:
			rhoV1   = value(i,j,k,2);  // rhoU
			rhoV2   = value(i,j,k,3);  // rhoV
			rhoV1V2 = value(i,j,k,Dimension+3);	// rhoUV
			break;

		case 3:
			rhoV1   = value(i,j,k,3);  // rhoV
			rhoV2   = value(i,j,k,3);  // rhoV
			rhoV1V2 = value(i,j,k,Dimension+4);	// rhoVV
			break;

		case 4:
			rhoV1   = value(i,j,k,2);  // rhoU
			rhoV2   = value(i,j,k,4);  // rhoW
			rhoV1V2 = value(i,j,k,Dimension+5);	// rhoUW
			break;

		case 5:
			rhoV1   = value(i,j,k,3);  // rhoV
			rhoV2   = value(i,j,k,4);  // rhoW
			rhoV1V2 = value(i,j,k,Dimension+6);	// rhoVW
			break;

		case 6:
			rhoV1   = value(i,j,k,4);  // rhoW
			rhoV2   = value(i,j,k,4);  // rhoW
			rhoV1V2 = value(i,j,k,Dimension+7);	// rhoWW
			break;
	};

	return ( ( rhoV1V2 - (rhoV1*rhoV2)/rho )/rho );

}

