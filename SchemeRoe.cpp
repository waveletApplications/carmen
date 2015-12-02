/***************************************************************************
                          SchemeRoe.cpp  -  description
                             -------------------
    begin                : jeu mar 11 2004
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

Returns the Euler flux using Roe's method an the Van Leer limiter (family of MUSCL schemes).
______________________________________________________________________________________________

*/

Vector SchemeRoe(const Cell& Cell1,const Cell& Cell2,const Cell& Cell3,const Cell& Cell4, int AxisNo)
{

		// Local variables

		Vector  V1, V2, V3, V4;                 // Velocities
		real 		rho1, rho2, rho3, rho4;					// Densities
		real		p1, p2, p3, p4;									// Pressures
		real		rhoE1, rhoE2, rhoE3, rhoE4;	    // Energies
		Vector  Result(QuantityNb);

		Vector	F1(QuantityNb);                 // Fluxes
		Vector  F2(QuantityNb);
		Vector  F3(QuantityNb);
		Vector  F4(QuantityNb);

		Vector Q1, Q2, Q3, Q4;									// Conservative quantities

		Vector  FL(QuantityNb),FR(QuantityNb);  // Left and right fluxex
		Vector  QL(QuantityNb),QR(QuantityNb);  // Left and right conservative quantities
		real    rhoL, rhoR;              				// Left and right densities
		real		pL, pR;                         // Left and right pressures
		real		rhoEL, rhoER;                   // Left and right energies
 		Vector	VL(Dimension), VR(Dimension);   // Left and right velocities

		Vector  One(QuantityNb);

		real		rho;														 // central density with Roe's average
		Vector  V(Dimension);										 // central velocity with Roe's average
    		real    H;                               // central enthalpy with Roe's average
    		real    c;                               // central speed of sound with Roe's average
		real		Roe;														 // Coefficient for Roe's average

		Matrix L, R;                              // left and right eigenmatrix
		Matrix Lambda(QuantityNb);								// diagonal matrix containing the eigenvalues
		Matrix A;																	// absolute value of the jacobian matrix

		Vector	 Lim;															// limiter (Van Leer)

	  int i;        														// coutner


		// vector one.

		for(i=1; i<=QuantityNb; i++ )
			One.setValue(i,1.);

 		// --- Get conservative quantities ---

		Q1 = Cell1.average();
		Q2 = Cell2.average();
		Q3 = Cell3.average();
		Q4 = Cell4.average();

		// --- Get primitive variables ---

		// density

    		rho1   = Cell1.density();
    		rho2   = Cell2.density();
    		rho3   = Cell3.density();
    		rho4   = Cell4.density();

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

    // --- Compute Euler fluxes ---

		F1.setValue(1,rho1*V1.value(AxisNo));
		F2.setValue(1,rho2*V2.value(AxisNo));
 		F3.setValue(1,rho3*V3.value(AxisNo));
		F4.setValue(1,rho4*V4.value(AxisNo));

		for(i=1; i<=Dimension; i++)
		{
			F1.setValue(i+1, rho1*V1.value(AxisNo)*V1.value(i) + ((AxisNo == i)? p1 : 0.));
			F2.setValue(i+1, rho2*V2.value(AxisNo)*V2.value(i) + ((AxisNo == i)? p2 : 0.));
 			F3.setValue(i+1, rho3*V3.value(AxisNo)*V3.value(i) + ((AxisNo == i)? p3 : 0.));
			F4.setValue(i+1, rho4*V4.value(AxisNo)*V4.value(i) + ((AxisNo == i)? p4 : 0.));
		}

		F1.setValue(QuantityNb,(rhoE1+p1)*V1.value(AxisNo));
		F2.setValue(QuantityNb,(rhoE2+p2)*V2.value(AxisNo));
 		F3.setValue(QuantityNb,(rhoE3+p3)*V3.value(AxisNo));
		F4.setValue(QuantityNb,(rhoE4+p4)*V4.value(AxisNo));

		// --- Van Leer limiter ---

	 	// Left

	 	Lim = Limiter(Q3-Q2, Q2-Q1);
		FL  = F2 + 0.5*(Lim|(F2-F1)) + 0.5*((One-Lim)|(F3-F2));
		QL  = Q2 + 0.5*(Lim|(Q2-Q1)) + 0.5*((One-Lim)|(Q3-Q2));

		// Right

		Lim = Limiter(Q3-Q2, Q4-Q3);
		FR  = F3 - 0.5*(Lim|(F4-F3)) - 0.5*((One-Lim)|(F3-F2));
		QR  = Q3 - 0.5*(Lim|(Q4-Q3)) - 0.5*((One-Lim)|(Q3-Q2));

/*
		FL = F2;
		FR = F3;
		QL = Q2;
		QR = Q3;
*/

		// --- Extract left and right primitive variables ---

		rhoL = QL.value(1);
		rhoR = QR.value(1);

		for (i=1; i<= Dimension; i++)
		{
			 VL.setValue(i,QL.value(i+1)/rhoL);
			 VR.setValue(i,QR.value(i+1)/rhoR);
		}

		rhoEL=QL.value(QuantityNb);
		rhoER=QR.value(QuantityNb);

		pL = (Gamma-1)*(rhoEL - .5*rhoL*(VL*VL));
		pR = (Gamma-1)*(rhoER - .5*rhoR*(VR*VR));

    		// --- Compute Roe's averages ---

		Roe = sqrt(rhoR/rhoL);

		rho = Roe*rhoL;
		V   = 1./(1.+Roe)*( Roe*VR + VL );
		H   = 1./(1.+Roe)*( Roe*(rhoER+pR)/rhoR + (rhoEL+pL)/rhoL );

		c   = sqrt ( (Gamma-1)*( H - 0.5*(V*V) ) );

    		// --- Compute diagonal matrix containing the absolute value of the eigenvalues ---

		for (i=1;i<=Dimension;i++)
			Lambda.setValue(i,i, fabs(V.value(AxisNo)));

		Lambda.setValue(Dimension+1, Dimension+1, fabs(V.value(AxisNo)+c));
		Lambda.setValue(Dimension+2, Dimension+2, fabs(V.value(AxisNo)-c));

		// --- Set left and right eigenmatrices ---

		L.setEigenMatrix(true, AxisNo, V, c);
		R.setEigenMatrix(false, AxisNo, V, c, H);

		// --- Compute absolute Jacobian matrix ---

  		A = R*Lambda*L;
	
    		// --- Compute Euler Flux ---

		Result = 0.5*(FL+FR) - 0.5*(A*(QR-QL));

		return Result;
}

