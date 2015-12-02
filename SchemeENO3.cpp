/***************************************************************************
                          SchemeENO3.cpp  -  description
                             -------------------
    begin                : ven mar 12 2004
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

Returns the Euler flux using an ENO3 reconstruction
______________________________________________________________________________________________

*/

Vector SchemeENO3(const Cell& Cell1,const Cell& Cell2,const Cell& Cell3,const Cell& Cell4, int AxisNo)
{

	// --- Local variables ---

		Vector  V1(Dimension), V2(Dimension);      //Velocities
		Vector  V3(Dimension), V4(Dimension);
		Vector  V5(Dimension), V6(Dimension);
		real 	rho1, rho2, rho3, rho4, rho5, rho6;		// Densities
		real	p1, p2, p3, p4, p5, p6;				// Pressures
		real	rhoE1, rhoE2, rhoE3, rhoE4, rhoE5, rhoE6; 	// Energies
		Vector  Result(QuantityNb);
		Vector  Flux1(QuantityNb), Flux2(QuantityNb), Flux6(QuantityNb);
		Vector  Flux3(QuantityNb), Flux4(QuantityNb), Flux5(QuantityNb);
    		real    cof3;                  // sqrt(Q[i+1]/Q[i])
    		real    usupcof3;      				 //1/(1+cof)

		Vector  U34(Dimension);    // vitesse demi-maille
    		real    h34;               // enthalpie demi-maille
    		real    c34;               // celerite demi-maille
		Vector  Vp3(QuantityNb);   // diagonal matrix
		Matrix  L34, R34;
		Vector		rm1f(6);
		Vector		Bilan(QuantityNb);
		Vector		dnd1(6), dnd2(5), dnd3(4);//Tables des divisions non divisees de Newton
		Matrix		ceno(3,4);  				// coefficients du schema ENO3
		real dx;
		int i, ig;


		// --- Get cell size ---

	  	dx = .5*(Cell2.size(AxisNo)+Cell3.size(AxisNo));

		// --- Initialisation des coefficients du schema ENO2

		ceno.setValue(1,1,  1./3.);
		ceno.setValue(2,1, -7./6.);
		ceno.setValue(3,1,  11./6.);

		ceno.setValue(1,2, -1./6.);
		ceno.setValue(2,2,  5./6.);
		ceno.setValue(3,2,  1./3.);

		ceno.setValue(1,3,  1./3.);
		ceno.setValue(2,3,  5./6.);
 		ceno.setValue(3,3, -1./6.);

		ceno.setValue(1,4,  11./6.);
		ceno.setValue(2,4, -7./6.);
 		ceno.setValue(3,4,  1./3.);


	// --- Get physical quantities ---

			// density
		rho1 = Cell1.density() - Cell1.gradient(AxisNo,1)*dx;
		rho2 = Cell1.density();
		rho3 = Cell2.density();
		rho4 = Cell3.density();
		rho5 = Cell4.density();
		rho6 = Cell4.density() - Cell4.gradient(AxisNo,1)*dx;

			// velocity
		for (i=1;i<=Dimension;i++)
 			V1.setValue(i, Cell1.velocity(i) - Cell1.gradient(AxisNo,i+1)*dx);

		V2 = Cell1.velocity();
		V3 = Cell2.velocity();
		V4 = Cell3.velocity();
		V5 = Cell4.velocity();

		for (i=1;i<=Dimension;i++)
 			V6.setValue(i, Cell4.velocity(i) - Cell4.gradient(AxisNo,i+1)*dx);

			// energy
		rhoE1 = Cell1.energy() - Cell1.gradient(AxisNo, QuantityNb)*dx;
		rhoE2 = Cell1.energy();
		rhoE3 = Cell2.energy();
		rhoE4 = Cell3.energy();
		rhoE5 = Cell4.energy();
		rhoE6 = Cell4.energy()- Cell4.gradient(AxisNo, QuantityNb)*dx;


			// pressure
		p1 = (Gamma-1.)*(rhoE1 - .5*rho1*(V1*V1));
		p2 = Cell1.pressure();
		p3 = Cell2.pressure();
		p4 = Cell3.pressure();
		p5 = Cell4.pressure();
		p6 = (Gamma-1.)*(rhoE6 - .5*rho6*(V6*V6));


			// fluxes
		Flux1.setValue(1, rho1*V1.value(AxisNo));
		Flux2.setValue(1, rho2*V2.value(AxisNo));
		Flux3.setValue(1, rho3*V3.value(AxisNo));
		Flux4.setValue(1, rho4*V4.value(AxisNo));
		Flux5.setValue(1, rho5*V5.value(AxisNo));
		Flux6.setValue(1, rho6*V6.value(AxisNo));


		for (i=1;i<=Dimension;i++)
		{
			Flux1.setValue(i+1, (AxisNo==i)? rho1*V1.value(AxisNo)*V1.value(i)+p1 : rho1*V1.value(AxisNo)*V1.value(i));
			Flux2.setValue(i+1, (AxisNo==i)? rho2*V2.value(AxisNo)*V2.value(i)+p2 : rho2*V2.value(AxisNo)*V2.value(i));
			Flux3.setValue(i+1, (AxisNo==i)? rho3*V3.value(AxisNo)*V3.value(i)+p3 : rho3*V3.value(AxisNo)*V3.value(i));
			Flux4.setValue(i+1, (AxisNo==i)? rho4*V4.value(AxisNo)*V4.value(i)+p4 : rho4*V4.value(AxisNo)*V4.value(i));
			Flux5.setValue(i+1, (AxisNo==i)? rho5*V5.value(AxisNo)*V5.value(i)+p5 : rho5*V5.value(AxisNo)*V5.value(i));
			Flux6.setValue(i+1, (AxisNo==i)? rho6*V6.value(AxisNo)*V6.value(i)+p6 : rho6*V6.value(AxisNo)*V6.value(i));
     }

 		Flux1.setValue(Dimension+2, (rhoE1+p1)*V1.value(AxisNo));
 		Flux2.setValue(Dimension+2, (rhoE2+p2)*V2.value(AxisNo));
 		Flux3.setValue(Dimension+2, (rhoE3+p3)*V3.value(AxisNo));
 		Flux4.setValue(Dimension+2, (rhoE4+p4)*V4.value(AxisNo));
 		Flux5.setValue(Dimension+2, (rhoE5+p5)*V5.value(AxisNo));
 		Flux6.setValue(Dimension+2, (rhoE6+p6)*V6.value(AxisNo));

     // --- compute Roe's averages ---

		cof3 = sqrt(rho4/rho3);
		usupcof3=1./(1.+ cof3);
		U34 =usupcof3*(cof3*V4 + V3);


		h34 = ((cof3*(rhoE4+p4)/rho4)+((rhoE3+p3)/rho3))*usupcof3;
		c34 = sqrt((Gamma-1.)*(h34 -0.5*(U34*U34)));


    // --- Eigenvalues in half-cells --- // Vp_cell_no[m=onde]

		for (i=1;i<=Dimension;i++)
			Vp3.setValue(i,U34.value(AxisNo));

		Vp3.setValue(Dimension+1,U34.value(AxisNo)+c34);
		Vp3.setValue(Dimension+2,U34.value(AxisNo)-c34);

		// --- Set left and right eigenmatrices ---

		L34.setEigenMatrix(true, AxisNo, U34, c34);
		R34.setEigenMatrix(false, AxisNo, U34, c34, h34);

		// --- Reconstruction des flux aux demi-mailles
 		// --- Projection des flux sur les valeurs propres


		for(int m=1; m<=QuantityNb; m++)
		{
			rm1f.setValue(1, L34.value(m,1)*Flux1.value(1)+L34.value(m,2)*Flux1.value(2)+L34.value(m,3)*Flux1.value(3));
			rm1f.setValue(2, L34.value(m,1)*Flux2.value(1)+L34.value(m,2)*Flux2.value(2)+L34.value(m,3)*Flux2.value(3));
			rm1f.setValue(3, L34.value(m,1)*Flux3.value(1)+L34.value(m,2)*Flux3.value(2)+L34.value(m,3)*Flux3.value(3));
			rm1f.setValue(4, L34.value(m,1)*Flux4.value(1)+L34.value(m,2)*Flux4.value(2)+L34.value(m,3)*Flux4.value(3));
			rm1f.setValue(5, L34.value(m,1)*Flux5.value(1)+L34.value(m,2)*Flux5.value(2)+L34.value(m,3)*Flux5.value(3));
			rm1f.setValue(6, L34.value(m,1)*Flux6.value(1)+L34.value(m,2)*Flux6.value(2)+L34.value(m,3)*Flux6.value(3));

			// --- Tables des differences non divisees

			dnd1.setValue(1, rm1f.value(1));
			dnd1.setValue(2, rm1f.value(2));
			dnd1.setValue(3, rm1f.value(3));
			dnd1.setValue(4, rm1f.value(4));
			dnd1.setValue(5, rm1f.value(5));
			dnd1.setValue(6, rm1f.value(6));

			dnd2.setValue(1, dnd1.value(2)-dnd1.value(1));
			dnd2.setValue(2, dnd1.value(3)-dnd1.value(2));
			dnd2.setValue(3, dnd1.value(4)-dnd1.value(3));
			dnd2.setValue(4, dnd1.value(5)-dnd1.value(4));
			dnd2.setValue(5, dnd1.value(6)-dnd1.value(5));

			dnd3.setValue(1, dnd2.value(2)-dnd2.value(1));
			dnd3.setValue(2, dnd2.value(3)-dnd2.value(2));
			dnd3.setValue(3, dnd2.value(4)-dnd2.value(3));
			dnd3.setValue(4, dnd2.value(5)-dnd2.value(4));

	// --- Construction du stencil

			ig=3;  
			if(Vp3.value(m) >= 0.)
			{
				if( fabs(dnd2.value(ig)) >= fabs(dnd2.value(ig-1)) )
					ig=ig-1;

				if( fabs(dnd3.value(ig)) >= fabs(dnd3.value(ig-1)) )
     			ig=ig-1;

			}
			else
			{
      		ig=4;
				if( fabs(dnd2.value(ig)) >= fabs(dnd2.value(ig-1)) )
					ig=ig-1;

				if( fabs(dnd3.value(ig)) >= fabs(dnd3.value(ig-1)) )
					ig=ig-1;
			}

	// --- Calcul de Bilan(m)=> reconstruction

			Bilan.setValue(m, 0.);     //initialisation
			for(int i=1; i<=3; i++)
				Bilan.setValue(m, Bilan.value(m) + ceno.value(i,ig) * rm1f.value(i+ig-1) );

		}	//end loop of m


		Result = R34*Bilan;

	return Result;

}

