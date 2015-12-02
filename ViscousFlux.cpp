/***************************************************************************
                          ViscousFlux.cpp  -  description
                             -------------------
    begin                : jeu fév 26 2004
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

Vector ViscousFlux(const Cell& Cell1, const Cell& Cell2, int AxisNo)
{
			Vector Result(QuantityNb);

			real   rho1=0., rho2=0., rho=0.;  	// Densities
			Vector V1, V2;  			// Velocties
			Vector V, tau(Dimension);
			Matrix GradV1(Dimension);   		// Velocity gradient
		  	Matrix GradV2(Dimension);
			Matrix GradV(Dimension);
			real   divV=0.;				// Velocity divergence
			real   rhoE1=0., rhoE2=0., rhoE=0.; 	// Energies
			real   T1=0., T2=0., T=0.;		// Temperature
			real   mu=0.;				// Viscosity (Sutherland's law)
			real   MuT=0.;				// Eddy viscosity (only if LES=true)
			real	PrT=0.6;			// Turbulent Prandtl number
			real	dx=0.;				// Size in direction AxisNo
			int  	i=0, j=0;
			real	InvDx=0.;			// Inverse of dx
			real	Y1=0.;				// Partial mass
			real	Y2=0.;	

			// --- Get conservative quantities ---

			// cell size

			dx = .5*(Cell1.size(AxisNo)+Cell2.size(AxisNo));

			// correct dx when close to walls
			
			if (UseBoundaryRegions && (BoundaryRegion(Cell1.center()) > 3 || BoundaryRegion(Cell2.center()) > 3))
				dx *= 0.5;
			
			InvDx = 1./dx;	
				
			// density

      			rho1   = Cell1.density();
      			rho2   = Cell2.density();

			rho 	 = .5*(rho1+rho2);

			// velocity and its gradient

			V1 = Cell1.velocity();
			V2 = Cell2.velocity();

      			V	 = .5*(V1+V2);
      
			// Partial mass
			
			if (ScalarEqNb == 1)
			{
				Y1 = Cell1.average(Dimension+3)/rho1;
				Y2 = Cell2.average(Dimension+3)/rho2;
			}
			
			for (i=1;i<=Dimension;i++)
			for (j=1;j<=Dimension;j++)
			{
				GradV1.setValue(i, j, Cell1.gradient(i,j+1));
				GradV2.setValue(i, j, Cell2.gradient(i,j+1));
      }

			GradV	= (GradV1+GradV2)*0.5;
			
			// If one of the cells is in the boundary, GradV is equal to the one in the fluid cell
			
			if (UseBoundaryRegions)
			{
				if (BoundaryRegion(Cell1.center()) != 0)
					GradV = GradV2;
				else if (BoundaryRegion(Cell2.center()) != 0)
					GradV = GradV1;
			}
					
			// energy

			rhoE1 = Cell1.energy();
			rhoE2 = Cell2.energy();

			rhoE	 = .5*(rhoE1+rhoE2);

			// --- Compute temperatures ---

			T1 = Cell1.temperature();
			T2 = Cell2.temperature();

			T = .5*(T1+T2);

			// --- Correct gradient in the direction AxisNo: grad = (V2-V1)/dx --

			for (j=1; j<=Dimension; j++)
				GradV.setValue( AxisNo, j, (V2.value(j)-V1.value(j))*InvDx );

			// --- Compute divergence of V ---

			divV = 0.;
			for (i=1; i<=Dimension; i++)
      				divV += GradV.value(i, i);

			// --- Compute components of the strain tensor ---

			mu = Sutherland(T);
			if (LES)
				MuT = 0.5*(Smagorinsky(Cell1)+Smagorinsky(Cell2));

			for (i=1;i<=Dimension;i++)
				tau.setValue( i, (i == AxisNo) ?
					 (mu+MuT)/Re*( GradV.value(AxisNo,i) + GradV.value(i, AxisNo)-2.*divV/3. )
					:(mu+MuT)/Re*( GradV.value(AxisNo,i) + GradV.value(i, AxisNo) )
							  );

			// --- Compute fluxes ---

			// Mass

			Result.setValue(1, 0. );

  			// Momentum

			for (i=1;i<=Dimension;i++)
				Result.setValue(i+1,-tau.value(i));

			// Energy
			
			// IMPORTANT : ThermalConduction = 1 / ((Gamma-1.)*Ma*Ma*Pr*Re)

			Result.setValue(Dimension+2, - V*tau - (mu + MuT*Pr/PrT)*ThermalConduction*(T2-T1)*InvDx);
			
			// Partial mass
			
			if (ScalarEqNb == 1)
				Result.setValue(Dimension+3, -0.5*(rho1+rho2)/(Re*Pr*Le)*(Y2-Y1)*InvDx);

	return Result;
}
