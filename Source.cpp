/***************************************************************************
                          Source.cpp  -  description
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

Vector Source(Cell& UserCell)
{
	// ---  Local variables ---

   	real 	s;			// Heat loss due to radiation
   	real 	T;			// Temperature
   	real 	Y;			// Concentration
   	real 	Omega;			// Arrhenius term for the reaction rate
	int  	AxisNo; 		// Loop on dimension

	Vector Force(Dimension);

   	Vector Result(QuantityNb);
	
	switch(EquationType)
	{
    		case 1:
		case 2:
		default:

	 		// --- LINEAR ADVECTION AND BURGERS ---
			Result.setZero();	
			return Result;
			break;

		case 3:
		case 4:
		case 5:

	 		// --- FLAME BALL, FLAME FRONT AND INTERACTION FLAME-CURL ---
	  		T = UserCell.average(1);
   			Y = UserCell.average(2);

   			// --- Compute reaction rate and heat loss due to radiation ---

   			Omega = ReactionRate(T,Y);
	 		s = RadiationRate(T);	

  			// --- Set source term ---

  			Result.setValue(1,Omega - s);
  			Result.setValue(2, -Omega);
			break;

		case 6:

			// --- NAVIER-STOKES ---

			for (AxisNo = 1; AxisNo <= Dimension; AxisNo++)
				Force.setValue(AxisNo,(AxisNo==1)? ForceX : ((AxisNo==2) ? ForceY : ForceZ));

			Result.setZero();
				
			for (AxisNo = 1; AxisNo <= Dimension; AxisNo++)
				Result.setValue(AxisNo+1, Force.value(AxisNo));
					
        		Result.setValue(Dimension+2, Force*UserCell.velocity());
			
			// Gravity forces (vertical dimension is y in 2D and z in 3D, default is down)
			
			if (Fr > 0. && Dimension > 1)
				Result.setValue(Dimension+1, Result.value(Dimension+1) -UserCell.density()/Fr);
			
			break;
		};

  return Result;
}
