/***************************************************************************
                          SchemeOsmp3Mc.cpp  -  description
                             -------------------
    begin                : lun dec 1 2003
    copyright            : (C) 2003 by Olivier Roussel and Alexei Tsigulin
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

Vector SchemeOsmp3Mc(const Cell& Cell1,const Cell& Cell2,const Cell& Cell3,const Cell& Cell4, int AxisNo)
{

		// Local variables

		Vector  V2(Dimension);
		Vector  V3(Dimension), V4(Dimension), V5(Dimension);
		Vector  V2Mc(Dimension);
		Vector  V3Mc(Dimension), V4Mc(Dimension), V5Mc(Dimension);
		real 	rho2, rho3, rho4, rho5;				// Densities
		real 	rho2Mc, rho3Mc, rho4Mc, rho5Mc;	
		real	p2, p3, p4, p5;					// Pressures
		real	p2Mc, p3Mc, p4Mc, p5Mc;
		real	rhoE2, rhoE3, rhoE4, rhoE5;	     		// Energies
		real	rhoE2Mc, rhoE3Mc, rhoE4Mc, rhoE5Mc;
		Vector  Result(QuantityNb);
		Vector  Flux3Mc(QuantityNb), Flux4Mc(QuantityNb);  
		Vector  One(QuantityNb);
    		real    cof2, cof3, cof4;                         	// sqrt(Q[i+1]/Q[i])
    		real    usupcof2, usupcof3, usupcof4;         		// 1/(1+cof)
		Vector  U12(Dimension);                                	// half-mesh velocity
		Vector  U23(Dimension);
		Vector  U34(Dimension);   
		Vector  U45(Dimension);
		real    h23, h34, h45;                             	// half-mesh enthalpy
		real    c23, c34, c45;                             	// half-mesh celerity
		
                // --- Diagonal matrix and eigenvectors ---
		
		Matrix  L23, L34, L45;
		Matrix  R34;
		
		Vector  Vp2(QuantityNb);
		Vector  Vp3(QuantityNb);  
    Vector  Vp4(QuantityNb);
		
 		// --- Riemann invariants ---
		
		Vector  DeltaU2(QuantityNb);
		Vector  DeltaU3(QuantityNb);  
		Vector  DeltaU4(QuantityNb);

		Vector   DeltaW2(QuantityNb);
		Vector   DeltaW3(QuantityNb);
		Vector   DeltaW4(QuantityNb);

		Vector   DeltaWNu1(QuantityNb);
		Vector   DeltaWNu2(QuantityNb);
		Vector   DeltaWNu3(QuantityNb);

		real     dx;		// Space step
		
		// --- Local CFL values  ---
		
		Vector   XNu2(QuantityNb);
		Vector   XNu3(QuantityNb);
		Vector   XNu4(QuantityNb);

		// --- Coefficients of the 3rd order scheme ---
		
		Vector   co3_1(QuantityNb);
		Vector   co3_2(QuantityNb);
		Vector   co3_3(QuantityNb); 

		// --- Sign of the accuracy function Phi ---
		
		real    SignPhiCellLeft;
		real	SignPhiPlus1CellLeft;
		real 	SignPhiCellRight;
		real	SignPhiPlus1CellRight;
		
		// --- Rate of the left and right slopes in the characteristic space
		
    		Vector   rRight(QuantityNb);
    		Vector   rLeft(QuantityNb);
		
		// --- Intermediate variables ---
		
    		Matrix   d2(2,QuantityNb);
		real   	 dmm, d4j;
		real   	 d4jp1, d4min;
		real   	 dm4, djp12;
		
		// --- Accuracy functions Phi ---
		
    		Vector   Phi0(QuantityNb); 		
    		real     PhiTVDMin, PhiTVDMax, Phi_TVD;
    		Vector   PhiLim(QuantityNb);
		
    		real     GammaMdPlus, GammaLcMinus, PhiLc, PhiMd;
    		real     PhiMin, PhiMax;

    		// --- Time step ---
		
		real     dt=TimeStep;
		
		// --- Order of the OS scheme ---
		
    		int      OsOrder=3;      
		
		// --- Counters and direction ---
		                 
	 	int 	 i, it, Direction;


		// --- One = (1,1,1) ---

		for(i=1; i<=QuantityNb; i++ )
			One.setValue(i,1.);

		// --- Get cell size ---

	  dx = .5*(Cell2.size(AxisNo)+Cell3.size(AxisNo));

 		// --- Get conservative quantities ---

		// density  at previous time


    rho2   = Cell1.tempDensity(); 
    rho3   = Cell2.tempDensity();
    rho4   = Cell3.tempDensity();
    rho5   = Cell4.tempDensity();

		// density  at intermediary time

    rho2Mc   = Cell1.density();  
    rho3Mc   = Cell2.density();
    rho4Mc   = Cell3.density();
    rho5Mc   = Cell4.density();

		// --- Get velocity  at previous time ---

		V2 = Cell1.tempVelocity();  
		V3 = Cell2.tempVelocity();
		V4 = Cell3.tempVelocity();
		V5 = Cell4.tempVelocity();

		// --- Get velocity  at intermediary time ---

 		V2Mc = Cell1.velocity();       
 		V3Mc = Cell2.velocity();
 		V4Mc = Cell3.velocity();
 		V5Mc = Cell4.velocity();

		// --- Get energy at previous time ---

		rhoE2 = Cell1.tempEnergy();
		rhoE3 = Cell2.tempEnergy();
		rhoE4 = Cell3.tempEnergy();
		rhoE5 = Cell4.tempEnergy();

		// --- Get energy at intermediary time ---

		rhoE2Mc = Cell1.energy();
		rhoE3Mc = Cell2.energy();
		rhoE4Mc = Cell3.energy();
		rhoE5Mc = Cell4.energy();

		// --- Get pressure at previous time ---

		p2 = Cell1.tempPressure(); 
		p3 = Cell2.tempPressure();
		p4 = Cell3.tempPressure();
		p5 = Cell4.tempPressure();

		// --- Get pressure at intermediary time ---

		p2Mc = Cell1.pressure(); 
		p3Mc = Cell2.pressure();
		p4Mc = Cell3.pressure();
		p5Mc = Cell4.pressure();

    		// --- Compute Euler fluxes at intermediary time ---

 		Flux3Mc.setValue(1,rho3Mc*V3Mc.value(AxisNo));   
		Flux4Mc.setValue(1,rho4Mc*V4Mc.value(AxisNo));

		for(i=1; i<=Dimension; i++)
		{
 			Flux3Mc.setValue(i+1,(AxisNo != i)? rho3Mc*V3Mc.value(AxisNo)*V3Mc.value(i) :rho3Mc*V3Mc.value(AxisNo)*V3Mc.value(i)+p3Mc);
			Flux4Mc.setValue(i+1,(AxisNo != i)? rho4Mc*V4Mc.value(AxisNo)*V4Mc.value(i) :rho4Mc*V4Mc.value(AxisNo)*V4Mc.value(i)+p4Mc);
		}

 		Flux3Mc.setValue(QuantityNb,(rhoE3Mc+p3Mc)*V3Mc.value(AxisNo));
		Flux4Mc.setValue(QuantityNb,(rhoE4Mc+p4Mc)*V4Mc.value(AxisNo));

		// --- Compute direction for predictor-corrector phases ---
		
    		it = IterationNo%(1<<Dimension);

		if (AxisNo==1)
			Direction = it%2;
			
		else if (AxisNo==2)
			Direction = (it/2)%2;

		else
			Direction = (it/2)/2;

		Direction = (Direction + (StepNo-1))%2;

		if (Direction == 1)
			Result = Flux3Mc;
	
		else
			Result = Flux4Mc;

		// --- PREDICTOR PHASE : return Flux3Mc --------------------------------------------------	
			
		if (StepNo == 1)
			return Result;

		// --- CORRECTOR PHASE -------------------------------------------------------------------
			
     		// --- Compute Roe's averages ---

		cof2 = sqrt(rho3/rho2);
		cof3 = sqrt(rho4/rho3);
		cof4 = sqrt(rho5/rho4);

		usupcof2=1./(1.+ cof2);
		usupcof3=1./(1.+ cof3);
		usupcof4=1./(1.+ cof4);

		// --- Compute half-mesh velocity, enthalpy and speed of sound ---
		
		U23 = usupcof2*(cof2*V3+V2);
		U34 = usupcof3*(cof3*V4+V3);
		U45 = usupcof4*(cof4*V5+V4);
	
		h23 = ((cof2*(rhoE3+p3)/rho3)+((rhoE2+p2)/rho2))*usupcof2;
		h34 = ((cof3*(rhoE4+p4)/rho4)+((rhoE3+p3)/rho3))*usupcof3;
		h45 = ((cof4*(rhoE5+p5)/rho5)+((rhoE4+p4)/rho4))*usupcof4;

		c23 = sqrt((Gamma-1.)*(h23 -0.5*(U23*U23)));
		c34 = sqrt((Gamma-1.)*(h34 -0.5*(U34*U34)));
		c45 = sqrt((Gamma-1.)*(h45 -0.5*(U45*U45)));


    		// --- Compute eigenvalues in half-cells ---

		for (i=1;i<=Dimension;i++)
		{
			Vp2.setValue(i,U23.value(AxisNo));
			Vp3.setValue(i,U34.value(AxisNo));
			Vp4.setValue(i,U45.value(AxisNo));
    }

		Vp2.setValue(Dimension+1,U23.value(AxisNo)+c23);
		Vp3.setValue(Dimension+1,U34.value(AxisNo)+c34);
		Vp4.setValue(Dimension+1,U45.value(AxisNo)+c45);

		Vp2.setValue(Dimension+2,U23.value(AxisNo)-c23);
		Vp3.setValue(Dimension+2,U34.value(AxisNo)-c34);
		Vp4.setValue(Dimension+2,U45.value(AxisNo)-c45);

		// --- Set left and right eigenmatrices ---

	  	L23.setEigenMatrix(true, AxisNo, U23, c23);
		L34.setEigenMatrix(true, AxisNo, U34, c34);
		L45.setEigenMatrix(true, AxisNo, U45, c45);

		R34.setEigenMatrix(false, AxisNo, U34, c34, h34);

		// --- Compute Riemann invariants --- 

		// Between 1 and 2
		
    		DeltaU2.setValue(1, rho3-rho2);
		
		for (i=1; i<=Dimension; i++)
    			DeltaU2.setValue(i+1, rho3*V3.value(i)-rho2*V2.value(i));
			
    		DeltaU2.setValue(QuantityNb, rhoE3-rhoE2);
		
		// Between 2 and 3

    		DeltaU3.setValue(1, rho4-rho3);
		
		for (i=1; i<=Dimension; i++)
    			DeltaU3.setValue(i+1, rho4*V4.value(i)-rho3*V3.value(i));

	  	DeltaU3.setValue(QuantityNb, rhoE4-rhoE3);
		
		// Between 3 and 4

    		DeltaU4.setValue(1, rho5-rho4);
		
		for (i=1; i<=Dimension; i++)
    			DeltaU4.setValue(i+1, rho5*V5.value(i)-rho4*V4.value(i));

	  	DeltaU4.setValue(QuantityNb, rhoE5-rhoE4);

    		// --- Compute local CFL ---

		XNu2= abs(Vp2)*dt/dx;
		XNu3= abs(Vp3)*dt/dx;
		XNu4= abs(Vp4)*dt/dx;

    		// --- Compute the projection of the Riemann invariants on the left eigenmatrix --

		DeltaW2=L23*DeltaU2;
		DeltaW3=L34*DeltaU3;
		DeltaW4=L45*DeltaU4;

    		// --- Compute DeltaWNu = (1- XNu)*DeltaW --- 

		// Term by term product (Warning: index moves from 1 step !)

    		DeltaWNu1= (One-XNu2) | DeltaW2;   
    		DeltaWNu2= (One-XNu3) | DeltaW3;
    		DeltaWNu3= (One-XNu4) | DeltaW4;

   		 // --- Correction term per wave ---

		// 3rd order correction coefficients
		// (Warning: index moves from 1 step !)

		co3_1= (One+XNu2)/3.;
		co3_2= (One+XNu3)/3.;
		co3_3= (One+XNu4)/3.;

    		// --- invariant divisions ---

		for (int m=1; m<=QuantityNb; m++)
		{
			// --- Compute slope rate for the Riemann invariants r+ ---
			
			SignPhiCellLeft      = Sign(DeltaWNu1.value(m));
			SignPhiPlus1CellLeft = Sign(DeltaWNu2.value(m));
			
     			rRight.setValue(m, SignPhiCellLeft*SignPhiPlus1CellLeft*(Abs(DeltaWNu1.value(m))+1.e-12)/(Abs(DeltaWNu2.value(m))+1.e-12));

			// --- Compute slope rate for the Riemann invariants r- ---   

			SignPhiCellRight      = Sign(DeltaWNu3.value(m));
			SignPhiPlus1CellRight = Sign(DeltaWNu2.value(m));
			
			rLeft.setValue(m, SignPhiCellRight*SignPhiPlus1CellRight*(Abs(DeltaWNu3.value(m))+1.e-12)/(Abs(DeltaWNu2.value(m))+1.e-12));

    			// --- Accuracy function --- per wave

			// --- 3rd order ---

			Phi0.setValue(m, 1.);

		  	switch(OsOrder)
		  	{
		    		case 3:
			    	default :
					if( Vp3.value(m) < 0.)	
	      					Phi0.setValue(m, Phi0.value(m) - (co3_2.value(m)- rLeft.value(m) * co3_3.value(m)));
          				 else
        					Phi0.setValue(m, Phi0.value(m) - (co3_2.value(m)- rRight.value(m) * co3_1.value(m)));
          
		     		 break;
			};

			// --- Compute curvature

			d2.setValue(1, m, Vp3.value(m)*DeltaW3.value(m) - Vp2.value(m)*DeltaW2.value(m));
			d2.setValue(2, m, Vp4.value(m)*DeltaW4.value(m) - Vp3.value(m)*DeltaW3.value(m));

			// --- Minmod between i and i+1

			if( d2.value(1,m)*d2.value(2,m) < 0.)  
				dmm = 0.;
			else
			{
				dmm = Min( Abs(d2.value(1,m)), Abs(d2.value(2,m)));
				
				if( d2.value(1,m) < 0. )
					dmm= -dmm;
			}

			// More severe

 			d4j= 4.*d2.value(1,m)- d2.value(2,m);
			d4jp1= 4.*d2.value(2,m)- d2.value(1,m);

			if( d4j*d4jp1 < 0.)
				d4min=0.;

			else
			{
				 d4min= Min(Abs(d4j),Abs(d4jp1)) ;
				 if( d4j<0.)  d4min= -d4min;
			}


			if( d4min*dmm < 0.)
				 dm4= 0.;

			else
			{
				 dm4= Min(Abs(d4min),Abs(dmm));
				 if(d4min <0.) dm4= -dm4;
			}
			djp12= dm4;
			
    			// --- TVD flux limiter

			if( Vp3.value(m) > 0.)
			  	PhiTVDMin=2.*rRight.value(m)/(XNu2.value(m)+1.e-12);			  	

			else		
			  	PhiTVDMin=2.*rLeft.value(m)/(XNu4.value(m)+1.e-12);			  	

			PhiTVDMax=2./(1.-XNu3.value(m)+1.e-12);			
			Phi_TVD= Max(0.,Min3(PhiTVDMin,Phi0.value(m),PhiTVDMax));

    			// --- MP flux limiter

			GammaMdPlus=.5*(1.- djp12/(DeltaWNu2.value(m)+1.e-12));
			GammaLcMinus=.5*(1.+ djp12/(DeltaWNu2.value(m)+1.e-12));
			PhiMd=2.*GammaMdPlus/(1.-XNu3.value(m)+1.e-12);
			
			if( Vp3.value(m) > 0.)
				PhiLc=2.*rRight.value(m)*GammaLcMinus*(1.-XNu2.value(m))/(XNu3.value(m)*(1.-XNu3.value(m))+1.e-12);				

			else
				PhiLc=2.*rLeft.value(m)*GammaLcMinus*(1.-XNu4.value(m))/(XNu3.value(m)*(1.-XNu3.value(m))+1.e-12);
				
			PhiMin=Max(Min(0., PhiMd), Min3(0., PhiTVDMin, PhiLc));
			PhiMax=Min(Max(PhiTVDMax, PhiMd), Max3(0., PhiTVDMin, PhiLc));

			// --- Depending on the limiter
			
			switch(LimiterNo)
			{
				default:
				case 1:
//					OSMP3
					PhiLim.setValue(m, Max(PhiMin,Min(Phi0.value(m), PhiMax)));
        			break;

				case 2:
//					Roe
					PhiLim.setValue(m, 0.);
        			break;

 				case 3:
//					Lax-Wendroff
					PhiLim.setValue(m, 1.);
				break;

 				case 4:
//					Beam-Warming
//					PhiLim.setValue(m, ri.value(m));
 				break;

				case 5:
//					OSTVD
					PhiLim.setValue(m, Phi_TVD);
 				break;

 				case 6:
//					OS
					PhiLim.setValue(m, Phi0.value(m));	
  				break;
			}
	}  //end loop on m

    // --- Compute right Euler Flux ---

	for(int m=1; m<=QuantityNb; m++)
		DeltaU3.setValue( m, Abs(Vp3.value(m))*(1.-PhiLim.value(m))*(1.-XNu3.value(m))*DeltaW3.value(m));		

		// --- Write into Result ---

	Result -= (R34*DeltaU3); 
	
	return Result;
}

