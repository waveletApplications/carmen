/***************************************************************************
                          SchemeOsmp5Mc  -  description
                             -------------------
    begin                : mer sep 21 2005
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

Vector SchemeOsmp5Mc(const Cell& Cell1,const Cell& Cell2,const Cell& Cell3,const Cell& Cell4, int AxisNo)
{

		// Local variables

		Vector  V1(Dimension), V2(Dimension), V3(Dimension);
		Vector  V4(Dimension), V5(Dimension), V6(Dimension);
		Vector  V1Mc(Dimension), V2Mc(Dimension), V3Mc(Dimension);
		Vector  V4Mc(Dimension), V5Mc(Dimension), V6Mc(Dimension);
		real 	  rho1, rho2, rho3, rho4, rho5, rho6;			 // Densities
		real 	  rho1Mc, rho2Mc, rho3Mc, rho4Mc, rho5Mc, rho6Mc;	
		real			p1, p2, p3, p4, p5, p6;						           // Pressures
		real			p1Mc, p2Mc, p3Mc, p4Mc, p5Mc, p6Mc;
		real			rhoE1, rhoE2, rhoE3, rhoE4, rhoE5, rhoE6;	     // Energies
		real			rhoE1Mc, rhoE2Mc, rhoE3Mc, rhoE4Mc, rhoE5Mc, rhoE6Mc;
		Vector  Result(QuantityNb);

		Vector  Flux3Mc(QuantityNb), Flux4Mc(QuantityNb);
		Vector  One(QuantityNb);
		Vector  Two(QuantityNb);
		Vector  Three(QuantityNb);

    real    cof1, cof2, cof3, cof4, cof5;                   // sqrt(Q[i+1]/Q[i])
    real    usupcof1, usupcof2, usupcof3, usupcof4, usupcof5; //1/(1+cof)

    Vector  U12(Dimension);                                // vitesse demi-maille
		Vector  U23(Dimension);
		Vector  U34(Dimension);   
    Vector  U45(Dimension);
    Vector  U56(Dimension);

    real    h12, h23, h34, h45, h56;        // enthalpie demi-maille
    real    c12, c23, c34, c45, c56;        // celerite demi-maille

    Vector  Vp1(QuantityNb);                // diagonal matrix
		Vector  Vp2(QuantityNb);
		Vector  Vp3(QuantityNb);  
    Vector  Vp4(QuantityNb);
    Vector  Vp5(QuantityNb);

		Matrix  L12, L23, L34, L45, L56;
		Matrix  R34;

    Vector   DeltaU1(QuantityNb);
		Vector   DeltaU2(QuantityNb);
		Vector   DeltaU3(QuantityNb);  //Invariants de Riemann
    Vector   DeltaU4(QuantityNb);
    Vector   DeltaU5(QuantityNb);

    Vector   DeltaW1(QuantityNb);
		Vector   DeltaW2(QuantityNb);
		Vector   DeltaW3(QuantityNb);
    Vector   DeltaW4(QuantityNb);
    Vector   DeltaW5(QuantityNb);

    Vector   DeltaWNu1(QuantityNb);
		Vector   DeltaWNu2(QuantityNb);
		Vector   DeltaWNu3(QuantityNb);
		Vector   DeltaWNu4(QuantityNb);
		Vector   DeltaWNu5(QuantityNb);

		real     dx;

    Vector   XNu1(QuantityNb);
		Vector   XNu2(QuantityNb);
		Vector   XNu3(QuantityNb);   // CFL 1er direction propre
    Vector   XNu4(QuantityNb);
    Vector   XNu5(QuantityNb);

//    Vector   co3_1(QuantityNb);
		Vector   co3_2(QuantityNb);
		Vector   co3_3(QuantityNb); // Coeff du schema O(3)
		Vector   co3_4(QuantityNb);
//		Vector   co3_5(QuantityNb);

    Vector   co4_1(QuantityNb);
		Vector   co4_2(QuantityNb);
		Vector   co4_3(QuantityNb); // Coeff du schema O(4)
		Vector   co4_4(QuantityNb);
		Vector   co4_5(QuantityNb);

    Vector   co5_1(QuantityNb);
		Vector   co5_2(QuantityNb);
		Vector   co5_3(QuantityNb); // Coeff du schema O(5)
		Vector   co5_4(QuantityNb);
		Vector   co5_5(QuantityNb);



    real     SignPhiCellLeft, SignPhiPlus1CellLeft, SignPhiCellRight, SignPhiPlus1CellRight;
    real     SignPhiCell, SignPhiPlus1Cell;
    Vector   rp_1(QuantityNb);
		Vector   rp_2(QuantityNb);
		Vector   rp_3(QuantityNb);
    Vector   rm_1(QuantityNb);
		Vector   rm_2(QuantityNb);
		Vector   rm_3(QuantityNb);
    Matrix   d2(4,QuantityNb);
    Vector   dmm(3), d4j(3);
		Vector   d4jp1(3), d4min(3);
		Vector   dm4(3), djp12(3);
    Vector   ri(QuantityNb);
		Vector   rip1(QuantityNb);
		Vector   rim1(QuantityNb);
    Vector   CellNo(QuantityNb);

    Vector   Phi0(QuantityNb); 								// Fonction de precision
    real     PhiTVDMin, PhiTVDMax, Phi_TVD;
    Vector   PhiLim(QuantityNb);
    real     GammaMdPlus, GammaLcMinus, PhiLc, PhiMd;
    real     PhiMin, PhiMax;

    real     dt=TimeStep;
    int      osordre=5;                       // Ordre=3,4 ou5 du schema osmp
	  int 			i;
		int 			it, Direction;


		// vecteur dont toutes les composantes sont =1.

		for(i=1; i<=QuantityNb; i++ )
			One.setValue(i,1.);

		// vecteur dont toutes les composantes sont =2.

		for(i=1; i<=QuantityNb; i++ )
			Two.setValue(i,2.);

		// vecteur dont toutes les composantes sont =3.

		for(i=1; i<=QuantityNb; i++ )
			Three.setValue(i,3.);

		// --- Get cell size ---

	  dx = .5*(Cell2.size(AxisNo)+Cell3.size(AxisNo));

 		// --- Get conservative quantities ---

		// density  at previous time

    rho1   = Cell1.tempDensity() - Cell1.tempGradient(AxisNo,1)*dx;
    rho2   = Cell1.tempDensity();
    rho3   = Cell2.tempDensity();
    rho4   = Cell3.tempDensity();
    rho5   = Cell4.tempDensity();
    rho6   = Cell4.tempDensity() - Cell4.tempGradient(AxisNo,1)*dx;

		// density  at Mc time

    rho1Mc   = Cell1.density() - Cell1.gradient(AxisNo,1)*dx;
    rho2Mc   = Cell1.density();
    rho3Mc   = Cell2.density();
    rho4Mc   = Cell3.density();
    rho5Mc   = Cell4.density();
    rho6Mc   = Cell4.density() - Cell4.gradient(AxisNo,1)*dx;

		// velocity  at previous time

		for (i=1;i<=Dimension;i++)
 			V1.setValue(i, Cell1.tempVelocity(i) - Cell1.tempGradient(AxisNo,i+1)*dx);

		V2 = Cell1.tempVelocity();
		V3 = Cell2.tempVelocity();
		V4 = Cell3.tempVelocity();
		V5 = Cell4.tempVelocity();

		for (i=1;i<=Dimension;i++)
 			V6.setValue(i, Cell4.tempVelocity(i) - Cell4.tempGradient(AxisNo,i+1)*dx);



		// velocity  at Mc time

		for (i=1;i<=Dimension;i++)
 			V1Mc.setValue(i, Cell1.velocity(i) - Cell1.gradient(AxisNo,i+1)*dx);

 		V2Mc = Cell1.velocity();
 		V3Mc = Cell2.velocity();
 		V4Mc = Cell3.velocity();
 		V5Mc = Cell4.velocity();

		for (i=1;i<=Dimension;i++)
 			V6Mc.setValue(i, Cell4.velocity(i) - Cell4.gradient(AxisNo,i+1)*dx);


		// energy at previous time

 		rhoE1 = Cell1.tempEnergy() - Cell1.tempGradient(AxisNo, QuantityNb)*dx;
		rhoE2 = Cell1.tempEnergy();
		rhoE3 = Cell2.tempEnergy();
		rhoE4 = Cell3.tempEnergy();
		rhoE5 = Cell4.tempEnergy();
		rhoE6 = Cell4.tempEnergy() - Cell4.tempGradient(AxisNo, QuantityNb)*dx;

		// energy   at Mc time

 		rhoE1Mc = Cell1.energy() - Cell1.gradient(AxisNo, QuantityNb)*dx;
		rhoE2Mc = Cell1.energy();
		rhoE3Mc = Cell2.energy();
		rhoE4Mc = Cell3.energy();
		rhoE5Mc = Cell4.energy();
 		rhoE6Mc = Cell4.energy() - Cell4.gradient(AxisNo, QuantityNb)*dx;

		// --- Compute pressures --- at previous time

		p1 = (Gamma-1.)*(rhoE1 - .5*rho1*(V1*V1));
		p2 = Cell1.tempPressure();
		p3 = Cell2.tempPressure();
		p4 = Cell3.tempPressure();
		p5 = Cell4.tempPressure();
		p6 = (Gamma-1.)*(rhoE6 - .5*rho6*(V6*V6));
          
		// --- Compute pressures --- at Mc time

		p1Mc = (Gamma-1.)*(rhoE1Mc - .5*rho1Mc*(V1Mc*V1Mc));
		p2Mc = Cell1.pressure();
		p3Mc = Cell2.pressure();
		p4Mc = Cell3.pressure();
		p5Mc = Cell4.pressure();
		p6Mc = (Gamma-1.)*(rhoE6Mc - .5*rho6Mc*(V6Mc*V6Mc));

			if(p1Mc < 0.)
			{
				cout << "DEBUG: p1Mc < 0. !!" << endl;
  	    		cout << p1Mc << endl;
				cout << "IterationNo Mc= " << IterationNo << endl;
				cout << "rhoE1Mc = " << rhoE1Mc << endl;
				cout << "rho1Mc = " << rho1Mc << endl;
				cout << "V1Mc = " << V1Mc << endl;
				exit(1);
			}

    // --- Compute Euler fluxes --- at Mc Time

 		Flux3Mc.setValue(1,rho3Mc*V3Mc.value(AxisNo));
		Flux4Mc.setValue(1,rho4Mc*V4Mc.value(AxisNo));

		for(i=1; i<=Dimension; i++)
		{
 			Flux3Mc.setValue(i+1,(AxisNo != i)? rho3Mc*V3Mc.value(AxisNo)*V3Mc.value(i) :rho3Mc*V3Mc.value(AxisNo)*V3Mc.value(i)+p3Mc);
			Flux4Mc.setValue(i+1,(AxisNo != i)? rho4Mc*V4Mc.value(AxisNo)*V4Mc.value(i) :rho4Mc*V4Mc.value(AxisNo)*V4Mc.value(i)+p4Mc);
		}

 		Flux3Mc.setValue(QuantityNb,(rhoE3Mc+p3Mc)*V3Mc.value(AxisNo));
		Flux4Mc.setValue(QuantityNb,(rhoE4Mc+p4Mc)*V4Mc.value(AxisNo));


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

	if (StepNo == 1)
			return Result;

     // --- compute Roe's averages ---

		cof1 = sqrt(rho2/rho1);
		cof2 = sqrt(rho3/rho2);
		cof3 = sqrt(rho4/rho3);
		cof4 = sqrt(rho5/rho4);
		cof5 = sqrt(rho6/rho5);

		usupcof1=1./(1.+ cof1);
		usupcof2=1./(1.+ cof2);
		usupcof3=1./(1.+ cof3);
		usupcof4=1./(1.+ cof4);
		usupcof5=1./(1.+ cof5);


		U12 =usupcof1*(cof1*V2+ V1);
		U23 =usupcof2*(cof2*V3+ V2);
		U34 =usupcof3*(cof3*V4+ V3);
		U45 =usupcof4*(cof4*V5+ V4);
		U56 =usupcof5*(cof5*V6+ V5);
    

		h12 = ((cof1*(rhoE2+p2)/rho2)+((rhoE1+p1)/rho1))*usupcof1;
		h23 = ((cof2*(rhoE3+p3)/rho3)+((rhoE2+p2)/rho2))*usupcof2;
		h34 = ((cof3*(rhoE4+p4)/rho4)+((rhoE3+p3)/rho3))*usupcof3;
		h45 = ((cof4*(rhoE5+p5)/rho5)+((rhoE4+p4)/rho4))*usupcof4;
		h56 = ((cof5*(rhoE6+p6)/rho6)+((rhoE5+p5)/rho5))*usupcof5;

		c12 = sqrt((Gamma-1.)*(h12 -0.5*(U12*U12)));
		c23 = sqrt((Gamma-1.)*(h23 -0.5*(U23*U23)));
		c34 = sqrt((Gamma-1.)*(h34 -0.5*(U34*U34)));
		c45 = sqrt((Gamma-1.)*(h45 -0.5*(U45*U45)));
		c56 = sqrt((Gamma-1.)*(h56 -0.5*(U56*U56)));

    // --- Eigenvalues in half-cells --- // Vp_cell_no[m=onde]

		for (i=1;i<=Dimension;i++)
		{
			Vp1.setValue(i,U12.value(AxisNo));
			Vp2.setValue(i,U23.value(AxisNo));
			Vp3.setValue(i,U34.value(AxisNo));
			Vp4.setValue(i,U45.value(AxisNo));
			Vp5.setValue(i,U56.value(AxisNo));
    }

		Vp1.setValue(Dimension+1,U12.value(AxisNo)+c12);
		Vp2.setValue(Dimension+1,U23.value(AxisNo)+c23);
		Vp3.setValue(Dimension+1,U34.value(AxisNo)+c34);
		Vp4.setValue(Dimension+1,U45.value(AxisNo)+c45);
		Vp5.setValue(Dimension+1,U56.value(AxisNo)+c56);

		Vp1.setValue(Dimension+2,U12.value(AxisNo)-c12);
		Vp2.setValue(Dimension+2,U23.value(AxisNo)-c23);
		Vp3.setValue(Dimension+2,U34.value(AxisNo)-c34);
		Vp4.setValue(Dimension+2,U45.value(AxisNo)-c45);
		Vp5.setValue(Dimension+2,U56.value(AxisNo)-c56);

		// --- Set left and right eigenmatrices ---

		L12.setEigenMatrix(true, AxisNo, U12, c12);
	  L23.setEigenMatrix(true, AxisNo, U23, c23);
		L34.setEigenMatrix(true, AxisNo, U34, c34);
		L45.setEigenMatrix(true, AxisNo, U45, c45);
		L56.setEigenMatrix(true, AxisNo, U56, c56);

		R34.setEigenMatrix(false, AxisNo, U34, c34, h34);

		// --- Riemann invariants --- delu_cell_no[m=onde]

	  DeltaU1.setValue(1, rho2-rho1);
		for (i=1; i<=Dimension; i++)
		{
	    DeltaU1.setValue(i+1, rho2*V2.value(i)-rho1*V1.value(i));
		}
	    DeltaU1.setValue(QuantityNb, rhoE2-rhoE1);
//
    DeltaU2.setValue(1, rho3-rho2);
		for (i=1; i<=Dimension; i++)
		{
    		DeltaU2.setValue(i+1, rho3*V3.value(i)-rho2*V2.value(i));
		}
    DeltaU2.setValue(QuantityNb, rhoE3-rhoE2);
//
    DeltaU3.setValue(1, rho4-rho3);
		for (i=1; i<=Dimension; i++)
		{
    		DeltaU3.setValue(i+1, rho4*V4.value(i)-rho3*V3.value(i));
		}
	  DeltaU3.setValue(QuantityNb, rhoE4-rhoE3);
//
    DeltaU4.setValue(1, rho5-rho4);
		for (i=1; i<=Dimension; i++)
		{
    		DeltaU4.setValue(i+1, rho5*V5.value(i)-rho4*V4.value(i));
		}
	  DeltaU4.setValue(QuantityNb, rhoE5-rhoE4);
//
    DeltaU5.setValue(1, rho6-rho5);
		for (i=1; i<=Dimension; i++)
		{
    		DeltaU5.setValue(i+1, rho6*V6.value(i)-rho5*V5.value(i));
		}
	  DeltaU5.setValue(QuantityNb, rhoE6-rhoE5);

    // --- Compute local CFL ---

		XNu1= abs(Vp1)*dt/dx;
		XNu2= abs(Vp2)*dt/dx;
		XNu3= abs(Vp3)*dt/dx;
		XNu4= abs(Vp4)*dt/dx;
		XNu5= abs(Vp5)*dt/dx;

    // --- Compute DeltaW --- delw_cell_no[m=onde]

		DeltaW1=L12*DeltaU1;
		DeltaW2=L23*DeltaU2;
		DeltaW3=L34*DeltaU3;
		DeltaW4=L45*DeltaU4;
		DeltaW5=L56*DeltaU5;

    // --- Compute DeltaWNu = (1- xnu)*delw --- delwnu_cell_no[m=onde]

//	Term by term product 

    DeltaWNu1= (One-XNu1) | DeltaW1;    //cell 12
    DeltaWNu2= (One-XNu2) | DeltaW2;    //cell 23
    DeltaWNu3= (One-XNu3) | DeltaW3;    //cell 34 
    DeltaWNu4= (One-XNu4) | DeltaW4;    //cell 45
    DeltaWNu5= (One-XNu5) | DeltaW5;    //cell 56

    // --- Correction term per wave ---

//		3rd order correction coefficients

//		co3_1= (One+XNu1)/3.;
		co3_2= (One+XNu2)/3.;
		co3_3= (One+XNu3)/3.;
		co3_4= (One+XNu4)/3.;
//		co3_5= (One+XNu5)/3.;

//		4th order correction coefficients

		co4_1= (One+XNu1)/3. | (XNu1 -Two)/4.;
		co4_2= (One+XNu2)/3. | (XNu2 -Two)/4.;
		co4_3= (One+XNu3)/3. | (XNu3 -Two)/4.;
		co4_4= (One+XNu4)/3. | (XNu4 -Two)/4.;
		co4_5= (One+XNu5)/3. | (XNu5 -Two)/4.;


//		5th order correction coefficients

		co5_1= (One+XNu1)/3. | (XNu1 -Two)/4. | (XNu1 -Three)/5.;
		co5_2= (One+XNu2)/3. | (XNu2 -Two)/4. | (XNu2 -Three)/5.;
		co5_3= (One+XNu3)/3. | (XNu3 -Two)/4. | (XNu3 -Three)/5.;
		co5_4= (One+XNu4)/3. | (XNu4 -Two)/4. | (XNu4 -Three)/5.;
		co5_5= (One+XNu5)/3. | (XNu5 -Two)/4. | (XNu5 -Three)/5.;

    // --- invariant divisions ---

// --- Compute r+   ---rp_cell[m=onde]

		for (i=1; i<=QuantityNb; i++)
		{
			SignPhiCellLeft= Sign(DeltaWNu1.value(i));
			SignPhiPlus1CellLeft= Sign(DeltaWNu2.value(i));

			SignPhiCell= Sign(DeltaWNu2.value(i));
			SignPhiPlus1Cell= Sign(DeltaWNu3.value(i));

			SignPhiCellRight= Sign(DeltaWNu3.value(i));
			SignPhiPlus1CellRight= Sign(DeltaWNu4.value(i));

     	rp_1.setValue(i, SignPhiCellLeft*SignPhiPlus1CellLeft*(fabs(DeltaWNu1.value(i))+1.e-12)/(fabs(DeltaWNu2.value(i))+1.e-12));
			rp_2.setValue(i, SignPhiCell*SignPhiPlus1Cell*(fabs(DeltaWNu2.value(i))+1.e-12)/(fabs(DeltaWNu3.value(i))+1.e-12));
			rp_3.setValue(i, SignPhiCellRight*SignPhiPlus1CellRight*(fabs(DeltaWNu3.value(i))+1.e-12)/(fabs(DeltaWNu4.value(i))+1.e-12));

		}

// --- Compute r-    ---rm_cell[m=onde]

		for (i=1; i<=QuantityNb; i++)
		{
			SignPhiCellLeft= Sign(DeltaWNu3.value(i));
			SignPhiPlus1CellLeft= Sign(DeltaWNu2.value(i));

			SignPhiCell= Sign(DeltaWNu4.value(i));
			SignPhiPlus1Cell= Sign(DeltaWNu3.value(i));

			SignPhiCellRight= Sign(DeltaWNu5.value(i));
			SignPhiPlus1CellRight= Sign(DeltaWNu4.value(i));


			rm_1.setValue(i, SignPhiCellLeft*SignPhiPlus1CellLeft*(fabs(DeltaWNu3.value(i))+1.e-12)/(fabs(DeltaWNu2.value(i))+1.e-12));
			rm_2.setValue(i, SignPhiCell*SignPhiPlus1Cell*(fabs(DeltaWNu4.value(i))+1.e-12)/(fabs(DeltaWNu3.value(i))+1.e-12));
			rm_3.setValue(i, SignPhiCellRight*SignPhiPlus1CellRight*(fabs(DeltaWNu5.value(i))+1.e-12)/(fabs(DeltaWNu4.value(i))+1.e-12));
		}

// --- Case Vp > 0

		for(int m=1; m<=QuantityNb; m++)
		{
			if( Vp3.value(m) > 0.)			
			{
			  rim1.setValue(m, rp_1.value(m));
			  ri.setValue(m, rp_2.value(m));
			  rip1.setValue(m, rp_3.value(m));
			}

// --- Case Vp <= 0

			else
			{    
			  rip1.setValue(m, rm_1.value(m));
			  ri.setValue(m, rm_2.value(m));
			  rim1.setValue(m, rm_3.value(m));
			}

    
    // --- Accuracy function --- per wave


// --- 3rd order

			Phi0.setValue(m, 1.);

		  	switch(osordre)
		  	{
		    	case 3:
					if( Vp3.value(m) > 0.)
					{
        			Phi0.setValue(m, Phi0.value(m) - (co3_3.value(m)- ri.value(m) * co3_2.value(m)));
           }
           else
					{
	      			Phi0.setValue(m, Phi0.value(m) - (co3_3.value(m)- ri.value(m) * co3_4.value(m)));
					}
		      break;

				case 4:
					if( Vp3.value(m) > 0.)
					{
        			Phi0.setValue(m, Phi0.value(m)+(co4_3.value(m)- 2.*ri.value(m)*co4_2.value(m)+ri.value(m)*rim1.value(m)*co4_1.value(m)));
          	}
					else
					{
	      			Phi0.setValue(m, Phi0.value(m)+(co4_3.value(m)- 2.*ri.value(m)*co4_4.value(m)+ri.value(m)*rim1.value(m)*co4_5.value(m)));
					}
          break;

				case 5:
					if( Vp3.value(m) > 0.)
					{
        			Phi0.setValue(m, Phi0.value(m)-(co5_4.value(m)/rip1.value(m)-3.*co5_3.value(m)+3.*ri.value(m)*co5_2.value(m)-ri.value(m)*rim1.value(m)*co5_1.value(m) ));
          	}
					else
					{
	      			Phi0.setValue(m, Phi0.value(m)-(co5_2.value(m)/rip1.value(m)-3.*co5_3.value(m)+3.*ri.value(m)*co5_4.value(m)-ri.value(m)*rim1.value(m)*co5_5.value(m) ));
					}
          break;

			  };



// --- Compute curvature

			d2.setValue(1, m, Vp2.value(m)*DeltaW2.value(m) - Vp1.value(m)*DeltaW1.value(m));
			d2.setValue(2, m, Vp3.value(m)*DeltaW3.value(m) - Vp2.value(m)*DeltaW2.value(m));
			d2.setValue(3, m, Vp4.value(m)*DeltaW4.value(m) - Vp3.value(m)*DeltaW3.value(m));
			d2.setValue(4, m, Vp5.value(m)*DeltaW5.value(m) - Vp4.value(m)*DeltaW4.value(m));

// --- Minmod between i and i+1

			for(int i=1; i<=3; i++)

			{
				if( d2.value(i,m)*d2.value(i+1,m) < 0.)  
				{
					dmm.setValue(i,0.);
				}
				else
				{
					dmm.setValue(i, Min( fabs(d2.value(i,m)), fabs(d2.value(i+1,m))));
					if( d2.value(i,m) < 0. ) dmm.setValue(i, -dmm.value(i));
				}
  
      }
// More severe

			for(int i=1; i<=3; i++)
			{
 				d4j.setValue(i,   4.*d2.value(i,m)- d2.value(i+1,m));
				d4jp1.setValue(i, 4.*d2.value(i+1,m)- d2.value(i,m));

				if( (d4j.value(i)*d4jp1.value(i) ) <0.)
					d4min.setValue(i, 0.);

				else
				{
				 	d4min.setValue(i, Min(fabs(d4j.value(i)),fabs(d4jp1.value(i)))) ;
				 	if( d4j.value(i)<0.)  d4min.setValue(i, -d4min.value(i)) ;
				}

				if( (d4min.value(i)*dmm.value(i)) <0.)
				 	dm4.setValue(i, 0.);

				else
				{
				 	dm4.setValue(i, Min(fabs(d4min.value(i)),fabs(dmm.value(i))));
				 	if(d4min.value(i) <0.) dm4.setValue(i, -dm4.value(i));
				}
				djp12.setValue(i, dm4.value(i));
			}

    // --- TVD flux limiter

			if( Vp3.value(m) > 0.)
			  	PhiTVDMin=2.*ri.value(m)/(XNu2.value(m)+1.e-12);			  	
			else
			  	PhiTVDMin=2.*ri.value(m)/(XNu4.value(m)+1.e-12);			  	

			PhiTVDMax=2./(1.-XNu3.value(m)+1.e-12);			
			Phi_TVD= Max(0.,Min3(PhiTVDMin,Phi0.value(m),PhiTVDMax)); 


    // --- MP flux limiter

			GammaMdPlus=.5*(1.- djp12.value(2)/(DeltaWNu3.value(m)+1.e-12));


			if( Vp3.value(m) > 0.)
				GammaLcMinus=.5*(1.+ djp12.value(1)/(DeltaWNu3.value(m)+1.e-12));
			else
				GammaLcMinus=.5*(1.+ djp12.value(3)/(DeltaWNu3.value(m)+1.e-12));

			PhiMd=2.*GammaMdPlus/(1.-XNu3.value(m)+1.e-12);

			
			if( Vp3.value(m) > 0.)
				PhiLc=2.*ri.value(m)*GammaLcMinus*(1.-XNu2.value(m))/(XNu3.value(m)*(1.-XNu3.value(m))+1.e-12);				
			else
				PhiLc=2.*ri.value(m)*GammaLcMinus*(1.-XNu4.value(m))/(XNu3.value(m)*(1.-XNu3.value(m))+1.e-12);
				
			PhiMin=Max(Min(0., PhiMd), Min3(0., PhiTVDMin, PhiLc));
			PhiMax=Min(Max(PhiTVDMax, PhiMd), Max3(0., PhiTVDMin, PhiLc));

//				Schema OSMP
			PhiLim.setValue(m, Max(PhiMin,Min(Phi0.value(m), PhiMax)));
////			Schema de Roe
//			PhiLim.setValue(m, 0.);
////			Schema de Lax-Wengroff
//			PhiLim.setValue(m, 1.);
////			Schema de Beam-Warming
//			PhiLim.setValue(m, ri.value(m));
////			Schema OSTVD 
//			PhiLim.setValue(m, Phi_TVD);
//				Schema OS 
//			PhiLim.setValue(m, Phi0.value(m));	//Lax instable OS5 et OS4		

		}  //end loop of m


    // --- Compute right Euler Flux

		for(int m=1; m<=QuantityNb; m++)
		{
			DeltaU3.setValue( m,(1.-PhiLim.value(m))*(1.-XNu3.value(m))*fabs(Vp3.value(m))*DeltaW3.value(m));
		}

		// --- Write into Result ---

		Result -= (R34*DeltaU3); 

		return Result;

}
