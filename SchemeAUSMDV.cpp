/***************************************************************************
                          SchemeAUSMDV  -  description
                             -------------------
    begin                : ven Feb 16 2011
    copyright            : (C) 2011 by  Ralf Deiterding, Margarete Domingues
    email                :
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

Returns the Euler flux using the AUSMDV method with min_mod Roe 1986 limiter (family of MUSCL schemes).

Solve Riemann problems for the 3D Euler equations using  an improved version of the Liou-Steffen Flux-Vector-Splitting

Yasuhiro Wada, Meng-Sing Liou "An accurate and robust flux  splitting scheme for shock and contact discontinuities", SIAM J. Sci. Comput., Vol. 18, No.2, pp 633-657, May 1997.


For 3D

Input: AxisNo = 1, 2 or 3
______________________________________________________________________________________________

*/

Vector SchemeAUSMDV(const Cell& Cell1,const Cell& Cell2,const Cell& Cell3,const Cell& Cell4, int AxisNo)
{

	if(Dimension<3){
           cout<<"Flux AUSMDV is just implement to 3D problems. Program exit."<<endl;
           exit(1);
        }


	// --- Local variables ---------------------------------------------------------------

	// General variables

	Vector	LeftAverage(QuantityNb);	//
	Vector	RightAverage(QuantityNb);	// Conservative quantities
	Vector	Result(QuantityNb);		// Euler flux

	// Variables for the AUSMDV scheme

	Vector	VL(Dimension), VR(Dimension);	// Left and right velocities
	real	rhoL=0., rhoR=0.;		// Left and right densities
	real	eL=0., eR=0.;			// Left and right energies per unit of mass
	real	HL=0., HR=0.;			// Left and right enthalpies
	//real    YL=0., YR=0.;     		// Left and right partial masses

	real	pL =0., pR =0.;			// Left, right  pressures


	// Variables for the limiter

	real   r, Limiter, LeftSlope = 0., RightSlope = 0.; // Left and right slopes
  	//real   DefaultLimiter = (LimiterNo >= 3)? 2.:1.;
	int i;

//	--- Limiter function ---------------------------------------------------------

for (i=1; i<=QuantityNb; i++)
	{
		// --- Compute left cell-average value ---
			
		if (Cell2.average(i) != Cell1.average(i))
		{
			RightSlope 	= Cell3.average(i)-Cell2.average(i);
			LeftSlope	= Cell2.average(i)-Cell1.average(i);
			r		= RightSlope/LeftSlope;
			Limiter 	=  (r > 0) ? (r*r+r)/(1+r*r) : 0.; // same as AUSM
                                             //max (0.0,min(1.0,r));
				            //( (r>=1.0) ? 1.0 : ( (r>=0.0)? r : 0.0) ) ;
			LeftAverage.setValue(i, Cell2.average(i) + double(0.5)*Limiter*LeftSlope);
		}
		else
			LeftAverage.setValue(i, Cell2.average(i));
				
		// --- Compute right cell-average value ---
			
		if (Cell3.average(i) != Cell2.average(i))
		{
			RightSlope = Cell4.average(i)-Cell3.average(i);
			LeftSlope    = Cell3.average(i)-Cell2.average(i);
			r		= RightSlope/LeftSlope;
			Limiter 	=  (r > 0) ? (r*r+r)/(1+r*r) : 0.; //same as AUSM
 						// max(0.0,min(1.0,r)) ;
				            ( (r>=1.0) ? 1.0 : ( (r>=0.0)? r : 0.0) ) ;
			RightAverage.setValue(i, Cell3.average(i) - 0.5*Limiter*LeftSlope);
		}
		else
			RightAverage.setValue(i, Cell3.average(i));
	}
	



	// --- Scheme -------------------------------------------------------------

	// --- Extract left and right natural variables ---

	// Left and right densities
	rhoL = LeftAverage.value(1);
	rhoR = RightAverage.value(1);

	// Left and right velocities
	for (i=1;i<=Dimension;i++)
	{
		VL.setValue( i, LeftAverage.value(i+1)/rhoL );
		VR.setValue( i, RightAverage.value(i+1)/rhoR );
	}

	// Left and right energies per unit of mass
	eL = LeftAverage.value(Dimension+2)/rhoL;
	eR = RightAverage.value(Dimension+2)/rhoR;
	
	
        pL = (Gamma -1.)*rhoL* ( eL - double(0.5)*(VL*VL) );
        pR = (Gamma -1.)*rhoR*( eR - double(0.5)*(VR*VR) );

	// Left and right enthalpies per unit of mass
	HL = eL + pL/rhoL;
	HR = eR + pR/rhoR;

	//Set mu to point to  the component of the system that corresponds to momentum in the direction of this slice, mv and mw to the orthogonal momentum:

	int mu,mv,mw;

//AxisNo=1 dimension=3, velocity positions 2,3,4

	switch (AxisNo){
	case 1:
       		mu = 1;  mv = 2;  mw = 3;
       		break;
	case 2:
		mu = 2 ; mv = 3;  mw = 1;
                break;
      	default:
          	mu = 3;	 mv = 1;  mw = 2;
 }

real uL,uR, vR,vL,wR,wL;
  
uL=VL.value(mu); 	uR=VR.value(mu);
vL=VL.value(mv); 	vR=VR.value(mv);
wL=VL.value(mw);	wR=VR.value(mw);


// -------------------------------------------------------------Compute momentum AUSMD pages 639-640,  eq 31

// ... Auxiliar variables

real	cL =0., cR =0., cMax;			// Left, right  speeds of sound


real aux=0., 
	pLrhoL = pL/rhoL, 
	pRrhoR = pR/rhoR;

// ....... Compute max sound speed  [ Eq. 26 c_m = max(c_R,c_L)]  based on cL, cR, left and right speeds of sound
cL = sqrt(Gamma*pL/rhoL);
cR = sqrt(Gamma*pR/rhoR);
cMax= ( (cL>=cR)?cL:cR ); 


real   alphaL, alphaR,
       uLplus, uRminus, 
       pLplus, pRminus;


// ....... Compute alpha_L and alpha_R Eq 25

aux= pLrhoL+pRrhoR;   
alphaL= double(2.0)* (pLrhoL) / aux;
alphaR= double(2.0)* (pRrhoR) / aux;

// Left-plus u and p. Here we are using aL instead of cm (for us aM).

// ....... Compute uL+, pL+, to avoid the if we do it in steps, first the otherwise

uLplus= double(0.5) * (uL+ fabs(uL));  // Eq 23, otherwise case
pLplus= pL * uLplus/uL;        // Eq 28, otherwise case
	 
if(fabs(uL)<=cMax)
{
    aux=  double(0.25)* (uL+cL)*(uL+cL) /cMax; //ralf use cL  instead of cMax
    uLplus= alphaL* aux + (double(1.0)-alphaL) * uLplus; //Eq. 23 WL
    pLplus= pL* aux/cMax * (double(2.0)-uL/cMax); //Eq. 28 WL
}
	//right-minus u and p. Here we are using aR instead of cm (for us aM).

uRminus= double(0.5) * (uR - fabs(uR)); // Eq 24, otherwise case 
pRminus= pR * uRminus/uR;// Eq 29, otherwise case

        if(fabs(uR)<=cMax)
	{
           aux= (uR-cR)*(uR-cR)/(double(4.0)*cMax);
           uRminus= -alphaR* aux + (double(1.0)-alphaR) * uRminus; //Eq 24 WL
           pRminus=  pR* aux/cMax * (double(2.0)+uR/cMax); //Eq. 29 WL
	}


real auxRhoR, auxRhoL, auxP12;

// MassFlux_Eq22, MassFlux_Eq_31_AUSMD, MassFlux_Eq_30_AUSMV, MassFlux_Eq_33_AUSMDV,Abs_MassFlux_Eq_33_AUSMDV ;

auxP12=pLplus +pRminus;

//MassFlux_Eq22= uLplus*rhoL + uRminus*rhoR;
//MassFlux_Eq_31_AUSMD= 0.5 * (MassFlux_Eq22 * (uL+uR) - fabs(MassFlux_Eq22) * (uR-uL));
//MassFlux_Eq_30_AUSMV= uLplus*rhoL*uL + uRminus*rhoR*uR;

	// -------- Blending between AUSMV and AUSMD  section 2.3 end of page 640, ref eq 30 and 31 -------------
	// sf=1 gives AUSMV, sf= -1 gives AUSMD

real s, Kfactor= double(10.0); //page 642 eq. 34, Kfactor is a forced parameter

aux= (pL<pR)?pL:pR; 
s = min(double(1.0), Kfactor*fabs(pR-pL)/aux); // 0<= sf<= 1/2

//MassFlux_Eq_33_AUSMDV= double(0.5) *( (1+s) * MassFlux_Eq_30_AUSMV + (1-s) * MassFlux_Eq_31_AUSMD);
//Abs_MassFlux_Eq_33_AUSMDV = fabs(MassFlux_Eq_33_AUSMDV);



auxRhoR= double(0.5)* (uRminus*rhoR - fabs(uRminus*rhoR)); 
auxRhoL= double(0.5)*  (uLplus*rhoL  + fabs(uLplus*rhoL)); 

aux = (auxRhoL+auxRhoR);
Result.setValue(1, aux); // for rho

aux = (auxRhoL*HL + auxRhoR*HR);
Result.setValue(Dimension+2, aux);//for rho*H


aux = double(0.5) *( (double(1.0)+s)*uLplus*rhoL*uL  
								+ (double(1.0)-s) *auxRhoL*uL ) 
           + auxP12 + 
           double(0.5) *( (double(1.0)+s)*uRminus*rhoR*uR 
                                                      + (double(1.0)-s) *auxRhoR*uR ) ;
Result.setValue(mu+1,aux); //for velocity component of the axis

aux =  (auxRhoL*vL+ auxRhoR*vR); // for rho
Result.setValue(mv+1,aux);//for velocity component perpendicular to the axis

aux = (auxRhoL*wL + auxRhoR*wR); // for rho
Result.setValue(mw+1,aux);//for velocity component perpendicular to the axis

// --- Return Euler flux ---------------------------------------------------------------

return Result;
}

