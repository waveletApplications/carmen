/***************************************************************************
                          Parameters.cpp  -  description
                             -------------------
    begin                : Thu Jun 7 2001, corrected in May 18 2011 by Margarete Domingues
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

User parameters
______________________________________________________________________________________________

*/
	real	CarmenVersion=1.541;		// Version release
/*
______________________________________________________________________________________________

*/

// 1) Time integration ----------------------------------------------

	int   	StepNb=2;			// Number of steps for the time integration
	int   	IterationNb=0;  		// Number of iterations
	int   	IterationNbRef=1000;  		// Number of iterations for the FV reference computation
	real  	CFL=0.5;			// Courant-Friedrich-Levy number
	real   	TimeStep=0.;			// Time step
	real   	PhysicalTime=0.;		// Physical elapsed time
	int	Refresh=0;			// CPU time refresh rate
	int	RefreshNb=200;			// Number of refreshments
	real  	PrintTime1=0.;			// Time for the first print
	real  	PrintTime2=0.;			// Time for the second print
	real  	PrintTime3=0.;			// Time for the third print
	real  	PrintTime4=0.;			// Time for the fourth print
	real  	PrintTime5=0.;			// Time for the fifth print
	real  	PrintTime6=0.;			// Time for the sixth print
	int	PrintIt1=0;			// Iteration for first print
	int	PrintIt2=0;			// Iteration for second print
	int	PrintIt3=0;			// Iteration for third print
	int	PrintIt4=0;			// Iteration for fourth print
	int	PrintIt5=0;			// Iteration for fifth print
	int	PrintIt6=0;			// Iteration for sixth print
	int 	PrintEvery=0;			// Print every PrintEvery iteration
	int 	ImageNb=0;			// Print ImageNb images
	bool    UseBackup = true;		// true = use Backup procedure
	bool    Recovery=false;			// true = restore data from previous computation
	bool 	TimeAveraging=false;		// true = use time-averaging grid (for turbulence only)
	real	StartTimeAveraging=0.;		// Time where the time-averaging procedure must start
	bool 	ComputeCPUTimeRef=false;	// The reference CPU is not being computed
	bool    centredIC=false;			// Values in initial condition  are put in a mean values of the neighbors, if false they are point in point values directly 
// 2) Solved equations ----------------------------------------------

	int   	EquationType=6;			// Type of equation

/*
---------------------------------------------------------------------

EquationType	= 1	Advection-diffusion
		= 2 	Burgers-diffusion
		= 3	Flame front
		= 4	Flame ball
		= 5	Interaction flame-curl
		= 6	Navier-Stokes

---------------------------------------------------------------------
*/
	int SchemeNb=1;		// Scheme number
	int LimiterNo=1;	// Number
	int ScalarEqNb=0;	// Number of additional scalar equations

/*
---------------------------------------------------------------------

SchemeNb	= 1	McCormack (default)
		= 2	Centered
		= 3	AUSM - TVD-MUSCL
		= 31	AUSMDV
		= 4	ENO 2rd order
		= 5     ENO 3rd order
		= 6 	Osmp3Mc

---------------------------------------------------------------------
*/

// 3) Geometry ------------------------------------------------------

	int	Dimension=1;  			// Dimension (1,2,3)
	int	Coordinate=1;			// 1 = Cartesian, 2 = Spherical in x (only in 1D)
	real 	XMin[4]={-1.,-1.,-1.,-1.};	// Minimal values of coordinates
	real	XMax[4]={1.,1.,1.,1.};		// Maximal values of coordinates
	int	CMin[4]={3,3,3,3};		// Min. boundary condition (1 = Dirichlet, 2 = Neuman, 3 = Periodic)
	int	CMax[4]={3,3,3,3};		// Max. boundary condition (1 = Dirichlet, 2 = Neuman, 3 = Periodic)
	bool    UseBoundaryRegions = false;	// true = use file carmen.bc

// 4) Multiresolution -----------------------------------------------

	bool	Multiresolution=true;		// true = Multiresolution (MR), false = Finite volume (FV) on fine mesh
	bool	TimeAdaptivity=false;		// true = use time adaptivity (only when Multiresolution = true, dummy else)
	int	TimeAdaptivityFactor = 1;	// The factor of time step between two scales is 2^TimeAdaptivityFactor
	int	ScaleNb=5;			// Maximal number of scales allowed
	int	ScaleNbRef=2;			// Number of scales for the FV reference computation
	int	PrintMoreScales=0;		// Scales to add for the printing
	real  	Tolerance=0.;			// Prescribed tolerance	(only for MR)
	real  	ToleranceScale=1.;		// Scale for tolerance	(only for MR when Tolerance = 0.)
	real 	PenalizeFactor=1.;		// Penalization factor (obsolete)
	bool  	ConstantTimeStep=true; 		// true = constant TimeStep
  	real  	ElapsedTime=0.;			// Elapsed time
	real  	ExpectedCompression = 0.;	// Expected compression (for variable tolerance)
	bool  	CVS=false;			// True = use Donoho thresholding in order to perform CVS
	bool	LES=false;			// True = use eddy-viscosity
	real	SmoothCoeff=0.;			// Smoothing coefficient (between 0 and 1) for CVS
	int	ThresholdNorm=0;		// Normalization of the wavelet basis for the threshold
						// (0 = default, 1 = L1, 2 = L2, 3 = L2 on derivatives)

// 5) Physics -------------------------------------------------------

	real	Celerity = 1.;          // Advection-diffusion celerity (0, 1, or -1)
	real 	Viscosity = 1.;		// 0 = no diffusion, 1 = Diffusion
	real	Re=100.;		// Reynolds number  (only for EquationType = 1, 2 or 6)
	real	Pr=0.71;		// Prandtl number   (only for EquationType = 6)
	real 	Gamma=1.4;		// Thermodynamic parameter Gamma   (only for EquationType = 6)
	real	Ma=0.3;			// Mach number  (only for EquationType = 6)
	real	Fr=0.;			// Froude number (only for EquationType = 6). 0 means no gravity
  	real	TRef=273.;		// Reference temperature for Sutherland's law  (only for EquationType = 6)
	real	Circulation=10.;	// Circulation of the vortex (only for EquationType = 5)
	real  	ForceX=0.;		// Volumic force in x-direction
	real  	ForceY=0.;		// Volumic force in y-direction
	real  	ForceZ=0.;		// Volumic force in z-direction
	real	ModelConstant=0.1;	// Constant used in the turbulence model
	real	ThermalConduction=0.;	// Dimensionless thermal conduction
	real	ConstantForce=true;	// False = adapt force to maintain constant energy (in periodic problems)


// 6) Chemistry -----------------------------------------------------

	real 	Alpha=0.64; 		// Temperature ratio		(only for EquationType = 3, 4 or 5)
	real 	Ze=10.;			// Zeldovic number  		(only for EquationType = 3, 4 or 5)
	real 	Le=1.;			// Lewis number      		(only for EquationType = 3, 4 or 5)
	real 	Sigma=5.E-02;		// Radiation coefficient	(only for EquationType = 3, 4 or 5)

// 7) Others --------------------------------------------------------

	bool	debug=false;	   	// true = check if tree is graded
	bool 	FluxCorrection=true;	// true = conservative flux correction
	int	PostProcessing=2;	// Post-processing format. Default = Gnuplot in 1D, DX elsewhere
	bool	DataIsBinary=true;	// Default = use binary, false = use ASCII
	bool 	ZipData=true;		// true = zip data files
	bool	WriteAsPoints=false;	// true = write data as point-values, false = write data as cell-averages

/*
______________________________________________________________________________________________

Internal global variables
______________________________________________________________________________________________

*/
	real	SpaceStep=0.;			// Minimal space step
	real	EigenvalueMax=0.;       	// Maximal eigenvalue
	Vector  QuantityMax;			// Vector containing the maximal quantities
	Vector  QuantityAverage;		// Vector containing the average quantities

	real	pi=acos(-1.);

	int   	StepNo=1;			// Current step number for the time integration
	int   	IterationNo=0;  		// Iteration number
	Timer	CPUTime;			// Timer for CPU time
	int 	Cluster=1;              	// 0 for local execution, 1 for cluster execution
	int   	QuantityNb=2;			// Number of conservative quantites
	double  FVTimeRef=0.;			// FV reference CPU time for 1 iteration

	int	CellNb=0;			// Number of cells
	int	LeafNb=0;			// Number of leaves
	real	TotalLeafNb=0.;			// Total number of leaves on all iterations
	real	TotalCellNb=0.;			// Total number of cells on all iterations

	real	ErrorMax=0.;			// Error Max on grid
	real	ErrorGlobalMax=0.;		// Error Max on grid and on time
	real  	ErrorMid=0.;			// Mean error on the grid
	real  	ErrorGlobalMid=0.;		// Global mean error on space and time
	real  	ErrorL2=0.;			// L2 error on the grid
	real  	ErrorGlobalL2=0.;		// Global L2 error on space and time
	int	ErrorNb=0;			// Number of points for the computation of the mean error
	int 	ErrorGlobalNb=0;		// Number of points for the computation of the global mean error

	// ---Parameters used in the Runge-Kutta-Fehlberg (RKF) method to adapt time step ---

	real 	RKFError=0.;			// Maximum of the relative errors between RK2 and RK3 (only if ConstantTimeStep=false)
	real	RKFAccuracyFactor=1.E-03;	// Desired value of RKFError (only if ConstantTimeStep=false)
	real	RKFSafetyFactor=0.01;		// Safety factor for the computation of the time step (only if ConstantTimeStep=false)

	real  	GlobalMomentum=0.;		// Global momentum
	real  	GlobalEnergy=0.;		// Global energy
	real  	GlobalEnstrophy=0.;		// Global enstrophy (only for EquationType = 6)
	real  	ExactMomentum=0.;		// Global momentum for the exact solution (only for EquationType = 1 or 2)
	real  	ExactEnergy=0.;			// Global energy for the exact solution (only for EquationType = 1 or 2)
	real	FlameVelocity=0.;		// Flame velocity
	real	GlobalMomentumOld=0.;  		// Old global momentum
	real	GlobalVolume=0.;		// Global volume

	real	GlobalReactionRate=0.;		// Global reaction rate
	real  	AverageRadius=0.;		// Average radius of the flame ball
	real  	PreviousAverageRadius=0.; 	// Previous average radius
	real  	PreviousAverageRadius2=0.; 	// Previous average radius
	real	XCenter[4]={0.,0.,0.,0.};	// Coordinate of the center of the flame ball
	real	ReactionRateMax=0.;		// Maximum of the reaction rate

	FILE*   GlobalFile;			// Global file
	int  	ChildNb=0;			// Number of children for a given parent (equal to 2**Dimension)

	real 	IntVorticity=0.;		// Integral of the modulus of the vorticity
	real 	IntDensity=0.;			// Integral of the density
	Vector 	IntMomentum;			//  Integral of the modulus of the momentum
	real 	IntEnergy=0.;			// Integral of the energy
	real 	BaroclinicEffect=0.;		// Intensity of the baroclinic effects


//parallel
 	int rank;
	int size;
	int CPUScales;
	real AllXMin[4],AllXMax[4];
  	int AllTaskScaleNb;
	int NeighbourNb;

	int coords[3];
	int CartDims[3];

	int CellElementsNb;
	int MaxCellElementsNb;

	int one_D,two_D;

	int MPISendType;
	int MPIRecvType;

  	Timer CommTimer;

  	int rank_il, rank_iu;
  	int rank_jl, rank_ju;
  	int rank_kl, rank_ku;

  	int WhatSend; //0 - Gradient, 1 - Average

  	// --- See manual for more information about this ---
  	int SendD     = 1 << 0;
  	int SendGrad  = 1 << 1;
  	int SendQ     = 1 << 2;
  	int SendQs    = 1 << 3;
  	int SendX     = 1 << 4;
  	int SenddX    = 1 << 5;

  	char BackupName[255];

#if defined PARMPI
	MPI_Status st[12];
	MPI_Request req[12];
	MPI_Datatype MPItypeSiL,MPItypeSiU,MPItypeRiL,MPItypeRiU;
	MPI_Datatype MPItypeSjL,MPItypeSjU,MPItypeRjL,MPItypeRjU;
	MPI_Datatype MPItypeSkL,MPItypeSkU,MPItypeRkL,MPItypeRkU;

	MPI_Aint *disp=NULL;
	int *blocklen=NULL;
	MPI_Comm comm_cart;

	real *MPIbuffer;
#endif

/*
______________________________________________________________________________________________

*/

/***************************************************************************
 * Init parameters from file <i>carmen.par<i>. If a parameter is not mentioned
 * in this file, the default value is used.
 *
 */
void InitParameters()
{
	// --- Local variables ---------------------------------------------------------------------

	int   i;    	// Counter

	// --- Set global variables from file "carmen.par" -----------------------------------------

	#include "carmen.par"

	// --- Adapt IterationNbRef to the dimension -----------------------------------------------

	IterationNbRef=(int)(exp((4.-Dimension)*log(10.)));

	// --- Compute the number of children of a given parent cell ---

	ChildNb = (1<<Dimension);

#if defined PARMPI

  	AllTaskScaleNb=ScaleNb;
  	for (i=0;i<4;i++)
	{
		AllXMax[i]=XMax[i];
  		AllXMin[i]=XMin[i];
	}

	//some combinations give deadlock...
	MPISendType = 10;  	//0 - Ibsend; 10 - Isend; 20 - Issend;
	MPIRecvType = 1;	//0 - Recv;  1 - Irecv;

	CPUScales=0;
	int tmp=size;
	while ((tmp=(tmp>>1))>0) CPUScales++;
	ScaleNb-=CPUScales/Dimension;

	one_D=1; two_D=1;
	if (Dimension >= 2) one_D=1<<ScaleNb;
	if (Dimension == 3) two_D=1<<ScaleNb;

//#if defined PARMPI

	NeighbourNb=2;
	MaxCellElementsNb=6;

  	// -- Create memory arrays thats are needs for the MPI Type creation ---
  	disp = new MPI_Aint[NeighbourNb*MaxCellElementsNb*one_D*two_D];
	blocklen = new int [NeighbourNb*MaxCellElementsNb*one_D*two_D];

  	// --- Allocate additional memory for MPI buffer send---
	Cell tc;
	int CellElNb,bufsize;
	CellElNb=tc.size().dimension()+tc.center().dimension()+tc.average().dimension()+tc.tempAverage().dimension()+tc.divergence().dimension();

  	if (EquationType==6)
		CellElNb += tc.gradient().lines()*tc.gradient().columns();

	bufsize=(CellElNb*one_D*two_D*NeighbourNb+MPI_BSEND_OVERHEAD)*2*Dimension+1024;
	MPIbuffer=new real[bufsize];
	MPI_Buffer_attach(MPIbuffer,bufsize*sizeof(real));

#else
	NeighbourNb=0;
#endif

#if defined PARMPI
	CreateMPITopology();

  // --- Compute domain coordinates for the processors ---
	XMin[1] = AllXMin[1] + coords[0]*(AllXMax[1]-AllXMin[1])/CartDims[0];
	XMax[1] = AllXMin[1] + (coords[0]+1)*(AllXMax[1]-AllXMin[1])/CartDims[0];

	if (Dimension >= 2)
	{
  		XMin[2] = AllXMin[2] + coords[1]*(AllXMax[2]-AllXMin[2])/CartDims[1];
		XMax[2] = AllXMin[2] + (coords[1]+1)*(AllXMax[2]-AllXMin[2])/CartDims[1];
	}

	if (Dimension == 3)
	{
		XMin[3] = AllXMin[3] + coords[2]*(AllXMax[3]-AllXMin[3])/CartDims[2];
		XMax[3] = AllXMin[3] + (coords[2]+1)*(AllXMax[3]-AllXMin[3])/CartDims[2];
	}

  	// --- Set the backup file name for the current processor
  	sprintf(BackupName,"%d_%d_%d_%s",coords[0],coords[1],coords[2],"carmen.bak");
#else
	sprintf(BackupName,"%s","carmen.bak");
#endif


	// --- Use CVS only if Dimension > 1 -------------------------------------------------------

	if (Dimension == 1)
		CVS = false;

	// --- TimeAveraging always false if not Navier-Stokes -------------------------------------

	if (EquationType != 6)
		TimeAveraging = false;

	// --- If there is no file "carmen.bak", set Recovery=false --------------------------------

	if (!fopen(BackupName,"r"))
  		Recovery = false;

	// --- If PrintMoreScales != 0 or 1 with Multiresolution = false, print error and stop ---

	if (!Multiresolution && (!(PrintMoreScales == 0 || PrintMoreScales == -1)) )
	{
		cout << "Parameters.cpp: In method `void InitParameters()':\n";
		cout << "Parameters.cpp: value of PrintMoreScales incompatible with FV computations\n";
		cout << "Parameters.cpp: must be 0 or -1\n";
		cout << "carmen: *** [Parameters.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}

  	// --- Compute global volume ---------------------------------------------------------------

	GlobalVolume = fabs(XMax[1]-XMin[1]);

	if (Dimension > 1)
		GlobalVolume *= fabs(XMax[2]-XMin[2]);

	if (Dimension > 2)
		GlobalVolume *= fabs(XMax[3]-XMin[3]);

	// --- Compute PostProcessing and DataIsBinary ---------------------------------------------

	// In 1D, use Gnuplot instead of Data Explorer

	if (Dimension == 1 && PostProcessing == 2)
		PostProcessing = 1;

	// In 2D-3D, use Data Explorer instead of Gnuplot

	if (Dimension != 1 && PostProcessing == 1)
		PostProcessing = 2;

	// --- Compute number of conservative quantities -------------------------------------------

	switch(EquationType)
	{
		// --- ADVECTION-DIFFUSION AND BURGERS ---

		case 1:
		case 2:
			QuantityNb = 1;
			break;

		// --- FLAME FRONT, FLAME BALL AND INTERACTION FLAME-CURL ---

		case 3:
		case 4:
			QuantityNb = 2;
			break;

		case 5:
			QuantityNb = 2;
			break;

		// --- NAVIER-STOKES ---

		case 6:
			QuantityNb = Dimension + 2 + ScalarEqNb;
			break;
	};

	// --- Set the dimension of QuantityMax to QuantityNb ---------------------------------------

  	QuantityMax.setDimension(QuantityNb);

  	// --- Set the dimension of QuantityAverage to 4 (pressure, vorticity, entropy, volume)

	QuantityAverage.setDimension(4);

  	// --- Set the dimension of IntMomentum to dimension ----------------------------------------

	IntMomentum.setDimension(Dimension);

	// --- Compute the dimensionless thermal conduction (Navier-Stokes only) --------------------

	if (EquationType == 6)
		ThermalConduction = 1./((Gamma-1.)*Ma*Ma*Pr*Re);

	// --- Compute minimal space step -----------------------------------------------------------

  	SpaceStep = fabs(XMax[1]-XMin[1]);

	for (i = 2; i <= Dimension; i++)
  		SpaceStep = Min(SpaceStep, fabs(XMax[i]-XMin[i]));

	SpaceStep /= (1<<ScaleNb);

	// --- Compute time step from CFL if TimeStep = 0 -------------------------------------------

	if (TimeStep == 0.)
	{

		// --- Compute time step with CFL ---

		switch (EquationType)
		{
			// ADVECTION-DIFFUSION
			case 1:
				if (Viscosity == 0.)
					TimeStep = CFL*SpaceStep;
				else if (Celerity == 0.)
					TimeStep = 0.5*CFL*Re*SpaceStep*SpaceStep;
				else
					TimeStep = CFL*(Re*SpaceStep*SpaceStep)/(Re*SpaceStep+4.);
				break;

			// BURGERS
			case 2:
				if (Viscosity == 0.)
          				TimeStep = CFL*SpaceStep;
				else
					TimeStep = CFL*(Re*SpaceStep*SpaceStep)/(Re*SpaceStep+4.);
				break;

			// FLAME FRONT OR FLAME BALL
			case 3:
			case 4:
			case 5:
				TimeStep = 0.5*CFL*SpaceStep*SpaceStep*Min(Le, 1.);
				break;

			// NAVIER-STOKES : Modified in InitTimeStep if TimeStep == 0.0
			case 6:
				if (Viscosity == 0.)
          				TimeStep =  0.0 ; //CFL*SpaceStep/(1.+1./Ma), commented by Margarete because we need the EigenvalueMax to compute TimeStep;
				else
					TimeStep = CFL*SpaceStep*Min(1./(1.+1./Ma),0.5*Re*SpaceStep);
				break;
		};
	}
	else
		ConstantTimeStep = true;

	// --- Define reference tolerance if this one is equal to zero ---
	// --- Only for convection-diffusion and Burgers ---

/*
	if (Tolerance == 0.)
	{
		if (EquationType == 1 || EquationType == 2)
                {
			Tolerance = 5.E+08*exp(-log(2.)*(StepNb+1)*ScaleNb) / ( Re + exp(log(2.)*(ScaleNb+2.)));
			printf("Tolerance= %e\n",Tolerance);
		}
	}
*/
}
