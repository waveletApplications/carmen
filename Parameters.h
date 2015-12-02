/***************************************************************************
                          Parameters.h  -  description
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
 ***************************************************************************
______________________________________________________________________________________________

DESCRIPTION

This header contains all parameters as global variables.
______________________________________________________________________________________________

PUBLIC VARIABLES

1/ Parameters for time integration and performance analysis
______________________________________________________________________________________________

*/		

#include <stdlib.h>

	extern int 		Cluster;					// 0 for local execution, 1 for cluster
	extern int   	StepNb;						// Number of steps for the time integration
	extern int   	StepNo;						// Current step number for the time integration
	extern int   	IterationNb;  		// Number of iterations
	extern int	IterationNbRef; 	// Number of iterations for the FV reference computation
	extern int	IterationNo;			// Current iteration number
	extern real  	CFL;							// Courant-Friedrich-Levy number
	extern real  	TimeStep;					// Time step
	extern real  	PhysicalTime;			// Physical elapsed time
	extern Timer  	CPUTime;					// Timer for CPU time
	extern double 	FVTimeRef;				// FV reference CPU time for 1 iteration
	extern int	Refresh;					// Refresh rate
	extern int	RefreshNb;				// Number of refreshments 
	extern real  	PrintTime1;				// Time for  print
	extern real  	PrintTime2;				//
	extern real  	PrintTime3;				//
	extern real  	PrintTime4;				//
	extern real  	PrintTime5;				//
	extern real  	PrintTime6;				//
	extern int		PrintIt1;					// Iteration for print
	extern int		PrintIt2;					//
	extern int		PrintIt3;					//
	extern int		PrintIt4;					//
	extern int		PrintIt5;					//
	extern int		PrintIt6;					//
	extern int 		PrintEvery;				//	Print every PrintEvery iteration
	extern int 		ImageNb;					//	Print ImageNb images
	extern bool 	TimeAveraging;		// true = use a time-averaging grid (only for turbulence)
	extern real   StartTimeAveraging; // Time where the time-averaging procedure must start
	extern bool   Recovery;					// true = restore data from previous computation
	extern bool   UseBackup;				// true = use Backup procedure.
	extern bool   ConstantTimeStep; // true = constant TimeStep
	extern real   ElapsedTime;			// ElapsedTime
	extern real   SpaceStep;				// Space step
	extern bool   ComputeCPUTimeRef; // True = the reference CPU time is being computed
	extern bool   centredIC;  // Values in initial condition  are put in a mean values of the neighbors, if false they are point in point values directly					


/*
______________________________________________________________________________________________

2/ Parameters for equations

* EquationType
  ------------
  	1 = advection-diffusion
  	2 = diffusive Burgers
   	3 = Flame front
  	4 = Flame ball
		5 = Flame ball-vortex
    6 = Navier-Stokes

______________________________________________________________________________________________

*/		
	extern int	SchemeNb;		// Scheme number
	extern int 	LimiterNo;		// Limiter number
	extern int   	QuantityNb;		// Number of conservative quantites
	extern int   	EquationType;		// Type of equation
	extern int 	ScalarEqNb;		// Number of additional scalar equations

/*
______________________________________________________________________________________________

3/ Geometry
______________________________________________________________________________________________

*/
	extern int	Dimension;  		// Dimension (1,2,3)
	extern int	ChildNb;		// Number of children for a given parent (equal to 2**Dimension)
	extern int	Coordinate;		// 1 = Cartesian, 2 = Spherical in x
	extern real 	XMin[4];		// Minimal values of coordinates;
	extern real	XMax[4];		// Maximal values of coordinates;
	extern int	CMin[4];		// Min. boundary condition (1 = Boundary, 2 = Symetric, 3 = Periodic)
	extern int	CMax[4];		// Max. boundary condition (1 = Boundary, 2 = Symetric, 3 = Periodic)
	extern bool   	UseBoundaryRegions;	// true = use file carmen.bc

/*
______________________________________________________________________________________________

4/ Parameters for multiresolution
______________________________________________________________________________________________

*/		
	extern bool	Multiresolution;	// true = Multiresolution, false = FV on fine mesh
	extern bool	TimeAdaptivity;		// true = use time adaptivity (only when Multiresolution = true, dummy else)
	extern int	TimeAdaptivityFactor;	// The factor of time step between two scales is 2^TimeAdaptivityFactor
	extern int 	CellNb;			// Number of cells	
	extern int	LeafNb;			// Number of leaves
	extern real	TotalCellNb;		// Total number of cells for all iterations
	extern real	TotalLeafNb;		// Total number of leaves for all iterations
	extern int	ScaleNb;		// Maximal number of scales allowed
	extern int	ScaleNbRef;		// Number of scales for the FV reference computation
	extern int	PrintMoreScales;	// More scales to print
	extern real  	Tolerance;		// Prescribed tolerance	
	extern real	ToleranceScale;		// Scale factor for tolerance (when Tolerance = 0).
	extern real  	ToleranceInit;		// Prescribed tolerance	for initial condition (obsolete)
	extern real	PenalizeFactor; 	// Factor of penalization (obsolete)
	extern real   	ExpectedCompression; 	// Expected memory compression
	extern bool   	CVS;			// True = use Donoho thresholding to perform CVS.
	extern bool   	LES;			// True = use eddy-viscosity.
	extern int	ThresholdNorm;		// Normalization of the wavelet basis for the threshold	
	extern real   	SmoothCoeff;		// Smoothing coefficient
/*
______________________________________________________________________________________________

5/ Physics
______________________________________________________________________________________________

*/		
	extern real	Celerity;		// Advection-diffusion celerity (0, 1, -1)
	extern real	Viscosity;		// 0 or 1. 0 means no viscosity.
	extern real	Re;			// Reynolds number
	extern real	Pr;			// Prandtl number
	extern real 	Ma;			// Mach number
	extern real	Fr;			// Froude number
	extern real 	TRef;			// Reference temperature for Sutherland's law
	extern real 	Gamma;			// Thermodynamic parameter Gamma (default=1.4)
	extern real 	Circulation;		// Circulation of the vortex (only for EquationType = 5)
	extern real 	ForceX;			// Volumic force in x-direction
	extern real 	ForceY;			// Volumic force in y-direction
	extern real 	ForceZ;			// Volumic force in z-direction
	extern real 	ModelConstant;		// Constant used in the turbulence model
	extern real	ThermalConduction;	// Dimensionless thermal conduction
	extern real	ConstantForce;		// False = adapt force to maintain constant energy
/*
______________________________________________________________________________________________

6/ Chemistry
______________________________________________________________________________________________

*/
	extern real 	Alpha; 			// Temperature ratio
	extern real 	Ze;			// Equivalent to Beta
	extern real 	Le;			// Lewis number
	extern real 	Sigma;			// Radiation coefficient			
/*
______________________________________________________________________________________________

Integral values
______________________________________________________________________________________________

*/
	extern real	ErrorMax;			// Error Max on the grid
	extern real	ErrorGlobalMax; 		// Global error max on space and time
	extern real 	ErrorMid;			// Mean error on the grid
	extern real 	ErrorGlobalMid;			// Global mean error on space and time
	extern real 	ErrorL2;			// L2 error on the grid
	extern real 	ErrorGlobalL2;			// Global L2 error on space and time
	extern int	ErrorNb;			// Number of points for the computation of the mean error
	extern int 	ErrorGlobalNb;			// Number of points for the computation of the global mean error

	extern real 	GlobalMomentum;			// Global momentum (only for EquationType = 1 or 2)
	extern real 	GlobalEnergy;			// Global energy (only for EquationType = 1 or 2)
	extern real 	GlobalEnstrophy;        	// Global enstrophy (only for EquationType = 6)
	extern real 	ExactMomentum;			// Global momentum for the exact solution (only for EquationType = 1 or 2)
	extern real 	ExactEnergy;			// Global energy for the exact solution (only for EquationType = 1 or 2)
	extern real 	GlobalMomentumOld; 		// Old global momentum (only for EquationType = 6)
	extern real	GlobalVolume;			// Global volume
	
	extern real	EigenvalueMax;      		// Maximal eigenvalue
	extern Vector 	QuantityMax;			// Vector containing the maximal quantities
	extern Vector 	QuantityAverage;		// Vector containing the average quantities
	
	// ---Parameters used in the Runge-Kutta-Fehlberg (RKF) method to adapt time step ---

	extern real 	RKFError;			// Maximum of the relative errors between RK2 and RK3 (only if TimeAdaptivity=true)
	extern real	RKFAccuracyFactor;		// Desired value of RKFError (only if TimeAdaptivity=true)
	extern real	RKFSafetyFactor;		// Safety factor for the computation of the time step (only if TimeAdaptivity=true)

	// --- Global variables for flame computations (EquationType = 3,4,5)
	
	extern real 	FlameVelocity;			// Flame velocity
	extern real	GlobalReactionRate; 		// Global reaction rate
	extern real 	AverageRadius;			// Average radius of the flame ball
	extern real 	PreviousAverageRadius; 		// Previous average radius
	extern real 	PreviousAverageRadius2; 	// Previous average radius
	extern real	XCenter[4];			// Coordinate of the center of the flame ball
	extern real	ReactionRateMax;		// Maximum of the reaction rate
		
	// ---Integral values for Navier-Stokes computations ---
	
	extern real 	IntVorticity;			// Integral of the modulus of the vorticity
	extern real 	IntDensity;			// Integral of the density
	extern Vector 	IntMomentum;			//  Integral of the modulus of the momentum
	extern real 	IntEnergy;			// Integral of the energy
	extern real 	BaroclinicEffect;		// Intensity of the baroclinic effects
			
/*
______________________________________________________________________________________________

Others
______________________________________________________________________________________________

*/
	extern real pi;
	extern bool debug;			// true = check if tree is graded
	extern bool FluxCorrection;		// true = conservative flux correction
	extern real CarmenVersion;		// Version release
	extern FILE* GlobalFile;		// Global file
	extern int PostProcessing;		// 1 = Gnuplot, 2 = Data Explorer, 3 = Tecplot
  	extern bool DataIsBinary;     		// true = write data in binary format, false = write data in ASCII format
  	extern bool ZipData;			// true = zip data files
	extern bool WriteAsPoints;		// true = write data as point-values, false = write data as cell-averages


/*
______________________________________________________________________________________________

parallel
______________________________________________________________________________________________

*/


extern int AllTaskScaleNb;  //Global Scale number
extern int CPUScales;       //CPUScales=log2(Number of processors)
extern int size;            //Number of processors
extern int rank;            //Current CPU
extern real AllXMin[4];     //Global domain parameters
extern real AllXMax[4];

extern int NeighbourNb;     //Important parameter: The deep of the inter-CPU domain overlapping

extern int coords[3];       //Current CPU coodinates in the virtual CPU processors cart
extern int CartDims[3];

extern int CellElementsNb;  
extern int MaxCellElementsNb;

extern int one_D;          
extern int two_D;

extern int MPISendType;    //Type of calling MPI communication functions. See Parameters.cpp for the more information
extern int MPIRecvType;

extern Timer CommTimer;     //Communication timer for perfomance analyse

            //ranks of the neighbours processors
extern int rank_il;        //axis X, direction - low
extern int rank_iu;        //axis X, direction - high
extern int rank_jl;        //axis Y, direction - low
extern int rank_ju;        //etc..
extern int rank_kl;
extern int rank_ku;

extern int WhatSend;       //see the manual for the more information about this parameters
extern int SendD;          
extern int SendGrad;
extern int SendQ;
extern int SendQs;
extern int SendX;
extern int SenddX;


#if defined PARMPI              //MPI Data Types
extern 	MPI_Datatype MPItypeSiL;
extern 	MPI_Datatype MPItypeSiU;
extern 	MPI_Datatype MPItypeRiL;
extern 	MPI_Datatype MPItypeRiU;
extern 	MPI_Datatype MPItypeSjL;
extern 	MPI_Datatype MPItypeSjU;
extern 	MPI_Datatype MPItypeRjL;
extern 	MPI_Datatype MPItypeRjU;
extern	MPI_Datatype MPItypeSkL;
extern  MPI_Datatype MPItypeSkU;
extern  MPI_Datatype MPItypeRkL;
extern  MPI_Datatype MPItypeRkU;

extern MPI_Aint *disp;         //MPI Types buffers
extern int *blocklen;

extern	MPI_Comm comm_cart;   //MPI Topology

extern real *MPIbuffer;       //Additional Buffer for the MPI

extern MPI_Status st[12];     //Status of MPI non-blocking operations
extern	MPI_Request req[12];
#endif

extern char BackupName[255];    //Name of the carmen back up files (cart coordinates of the current CPU + standart name "carmen.bak")

// End file
