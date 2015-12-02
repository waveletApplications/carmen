/***************************************************************************
                          Carmen.h  -  description
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



#include "carmen.pre"


#if defined PARMPI
#include <mpi.h>
#endif

#include "PreProcessor.h"


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;


#include "Vector.h"
#include "Matrix.h"
#include "Timer.h"
#include "Parameters.h"
#include "PrintGrid.h"
#include "TimeAverageGrid.h"
#include "Cell.h"
#include "Node.h"
#include "FineMesh.h"



#define Max(x,y) (((x) > (y)) ? (x):(y))
#define Max3(x,y,z) (Max((x),Max((y),(z))))
#define Min(x,y) (((x) < (y)) ? (x):(y))
#define Min3(x,y,z) (Min((x),Min((y),(z))))
#define power2(x)  ((x)*(x))
#define power3(x)  ((x)*(x)*(x))
#define Abs(x) ( ((x) < 0)? -(x):(x) ) 

/*
______________________________________________________________________________________________

External functions
______________________________________________________________________________________________

*/

/***************************************************************************
 * Adapts time step when required.
 *
 */
void AdaptTimeStep();

/***************************************************************************
 * Adds a L1-error <i>NewError</i> to <i>N-1</i> L1-errors stored in <i>TotalError</i>.
 *
 */

real AddErrorL1(int N, real TotalError, real NewError);


/***************************************************************************
 * Adds a L2-error <i>NewError</i> to <i>N-1</i> L2-errors stored in <i>TotalError</i>.
 *
 */

real AddErrorL2(int N, real TotalError, real NewError);

/***************************************************************************
 * Returns the analytical solution of the 1D convection-diffusion equation
 * for an initial step.
 *
 */
real Advection(real x, real t);

/***************************************************************************
 * Returns the cell-average analytical solution of a 1D equation.
 *
 */
real AnalyticAverage(real x, real dx, real t);

/***************************************************************************
 * Returns the cell-average analytical solution of a 2D equation.
 *
 */
real AnalyticAverage(real x, real dx, real y, real dy, real t);

/***************************************************************************
 * Returns the cell-average analytical solution of a 3D equation.
 *
 */
real AnalyticAverage(real x, real dx, real y, real dy, real z, real dz, real t);


/***************************************************************************
 * Returns the position of <I>i</I> taking into account the boundary conditions
 * in the direction <I>AxisNo</I>. The number of points in this direction is <I>N</I>. <BR>
 * Example: for AxisNo=1 and for N=256, i must be between 0 and 255.
 * If i=-1, the function returns 255 for periodic boundary conditions and 0 for Neumann
 * boundary conditions.
 *
 */
int BC(int i, int AxisNo, int N=(1<<ScaleNb));

/***************************************************************************
 * Returns the boundary region type at the position <I>X=(x,y,z)</I>.<BR>
 * The returned value correspond to:
 * 0 = Fluid (not in the boundary)
 * 1 = Inflow
 * 2 = Outflow
 * 3 = Solid with free-slip walls
 * 4 = Solid with adiabatic walls
 * 5 = Solid with isothermal walls
 */
int BoundaryRegion(const Vector& X);

/***************************************************************************
 * Stores the tree structure and data in order to restart a multiresolution computation.
 + <i>Root</i> denotes the pointer to the first node of the tree structure.
 *
 */
void Backup(Node* Root);

/***************************************************************************
 * Stores the data countained in a regular mesh <i>Root</i> in order to restart a
 * finite volume computation.
 *
 */
void Backup(FineMesh* Root);

/***************************************************************************
 * Returns the analytical solution of the 1D diffusive Burgers equation
 * for an initial step.
 *
 */
real Burgers(real x, real t);

/***************************************************************************
 * Returns the time required by a finite volume computation using <i>iterations</i>
 * iterations and <i>scales</i> scales. It is use to estimate the CPU time compression.
 *
 */
double CPUTimeRef(int iterations, int scales);

/***************************************************************************
 * Returns the computed tolerance at the scale <i>ScaleNo</i>, either using
 * Harten or Donoho thresholding (if <i>CVS</i>=true).
 *
 */
real ComputedTolerance(const int ScaleNo);


/***************************************************************************
 * Returns the analytical solution of the 2D Navier-Stokes equations
 * for a Hamel-Oseen vortex flow.
 *
 */
real CurlVelocity(real x, real y, real t, int AxisNo=1);

/***************************************************************************
 * Returns the number of digits of the integer <i>arg</i>.
 *
 */
int DigitNumber(int arg);

/***************************************************************************
 * Writes in binary or ASCII mode the real number <i>arg</i> into the file <i>f</i>
 * with the format <i>format</i>.
 * The global parameter <i>DataIsBinary</i> determines this choice.
 *
 */
int FileWrite(FILE *f, const char* format, real arg);

/***************************************************************************
 * Returns the flux at the interface between <i>Cell2</i> and <i>Cell3</i>.
 * Here a 4-point space scheme is used.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector Flux(Cell& Cell1, Cell& Cell2, Cell& Cell3, Cell& Cell4, int AxisNo);

/***************************************************************************
 * Transform the 4 cells of the flux <i>Cell1</i>, <i>Cell2</i>, <i>Cell3</i>, <i>Cell4</i> into
 * <i>C1</i>, <i>C2</i>, <i>C3</i>, <i>C4</i>, to take into account boundary conditions (used in Flux.cpp).
 */
void GetBoundaryCells(Cell& Cell1, Cell& Cell2, Cell& Cell3, Cell& Cell4,
			Cell& C1, Cell& C2, Cell& C3, Cell& C4, int AxisNo);

/***************************************************************************
 * Returns the initial condition in (<i>x, y, z</i>) form the one defined in <i>carmen.ini</i>.
 *
 */
Vector InitAverage(real x, real y=0., real z=0.);

/***************************************************************************
 * Inits parameters from file <i>carmen.par</i>. If a parameter is not mentioned
 * in this file, the default value is used.
 *
 */
void InitParameters();

/***************************************************************************
 * Inits time step and all the parameters which depend on it.
 *
 */
void InitTimeStep();

/***************************************************************************
 * Inits tree structure from initial condition, starting form the node <i>Root</i>.
 * Only for multiresolution computations.
 *
 */
void InitTree(Node* Root);

/***************************************************************************
 * Returns the value of the slope limiter between the slopes <i>u</i> and <i>v</i>.
 *
 */
Vector Limiter(const Vector u, const Vector v);
real Limiter(const real x);

/***************************************************************************
 * Returns the minimum in module of <i>a</i> and <i>b</i>.
 *
 */
real MinAbs(const real a, const real b);

/***************************************************************************
 * Returns the analytical solution of the 1D convection-diffusion equation
 * for an initial Gaussian bump.
 *
 */
real MoveGauss(const real x, const real t);

/***************************************************************************
 * Returns the analytical solution of the 2D convection-diffusion equation
 * for an initial Gaussian bump.
 *
 */
real MoveGauss(const real x, const real y, const real t);

/***************************************************************************
 * Returns the analytical solution of the 3D convection-diffusion equation
 * for an initial Gaussian bump.
 *
 */
real MoveGauss(const real x, const real y, const real z, const real t);

/***************************************************************************
 * Returns the Max-norm of the vector where every quantity is divided by its
 * characteristic value.
 *
 */
real NormMaxQuantities(const Vector& V);

/***************************************************************************
 * Computes the performance of the computation and, for cluster computations, write it
 * into file <i>FileName</i>.
 *
 */
void Performance(const char* FileName);

/***************************************************************************
 * Writes the integral values, like e.g flame velocity, global error, into file <i>FileName</i>.
 *
 */
void PrintIntegral(const char* FileName);

/***************************************************************************
 * Returns the dimensionless heat loss due to radiation in function of
 * the dimensionless temperature <i>T</i> for a black-body radiation model.
 *
 */
real RadiationRate(const real T);

/***************************************************************************
 * Returns the dimensionless reaction rate for a one-step Arrhenius kinetics
 * in function of the dimensionless temperature <i>T</i> and the partial mass
 * of the limiting reactant <i>Y</i>.
 *
 */
real ReactionRate(const real T, const real Y);

/***************************************************************************
 * Refresh the tree structure, i.e. compute the cell-averages of the
 * internal nodes by projection and those of the virtual leaves by prediction.
 * The root node is <i>Root</i>. Only for multiresolution computations.
 *
 */
void RefreshTree(Node* Root);

/***************************************************************************
 * Remesh the tree structure after a time evolution. The root node is <i>Root</i>.
 * Only for multiresolution computations.
 *
 */
void Remesh(Node* Root);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using the AUSM+ scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeAUSM(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, const int AxisNo);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using a 4th order centered scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */

Vector SchemeAUSMDV(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, const int AxisNo);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using a 4th order centered scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */


Vector SchemeCentered(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, const int AxisNo);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using a second-order ENO scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeENO2(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, const int AxisNo);


/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using a thrid-order ENO scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeENO3(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, const int AxisNo);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using a 2-4 McCormack scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeMcCormack(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, const int AxisNo);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using an OSMP3 scheme with McCormack forward-backward expression.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeOsmp3Mc(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, int AxisNo);


/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using an OSMP5 scheme with McCormack forward-backward expression.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeOsmp5Mc(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, int AxisNo);

/***************************************************************************
 * Returns the inviscid flux for the Navier-Stokes equation using Roe's scheme.
 * The scheme uses four cells to estimate the flux at the interface.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 * <i>Cell1</i> and <i>Cell4</i> are the second neighbours on the left and right sides.
 *
 */
Vector SchemeRoe(const Cell& Cell1, const Cell& Cell2, const Cell& Cell3, const Cell& Cell4, int AxisNo);

/***************************************************************************
 * Writes on screen the estimation of total and remaining CPU times.
 * These informations are stored in the timer <i>arg</i>.
 *
 */
void ShowTime(Timer arg);

/***************************************************************************
 * Returns 1 if <i>a</i> is non-negative, -1 elsewhere.
 *
 */
int Sign(const real a);

/***************************************************************************
 * Returns the eddy-viscosity in the cell <i>UserCell</i> using Smagorinsky's model.
 *
 */
real Smagorinsky(const Cell& UserCell);

/***************************************************************************
 * Returns the source term in the cell <i>UserCell</i>.
 *
 */
Vector Source(Cell& UserCell);

/***************************************************************************
 * Returns a step (1 if x <0, 0 if x >0, 0.5 if x=0)
 *
 */
real Step(real x);

/***************************************************************************
 * Returns the dimensionless viscosity in function of the dimensionless temperature
 * <i>T</i> following Sutherland's law
 *
 */
real Sutherland(real T);

/***************************************************************************
 * Computes a time evolution on the regular fine mesh <i>Root</i>.
 * Only for finite volume computations.
 *
 */
void TimeEvolution(FineMesh* Root);

/***************************************************************************
 * Computes a time evolution on the tree structure, the root node being <i>Root</i>.
 * Only for multiresolution computations.
 *
 */
void TimeEvolution(Node* Root);

/***************************************************************************
 * Writes the current cel--averages of the fine mesh <i>Root</i> into file <i>AverageFileName</i>.
 * Only for finite volume computations.
 *
 */
void View(FineMesh* Root, const char* AverageFileName);

/***************************************************************************
 * Writes the data of the tree structure into files <i>TreeFileName</i> (tree structure),
 * <i>MeshFileName</i> (mesh) and <i>AverageFileName</i> (cell-averages).
 * The root node is <i>Root</i>. Only for multiresolution computations.
 *
 */
void View(Node* Root, const char* TreeFileName, const char* MeshFileName, const char* AverageFileName);

/***************************************************************************
 * Same as previous for a fine mesh <i>Root</i>. Only for finite volume
 * computations.
 *
 */
void ViewEvery(FineMesh* Root, int arg);

/***************************************************************************
 * Writes into file the data of the tree structure at iteration <i>arg</i>.
 * The output file names are <i>AverageNNN.dat</i> and <i>MeshNNN.dat</i>, <i>NNN</i> being the
 * iteration in an accurate format.
 * The root node is <i>Root</i>. Only for multiresolution computations.
 *
 */
void ViewEvery(Node* Root, int arg);

/***************************************************************************
 * Same as previous for a fine mesh <i>Root</i>. Only for finite volume
 * computations.
 *
 */
void ViewIteration(FineMesh* Root);

/***************************************************************************
 * Writes into file the data of the tree structure from physical time <i>PrintTime1</i>
 * to physical time <i>PrintTime6</i>.
 * The output file names are <i>Average_N.dat</i> and <i>Mesh_N.dat</i>, <i>N</i> being between 1 and 6.
 * The root node is <i>Root</i>. Only for multiresolution computations.
 *
 */
void ViewIteration(Node* Root);

/***************************************************************************
 * Returns the viscous flux for the Navier-Stokes equation. This flux is
 * computed using a second-order centered scheme.
 * <i>Cell2</i> and <i>Cell3</i> are the first neighbours on the left and right sides.
 *
 */

 Vector ViscousFlux(const Cell& Cell1, const Cell& Cell2, int AxisNo);



//Parallel
void CreateMPIType(FineMesh *Root);
void CPUExchange(FineMesh *Root,int);
void FreeMPIType();
void CreateMPITopology();
void CreateMPILinks();
void ReduceIntegralValues();




