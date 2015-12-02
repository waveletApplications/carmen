/***************************************************************************
                          FineMesh.h  -  description
                             -------------------
    begin                : Wed Jun 13 2001
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

/***************************************************************************
 * An object FineMesh is a regular fine mesh, used for
 * finite volume computations.
 * It is not used for multiresolution computations.
 *
 * It contains an array of cells <i>*MeshCell</i>.
 *
 * NOTE: for reasons of simplicity, only periodic and Neuman boundary conditions have been
 * implemented.
 */
                    
class FineMesh
{
/*
______________________________________________________________________________________________

PUBLIC METHODS

	Constructor and distructor
______________________________________________________________________________________________

*/
public:

/***************************************************************************
 * Generates a regular fine mesh containing <i>2**(Dimension*ScaleNb)</i>
 * cells. The cell-averages are initialized from the initial condition
 * contained in file <i>carmen.ini</i>.
 *
 */
	FineMesh();

/***************************************************************************
 * Destroys the regular fine mesh.
 *
 */
	~FineMesh();
/*
______________________________________________________________________________________________

	Time evolution procedures
______________________________________________________________________________________________

*/

/***************************************************************************
 * Stores cell-average values into temporary cell-average values.
 *
 */
	void store();                 			

/***************************************************************************
 * Stores gradient values into temporary gradient values.
 *
 */	
	void storeGrad();      

	
/***************************************************************************
 * Computes the divergence vector with the space discretization scheme.
 *
 */
	void computeDivergence(int);


//Computes one Cell Divirgence
	void computeDivergence_cell(int);


/***************************************************************************
 * Computes one Runge-Kutta step.
 *
 */
	void RungeKutta_cell(int);

	//each cells
	void RungeKutta(int);

  

/***************************************************************************
 * Computes integral values like e.g. flame velocity, global error, etc.
 *
 */
	void computeIntegral();


/***************************************************************************
 * Computes velocity gradient (only for Navier-Stokes).
 *
 */
  //one cell
  void computeGradient(int);

  //each cells
  void computeGradient_cell(int);
 
/***************************************************************************
 * Computes the time-average value in every cell.
 */
	void computeTimeAverage();

/***************************************************************************
 * Checks if the computation is numerically unstable, i.e. if one of the
 * cell-averages is overflow. In case of numerical instability, the computation is
 * stopped and a message appears.
 *
 */
	void checkStability() const;

/*
______________________________________________________________________________________________

	Output procedures
______________________________________________________________________________________________

*/

/***************************************************************************
 * Write header for Data Explorer into file <i>FileName</i>.
 *
 */
	void writeHeader(const char* FileName) const;

/***************************************************************************
 * Write cell-averages for Data Explorer into file <i>FileName</i>.
 *
 */
	void writeAverage(const char* FileName);

/***************************************************************************
 * Write time-averages into file <i>FileName</i>.
 *
 */
	void writeTimeAverage(const char* FileName) const;
/*
______________________________________________________________________________________________

	Backup-restore procedure (to restart a computation)
______________________________________________________________________________________________

*/
/***************************************************************************
 * Backs up the tree structure and the cell-averages into a file <i>carmen.bak</i>.
 * In further computations, the data can be recovered using <b>Restore()</b>.
 *
 */
	void backup();

/***************************************************************************
 * Restores the tree structure and the cell-averages from the file <i>carmen.bak</i>.
 * This file was created by the method <b>Backup()</b>.
 *
 */
	void restore();


//Parallel
	Cell ***Neighbour_iL,***Neighbour_iU,***Neighbour_jL,***Neighbour_jU,***Neighbour_kL,***Neighbour_kU;
  
/*
______________________________________________________________________________________________

PRIVATE METHODS AND VARIABLES
______________________________________________________________________________________________

*/		

//parallel modification

/***************************************************************************
 * Returns pointer to the cell located at <i>i, j, k</i>.
 *
 */
public:
  Cell *MeshCell;

private:
	Cell* cell(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns pointer to the <i>n</i>-th cell.
 *
 */
	inline Cell* cell(const int n) const;

	
/***************************************************************************
 * Compute coarse grid data and store fine grid data in temporary
 *
 */
	void coarsen();
			

/***************************************************************************
 * Restore fine grid data from temporary. Inverse of : coarsen()
 *
 */
	void refine();
		

/***************************************************************************
 * Array of cells
 *
 */
	TimeAverageGrid *MyTimeAverageGrid;
		
};


/*
______________________________________________________________________________________________

INLINE FUNCTIONS
______________________________________________________________________________________________

*/

inline Cell* FineMesh::cell(const int n) const
{
	return (MeshCell +n);
}



