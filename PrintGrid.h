/***************************************************************************
                          PrintGrid.h  -  description
                             -------------------
    begin                : Tue Mar 19 2002
    copyright            : (C) 2002 by Olivier Roussel & Alexei Tsigulin
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
 * An object PrintGrid is a special regular grid created to write
 * tree-structured data into an output file.
 *
 * It contains the following data: <br>
 * * the scale number of the grid <i>LocalScaleNb</i> ;<br>
 * * the current number of elements used in the grid <i>ElementNb</i> ;<br>
 * * the array of cell-average values <i>*Q</i> ;<br>
 * * an array of temporary cell-average values <i>*Qt</i>.<br>
 *
 * To write tree-structured data into a regular fine grid, one starts
 * with the grid of level 0 and one stops at the level <i>L</i>.
 * For a given grid of level <i>l</i> and a given element <i>i, j, k</i> of this grid, if the node
 * of the tree corresponding to the element is
 * a leaf, the value is replaced by the one of the node, else it is predicted from parents.
 *
 * Such grid does not contain any cell information, in order to reduce memory requirements.
 *
 * For Navier-Stokes equations, derived quantities can be obtained, such as pressure, temperature,
 * vorticity, etc.
 *
 */

class PrintGrid
{
/*
______________________________________________________________________________________________

DESCRIPTION

	This is a grid of level L used to print data into file.

______________________________________________________________________________________________

PUBLIC METHODS

	Constructor and distructor
______________________________________________________________________________________________

*/
/***************************************************************************
 * Generates a grid of 2**(<i>Dimension*L</i>) elements.
 */

public:
	PrintGrid(int L);

/***************************************************************************
 * Removes the grid.
 *
 */
	~PrintGrid();

/*
______________________________________________________________________________________________

	Set and get procedures
______________________________________________________________________________________________

*/
/***************************************************************************
 * Sets the cell-average vector located at <i>i, j, k</i> to <i>UserAverage</i>.
 *
 */
	void   setValue(const int i, const int j, const int k, const Vector& UserAverage);

/***************************************************************************
 * Returns the cell-average vector located at <i>i, j, k</i>.
 *
 */
	Vector value(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns the quantity <i>QuantityNo</i> of the cell-average vector located at <i>i, j, k</i>.
 *
 */
	real value(const int i, const int j, const int k, const int QuantityNo) const;

/***************************************************************************
 * Returns the quantity <i>QuantityNo</i> of the cell-average vector located at <i>i, j, k</i>,
 * taking into account boundary conditions.
 *
 */
	real cellValue(const int i, const int j, const int k, const int QuantityNo) const;

/***************************************************************************
 * Returns the cell-average density located at <i>i, j, k</i>.
 *
 */
	inline real density(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns the cell-average pressure located at <i>i, j, k</i>.
 *
 */
	real pressure(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns the cell-average temperature located at <i>i, j, k</i>.
 *
 */
	real temperature(const int i, const int j, const int k) const;


/***************************************************************************
 * Returns the cell-average concentration of the limiting reactant, located at <i>i, j, k</i>.
 *
 */
	real concentration(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns the cell-average energy per unit of volume located at <i>i, j, k</i>.
 *
 */
	inline real energy(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns the <i>AxisNo</i>-th component of the cell-average velocity located at <i>i, j, k</i>.
 *
 */
	inline real velocity(const int i, const int j, const int k, const int AxisNo) const;

/***************************************************************************
 * Returns the cell-average velocity located at <i>i, j, k</i>.
 *
 */
	Vector	velocity(const int i, const int j, const int k) const;

/***************************************************************************
 * Returns 0 in 1D, the scalar vorticity in 2D, the norm of the cell-average vorticity in 3D, located at <i>i, j, k</i>.
 *
 */
	real vorticity(const int i, const int j, const int k) const;
/*
______________________________________________________________________________________________

	Refresh and predict procedure
______________________________________________________________________________________________

*/

/***************************************************************************
 * Stores the cell-average values of the current grid into temporary values,
 * in order to compute cell-averages in the next finer grid.
 *
 */
	void refresh();

/***************************************************************************
 * Predicts the cell-average values of the current grid from the values stored in the
 * temporary ones.
 *
 */
	void predict(const int l, const int i, const int j, const int k);

/***************************************************************************
 * Transform cell-average values into point values.
 *
 */
	void computePointValue();
/*
______________________________________________________________________________________________

PRIVATE VARIABLES AND METHODS
______________________________________________________________________________________________

*/
		
/***************************************************************************
 * Returns temporary cell-average values at the level <i>l</i> and at
 * the position <i>i, j, k</i>.
 *
 */
private:
  Vector tempValue(const int l, const int i, const int j, const int k) const;

	int localScaleNb;	// Scale number for the printed grid
	int elementNb;		// Number of elements
	Vector *Q;		// Quantities vector
	Vector *Qt;		// Temporary quantities vector
};


/*
______________________________________________________________________________________________

INLINE FUNCTIONS
______________________________________________________________________________________________

*/

inline real PrintGrid::density(const int i, const int j, const int k) const
{
	return value(i,j,k,1);
}

/*
______________________________________________________________________________________________

*/

inline real PrintGrid::energy(const int i, const int j, const int k) const
{
	return value(i,j,k,Dimension+2);
}

/*
______________________________________________________________________________________________

*/

inline real PrintGrid::velocity(const int i, const int j, const int k, const int AxisNo) const
{
	return cellValue(i,j,k,AxisNo+1)/cellValue(i,j,k,1);
}

