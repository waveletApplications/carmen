/***************************************************************************
                          Cell.h  -  description
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

/***************************************************************************
 * An object Cell contains all the informations of a cell,
 * for both multiresolution and finite volume computations.
 *
 * These data are:
 *
 * * the coordinates of the cell-center <i>X</i> ;<br>
 * * the cell size <i>dX</i> ;<br>
 * * the cell-average quantities <i>Q</i> ;<br>
 * * the temporary cell-averages <i>Qs</i> ;<br>
 * * the divergence vector <i>D</i> ;<br>
 * * the velocity gradient <i>Grad</i> (only for Navier-Stokes).
 *
 * For the Navier-Stokes equations, some methods enable to compute derived quantities:
 * pressure, temperature, velocity, vorticity, and energy.
 *
 */

class Cell
{
/*
______________________________________________________________________________________________

PUBLIC METHODS

	Constructor and distructor
______________________________________________________________________________________________

*/
public:

/***************************************************************************
 * Allocate memory for the vectors in function of
 * the dimension of the problem <i>Dimension</i> and the number of
 * physical quantities <i>QuantityNb</i>.
 *
 * The position of the cell-center <i>X</i> and the cell size <i>dX</i> have the
 * dimension <i>Dimension</i>.
 * The velocity gradient <i>Grad</i> has the dimension <i>Dimension*Dimension</i>. The vectors
 * containing the cell-averages <i>Q</i>, the temporary cell-averages <i>Qs</i> and
 * the divergence <i>D</i> have the dimension <i>QuantityNb</i>.
 *
 */
	Cell ();

/***************************************************************************
 * Deallocate memory.
 *
 */
	~Cell ();
/*
______________________________________________________________________________________________

	Operators
______________________________________________________________________________________________

*/
void operator=(const Cell& C);

/*
______________________________________________________________________________________________

	Set procedures
______________________________________________________________________________________________

*/
/***************************************************************************
 * Sets the size of the cell in the direction <i>AxisNo</i> to <i>UserSize</i>.
 *
 * Example:
 *
 * <tt>#include "Carmen.h"
 *
 * <tt>Dimension = 2;<br>
 * real dx=1., dy=2.;<br>
 * Cell c;           <br>
 * c.setSize(1, dx); <br>
 * c.setSize(2, dy);
 *
 */
	inline void setSize(const int AxisNo, const real UserSize);

/***************************************************************************
 * Sets the size of the cell in every direction to the vector <i>UserSize</i> .
 *
 * Example:
 *
 * <tt>
 * #include "Carmen.h"
 *
 * <tt>Dimension = 2;<br>
 * Vector W(1.,2.);<br>
 * Cell c;<br>
 * c.setSize(W);<br>
 * </tt>
 *
 */
	inline void setSize(const Vector& UserSize);

/***************************************************************************
 * Sets the coordinate of the cell-center in the direction <i>AxisNo</i> to
 * <i>UserX</i>.
 *
 * Example:
 *
 *  <tt>
 * #include "Carmen.h"
 *
 * <tt>Dimension = 3;<br>
 * real x=1., y=0., z=0.;<br>
 * Cell c;<br>
 * c.setCenter(1,x);<br>
 * c.setCenter(2,y);<br>
 * c.setCenter(3,z);<br>
 * </tt>
 *
 */
	inline void setCenter(const int AxisNo, const real UserX);

/***************************************************************************
 * Sets the position of the cell-center to the vector <i>UserX</i>.
 *
 * Example:
 *
 *  <tt>
 * #include "Carmen.h"
 *
 * <tt>Dimension = 3;<br>
 * Vector V(1.,0.,0.);<br>
 * Cell c;<br>
 * c.setCenter(V);<br>
 * </tt>
 *
 */
	inline void setCenter(const Vector& UserX);

/***************************************************************************
 * Sets the cell-average of the quantity <i>QuantityNo</i> to <i>UserAverage</i>.
 *
 * Example:
 *
 * <tt>
 * #include "Carmen.h"
 *
 * <tt>QuantityNb = 2;<br>
 * real T=0.,Y=1.;<br>
 * Cell c;<br>
 * c.setAverage(1,T);<br>
 * c.setAverage(2,Y);<br>
 * </tt>
 *
 */
	inline void setAverage(const int QuantityNo, const real UserAverage);

/***************************************************************************
 * Sets all the cell-average quantities to the vector <i>UserAverage</i>.
 *
 * Example:
 *
 *  <tt>
 * #include "Carmen.h"
 *
 * <tt>QuantityNb = 2;<br>
 * Vector Q(1.,0.);<br>
 * Cell c;<br>
 * c.setAverage(Q);<br>
 * </tt>
 *
 */
	inline void setAverage(const Vector& UserAverage);

/***************************************************************************
 * Sets all the cell-average values to zero.
 *
 */
	inline void setAverageZero();

 /***************************************************************************
 * Identical to <a href="#ref7">setAverage</a> (int QuantityNo, real UserAverage),
 * but for the vector of the temporary cell-average values.
 *
 */
	inline void setTempAverage(const int QuantityNo, const real UserAverage);

 /***************************************************************************
 * Identical to
 * void <a href="#ref8">setAverage</a> (const Vector& UserAverage),
 * but for the vector of the temporary cell-average values.
 *
 */
	inline void setTempAverage(const Vector& UserAverage);
	

/***************************************************************************
 * Sets all the temmporary cell-average values to zero.
 *
 */
	inline void setTempAverageZero();

 /***************************************************************************
 * Identical to <a href="#ref7">setAverage</a> (int QuantityNo, real UserAverage),
 * but for the vector of the cell-average values with low precision in the 
 * Runge-Kutta-Fehlberg method.
 *
 */
	inline void setLowAverage(const int QuantityNo, const real UserAverage);

 /***************************************************************************
 * Identical to
 * void <a href="#ref8">setAverage</a> (const Vector& UserAverage),
 * but for the vector of the cell-average values with low precision in the 
 * Runge-Kutta-Fehlberg method.
 *
 */
	inline void setLowAverage(const Vector& UserAverage);
	

 /***************************************************************************
 * Identical to <a href="#ref7">setAverage</a> (int QuantityNo, real UserAverage),
 * but for the vector of the old cell-average values.
 *
 */
	inline void setOldAverage(const int QuantityNo, const real UserAverage);

 /***************************************************************************
 * Identical to
 * void <a href="#ref8">setAverage</a> (const Vector& UserAverage),
 * but for the vector of the cell-average values.
 *
 */
	inline void setOldAverage(const Vector& UserAverage);
	
 /***************************************************************************
 * Identical to
 * void <a href="#ref7">setAverage</a> (int QuantityNo, real UserAverage),
 * but for the divergence vector.
 *
 */
	inline void setDivergence(const int QuantityNo, const real UserAverage);

 /***************************************************************************
 * Identical to
 * void <a href="#ref8">setAverage</a> (const Vector& UserAverage),
 * but for the divergence vector.
 */
	inline void setDivergence(const Vector& UserAverage);

/***************************************************************************
 * Sets all the components of the divergence vector to zero.
 *
 */
	inline void setDivergenceZero();

/***************************************************************************
 * Sets the component no. <i>i, j</i> of the quantity gradient to <i>UserAverage</i>.
 *
 * Example:
 *
 * <tt>
 * #include "Carmen.h"
 *
 * <tt>QuantityNb = 2;<br>
 * real Gxx=0., Gxy=1., Gyx=1., Gyy=0.;<br>
 * Cell c;<br>
 * c.setGradient(1,1,Gxx);<br>
 * c.setGradient(1,2,Gxy);<br>
 * c.setGradient(2,1,Gyx);<br>
 * c.setGradient(2,2,Gyy);<br>
 * </tt>
 *
 */
	inline void setGradient(const int i, const int j, const real UserAverage);
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	
	inline void setTempGradient(const int i, const int j, const real UserAverage);

/***************************************************************************
 * Sets the quantity gradient to the matrix <i>UserAverage</i>.
 *
 * Example:
 *
 * <tt>
 * #include "Carmen.h"
 *
 * <tt>QuantityNb = 5;<br>
 * Dimension = 3;<br>
 * Matrix G(3,5);<br>
 + Cell c<br>
 + ...<br>
 * c.SetGradient(G);<br>
 * </tt>
 *
 */
	inline void setGradient(const Matrix& UserAverage);
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	inline void setTempGradient(const Matrix& UserAverage);

/***************************************************************************
 * Sets all the components of the velocity gradient to zero.
 *
 */
	inline void setGradientZero();


/***************************************************************************
 * Returns true if the cell is inside the boundary
 *
 */
	bool isInsideBoundary() const;


/***************************************************************************
 * Returns true if the cell is inside the fluid
 *
 */
	bool isInFluid() const;
  
/*
______________________________________________________________________________________________

	Get procedures
______________________________________________________________________________________________

*/
	
/***************************************************************************
 * Returns the cell size in the direction <i>AxisNo</i>.
 *
 */
	inline real  size(const int AxisNo) const;

/***************************************************************************
 * Returns the vector cointaining the cell size in every direction.
 *
 */
  	inline Vector size() const;

/***************************************************************************
 * Returns the component no. <i>AxisNo</i> of the cell-center position
 *
 */
	inline real  center(const int AxisNo) const;	

/***************************************************************************
 * Returns the cell-center position vector.
 *
 */
	inline Vector center() const;
	
/***************************************************************************
 * Returns the component no. <i>QuantityNo</i> of the cell-average value.
 *
 */
	inline real  average(const int QuantityNo) const;

/***************************************************************************
 * Returns the cell-average value vector.
 *
 */
	inline Vector average() const;

/***************************************************************************
 * Returns the component no. <i>QuantityNo</i> of the temporary cell-average value.
 *
 */
	inline real  tempAverage(const int QuantityNo) const;

/***************************************************************************
 * Returns the temporary cell-average value vector.
 *
 */
	inline Vector tempAverage() const;		


/***************************************************************************
 * Returns the component no. <i>QuantityNo</i> of the cell-average value with low precision
 * in the Runge-Kutta-Fehlberg method.
 *
 */
	inline real  lowAverage(const int QuantityNo) const;

/***************************************************************************
 * Returns the cell-average value with low precision
 * in the Runge-Kutta-Fehlberg method.
 *
 */
	inline Vector lowAverage() const;		


/***************************************************************************
 * Returns the component no. <i>QuantityNo</i> of the old cell-average values.
 *
 */
	inline real  oldAverage(const int QuantityNo) const;

/***************************************************************************
 * Returns the old cell-average values.
 *
 */
	inline Vector oldAverage() const;		


/***************************************************************************
 * Returns the component no. <i>QuantityNo</i> of the divergence vector.
 *
 */
	inline real  divergence(const int QuantityNo) const;

/***************************************************************************
 * Returns the divergence vector.
 *
 */
	inline Vector divergence() const;		

/***************************************************************************
 * Returns the component no. <i>i, j</i> of the velocity gradient.
 *
 */
	inline real 	gradient(const int i, const int j) const;
 
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	inline real 	tempGradient(const int i, const int j) const;

/***************************************************************************
 * Returns the velocity gradient in matrix form.
 *
 */
  	inline Matrix gradient() const;
  
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
  	inline Matrix tempGradient() const;

/***************************************************************************
 * Returns the component no. <i>AxisNo</i> of the cell-average vorticity.
 *
 */
	real  vorticity(const int AxisNo) const;

/***************************************************************************
 * Returns the cell-average vorticity vector.
 *
 */
	Vector vorticity() const;

/***************************************************************************
 * Returns the cell-average density.
 *
 */
	inline real	density() const;
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	inline real	tempDensity() const;

/***************************************************************************
 * Returns the cell-average pressure.
 *
 */
	real	pressure() const;
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	real	tempPressure() const;

	
 /***************************************************************************
 * Identical to the previous one for the  values at the instant <i>n-1</i>.
 *
 */
	real	oldPressure() const;

/***************************************************************************
 * Returns the cell-average temperature.
 *
 */
	real  temperature() const;
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	real  tempTemperature() const;

/***************************************************************************
 * Returns the cell-average concentration of the limiting reactant.
 *
 */
	inline real  concentration() const;

/***************************************************************************
 * Returns the cell-average energy per unit of volume.
 *
 */
	inline real	energy() const;
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	inline real tempEnergy() const;

/***************************************************************************
 * Returns the component no. <i>AxisNo</i> of the cell-average velocity.
 *
 */
  	inline real  velocity(const int AxisNo) const;
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
 	inline real tempVelocity(const int AxisNo) const;

/***************************************************************************
 * Returns the cell-average velocity vector.
 *
 */
	Vector velocity() const;
	
 /***************************************************************************
 * Identical to the previous one for the temporary values.
 *
 */
	Vector tempVelocity() const;

/***************************************************************************
 * Returns the cell-average speed of sound.
 *
 */
	inline real	speedOfSound() const;


/***************************************************************************
 * Returns the entropy (p/rho^Gamma).
 *
 */
	inline real	entropy() const;

/***************************************************************************
 * Returns the volume of the cell (length in 1D, area in 2D, volume in 3D).
 *
 */
	real	volume() const;

/*
______________________________________________________________________________________________

	Test procedures
______________________________________________________________________________________________

*/

/***************************************************************************
 * Return true if one of the cell-average quantities is greater than the maximum.
 * This usually means the computation is numerically unstable.
 *
 */
	bool isOverflow() const;	
/*
______________________________________________________________________________________________

PRIVATE VARIABLES
______________________________________________________________________________________________

*/
//parallel modification		
//private:
public:

/***************************************************************************
 * Position of the cell center. Its dimension is <i>Dimension</i>.
 *
 */
Vector X;

/***************************************************************************
 * Cell size in each direction. Its dimension is <i>Dimension</i>.
 *
 */
Vector dX;

/***************************************************************************
 * Vector containing the cell-average values. Its dimension is <i>QuantityNb</i>.
 *
 */
Vector Q;		

/***************************************************************************
 * Temporary cell-average values. This vector is used to store the
 * intermediary value in a multi-step Runge-Kutta or McCormack time integration.
 * Its dimension is <i>QuantityNb</i>.
 *
 */
Vector Qs;


/***************************************************************************
 * Cell-average values with low precision. This vector is used to store the
 * cell-averages computed with the N-stage Runge-Kutta-Fehlberg N(N+1) 
 * method.
 * Its dimension is <i>QuantityNb</i>.
 *
 */
Vector Qlow;


/***************************************************************************
 * Cell-average values at the instant <i>n-1</i>. 
 * Its dimension is <i>QuantityNb</i>.
 *
 */
Vector Qold;

/***************************************************************************
 * Divergence vector. Its dimension is <i>QuantityNb</i>.
 *
 */
Vector D;			// Divergence vector

/***************************************************************************
 * Quantity gradient. Only necessary for a Navier-Stokes computation.
 * Its dimension is <i>Dimension*QuantityNb</i>.
 *
 */
Matrix Grad; 

/***************************************************************************
 * Temporary quantity gradient. Only necessary for a Navier-Stokes computation.
 * Its dimension is <i>Dimension*QuantityNb</i>.
 *
 */
Matrix Grads; 

};

/*
______________________________________________________________________________________________

INLINE METHODS
______________________________________________________________________________________________

*/
inline void Cell::setSize(const int AxisNo, const real UserSize)
{
	dX.setValue(AxisNo, UserSize);
}
/*
______________________________________________________________________________________________

*/

inline void Cell::setSize(const Vector& UserSize)
{
	dX = UserSize;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setCenter(const int AxisNo, const real UserX)
{
	X.setValue(AxisNo, UserX);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setCenter(const Vector& UserX)
{
	X = UserX;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setAverage(const int QuantityNo, const real UserAverage)
{
	Q.setValue(QuantityNo, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setAverage(const Vector& UserAverage)
{
	Q = UserAverage;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setAverageZero()
{
	Q.setZero();
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setTempAverage(const int QuantityNo, const real UserAverage)
{
	Qs.setValue(QuantityNo, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setTempAverage(const Vector& UserAverage)
{
	Qs = UserAverage;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setTempAverageZero()
{
	Qs.setZero();
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setLowAverage(const int QuantityNo, const real UserAverage)
{
	Qlow.setValue(QuantityNo, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setLowAverage(const Vector& UserAverage)
{
	Qlow = UserAverage;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setOldAverage(const int QuantityNo, const real UserAverage)
{
	Qold.setValue(QuantityNo, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setOldAverage(const Vector& UserAverage)
{
	Qold = UserAverage;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setDivergence(const int QuantityNo, const real UserAverage)
{
	D.setValue(QuantityNo, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setDivergence(const Vector& UserAverage)
{
	D = UserAverage;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setDivergenceZero()
{
	D.setZero();
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setGradient(const int i, const int j, const real UserAverage)
{
	Grad.setValue(i, j, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setTempGradient(const int i, const int j, const real UserAverage)
{
	Grads.setValue(i, j, UserAverage);
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setGradient(const Matrix& UserAverage)
{
	Grad = UserAverage;
}
/*
______________________________________________________________________________________________

*/
inline void Cell::setTempGradient(const Matrix& UserAverage)
{
	Grads = UserAverage;
}

/*
______________________________________________________________________________________________

*/
inline void Cell::setGradientZero()
{
	Grad.setZero();
}
/*
______________________________________________________________________________________________

	Get procedures
______________________________________________________________________________________________

*/
inline real Cell::size(const int AxisNo) const
{
	return dX.value(AxisNo);
}
/*
______________________________________________________________________________________________

*/

inline Vector Cell::size() const
{
	return dX;
}
/*
______________________________________________________________________________________________

*/
inline real Cell::center(const int AxisNo) const
{
	return X.value(AxisNo);
}
/*
______________________________________________________________________________________________

*/
inline Vector Cell::center() const
{
	return X;
}
/*
______________________________________________________________________________________________

*/
inline real Cell::average(const int QuantityNo) const
{
	return Q.value(QuantityNo);
}
/*
______________________________________________________________________________________________

*/
inline Vector Cell::average() const
{
	return Q;
}
/*
______________________________________________________________________________________________

*/
inline real Cell::tempAverage(const int QuantityNo) const
{
	return Qs.value(QuantityNo);
}
/*
______________________________________________________________________________________________

*/
inline Vector Cell::tempAverage() const
{
	return Qs;
}
/*
______________________________________________________________________________________________

*/
inline real Cell::lowAverage(const int QuantityNo) const
{
	return Qlow.value(QuantityNo);
}
/*
______________________________________________________________________________________________

*/
inline Vector Cell::lowAverage() const
{
	return Qlow;
}
/*
______________________________________________________________________________________________

*/
inline real Cell::oldAverage(const int QuantityNo) const
{
	return Qold.value(QuantityNo);
}
/*
______________________________________________________________________________________________

*/
inline Vector Cell::oldAverage() const
{
	return Qold;
}
/*
______________________________________________________________________________________________

*/
inline real Cell::divergence(const int QuantityNo) const
{
	return D.value(QuantityNo);
}
/*
______________________________________________________________________________________________

*/
inline Vector Cell::divergence() const
{
	return D;
}
/*
______________________________________________________________________________________________

*/

inline real Cell::gradient(const int i, const int j) const
{
	return Grad.value(i,j);
}

/*
______________________________________________________________________________________________

*/

inline real Cell::tempGradient(const int i, const int j) const
{
	return Grads.value(i,j);
}

/*
______________________________________________________________________________________________

*/

inline Matrix Cell::gradient() const
{
	return Grad;
}

/*
______________________________________________________________________________________________

*/
inline Matrix Cell::tempGradient() const
{
	return Grads;
}
/*
______________________________________________________________________________________________

*/
inline real	Cell::density() const
{
	return	Q.value(1);
}

/*
______________________________________________________________________________________________

*/
inline real	Cell::tempDensity() const
{
	return	Qs.value(1);
}


/*
______________________________________________________________________________________________

*/
inline real	Cell::energy() const
{
	return	Q.value(Dimension+2);
}

/*
______________________________________________________________________________________________

*/
inline real	Cell::tempEnergy() const
{
	return	Qs.value(Dimension+2);
}

/*
______________________________________________________________________________________________

*/
inline real	Cell::concentration() const
{
	// Warning: only for flame ball, flame front and flame-vortex equations
	
	return	Q.value(2);
}

/*
______________________________________________________________________________________________

*/
inline real  Cell::velocity(const int AxisNo) const
{
	return	Q.value(1+AxisNo)/Q.value(1);
}

/*
______________________________________________________________________________________________

*/
inline real  Cell::tempVelocity(const int AxisNo) const
{
	return	Qs.value(1+AxisNo)/Qs.value(1);
}
/*
______________________________________________________________________________________________

*/
inline real	Cell::speedOfSound() const
{
	return sqrt(Gamma*pressure()/density());
}
/*
______________________________________________________________________________________________

*/
inline real	Cell::entropy() const
{
	return pressure()*exp(-Gamma*log(density()));
}
