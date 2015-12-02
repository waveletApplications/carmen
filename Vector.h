/***************************************************************************
                          Vector.h  -  description
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
 * Standard class for every vector in Carmen.
 *
 * It contains the following data:<br>
 * * the dimension of the vector <i>Columns</i> ;<br>
 * * the array of reals <i>*U</i>.
 *
 */

class Vector
{

/*
______________________________________________________________________________________________

	Constructor and distructor
______________________________________________________________________________________________

*/

/***************************************************************************
 * Generates a 1D vector equal to zero.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V;
 *
 */
public:
	Vector();


/***************************************************************************
 * Generates a vector of dimension <i>n</i>, each component being equal to zero.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(4);
 *
 */
	Vector(const int n);

/***************************************************************************
 * Generates the 2D vector (<i>x,y</i>).
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(0.,1.);
 *
 */
	Vector(const real x, const real y);

/***************************************************************************
 * Generates the 3D vector (<i>x,y,z</i>).
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(0.,1.,0.);
 *
 */
	Vector(const real x, const real y, const real z);

/***************************************************************************
 * Generates a copy of the vector <i>V</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(0.,1.,0.);<br>
 * Vector W(V);
 *
 */
	Vector(const Vector& V);

/***************************************************************************
 * Deallocate memory of the vector.
 *
 */
	~Vector();
/*
______________________________________________________________________________________________

	Set and get
______________________________________________________________________________________________

*/

/***************************************************************************
 * Sets the component <i>n</i> to value <i>a</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(2);<br>
 * real x = 3.;<br>
 * real y = 1.;<br>
 * V.setValue(1,x);<br>
 * V.setValue(2,y);
 *
 */
	inline void setValue(const int n, const real a);

/***************************************************************************
 * Sets all the components to zero.
 *
 */
	void setZero();

/***************************************************************************
 * Sets the dimension of the vector to <i>n</i> and reset values to zero.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V;<br>
 * ...    <br>
 * V.setDimension(3); <br>
 */
	void setDimension(const int n);

/***************************************************************************
 * Returns the value of the component <i>n</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(2);<br>
 * real x;<br>
 * real y;<br>
 * ...    <br>
 * x = V.value(1); <br>
 * y = V.value(2); <br>
 *
 */
	inline real value(const int n) const;

/***************************************************************************
 * Returns the dimension of the vector.
 *
 */
	inline int dimension() const;


/*
______________________________________________________________________________________________

	Operators
______________________________________________________________________________________________

*/
/***************************************************************************
 * Compares the current vector to a vector <i>V</i> and returns true if they are equal.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(2);<br>
 * Vector W(2);<br>
 * real x; <br>
 * ...    <br>
 * if (V == W) <br>
 *    x = V.value(1); <br>
 */
	bool operator== (const Vector& V ) const;

/***************************************************************************
 * Set the current vector to the dimension and the value of <i>V</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W; <br>
 * ...    <br>
 * W = V; <br>
 *
 */
	void operator= (const Vector& V );

/***************************************************************************
 * Adds <i>V</i> to the current vector.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W(0.,-1.,2.); <br>
 * ...    <br>
 * W += V; <br>
 *
 */
	void operator+=(const Vector& V );

/***************************************************************************
 * Returns the addition of the current vector and <i>V</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W(0.,-1.,2.); <br>
 * Vector U; <br>
 * ...    <br>
 * U = V + W;
 *
 */
	Vector operator+ (const Vector& V ) const;


/***************************************************************************
 * Subtracts <i>V</i> to the current vector.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W(0.,-1.,2.); <br>
 * ...    <br>
 * W -= V; <br>
 *
 */
	void operator-=(const Vector& V );

/***************************************************************************
 * Returns the difference between the current vector and <i>V</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W(0.,-1.,2.); <br>
 * Vector U; <br>
 * ...    <br>
 * U = V - W; <br>
 *
 */
	Vector operator- (const Vector& V ) const;

/***************************************************************************
 * Returns the opposite of the current vector.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W; <br>
 * ...    <br>
 * W = -V; <br>
 *
 */
   Vector operator- () const;

/***************************************************************************
 * Multiplies the current vector by a real <i>a</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * real x = 2.; <br>
 * ...    <br>
 * V *= x; <br>
 *
 */
	void operator*=(const real a );

/***************************************************************************
 * Returns the product of the current vector and a real <i>a</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W;  <br>
 * real x = 2.; <br>
 * ...    <br>
 * W = V*x; <br>
 *
 * The operation W = x*V can also be done. See <b>operator*(const real a, const Vector& V)</b>.
 *
 */
	Vector operator* (const real a ) const;

/***************************************************************************
 * Divides the current vector by a real <i>a</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * real x = 2.; <br>
 * ...    <br>
 * V /= x; <br>
 *
 */
	void operator/=(const real a );

/***************************************************************************
 * Returns the division of the current vector by a real <i>a</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W; <br>
 * real x = 2.; <br>
 * ...    <br>
 * W = V / x; <br>
 *
 */
	Vector operator/ (const real a ) const;

/***************************************************************************
 * Returns the dot product of the current vector and <i>V</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W(1., 2., 1.);<br>
 * real x; <br>
 * ...    <br>
 * x = V*W; <br>
 *
 */
	real operator* (const Vector& V ) const;
	
/***************************************************************************
 * Returns the term-by-term product of the current vector and <i>V</i>.
 *
 */
	Vector operator| (const Vector& V) const;
	
	
/***************************************************************************
 * Returns the vectorial product of the current vector and <i>V</i>.
 *
 */
 	Vector operator^ (const Vector& V) const;
	
	
/***************************************************************************
 * Returns true if one of the components of the current vector is not a number.
 *
 */

	bool isNaN() const;

/*
______________________________________________________________________________________________

PRIVATE VARIABLES
______________________________________________________________________________________________

*/
//parallel modification
//private:
public:
	int  Columns; 	// Number of columns
	/* real *U;	// Components */
	real U[5];	// Components
};
/*
______________________________________________________________________________________________

EXTERNAL FUNCTIONS
______________________________________________________________________________________________

*/
/***************************************************************************
 * Returns the product of the current vector and a real <i>a</i>.
 *
 * Example :
 *
 * <tt>#include "Vector.h"
 *
 * <tt>Vector V(1.,0.,0.);<br>
 * Vector W;  <br>
 * real x = 2.; <br>
 * ...    <br>
 * W = x*V; <br>
 *
 * The operation W = V*x can also be done. See <b>Vector Vector::operator*(const real a) const</b>.
 *
 */
Vector operator* (const real a, const Vector& V);

/***************************************************************************
 * Returns the absolute value term by term of the vector.
 *
 */
Vector abs (const Vector& V);

/***************************************************************************
 * Returns the dimension of the vector. Similar to <b>int Vector::dimension()</b>.
 *
 */
int dim (const Vector& V);

/***************************************************************************
 * Returns the L1-norm of the vector.
 *
 */
real N1(const Vector& V);

/***************************************************************************
 * Returns the L2-norm of the vector.
 *
 */
real N2(const Vector& V);

/***************************************************************************
 * Returns the Max-norm of the vector.
 *
 */
real NMax(const Vector& V);

/***************************************************************************
 * Writes the components of the vector <i>V</i> on screen.
 *
 */

ostream& operator<<(ostream& out, const Vector& V);

/*
______________________________________________________________________________________________

INLINE FUNCTIONS
______________________________________________________________________________________________

*/
inline int Vector::dimension() const
{
	return Columns;
}

/*
______________________________________________________________________________________________

*/

inline void Vector::setValue(const int n, const real a)
{

#ifdef DEBUG
	if ( n <= 0 || n > Columns)
	{
		cout << "Vector.cpp: In method `void Vector::setValue(int, real)':\n";
		cout << "Vector.cpp: first argument out of range\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

		*(U+n-1) = a;
}
/*
______________________________________________________________________________________________

*/
inline real Vector::value(const int n) const
{

#ifdef DEBUG

	if ( n <= 0 || n > Columns)
	{
		cout << "Vector.cpp: In method `void Vector::value(int)':\n";
		cout << "Vector.cpp: argument out of range\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	return *(U+n-1);
}

