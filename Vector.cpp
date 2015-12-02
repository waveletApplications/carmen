/***************************************************************************
                          Vector.cpp  -  description
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

#include <math.h>
#include <iostream>
#include "PreProcessor.h"

using namespace std;

#include "Vector.h"


/*
______________________________________________________________________________________________

Constructors
______________________________________________________________________________________________

*/

Vector::Vector()
{
	Columns = 1;
	*U = 0.;
}
/*
______________________________________________________________________________________________

*/
Vector::Vector(const int n)
{
	Columns = n;
	
	// If the size of the vector is 0, do not allocate memory
	if (Columns == 0) return;
	
	if (Columns==5) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	  *(U+3)=0.;
	  *(U+4)=0.;
	}
	else if (Columns==4) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	  *(U+3)=0.;
	}
	else if (Columns==3) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	}
	else if (Columns==2) {
	  *U=0.;
	  *(U+1)=0.;
	}
	else if (Columns==1) {
	  *U=0.;
	}
}
/*
______________________________________________________________________________________________

*/
Vector::Vector(const real x, const real y)
{
	Columns = 2;
	*U        = x;
	*(U+1)    = y;
}
/*
______________________________________________________________________________________________

*/
Vector::Vector(const real x, const real y, const real z)
{
	Columns = 3;
	*U        = x;
	*(U+1)    = y;
	*(U+2)    = z;
}
/*
______________________________________________________________________________________________

*/
Vector::Vector(const Vector& V)
{
	Columns = V.dimension();
	if (Columns==5) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	  *(U+2) = V.value(3);
	  *(U+3) = V.value(4);
	  *(U+4) = V.value(5);
	}
	else if (Columns==4) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	  *(U+2) = V.value(3);
	  *(U+3) = V.value(4);
	}
	else if (Columns==3) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	  *(U+2) = V.value(3);
	}
	else if (Columns==2) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	}
	else if (Columns==1) {
	  *U = V.value(1);
	}
}
/*
______________________________________________________________________________________________

Distructor
______________________________________________________________________________________________

*/

Vector::~Vector()
{
	// If the size of the vector is equal to zero, do not deallocate memory
	if (Columns == 0) return;
}
/*
______________________________________________________________________________________________

Set and Get
______________________________________________________________________________________________

*/
void Vector::setZero()
{
	if (Columns==5) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	  *(U+3)=0.;
	  *(U+4)=0.;
	}
	else if (Columns==4) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	  *(U+3)=0.;
	}
	else if (Columns==3) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	}
	else if (Columns==2) {
	  *U=0.;
	  *(U+1)=0.;
	}
	else if (Columns==1) {
	  *U=0.;
	}
}
/*
______________________________________________________________________________________________

*/
void Vector::setDimension(const int n)
{
	Columns = n;
	if (Columns==5) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	  *(U+3)=0.;
	  *(U+4)=0.;
	}
	else if (Columns==4) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	  *(U+3)=0.;
	}
	else if (Columns==3) {
	  *U=0.;
	  *(U+1)=0.;
	  *(U+2)=0.;
	}
	else if (Columns==2) {
	  *U=0.;
	  *(U+1)=0.;
	}
	else if (Columns==1) {
	  *U=0.;
	}
}
/*
______________________________________________________________________________________________

Comparison
______________________________________________________________________________________________

*/

bool Vector::operator== (const Vector& V ) const
{
#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `bool Vector::operator==(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  if ( *U != V.value(1) ) return false;
	  if ( *(U+1) != V.value(2) ) return false;
	  if ( *(U+2) != V.value(3) ) return false;
	  if ( *(U+3) != V.value(4) ) return false;
	  if ( *(U+4) != V.value(5) ) return false;
	}
	else if (Columns==4) {
	  if ( *U != V.value(1) ) return false;
	  if ( *(U+1) != V.value(2) ) return false;
	  if ( *(U+2) != V.value(3) ) return false;
	  if ( *(U+3) != V.value(4) ) return false;
	}
	else if (Columns==3) {
	  if ( *U != V.value(1) ) return false;
	  if ( *(U+1) != V.value(2) ) return false;
	  if ( *(U+2) != V.value(3) ) return false;
	}
	else if (Columns==2) {
	  if ( *U != V.value(1) ) return false;
	  if ( *(U+1) != V.value(2) ) return false;
	}
	else if (Columns==1) {
	  if ( *U != V.value(1) ) return false;
	}
	return true;

}

/*
______________________________________________________________________________________________

Equality
______________________________________________________________________________________________

*/

void Vector::operator= ( const Vector& V )
{
	if (V.dimension() != Columns)
	{
	  Columns = V.dimension();
	}

	if (Columns==5) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	  *(U+2) = V.value(3);
	  *(U+3) = V.value(4);
	  *(U+4) = V.value(5);
	}
	else if (Columns==4) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	  *(U+2) = V.value(3);
	  *(U+3) = V.value(4);
	}
	else if (Columns==3) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	  *(U+2) = V.value(3);
	}
	else if (Columns==2) {
	  *U = V.value(1);
	  *(U+1) = V.value(2);
	}
	else if (Columns==1) {
	  *U = V.value(1);
	}
}

/*
______________________________________________________________________________________________

Sum
______________________________________________________________________________________________

*/

void Vector::operator+=( const Vector& V )
{
#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `void Vector::operator+=(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  *U += V.value(1);
	  *(U+1) += V.value(2);
	  *(U+2) += V.value(3);
	  *(U+3) += V.value(4);
	  *(U+4) += V.value(5);
	}
	else if (Columns==4) {
	  *U += V.value(1);
	  *(U+1) += V.value(2);
	  *(U+2) += V.value(3);
	  *(U+3) += V.value(4);
	}
	else if (Columns==3) {
	  *U += V.value(1);
	  *(U+1) += V.value(2);
	  *(U+2) += V.value(3);
	}
	else if (Columns==2) {
	  *U += V.value(1);
	  *(U+1) += V.value(2);
	}
	else if (Columns==1) {
	  *U += V.value(1);
	}
}

/*
______________________________________________________________________________________________

*/

Vector Vector::operator+ ( const Vector& V ) const
{
	Vector result(Columns);

#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `Vector& Vector::operator+(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  result.setValue(1, value(1) + V.value(1));
	  result.setValue(2, value(2) + V.value(2));
	  result.setValue(3, value(3) + V.value(3));
	  result.setValue(4, value(4) + V.value(4));
	  result.setValue(5, value(5) + V.value(5));
	}
	else if (Columns==4) {
	  result.setValue(1, value(1) + V.value(1));
	  result.setValue(2, value(2) + V.value(2));
	  result.setValue(3, value(3) + V.value(3));
	  result.setValue(4, value(4) + V.value(4));
	}
	else if (Columns==3) {
	  result.setValue(1, value(1) + V.value(1));
	  result.setValue(2, value(2) + V.value(2));
	  result.setValue(3, value(3) + V.value(3));
	}
	else if (Columns==2) {
	  result.setValue(1, value(1) + V.value(1));
	  result.setValue(2, value(2) + V.value(2));
	}
	else if (Columns==1) {
	  result.setValue(1, value(1) + V.value(1));
	}
	return result;

}

/*
______________________________________________________________________________________________

Difference
______________________________________________________________________________________________

*/

void Vector::operator-= ( const Vector& V )
{
#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `void Vector::operator-=(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  *U -= V.value(1);
	  *(U+1) -= V.value(2);
	  *(U+2) -= V.value(3);
	  *(U+3) -= V.value(4);
	  *(U+4) -= V.value(5);
	}
	else if (Columns==4) {
	  *U -= V.value(1);
	  *(U+1) -= V.value(2);
	  *(U+2) -= V.value(3);
	  *(U+3) -= V.value(4);
	}
	else if (Columns==3) {
	  *U -= V.value(1);
	  *(U+1) -= V.value(2);
	  *(U+2) -= V.value(3);
	}
	else if (Columns==2) {
	  *U -= V.value(1);
	  *(U+1) -= V.value(2);
	}
	else if (Columns==1) {
	  *U -= V.value(1);
	}
}

/*
______________________________________________________________________________________________

*/

Vector Vector::operator- ( const Vector& V ) const
{
	Vector result(Columns);

#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `Vector& Vector::operator-(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  result.setValue(1, value(1) - V.value(1));
	  result.setValue(2, value(2) - V.value(2));
	  result.setValue(3, value(3) - V.value(3));
	  result.setValue(4, value(4) - V.value(4));
	  result.setValue(5, value(5) - V.value(5));
	}
	else if (Columns==4) {
	  result.setValue(1, value(1) - V.value(1));
	  result.setValue(2, value(2) - V.value(2));
	  result.setValue(3, value(3) - V.value(3));
	  result.setValue(4, value(4) - V.value(4));
	}
	else if (Columns==3) {
	  result.setValue(1, value(1) - V.value(1));
	  result.setValue(2, value(2) - V.value(2));
	  result.setValue(3, value(3) - V.value(3));
	}
	else if (Columns==2) {
	  result.setValue(1, value(1) - V.value(1));
	  result.setValue(2, value(2) - V.value(2));
	}
	else if (Columns==1) {
	  result.setValue(1, value(1) - V.value(1));
	}
	return result;

}

/*
______________________________________________________________________________________________

*/

Vector Vector::operator- () const
{
	Vector result(Columns);

	if (Columns==5) {
	  result.setValue(1,-*U);
	  result.setValue(2,-*(U+1));
	  result.setValue(3,-*(U+2));
	  result.setValue(4,-*(U+3));
	  result.setValue(5,-*(U+4));
	}
	else if (Columns==4) {
	  result.setValue(1,-*U);
	  result.setValue(2,-*(U+1));
	  result.setValue(3,-*(U+2));
	  result.setValue(4,-*(U+3));
	}
	else if (Columns==3) {
	  result.setValue(1,-*U);
	  result.setValue(2,-*(U+1));
	  result.setValue(3,-*(U+2));
	}
	else if (Columns==2) {
	  result.setValue(1,-*U);
	  result.setValue(2,-*(U+1));
	}
	else if (Columns==1) {
	  result.setValue(1,-*U);
	}
	return result;
}
/*
______________________________________________________________________________________________

External product
______________________________________________________________________________________________

*/

void Vector::operator*= ( const real a )
{
	if (Columns==5) {
	  *U *= a;
	  *(U+1) *= a;
	  *(U+2) *= a;
	  *(U+3) *= a;
	  *(U+4) *= a;
	}
	else if (Columns==4) {
	  *U *= a;
	  *(U+1) *= a;
	  *(U+2) *= a;
	  *(U+3) *= a;
	}
	else if (Columns==3) {
	  *U *= a;
	  *(U+1) *= a;
	  *(U+2) *= a;
	}
	else if (Columns==2) {
	  *U *= a;
	  *(U+1) *= a;
	}
	else if (Columns==1) {
	  *U *= a;
	}
}
/*
______________________________________________________________________________________________

*/
Vector Vector::operator* ( const real a ) const
{
	Vector result(Columns);

	if (Columns==5) {
	  result.setValue(1, value(1) *a);
	  result.setValue(2, value(2) *a);
	  result.setValue(3, value(3) *a);
	  result.setValue(4, value(4) *a);
	  result.setValue(5, value(5) *a);
	}
	else if (Columns==4) {
	  result.setValue(1, value(1) *a);
	  result.setValue(2, value(2) *a);
	  result.setValue(3, value(3) *a);
	  result.setValue(4, value(4) *a);
	}
	else if (Columns==3) {
	  result.setValue(1, value(1) *a);
	  result.setValue(2, value(2) *a);
	  result.setValue(3, value(3) *a);
	}
	else if (Columns==2) {
	  result.setValue(1, value(1) *a);
	  result.setValue(2, value(2) *a);
	}
	else if (Columns==1) {
	  result.setValue(1, value(1) *a);
	}
	return result;
}
/*
______________________________________________________________________________________________

External division
______________________________________________________________________________________________

*/

void Vector::operator/= ( const real a )
{
#ifdef DEBUG
	if ( a == 0.)
	{
		cout << "Vector.cpp: In method `void Vector::operator/=(real)':\n";
		cout << "Vector.cpp: division by zero\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  *U /= a;
	  *(U+1) /= a;
	  *(U+2) /= a;
	  *(U+3) /= a;
	  *(U+4) /= a;
	}
	else if (Columns==4) {
	  *U /= a;
	  *(U+1) /= a;
	  *(U+2) /= a;
	  *(U+3) /= a;
	}
	else if (Columns==3) {
	  *U /= a;
	  *(U+1) /= a;
	  *(U+2) /= a;
	}
	else if (Columns==2) {
	  *U /= a;
	  *(U+1) /= a;
	}
	else if (Columns==1) {
	  *U /= a;
	}
}

/*
______________________________________________________________________________________________

*/

Vector Vector::operator/ ( const real a ) const
{
	Vector result(Columns);

#ifdef DEBUG
	if ( a == 0.)
	{
		cout << "Vector.cpp: In method `void Vector::operator/(real)':\n";
		cout << "Vector.cpp: division by zero\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  result.setValue(1, value(1) /a);
	  result.setValue(2, value(2) /a);
	  result.setValue(3, value(3) /a);
	  result.setValue(4, value(4) /a);
	  result.setValue(5, value(5) /a);
	}
	else if (Columns==4) {
	  result.setValue(1, value(1) /a);
	  result.setValue(2, value(2) /a);
	  result.setValue(3, value(3) /a);
	  result.setValue(4, value(4) /a);
	}
	else if (Columns==3) {
	  result.setValue(1, value(1) /a);
	  result.setValue(2, value(2) /a);
	  result.setValue(3, value(3) /a);
	}
	else if (Columns==2) {
	  result.setValue(1, value(1) /a);
	  result.setValue(2, value(2) /a);
	}
	else if (Columns==1) {
	  result.setValue(1, value(1) /a);
	}
	return result;
}
/*
______________________________________________________________________________________________

Dot product
______________________________________________________________________________________________

*/

real Vector::operator* ( const Vector& V ) const
{
	real result = 0.;

#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `real Vector::operator*(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	if (Columns==5) {
	  result += *U * V.value(1);
	  result += *(U+1) * V.value(2);
	  result += *(U+2) * V.value(3);
	  result += *(U+3) * V.value(4);
	  result += *(U+4) * V.value(5);
	}
	else if (Columns==4) {
	  result += *U * V.value(1);
	  result += *(U+1) * V.value(2);
	  result += *(U+2) * V.value(3);
	  result += *(U+3) * V.value(4);
	}
	else if (Columns==3) {
	  result += *U * V.value(1);
	  result += *(U+1) * V.value(2);
	  result += *(U+2) * V.value(3);
	}
	else if (Columns==2) {
	  result += *U * V.value(1);
	  result += *(U+1) * V.value(2);
	}
	else if (Columns==1) {
	  result += *U * V.value(1);
	}

	return result;

}
/*
______________________________________________________________________________________________

Term by term product
______________________________________________________________________________________________

*/


Vector Vector::operator| ( const Vector& V ) const
{
	int n;
	Vector result(Columns);

#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `real Vector::operator|(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif

	for (n=1; n<=Columns; n++)
		result.setValue(n, value(n) * V.value(n));

	return result;

}

/*
______________________________________________________________________________________________

Cross product
______________________________________________________________________________________________

*/

Vector Vector::operator^ ( const Vector& V ) const
{
	Vector result(3);

#ifdef DEBUG
	if ( Columns != V.dimension())
	{
		cout << "Vector.cpp: In method `real Vector::operator^(Vector&)':\n";
		cout << "Vector.cpp: vectors have different dimensions\n";
		cout << "carmen: *** [Vector.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
#endif
	
	if (Columns > 1)
	{
		result.setValue(3, value(1)*V.value(2)-value(2)*V.value(1));
		if (Columns > 2)
		{
			result.setValue(1, value(2)*V.value(3)-value(3)*V.value(2));
			result.setValue(2, value(3)*V.value(1)-value(1)*V.value(3));
		}
	}				

	return result;

}
/*
______________________________________________________________________________________________

Returns true if one of the components of the vector is not a number
______________________________________________________________________________________________

*/

bool Vector::isNaN() const
{
	int n;

  for (n = 1; n <= Columns; n++)

//On some machines should be math::isnan
#if defined SP_SMP
   if (::isnan(value(n)) )
#else
   if ( isnan(value(n)) )
#endif
    return true;
		
	return false;
}

/*
_____________________________________________________________________________________________

EXTERNAL FUNCTIONS
	External product (Vector = real * Vector)
______________________________________________________________________________________________

*/
Vector operator* (const real a, const Vector& V)
{
	return V*a;
}
/*
______________________________________________________________________________________________

Dimension
______________________________________________________________________________________________

*/
int dim (const Vector& V)
{
	return V.dimension();
}
/*
______________________________________________________________________________________________

Abs function
______________________________________________________________________________________________

*/
Vector abs(const Vector& V)
{
	int n;
	real a;
	Vector result( dim(V) );
	
	for (n = 1; n <= dim(V); n++)
	{	
		a = V.value(n);
		if (a < 0.)
			result.setValue(n, -a);
		else
			result.setValue(n, a);
	}
	return result;
}
/*
______________________________________________________________________________________________

L1 Norm
______________________________________________________________________________________________

*/
real N1(const Vector& V)
{
	int n;
	real result;
	
	result = 0.;
	for (n = 1; n <= dim(V); n++)
		result += fabs(V.value(n));

	return result;
}
/*
______________________________________________________________________________________________

L2 Norm
______________________________________________________________________________________________

*/
real N2(const Vector& V)
{
	int n;
	real result;
	
	result = 0.;
	for (n = 1; n <= dim(V); n++)
		result += V.value(n)*V.value(n);

	return sqrt(result);
}
/*
______________________________________________________________________________________________

Max Norm
______________________________________________________________________________________________

*/
real NMax(const Vector& V)
{
	int n;
	real result;
	
	result = 0.;
	for (n = 1; n <= dim(V); n++)
		if (result < fabs(V.value(n))) result = fabs(V.value(n));

	return result;
}
/*
______________________________________________________________________________________________

Writes the components of the vector V on screen
______________________________________________________________________________________________

*/

ostream& operator<<(ostream& out, const Vector& V)
{
	int n;

	for (n = 1; n <= dim(V); n++)
	{
		out<<n<<": "<<V.value(n)<<endl;
   }
	return out;
}

