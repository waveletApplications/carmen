/***************************************************************************
                          Matrix.cpp  -  description
                             -------------------
    begin                : mer fév 11 2004
    copyright            : (C) 2004 by Olivier Roussel and Alexei Tsigulin
    email                : roussel@ict.uni-karlsruhe.de; lpsoft@mail.ru
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

Constructors
______________________________________________________________________________________________

*/

Matrix::Matrix()
{
	Lines = Columns = 1;
	U  = new real;
	*U = 0.;
}
/*
______________________________________________________________________________________________

*/
Matrix::Matrix(const int i, const int j)
{
	int n; // Counter

	Lines = i;
  	Columns = j;

	// If the size of the matrix is equal to zero, do not allocate memory
	if (Lines==0 && Columns==0) return;
	
	U = new real[Lines*Columns];

	for (n = 1; n <= Lines*Columns; n++)
		*(U+n-1) = 0.;
}
/*
______________________________________________________________________________________________

*/
Matrix::Matrix(const int i)
{
	int n; // Counter

	Lines = i;
  Columns = i;

	U = new real[Lines*Columns];

	for (n = 1; n <= Lines*Columns; n++)
		*(U+n-1) = 0.;
}
/*
______________________________________________________________________________________________

*/
Matrix::Matrix(const Matrix& M)
{
	int i,j;

	Lines = M.lines();
  Columns = M.columns();

	U = new real[Lines*Columns];

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		setValue(i, j, M.value(i,j));
}
/*
______________________________________________________________________________________________

*/
Matrix::Matrix(const Vector& V)
{
	int i;

	Lines = V.dimension();
  Columns = 1;

	U = new real[Lines*Columns];

	for (i = 1; i <= Lines; i++)
		setValue(i, 1, V.value(i));
}
/*
______________________________________________________________________________________________

Distructor
______________________________________________________________________________________________

*/

Matrix::~Matrix()
{
	// If the size of the matrix is equal to zero, do not deallocate memory
	if (Lines==0 && Columns==0) return;
	
	delete[] U;
}
/*
______________________________________________________________________________________________

Set and Get
______________________________________________________________________________________________

*/

void Matrix::setZero()
{
	int n;

	for (n = 1; n <= Lines*Columns; n++)
		*(U+n-1) = 0.;
}

/*
______________________________________________________________________________________________

Comparison
______________________________________________________________________________________________

*/

bool Matrix::operator== (const Matrix& M ) const
{
	int i, j;

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		if ( value(i,j) != M.value(i,j) ) return false;

	return true;
}
/*
______________________________________________________________________________________________

Equality
______________________________________________________________________________________________

*/

void Matrix::operator= ( const Matrix& M )
{
	int i, j;

	if (M.lines() != Lines || M.columns() != Columns)
	{
		delete[] U;

		Lines   = M.lines();
  		Columns = M.columns();

		U = new real[Lines*Columns];
	}

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		setValue(i, j, M.value(i,j));
}
/*
______________________________________________________________________________________________

Sum
______________________________________________________________________________________________

*/

void Matrix::operator+=( const Matrix& M )
{
	int i, j;

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		setValue(i, j, value(i,j) + M.value(i,j));
}

/*
______________________________________________________________________________________________

*/

Matrix Matrix::operator+ ( const Matrix& M ) const
{
	int i, j;
	Matrix result(Lines, Columns);

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		result.setValue(i, j, value(i,j) + M.value(i,j));

	return result;
}
/*
______________________________________________________________________________________________

Difference
______________________________________________________________________________________________

*/

void Matrix::operator-= ( const Matrix& M )
{
	int i, j;

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		setValue(i, j, value(i,j) - M.value(i,j));
}

/*
______________________________________________________________________________________________

*/

Matrix Matrix::operator- ( const Matrix& M ) const
{
	int i, j;
	Matrix result(Lines, Columns);

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		result.setValue(i, j, value(i,j) - M.value(i,j));

	return result;
}
/*
______________________________________________________________________________________________

*/
Matrix Matrix::operator- () const
{
	int i, j;
	Matrix result(Lines, Columns);

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		result.setValue(i, j, -value(i,j));

	return result;
}
/*
______________________________________________________________________________________________

External product
______________________________________________________________________________________________

*/

void Matrix::operator*= ( const real a )
{
	int n;

	for (n=1; n<=Lines*Columns; n++)
		*(U+n-1) *= a;
}
/*
______________________________________________________________________________________________

*/

Matrix Matrix::operator* ( const real a ) const
{
	int i,j;
	Matrix result(Lines, Columns);

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		result.setValue(i, j, value(i,j)*a);

	return result;
}
/*
______________________________________________________________________________________________

External division
______________________________________________________________________________________________

*/

void Matrix::operator/= ( const real a )
{
	int n;

	for (n=1; n<=Lines*Columns; n++)
		*(U+n-1) /= a;
}

/*
______________________________________________________________________________________________

*/

Matrix Matrix::operator/ ( const real a ) const
{
	int i,j;
	Matrix result(Lines, Columns);

	for (i = 1; i <= Lines; i++)
	for (j = 1; j <= Columns; j++)
		result.setValue(i, j, value(i,j)/a);

	return result;
}
/*
______________________________________________________________________________________________

Matrix product
______________________________________________________________________________________________

*/

Matrix Matrix::operator* ( const Matrix& M ) const
{
	int i, j, k;
	Matrix result(lines(), M.columns());

	for (i = 1; i <= lines(); i++)
	for (j = 1; j <= M.columns(); j++)
	for (k = 1; k <= columns(); k++)
		result.setValue(i,j, result.value(i,j) + value(i,k)*M.value(k,j));

	return result;
}
/*
______________________________________________________________________________________________

Matrix-vector product
______________________________________________________________________________________________

*/
Vector Matrix::operator* (const Vector& V) const
{
	int i, k;

	Vector result(lines());

	for (i = 1; i <= lines(); i++)
	for (k = 1; k <= columns(); k++)
		result.setValue(i, result.value(i) + value(i,k)*V.value(k));

	return result;

}
/*
______________________________________________________________________________________________

Set eigenmatrix
______________________________________________________________________________________________

*/
void Matrix::setEigenMatrix(const bool isLeft, const int AxisNo, const Vector V, const real c, const real h)
{

	// --- Local variables ---

	int n; // Counter
	real u = V.value(1);
	real v = (Dimension > 1) ? V.value(2): 0.;
	real w = (Dimension > 2) ? V.value(3): 0.;

	 // --- Allocate memory and set to zero ---

	delete[] U;

	Lines = Columns = QuantityNb;

	U = new real[Lines*Columns];

	for (n = 1; n <= Lines*Columns; n++)
		*(U+n-1) = 0.;

	// --- Fill values ---

	switch (Dimension)
	{

		// --- 1D CASE --------------------------------------------------------------------------------

		case 1:
			if (isLeft)
			{
				// --- Left eigenmatrix ---

//  		  setValue(1,1, (c*c)*(1.- ((Gamma-1.)*u*u)/(2.*c*c)));
//  		  setValue(2,1, (c*c)*(.5*(((Gamma-1.)*u*u)/(2.*c*c) - (u/c))));
//  		  setValue(3,1, (c*c)*(.5*(((Gamma-1.)*u*u)/(2.*c*c) + (u/c))));
//
//  		  setValue(1,2, (c*c)*((Gamma-1.)*u/(c*c)));
//  		  setValue(2,2, (c*c)*(((Gamma-1.)/(2.*c*c))*(c/(Gamma-1.)-u)));
//  		  setValue(3,2, (c*c)*(((Gamma-1.)/(2.*c*c))*(-c/(Gamma-1.)-u)));
//
//  		  setValue(1,3, (c*c)*(-(Gamma-1.)/(c*c)));
//  		  setValue(2,3, (c*c)*((Gamma-1.)/(2.*c*c)));
//  		  setValue(3,3, (c*c)*((Gamma-1.)/(2.*c*c)));

//			matrices de Vinokur

   		setValue(1,1, (1.- ((Gamma-1.)*u*u)/(2.*c*c)));
  		  setValue(2,1, (.5*(((Gamma-1.)*u*u)/(2.*c*c) - (u/c))));
  		  setValue(3,1, (.5*(((Gamma-1.)*u*u)/(2.*c*c) + (u/c))));

  		  setValue(1,2, ((Gamma-1.)*u/(c*c)));
  		  setValue(2,2, (((Gamma-1.)/(2.*c*c))*(c/(Gamma-1.)-u)));
  		  setValue(3,2, (((Gamma-1.)/(2.*c*c))*(-c/(Gamma-1.)-u)));

  		  setValue(1,3, (-(Gamma-1.)/(c*c)));
  		  setValue(2,3, ((Gamma-1.)/(2.*c*c)));
  		  setValue(3,3, ((Gamma-1.)/(2.*c*c)));

			}
			else
			{
				// --- Right eigenmatrix ---

//  			setValue(1,1, 1./(c*c));
//  			setValue(2,1, u/(c*c));
//  			setValue(3,1, .5*(u*u)/(c*c));
//
//  			setValue(1,2, 1./(c*c));
//  			setValue(2,2, (u+c)/(c*c));
//  			setValue(3,2, (h+ c*u)/(c*c));
//
//  			setValue(1,3, 1./(c*c));
//  			setValue(2,3, (u-c)/(c*c));
//  			setValue(3,3, (h-u*c)/(c*c));


//			matrices de Vinokur

  			setValue(1,1, 1.);
  			setValue(2,1, u);
  			setValue(3,1, .5*(u*u));

  			setValue(1,2, 1.);
  			setValue(2,2, (u+c));
  			setValue(3,2, (h+ c*u));

  			setValue(1,3, 1.);
  			setValue(2,3, (u-c));
  			setValue(3,3, (h-u*c));


			}
			break;


		// --- 2D CASE --------------------------------------------------------------------------------

		case 2:
			if (isLeft)
			{
				// --- Left eigenmatrix ---

				switch(AxisNo)
				{
					case 1:
//  					setValue(1,1, (c*c)*(1.-((Gamma-1.)*(u*u+v*v)/(2.*c*c))));
//  					setValue(2,1, (c*c)*(-v));
//  					setValue(3,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) - u/c )));
//  					setValue(4,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) + u/c )));
//
//  					setValue(1,2, (c*c)*((Gamma-1.)*u/(c*c)));
//  					setValue(2,2, 0.);
//  					setValue(3,2, (c*c)*((c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
//  					setValue(4,2, (c*c)*((-c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
//
//  					setValue(1,3, (c*c)*((Gamma-1.)*v/(c*c)));
//  					setValue(2,3, (c*c)*1.);
//  					setValue(3,3, (c*c)*(-(Gamma-1.)*v/(2.*c*c)));
//  					setValue(4,3, (c*c)*(-(Gamma-1.)*v/(2.*c*c)));
//
//  					setValue(1,4, (c*c)*(-(Gamma-1.)/(c*c)));
//  					setValue(2,4, 0.);
//  					setValue(3,4, (c*c)*((Gamma-1.)/(2.*c*c)));
//  					setValue(4,4, (c*c)*((Gamma-1.)/(2.*c*c)));


//			matrix vinokur

  					setValue(1,1, (1.-((Gamma-1.)*(u*u+v*v)/(2.*c*c))));
  					setValue(2,1, (-v));
  					setValue(3,1, (.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) - u/c )));
  					setValue(4,1, (.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) + u/c )));

  					setValue(1,2, ((Gamma-1.)*u/(c*c)));
  					setValue(2,2, 0.);
  					setValue(3,2, ((c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
  					setValue(4,2, ((-c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));

  					setValue(1,3, ((Gamma-1.)*v/(c*c)));
  					setValue(2,3, 1.);
  					setValue(3,3, (-(Gamma-1.)*v/(2.*c*c)));
  					setValue(4,3, (-(Gamma-1.)*v/(2.*c*c)));

  					setValue(1,4, (-(Gamma-1.)/(c*c)));
  					setValue(2,4, 0.);
  					setValue(3,4, ((Gamma-1.)/(2.*c*c)));
  					setValue(4,4, ((Gamma-1.)/(2.*c*c)));


						break;

					case 2:
//					setValue(1,1, (c*c)*(1.-((Gamma-1.)*(u*u+v*v)/(2.*c*c))));
//  					setValue(2,1, (c*c)*(-v));
//  					setValue(3,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) - u/c )));
//  					setValue(4,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) + u/c )));
//
//  					setValue(1,2, (c*c)*((Gamma-1.)*u/(c*c)));
//  					setValue(2,2, 0.);
//  					setValue(3,2, (c*c)*((c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
//  					setValue(4,2, (c*c)*((-c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
//
//  					setValue(1,3, (c*c)*((Gamma-1.)*v/(c*c)));
//  					setValue(2,3, (c*c)*1.);
//  					setValue(3,3, (c*c)*(-(Gamma-1.)*v/(2.*c*c)));
//  					setValue(4,3, (c*c)*(-(Gamma-1.)*v/(2.*c*c)));
//
//  					setValue(1,4, (c*c)*(-(Gamma-1.)/(c*c)));
//  					setValue(2,4, 0.);
//  					setValue(3,4, (c*c)*((Gamma-1.)/(2.*c*c)));
//  					setValue(4,4, (c*c)*((Gamma-1.)/(2.*c*c)));


//			matrix vinokur
 					setValue(1,1, (1.-((Gamma-1.)*(u*u+v*v)/(2.*c*c))));
  					setValue(2,1, (-v));
  					setValue(3,1, (.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) - u/c )));
  					setValue(4,1, (.5*( (Gamma-1.)*(u*u+v*v)/(2.*c*c) + u/c )));

  					setValue(1,2, ((Gamma-1.)*u/(c*c)));
  					setValue(2,2, 0.);
  					setValue(3,2, ((c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
  					setValue(4,2, ((-c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));

  					setValue(1,3, ((Gamma-1.)*v/(c*c)));
  					setValue(2,3, 1.);
  					setValue(3,3, (-(Gamma-1.)*v/(2.*c*c)));
  					setValue(4,3, (-(Gamma-1.)*v/(2.*c*c)));

  					setValue(1,4, (-(Gamma-1.)/(c*c)));
  					setValue(2,4, 0.);
  					setValue(3,4, ((Gamma-1.)/(2.*c*c)));
  					setValue(4,4, ((Gamma-1.)/(2.*c*c)));


						break;
				};
			}
			else
			{
				// --- Right eigenmatrix ---

				switch(AxisNo)
				{
					case 1:
//  				  setValue(1,1, 1./(c*c));
//  				  setValue(2,1, u/(c*c));
//  				  setValue(3,1, v/(c*c));
//  				  setValue(4,1, .5*(u*u+v*v)/(c*c));
//
//  				  setValue(1,2, 0.);
//  				  setValue(2,2, 0.);
//  				  setValue(3,2, 1./(c*c));
//  				  setValue(4,2, v/(c*c));
//
//  				  setValue(1,3, 1./(c*c));
//  				  setValue(2,3, (u+c)/(c*c));
//  				  setValue(3,3, v/(c*c));
//  				  setValue(4,3, (h+ c*u)/(c*c));
//
//  				  setValue(1,4, 1./(c*c));
//  				  setValue(2,4, (u-c)/(c*c));
//  				  setValue(3,4, v/(c*c));
//  				  setValue(4,4, (h- c*u)/(c*c));


//			matrix vinokur
  				  setValue(1,1, 1.);
  				  setValue(2,1, u);
  				  setValue(3,1, v);
  				  setValue(4,1, .5*(u*u+v*v));

  				  setValue(1,2, 0.);
  				  setValue(2,2, 0.);
  				  setValue(3,2, 1.);
  				  setValue(4,2, v);

  				  setValue(1,3, 1.);
  				  setValue(2,3, (u+c));
  				  setValue(3,3, v);
  				  setValue(4,3, (h+ c*u));

  				  setValue(1,4, 1.);
  				  setValue(2,4, (u-c));
  				  setValue(3,4, v);
  				  setValue(4,4, (h- c*u));

						break;

					case 2:
//  				  setValue(1,1, 0.);
//  				  setValue(2,1, 0.);
//  				  setValue(3,1, 1./(c*c));
//  				  setValue(4,1, u/(c*c));
//
//  				  setValue(1,2, 1./(c*c));
//  				  setValue(2,2, u/(c*c));
//  				  setValue(3,2, v/(c*c));
//  				  setValue(4,2, .5*(u*u+v*v)/(c*c));
//
//  				  setValue(1,3, 1./(c*c));
//  				  setValue(2,3, u/(c*c));
//  				  setValue(3,3, (v+c)/(c*c));
//  				  setValue(4,3, (h+ c*v)/(c*c));
//
//  				  setValue(1,4, 1./(c*c));
//  				  setValue(2,4, u/(c*c));
//  				  setValue(3,4, (v-c)/(c*c));
//  				  setValue(4,4, (h- c*v)/(c*c));


//			matrix vinokur

   				setValue(1,1, 0.);
  				  setValue(2,1, 0.);
  				  setValue(3,1, 1.);
  				  setValue(4,1, u);

  				  setValue(1,2, 1.);
  				  setValue(2,2, u);
  				  setValue(3,2, v);
  				  setValue(4,2, .5*(u*u+v*v));

  				  setValue(1,3, 1.);
  				  setValue(2,3, u);
  				  setValue(3,3, (v+c));
  				  setValue(4,3, (h+ c*v));

  				  setValue(1,4, 1.);
  				  setValue(2,4, u);
  				  setValue(3,4, (v-c));
  				  setValue(4,4, (h- c*v));

						break;
				};
			}
			break;

		// --- 3D CASE --------------------------------------------------------------------------------

		case 3:
			if (isLeft)
			{
				// --- Left eigenmatrix ---

				switch(AxisNo)
				{
//					case 1:
//  					setValue(1,1, (c*c)*(1.-( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c))));
//  					setValue(2,1, (c*c)*(-v));
//  					setValue(3,1, (c*c)*(-w));
// 	 				setValue(4,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) - u/c )));
//  					setValue(5,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) + u/c )));
//
//  					setValue(1,2, (c*c)*((Gamma-1.)*u/(c*c)));
//  					setValue(2,2, 0.);
//  					setValue(3,2, 0.);
//  					setValue(4,2, (c*c)*((c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
//  					setValue(5,2, (c*c)*((-c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
//
//  					setValue(1,3, (c*c)*((Gamma-1.)*v/(c*c)));
//  					setValue(2,3, (c*c)*(1.));
//  					setValue(3,3, 0.);
//  					setValue(4,3, (c*c)*(-v*(Gamma-1.)/(2.*c*c)));
//  					setValue(5,3, (c*c)*(-v*(Gamma-1.)/(2.*c*c)));
//
//  					setValue(1,4, (c*c)*((Gamma-1.)*w/(c*c)));
//  					setValue(2,4, 0.);
//  					setValue(3,4, (c*c)*(1.));
//  					setValue(4,4, (c*c)*(-w*(Gamma-1.)/(2.*c*c)));
//  					setValue(5,4, (c*c)*(-w*(Gamma-1.)/(2.*c*c)));
//
//  					setValue(1,5, (c*c)*(-(Gamma-1.)/(c*c)));
//  					setValue(2,5, 0.);
//  					setValue(3,5, 0.);
//  					setValue(4,5, (c*c)*((Gamma-1.)/(2.*c*c)));
//  					setValue(5,5, (c*c)*((Gamma-1.)/(2.*c*c)));
//						break;

//			matrix vinokur
					case 1:
  					setValue(1,1, (1.-( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c))));
  					setValue(2,1, (-v));
  					setValue(3,1, (-w));
 	 				setValue(4,1, (.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) - u/c )));
  					setValue(5,1, (.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) + u/c )));

  					setValue(1,2, (Gamma-1.)*u/(c*c));
  					setValue(2,2, 0.);
  					setValue(3,2, 0.);
  					setValue(4,2, ((c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));
  					setValue(5,2, ((-c/(Gamma-1.)-u)*((Gamma-1.)/(2.*c*c))));

  					setValue(1,3, ((Gamma-1.)*v/(c*c)));
  					setValue(2,3, (1.));
  					setValue(3,3, 0.);
  					setValue(4,3, -v*(Gamma-1.)/(2.*c*c));
  					setValue(5,3, -v*(Gamma-1.)/(2.*c*c));

  					setValue(1,4, (Gamma-1.)*w/(c*c));
  					setValue(2,4, 0.);
  					setValue(3,4, (1.));
  					setValue(4,4, -w*(Gamma-1.)/(2.*c*c));
  					setValue(5,4, -w*(Gamma-1.)/(2.*c*c));

  					setValue(1,5, (-(Gamma-1.)/(c*c)));
  					setValue(2,5, 0.);
  					setValue(3,5, 0.);
  					setValue(4,5, (Gamma-1.)/(2.*c*c));
  					setValue(5,5, (Gamma-1.)/(2.*c*c));
						break;








//					case 2:
//					  setValue(1,1,(c*c)*( -u));
//					  setValue(2,1, (c*c)*(1.-( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c))));
//					  setValue(3,1, (c*c)*(-w));
//					  setValue(4,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) - v/c )));
//					  setValue(5,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) + v/c )));
//
//					  setValue(1,2, 0.);
//					  setValue(2,2, (c*c)*((Gamma-1.)*u/(c*c)));
//					  setValue(3,2, 0.);
//					  setValue(4,2, (c*c)*(-u*(Gamma-1.)/(2.*c*c)));
//					  setValue(5,2, (c*c)*(-u*(Gamma-1.)/(2.*c*c)));
//
//					  setValue(1,3, (c*c)*(1.));
//					  setValue(2,3, (c*c)*((Gamma-1.)*v/(c*c)));
//					  setValue(3,3, 0.);
//					  setValue(4,3, (c*c)*((c/(Gamma-1.)-v)*((Gamma-1.)/(2.*c*c))));
//					  setValue(5,3, (c*c)*((-c/(Gamma-1.)-v)*((Gamma-1.)/(2.*c*c))));
//
//					  setValue(1,4, 0.);
//					  setValue(2,4, (c*c)*((Gamma-1.)*w/(c*c)));
//					  setValue(3,4, (c*c)*(1.));
//					  setValue(4,4, (c*c)*(-w*(Gamma-1.)/(2.*c*c)));
//					  setValue(5,4, (c*c)*(-w*(Gamma-1.)/(2.*c*c)));
//
//					  setValue(1,5, 0.);
//					  setValue(2,5, (c*c)*(-(Gamma-1.)/(c*c)));
//					  setValue(3,5, 0.);
//					  setValue(4,5, (c*c)*((Gamma-1.)/(2.*c*c)));
//					  setValue(5,5, (c*c)*((Gamma-1.)/(2.*c*c)));
//						break;

//			matrix vinokur
					case 2:
					  setValue(1,1, -u);
					  setValue(2,1, (1.-( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c))));
					  setValue(3,1, -w);
					  setValue(4,1, (.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) - v/c )));
					  setValue(5,1, (.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) + v/c )));

					  setValue(1,2, 0.);
					  setValue(2,2, (Gamma-1.)*u/(c*c));
					  setValue(3,2, 0.);
					  setValue(4,2, (-u*(Gamma-1.)/(2.*c*c)));
					  setValue(5,2, (-u*(Gamma-1.)/(2.*c*c)));

					  setValue(1,3, (c*c)*(1.));
					  setValue(2,3, (c*c)*((Gamma-1.)*v/(c*c)));
					  setValue(3,3, 0.);
					  setValue(4,3, (c*c)*((c/(Gamma-1.)-v)*((Gamma-1.)/(2.*c*c))));
					  setValue(5,3, (c*c)*((-c/(Gamma-1.)-v)*((Gamma-1.)/(2.*c*c))));

					  setValue(1,4, 0.);
					  setValue(2,4, ((Gamma-1.)*w/(c*c)));
					  setValue(3,4, 1.);
					  setValue(4,4, -w*(Gamma-1.)/(2.*c*c));
					  setValue(5,4, -w*(Gamma-1.)/(2.*c*c));

					  setValue(1,5, 0.);
					  setValue(2,5, -(Gamma-1.)/(c*c));
					  setValue(3,5, 0.);
					  setValue(4,5, (Gamma-1.)/(2.*c*c));
					  setValue(5,5, (Gamma-1.)/(2.*c*c));
						break;



//					case 3:
//					  setValue(1,1, (c*c)*(-u));
//					  setValue(2,1, (c*c)*(-v));
//					  setValue(3,1, (c*c)*(1.-( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c))));
//					  setValue(4,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) - w/c )));
//					  setValue(5,1, (c*c)*(.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) + w/c )));
//
//					  setValue(1,2, (c*c)*(1.));
//					  setValue(2,2, 0.);
//					  setValue(3,2, (c*c)*((Gamma-1.)*u/(c*c)));
//					  setValue(4,2, (c*c)*(-(Gamma-1.)*u/(2.*c*c)));
//					  setValue(5,2, (c*c)*(-(Gamma-1.)*u/(2.*c*c)));
//
//					  setValue(1,3, 0.);
//					  setValue(2,3, (c*c)*(1.));
//					  setValue(3,3, (c*c)*((Gamma-1.)*v/(c*c)));
//					  setValue(4,3, (c*c)*(-(Gamma-1.)*v/(2.*c*c)));
//					  setValue(5,3, (c*c)*(-(Gamma-1.)*v/(2.*c*c)));
//
//					  setValue(1,4, 0.);
//					  setValue(2,4, 0.);
//					  setValue(3,4, (c*c)*((Gamma-1.)*w/(c*c)));
//					  setValue(4,4, (c*c)*((c/(Gamma-1.)-w)*((Gamma-1.)/(2.*c*c))));
//					  setValue(5,4, (c*c)*((-c/(Gamma-1.)-w)*((Gamma-1.)/(2.*c*c))));
//
//					  setValue(1,5, 0.);
//					  setValue(2,5, 0.);
//					  setValue(3,5, (c*c)*(-(Gamma-1.)/(c*c)));
//					  setValue(4,5, (c*c)*((Gamma-1.)/(2.*c*c)));
//					  setValue(5,5, (c*c)*((Gamma-1.)/(2.*c*c)));
//						break;


//			matrix vinokur
					case 3:
					  setValue(1,1, -u);
					  setValue(2,1, -v);
					  setValue(3,1, (1.-( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c))));
					  setValue(4,1, (.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) - w/c )));
					  setValue(5,1, (.5*( (Gamma-1.)*(u*u+v*v+w*w)/(2.*c*c) + w/c )));

					  setValue(1,2, 1.);
					  setValue(2,2, 0.);
					  setValue(3,2, (Gamma-1.)*u/(c*c));
					  setValue(4,2, -(Gamma-1.)*u/(2.*c*c));
					  setValue(5,2, -(Gamma-1.)*u/(2.*c*c));

					  setValue(1,3, 0.);
					  setValue(2,3, 1.);
					  setValue(3,3, (Gamma-1.)*v/(c*c));
					  setValue(4,3, -(Gamma-1.)*v/(2.*c*c));
					  setValue(5,3, -(Gamma-1.)*v/(2.*c*c));

					  setValue(1,4, 0.);
					  setValue(2,4, 0.);
					  setValue(3,4, (Gamma-1.)*w/(c*c));
					  setValue(4,4, ((c/(Gamma-1.)-w)*((Gamma-1.)/(2.*c*c))));
					  setValue(5,4, ((-c/(Gamma-1.)-w)*((Gamma-1.)/(2.*c*c))));

					  setValue(1,5, 0.);
					  setValue(2,5, 0.);
					  setValue(3,5, -(Gamma-1.)/(c*c));
					  setValue(4,5, (Gamma-1.)/(2.*c*c));
					  setValue(5,5, (Gamma-1.)/(2.*c*c));
						break;



				};
			}
			else
			{
				// --- Right eigenmatrix ---

				switch(AxisNo)
				{
//					case 1:
//					  setValue(1,1, 1./(c*c));
//					  setValue(2,1, u/(c*c));
//					  setValue(3,1, v/(c*c));
//					  setValue(4,1, w/(c*c));
//					  setValue(5,1, .5*(u*u+v*v+w*w)/(c*c));
//
//					  setValue(1,2, 0.);
//					  setValue(2,2, 0.);
//					  setValue(3,2, 1./(c*c));
//					  setValue(4,2, 0.);
//					  setValue(5,2, v/(c*c));
//
//					  setValue(1,3, 0.);
//					  setValue(2,3, 0.);
//					  setValue(3,3, 0.);
//					  setValue(4,3, 1./(c*c));
//					  setValue(5,3, w/(c*c));
//
//					  setValue(1,4, 1./(c*c));
//					  setValue(2,4, (u+c)/(c*c));
//					  setValue(3,4, v/(c*c));
//					  setValue(4,4, w/(c*c));
//					  setValue(5,4, (h+ c*u)/(c*c));
//
//					  setValue(1,5, 1.);
//					  setValue(2,5, (u-c)/(c*c));
//					  setValue(3,5, v/(c*c));
//					  setValue(4,5, w/(c*c));
//					  setValue(5,5, (h- c*u)/(c*c));
//						break;

//			matrix vinokur
					case 1:
					  setValue(1,1, 1.);
					  setValue(2,1, u);
					  setValue(3,1, v);
					  setValue(4,1, w);
					  setValue(5,1, .5*(u*u+v*v+w*w));

					  setValue(1,2, 0.);
					  setValue(2,2, 0.);
					  setValue(3,2, 1.);
					  setValue(4,2, 0.);
					  setValue(5,2, v);

					  setValue(1,3, 0.);
					  setValue(2,3, 0.);
					  setValue(3,3, 0.);
					  setValue(4,3, 1.);
					  setValue(5,3, w);

					  setValue(1,4, 1.);
					  setValue(2,4, u+c);
					  setValue(3,4, v);
					  setValue(4,4, w);
					  setValue(5,4, h+ c*u);

					  setValue(1,5, 1.);
					  setValue(2,5, u-c);
					  setValue(3,5, v);
					  setValue(4,5, w);
					  setValue(5,5, h- c*u);
						break;


//			matrix vinokur
					case 2:
					  setValue(1,1, 0.);
					  setValue(2,1, 0.);
					  setValue(3,1, 1.);
					  setValue(4,1, 0.);
					  setValue(5,1, u);

					  setValue(1,2, 1.);
					  setValue(2,2, u);
					  setValue(3,2, v);
					  setValue(4,2, w);
					  setValue(5,2, .5*(u*u+v*v+w*w));

					  setValue(1,3, 0.);
					  setValue(2,3, 0.);
					  setValue(3,3, 0.);
					  setValue(4,3, 1.);
					  setValue(5,3, w);

					  setValue(1,4, 1.);
					  setValue(2,4, u);
					  setValue(3,4, v+c);
					  setValue(4,4, w);
					  setValue(5,4, h+ c*v);

					  setValue(1,5, 1.);
					  setValue(2,5, u);
					  setValue(3,5, v-c);
					  setValue(4,5, w);
					  setValue(5,5, h- c*v);
						break;

//			matrix vinokur

					case 3:
					  setValue(1,1, 0.);
					  setValue(2,1, 1.);
					  setValue(3,1, 0.);
					  setValue(4,1, 0.);
					  setValue(5,1, u);

					  setValue(1,2, 0.);
					  setValue(2,2, 0.);
					  setValue(3,2, 1.);
					  setValue(4,2, 0.);
					  setValue(5,2, v);

					  setValue(1,3, 1.);
					  setValue(2,3, u);
					  setValue(3,3, v);
					  setValue(4,3, w);
					  setValue(5,3, .5*(u*u+v*v+w*w));

					  setValue(1,4, 1.);
					  setValue(2,4, u);
					  setValue(3,4, v);
					  setValue(4,4, w+c);
					  setValue(5,4, h+ c*w);

					  setValue(1,5, 1.);
					  setValue(2,5, u);
					  setValue(3,5, v);
					  setValue(4,5, w-c);
					  setValue(5,5, h- c*w);
						break;
				};
			}
			break;
	};

}

/*
_____________________________________________________________________________________________

EXTERNAL FUNCTIONS
	External product (Matrix = real * Matrix)
______________________________________________________________________________________________

*/
Matrix operator* (const real a, const Matrix& M)
{
	return M*a;
}

/*           
______________________________________________________________________________________________

Writes the components of the matrix M on screen.
______________________________________________________________________________________________

*/

ostream& operator<<(ostream& out, const  Matrix& M)
{
	int n;
	int m;

	for (n = 1; n <= M.lines(); n++)
	{
	  for (m = 1; m <= M.columns(); m++)
	  {
		  out<<n<<","<<m<<": "<<M.value(n,m)<<endl;
    }
  }
	return out;
}
