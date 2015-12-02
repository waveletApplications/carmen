/***************************************************************************
                          TimeAverageGrid.h  -  description
                             -------------------
    begin                : ven déc 3 2004
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

class TimeAverageGrid
{
/*
______________________________________________________________________________________________

PUBLIC METHODS
______________________________________________________________________________________________

*/
	public:
		TimeAverageGrid(const int UserScaleNb, const int UserQuantityNb);
		TimeAverageGrid(const int UserScaleNb);
		~TimeAverageGrid();

		void updateValue(const int ElementNo, const int QuantityNo, const real UserValue);
		void updateValue(const int i, const int j, const int k, const int QuantityNo, const real UserValue);
		void updateValue(const int i, const int j, const int k, const Vector arg);

		inline void updateSample();

		inline real value(const int ElementNo, const int QuantityNo) const;
		inline real value(const int i, const int j, const int k, const int QuantityNo) const;

		inline real density(const int i, const int j, const int k) const;
		inline real velocity(const int i, const int j, const int k, const int AxisNo) const;
		real stress(const int i, const int j, const int k, const int No) const;

/*
______________________________________________________________________________________________

PRIVATE VARIABLES
______________________________________________________________________________________________

*/
	private:
		int	LocalQuantityNb;		// Number of time-average quantities (ex: u, v, w, u'u', u'v', etc ...)
		int LocalScaleNb; 			// Number of scales for this grid
		int ElementNb;					// 2^(Dimension*LocalScaleNb)
		int SampleNb;						// Number of samples for the time-averaging procedure
		Vector *Q;    					// Vector containing the time-average quantities
};

/*
______________________________________________________________________________________________

INLINE FUNCTIONS
______________________________________________________________________________________________

Get density at the position i,j,k
______________________________________________________________________________________________

*/
inline real TimeAverageGrid::density(const int i, const int j, const int k) const
{
	return value(i,j,k,1);
}

/*
______________________________________________________________________________________________

Update number of samples
______________________________________________________________________________________________

*/
inline void TimeAverageGrid::updateSample()
{
	SampleNb++;
}

/*
______________________________________________________________________________________________

Get velocity at the position i,j,k
______________________________________________________________________________________________

*/
inline real TimeAverageGrid::velocity(const int i, const int j, const int k, const int AxisNo) const
{
	return value(i,j,k,AxisNo+1)/value(i,j,k,1);
}

/*
______________________________________________________________________________________________

Get value at the position ElementNo
______________________________________________________________________________________________

*/

inline real TimeAverageGrid::value(const int ElementNo, const int QuantityNo) const
{
	return (Q+ElementNo)->value(QuantityNo);
}

/*
______________________________________________________________________________________________

Get value at the position i,j,k
______________________________________________________________________________________________

*/
inline real TimeAverageGrid::value(const int i, const int j, const int k, const int QuantityNo) const
{
	return value(i + (1<<LocalScaleNb)*(j +(1<<LocalScaleNb)*k), QuantityNo);
}


