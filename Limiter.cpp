/***************************************************************************
                          Limiter.cpp  -  description
                             -------------------
    begin                : mar mar 9 2004
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

real Limiter(const real r)
{
  real Result = 0.;
  
  switch (LimiterNo)
  {
    case 0: // No limiter
      Result = 0.;
      break;
      
    case 1: // Min-Mod
      if (isnan(r))
        Result = 1.;
      else
        Result = Max(0., Min(1.,r));
      break;  
      
    case 2: // Van Albada  
      if (isnan(r))
        Result = 1.;
      else
        Result = (r<=0.)? 0.: (r*r+r)/(r*r+1.);
      break;  

      case 3: // Van Leer  
        if (isnan(r))
          Result = 2.;
        else
          Result = (r<=0.) ? 0.: (r+Abs(r))/(1.+Abs(r));
        break;  

      case 4: // Superbee 
          if (isnan(r))
            Result = 2.;
          else
            Result = (r<=0.) ? 0.: Max(0.,Max(Min(2.*r,1.),Min(r,2.)));
          break;  
      
  };
  
		return Result;
}





Vector Limiter(const Vector u, const Vector v)
{
		// Min Mod limiter 

		int LimiterNo = 3;

		Vector Result(u.dimension());
		int i;
		real x, y; 	// slopes

    for (i=1; i<=u.dimension(); i++)
		{
			x = u.value(i);
			y = v.value(i);

			switch(LimiterNo)
			{
				// MIN-MOD
				case 1:
					if (x == y)
						Result.setValue(i, 0.);
					else
						Result.setValue(i,Min(1., fabs(x)/fabs(x-y)));
					break;

        // VAN LEER
				case 3:
				default:
					if ((fabs(x) + fabs(y)) == 0.)
						Result.setValue(i, 0.);
					else
						Result.setValue(i,fabs(x)/(fabs(x)+fabs(y)));
					break;
 			};
		}

		return Result;
}
