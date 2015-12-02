/***************************************************************************
                          ViewIteration.cpp  -  description
                             -------------------
    begin                : Fri Jun 14 2002
    copyright            : (C) 2002 by Olivier Roussel and Alexei Tsigulin
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
__________________________________________________________________________

Print solution if IterationNo = PrintIt1 to PrintIt6
__________________________________________________________________________

Multiresolution solver
--------------------------------------------------------------------------
*/

void ViewIteration(Node* Root)
{
		if (IterationNo == PrintIt1)
				View(Root, "Tree.dat", "Mesh_1.dat", "Average_1.dat");

		if (IterationNo == PrintIt2)
				View(Root, "Tree.dat", "Mesh_2.dat", "Average_2.dat");

		if (IterationNo == PrintIt3)
				View(Root, "Tree.dat", "Mesh_3.dat", "Average_3.dat");

		if (IterationNo == PrintIt4)
				View(Root, "Tree.dat", "Mesh_4.dat", "Average_4.dat");

		if (IterationNo == PrintIt5)
				View(Root, "Tree.dat", "Mesh_5.dat", "Average_5.dat");

		if (IterationNo == PrintIt6)
				View(Root, "Tree.dat", "Mesh_6.dat", "Average_6.dat");
}

/*
--------------------------------------------------------------------------
Finite volume solver
--------------------------------------------------------------------------
*/

void ViewIteration(FineMesh* Root)
{
		if (IterationNo == PrintIt1)
				View(Root, "Average_1.dat");

		if (IterationNo == PrintIt2)
				View(Root, "Average_2.dat");

		if (IterationNo == PrintIt3)
				View(Root, "Average_3.dat");

		if (IterationNo == PrintIt4)
				View(Root, "Average_4.dat");

		if (IterationNo == PrintIt5)
				View(Root, "Average_5.dat");

		if (IterationNo == PrintIt6)
				View(Root, "Average_6.dat");
}
