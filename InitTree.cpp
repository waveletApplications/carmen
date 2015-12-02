/***************************************************************************
                          InitTree.cpp  -  description
                             -------------------
    begin                : mar fév 3 2004
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

Init graded tree (only for multiresolution solver)
______________________________________________________________________________________________

*/
void InitTree(Node* Root)
{
	// --- Local variables ---

	int l;	// Counter on levels

	// --- Init cell-average value in root and split it ---

	if (Recovery && UseBackup)

  		 Root->restore();
		 
	else
	{
		Root->initValue();

		// --- Add and init nodes in different levels, when necessary ---

		for (l=1; l <= ScaleNb; l++)
			Root->addLevel();
	}

	// -- Check if tree is graded ---

	if (debug) Root->checkGradedTree();

}
