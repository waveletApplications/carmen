/***************************************************************************
                          Remesh.cpp  -  description
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
#include "Carmen.h"

void Remesh(Node* Root)
{
	// --- Refresh tree structure ---
//	RefreshTree(Root);

	// --- Check if tree is graded ---
	if (debug) Root->checkGradedTree();

	// --- Adapt : depending on details, refine or combine ---
	Root->adapt();

	// --- Check if tree is graded ---
	if (debug) Root->checkGradedTree();
}
