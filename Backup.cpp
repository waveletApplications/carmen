/***************************************************************************
                          Backup.cpp  -  description
                             -------------------
    begin                : jeu fév 12 2004
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

Multiresolution backup
______________________________________________________________________________________________

*/
void Backup(Node* Root)
{
	Root->backup();
}

/*
______________________________________________________________________________________________

Fnite volume backup
______________________________________________________________________________________________

*/
void Backup(FineMesh* Root)
{
	Root->backup();
}
