/***************************************************************************
                          Node.cpp  -  description
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

/*
______________________________________________________________________________________________

Create static elements
______________________________________________________________________________________________

*/
Node* Node::Root;
//Node  Node::Boundary;
/*
______________________________________________________________________________________________

Constructor
______________________________________________________________________________________________

*/
Node::Node (const int l, const int i, const int j, const int k)
{
	// Set Nl, Ni, Nj, Nk
	Nl = l;
	Ni = i;
	Nj = j;
	Nk = k;

	// --- If l = 0, then node is root ---

	if (Nl == 0)
	{
		Root = this;
   		CellNb = 1;
		LeafNb = 1;
	}

	// --- Increase the total number of cells  ---

	CellNb ++;

	// --- Allocate array of pointers to children ---

	Child = new Node* [ChildNb];

	// --- A new node is a simple leaf ---

	setSimpleLeaf();

	// --- Set the coordinates and the size of the cell ---

  // x-direction
	ThisCell.setSize( 1, (XMax[1]-XMin[1])/(1<<Nl) );
	ThisCell.setCenter( 1, XMin[1] + (Ni + .5)*ThisCell.size(1) );

  // y-direction
	if (Dimension > 1)
	{
		ThisCell.setSize( 2, (XMax[2]-XMin[2])/(1<<Nl) );
		ThisCell.setCenter( 2, XMin[2] + (Nj + .5)*ThisCell.size(2) );
	}

  // z-direction
	if (Dimension > 2)
	{
		ThisCell.setSize( 3, (XMax[3]-XMin[3])/(1<<Nl) );
		ThisCell.setCenter( 3, XMin[3] + (Nk + .5)*ThisCell.size(3) );
	}
}
/*
______________________________________________________________________________________________

Distructor
______________________________________________________________________________________________

*/
Node::~Node()
{
	// --- Local variables --------------------------------------------------------------------

	int n=0;  // Counter on children

	// --- Distructor procedure ---------------------------------------------------------------

	// --- Decrease the total number of cells  ---

	CellNb --;

	// --- If the node has children, delete them ---

	if (hasChildren())
		for (n = 0; n < ChildNb; n++) delete Child[n];

	// --- Delete array of pointers to children ---

	delete[] Child;
}

/*
______________________________________________________________________________________________

Init cell-average value in node with initial condition
______________________________________________________________________________________________

*/
void Node::initValue()
{
	int i,j,k; // Counters on directions

	ThisCell.setAverageZero();

	if (UseBoundaryRegions && isInsideBoundary())
	{
		ThisCell.setAverage(InitAverage(	ThisCell.center(1),
							(Dimension >1)? ThisCell.center(2):0.,
							(Dimension > 2)? ThisCell.center(3):0. ) );
	}
	else
	{
	if(centredIC){
		switch (Dimension)
		{
			case 1:
				for (i=0;i<=1;i++)
				ThisCell.setAverage( ThisCell.average()+.5* InitAverage(
					ThisCell.center(1)+(i-0.5)*ThisCell.size(1)
				) );
				break;

			case 2:
				for (i=0;i<=1;i++)
				for (j=0;j<=1;j++)
				ThisCell.setAverage( ThisCell.average()+.25* InitAverage(
					ThisCell.center(1)+(i-0.5)*ThisCell.size(1),
					ThisCell.center(2)+(j-0.5)*ThisCell.size(2)
				) );
				break;

			case 3:
				for (i=0;i<=1;i++)
				for (j=0;j<=1;j++)
				for (k=0;k<=1;k++)
				ThisCell.setAverage( ThisCell.average()+.125* InitAverage(
					ThisCell.center(1)+(i-0.5)*ThisCell.size(1),
					ThisCell.center(2)+(j-0.5)*ThisCell.size(2),
					ThisCell.center(3)+(k-0.5)*ThisCell.size(3)
				) );
				break;
		};
            }
            else{
		switch (Dimension)
		{
			case 1:
				ThisCell.setAverage(InitAverage(
					ThisCell.center(1)) );
				break;

			case 2:
				ThisCell.setAverage(InitAverage(
					ThisCell.center(1),ThisCell.center(2) ) );
				break;

			case 3:
				ThisCell.setAverage(InitAverage(
					ThisCell.center(1), ThisCell.center(2),
					ThisCell.center(3)) );
				break;
		};

	}


	}
}
/*
______________________________________________________________________________________________

Add one more level when necessary
______________________________________________________________________________________________

*/
void Node::addLevel()
{
	// --- Local variables --------------------------------------------------------------------

	int  n=0; 	// Counter on children

	// --- If level higher or equal to maximum scale number allowed, do not split

	if ( Nl >= ScaleNb ) return;

	// --- If node is not a leaf, recurse on children

	if (isInternalNode())
	{
		for (n = 0 ; n < ChildNb; n++)
			Child[n]->addLevel();
	}
  	else
	{
		// If it is a virtual leaf, no splitting
		if (isVirtualLeaf()) return;

		// If it is on the wall, always split
		if (UseBoundaryRegions && isOnBoundary())
			split(true);

		// Test on prediction error (always split on levels 0 and 1)
		if (!detailIsSmall() || Nl <= 1)
			split(true);

   	}
}
/*
______________________________________________________________________________________________

Split
______________________________________________________________________________________________

*/
void Node::split(const bool init)
{
	// --- Local variables --------------------------------------------------------------------

	int  n; 				// Counter on children
	int  i,j,k,l;  	// Children position
	int li, lj, lk;	// Idem, but taking into account boundary conditions

	int ei = 1;
	int ej = (Dimension > 1)? 1:0;
	int ek = (Dimension > 2)? 1:0;

	Node *t;	// temporary pointer to node

	// --- If level higher or equal to maximum scale number allowed, do not split

	if ( Nl >= ScaleNb ) return;
	if (!isLeaf()) return;

	// --- Graded tree: if the neighbours do not exist or are virtual, create them by splitting the father of the neighbour --

  	for (i=-ei; i<=ei; i++)
  	for (j=-ej; j<=ej; j++)
  	for (k=-ek; k<=ek; k++)
  	{
		if (i!=0 || j!=0 || k!=0)
		{
			li = (CMin[1]==3) ? (Ni+i+(1<<Nl))%(1<<Nl) : Ni+i;
			lj = (CMin[2]==3) ? (Nj+j+(1<<Nl))%(1<<Nl) : Nj+j;
			lk = (CMin[3]==3) ? (Nk+k+(1<<Nl))%(1<<Nl) : Nk+k;

			t= node(Nl,li,lj,lk);
			if (t != 0 && t->Flag==3)
				node(Nl-1, (li+4)/2-2, (lj+4)/2-2, (lk+4)/2-2)->split(init);
		}
  }

	// --- If a leaf has neighbours, create virtual children in neighbours (when not existing) ---

 	for (i=-ei; i<=ei; i++)
  	for (j=-ej; j<=ej; j++)
  	for (k=-ek; k<=ek; k++)
  	{
		if (!(i==0 && j==0 && k==0))
		{
			li = (CMin[1]==3) ? (Ni+i+(1<<Nl))%(1<<Nl) : Ni+i;
			lj = (CMin[2]==3) ? (Nj+j+(1<<Nl))%(1<<Nl) : Nj+j;
			lk = (CMin[3]==3) ? (Nk+k+(1<<Nl))%(1<<Nl) : Nk+k;

			t= node(Nl,li,lj,lk);
			if (t != 0)
				t->makeVirtualChildren(init);
		}
  }

	// --- Splitting procedure ----------------------------------------------------------------

	switch(Flag)
	{
		// --- If the node is not a leaf or a virtual leaf, do not split ---
		case 0:
		case 3:
			break;

		// --- If the node is a simple leaf, allocate the children as leaves and compute values ---
		case 1:
			LeafNb += ChildNb-1;
			for (n = 0; n < ChildNb; n++)
			{
				// Compute (l,i,j,k) of children
				l = Nl + 1;
				i = 2*Ni + n%2;
				j = 2*Nj + (n/2) % 2;
				k = 2*Nk + n/4;
				Child[n] = new Node(l,i,j,k);

				// Compute value of children (if init, computed from initial condition, else predicted)
				if (init)
					Child[n]->initValue();
				else
					Child[n]->ThisCell.setAverage(Child[n]->predict());
			}

			break;

		// --- If the node is a leaf with virtual children, virtual children become simple leaves ---

		case 2:
			LeafNb += ChildNb-1;
			for (n = 0; n < ChildNb; n++)
			{
				// Children become simple leaves
				Child[n]->setSimpleLeaf();
			}

			break;
	};

	// Now the node is no more a leaf
	setInternalNode();
}
/*
______________________________________________________________________________________________

Split all the tree
______________________________________________________________________________________________

*/

void Node::splitAll()
{
	// --- Local variables --------------------------------------------------------------------

	int n=0; 			// Counter on children

	// --- If the smallest scale is reached, return ---

	if (Nl >= ScaleNb) return;

	// --- Split node ---

	split(true);

	// --- Recurse on children ---

	for (n = 0; n < ChildNb; n++)
		Child[n]->splitAll();

}

/*
______________________________________________________________________________________________

Combine
______________________________________________________________________________________________

*/
void Node::combine()
{
	// --- Local variables --------------------------------------------------------------------

	int n=0; 			// Counter on children
	int i=0,  j=0,  k=0;		// Counters in each direction
	int li=0, lj=0, lk=0;		// Idem, but taking into account boundary conditions

	int ei = 1;
	int ej = (Dimension > 1)? 1:0;
	int ek = (Dimension > 2)? 1:0;

	Node *t;	// temporary pointer to node

	// --- Graded tree: check that this node is not a leaf and that its children are simple leaves

	// If the node is on a wall, do not combine
	if (UseBoundaryRegions && isOnBoundary()) return;

	// If the node is not an internal node, do not combine
	if (!isInternalNode()) return;

	// If one of the children is not a simple leaf, do not combine
	for ( n = 0; n < ChildNb; n++ )
		if (! (Child[n]->isSimpleLeaf()) ) return;

	// --- Combine procedure ------------------------------------------------------------------

	// --- Always create leaf with virtual children ---

	// Node becomes a leaf with virtual children
	setLeafWithVirtualChildren();

	// Real children become virtual
	for ( n = 0; n < ChildNb; n++ )
		Child[n]->setVirtualLeaf();

	// If the virtual children of this leaf and its neighbours are not needed, delete them.

	// 1 - For the current node

	deleteVirtualChildren();

	// 2 - For the cousins

  	for (i=-ei; i<=ei; i++)
  	for (j=-ej; j<=ej; j++)
  	for (k=-ek; k<=ek; k++)
  	{
		li = (CMin[1]==3) ? (Ni+i+(1<<Nl))%(1<<Nl) : Ni+i;
		lj = (CMin[2]==3) ? (Nj+j+(1<<Nl))%(1<<Nl) : Nj+j;
		lk = (CMin[3]==3) ? (Nk+k+(1<<Nl))%(1<<Nl) : Nk+k;

		t= node(Nl,li,lj,lk);
		if (t != 0 && t->isLeafWithVirtualChildren())
			t->deleteVirtualChildren();
  }

	// Reduce number of leaves
	LeafNb -= ChildNb-1;
}
/*
______________________________________________________________________________________________

Get pointer to cell
______________________________________________________________________________________________

*/
Cell* Node::cell(int l, int i, int j, int k) const
{
	// --- Local variables --------------------------------------------------------------------

	int  li, lj, lk;	// i,j,k used to search node
	int  n;			// 2^l

	Node *t;  		// Temporary pointer to node

	// --- Take into account boundary conditions -----------------------------------------------

	// Init li, lj,lk, n

	n = (1<<l);

	li = BC(i,1,n);
	lj = BC(j,2,n);
	lk = BC(k,3,n);

	// --- Find the node (l, li, lj, lk) ------------------------------------------------------

	t = node(l, li, lj, lk);

	// --- Return cell ------------------------------------------------------------------------

/*	if ( (i =! li && CMin[1] == 1) || (j =! lj && CMin[2] == 1) || (k =! lk && CMin[3] == 1) )
	{
		// Init boundary value
		Boundary.ThisCell.setAverage(t->ThisCell.average());

		// Take into account the boundary value in Dirichlet case

		// x-direction
		if (CMin[1] == 1 && li == 0 && i != li)
			Boundary.ThisCell.setAverage( Boundary.ThisCell.average() - 2*(t->ThisCell.average()) );
		if (CMax[1] == 1 && li == (1<<l)-1 && i != li)
			Boundary.ThisCell.setAverage( Boundary.ThisCell.average() - 2*(t->ThisCell.average()) );

		// y-direction
		if (CMin[2] == 1 && lj == 0 && j != lj)
			Boundary.ThisCell.setAverage( Boundary.ThisCell.average() - 2*(t->ThisCell.average()) );
		if (CMax[2] == 1 && lj == (1<<l)-1 && j != lj)
			Boundary.ThisCell.setAverage( Boundary.ThisCell.average() - 2*(t->ThisCell.average()) );

		// z-direction
		if (CMin[3] == 1 && lk == 0 && k != lk)
			Boundary.ThisCell.setAverage( Boundary.ThisCell.average() - 2*(t->ThisCell.average()) );
		if (CMax[3] == 1 && lk == (1<<l)-1 && k != lk)
			Boundary.ThisCell.setAverage( Boundary.ThisCell.average() - 2*(t->ThisCell.average()) );

		return &(Boundary.ThisCell);
	}
  else
	{*/
		if (t!=0)
			return &(t->ThisCell);
		else
			return 0;
//	}
}
/*
______________________________________________________________________________________________

Get pointer to node
______________________________________________________________________________________________

*/
Node* Node::node(int l, int i, int j, int k) const
{
	// --- Local variables --------------------------------------------------------------------

	Node *t;  // Temporary pointer to node
   	int  n=0;   // Counter on children

	// Init on root node
	t = Root;

	// Init i,j,k for eventually periodic case
	if (CMin[1]==3) 			i = (i+(1<<l))%(1<<l);
	if (Dimension > 1 && CMin[2]==3) 	j = (j+(1<<l))%(1<<l);
	if (Dimension > 2 && CMin[3]==3) 	k = (k+(1<<l))%(1<<l);

	// If (l,i,j,k) is out of the mesh, return 0
	if ( i < 0 || i >= (1<<l)) return 0;
	if ( j < 0 || j >= (1<<l)) return 0;
	if ( k < 0 || k >= (1<<l)) return 0;

	// Iterative procedure: find the path
	while (l > 0)
	{
		// if i is on the first half of the tree, choose the first child in the x-direction
		// if i is on the second half of the tree, choose the second child in the x-direction
		if ( i >= (1 << (l-1)) )
			n = 1;
		else
			n = 0;

		// idem for j in the y-direction
		if ( j >= (1 << (l-1)) )
			n += 2;

		// idem for k in the z-direction
		if ( k >= (1 << (l-1)) )
			n += 4;

		// if the node is already a leaf or a virtual leaf and at least one more step
		// is needed, this means that the researched node does not exist: return 0
		if ((t->isSimpleLeaf()) || (t->isVirtualLeaf()))
			return 0;

		// Recursion
		t = (Node*)(t->Child[n]);

		// Set (l,i,j,k) for the chosen child
		l--;
		i %= (1 << l);
		j %= (1 << l);
		k %= (1 << l);
	}

	// node found ! => return pointer to it
	return t;
}
/*
______________________________________________________________________________________________

Test and eventually delete virtual children
______________________________________________________________________________________________

*/
void Node::deleteVirtualChildren()
{
	int  n=0; 				// Counter on children
	int i=0,j=0,k=0;			// Counters in each direction
	int li=0, lj=0, lk=0;			// Idem, but taking into account boundary conditions

	int ei = 1;
	int ej = (Dimension > 1)? 1:0;
	int ek = (Dimension > 2)? 1:0;

	Node *t;				// temporary pointer to node

	// --- If at least one neighbour is not a leaf, keep the virtual children -----------------

  for (i=-ei; i<=ei; i++)
  for (j=-ej; j<=ej; j++)
  for (k=-ek; k<=ek; k++)
  {
		li = (CMin[1]==3) ? (Ni+i+(1<<Nl))%(1<<Nl) : Ni+i;
		lj = (CMin[2]==3) ? (Nj+j+(1<<Nl))%(1<<Nl) : Nj+j;
		lk = (CMin[3]==3) ? (Nk+k+(1<<Nl))%(1<<Nl) : Nk+k;

		t= node(Nl,li,lj,lk);
		if (t != 0 && t->isInternalNode()) return;
	}

	// --- Delete virtual children ------------------------------------------------------------

	for (n = 0; n < ChildNb; n++)
		delete Child[n];

	// Now node is a simple leaf
	setSimpleLeaf();
}
/*
______________________________________________________________________________________________

	Make virtual children
______________________________________________________________________________________________

*/
void Node::makeVirtualChildren(bool init)
{
	// --- Local variables ---

	int n=0; 				// Counter on children
	int l,i,j,k;		// Child position

	// --- Test on flag value ---

	if( isSimpleLeaf())
	{
		// Change flag
		setLeafWithVirtualChildren();

		// Create virtual children
		for (n = 0; n < ChildNb; n++)
		{
			// Compute (l,i,j,k) of child
			l = Nl + 1;
			i = 2*Ni + n%2;
			j = 2*Nj + (n/2) % 2;
			k = 2*Nk + n/4;

			// Create virtual child and compute its value
			Child[n] = new Node(l,i,j,k);
    	Child[n]->setVirtualLeaf();
			// Compute value of children (if init, computed from initial condition, else predicted)
			if (init)
				Child[n]->initValue();
			else
			{
				Child[n]->ThisCell.setAverage(Child[n]->predict());

				// If time adaptivity, predict also temporary cell-average values
				if (TimeAdaptivity && !init)
					Child[n]->ThisCell.setTempAverage(Child[n]->predictTempAverage());
			}
		}
	}
	// Elsewhere, do nothing !
}
/*
______________________________________________________________________________________________

	Fill virtual children values
______________________________________________________________________________________________

*/
void Node::fillVirtualChildren()
{
	// --- Local variables ---

	int n=0; // Counter on children

	// --- Recursion ---

	switch(Flag)
	{
  	// If node is not a leaf or leaf with virtual children, recurse on children
		case 0:
		case 2:
			for (n = 0; n < ChildNb; n++)
				Child[n]->fillVirtualChildren();
			break;

		// If node is a simple leaf, stop procedure
		case 1:
			return;
			break;

		// If node is a virtual leaf, compute value with prediction
		case 3:
			ThisCell.setAverage(predict());
			if (EquationType==6) ThisCell.setGradient(parentCell()->gradient());
			break;

	};
}
/*
______________________________________________________________________________________________

	Project : compute the mean value
______________________________________________________________________________________________

*/
Cell* Node::project()
{
	// --- Local variables ---

	int n=0;	// Counter on children

	// --- If cell is not a leaf, compute projection ie mean value of children ---

	if (isInternalNode())
	{
		// Set value to zero
		ThisCell.setAverageZero();

		// Compute the mean value
		for (n = 0; n < ChildNb; n++)
			ThisCell.setAverage( ThisCell.average() + Child[n]->project()->average() );

		ThisCell.setAverage( ThisCell.average() / ChildNb );

	}

	return &ThisCell;
}
/*
______________________________________________________________________________________________

	Adapt procedure
______________________________________________________________________________________________

*/
int Node::adapt()
{
	// --- Local variables ---

	int 	n;		// Counter on children
	int 	isDeletable;  	// Test if children are deletable (0 = all children are deletable)

	// --- In case of time adaptivity, only remesh when a complete time evolution has been done ---

	// --- Init ---

	isDeletable = 0;

	// --- Test to stop recursion ---

	// If the node is not an internal node, this node can be deleted
	if (!isInternalNode() || Nl >= ScaleNb) return 0;

	// --- Recursion ---

	// Test if children are deletable
	for (n = 0; n < ChildNb; n++)
		isDeletable += Child[n]->adapt();

	// If all children are deletable, test if this node is also deletable

  	if (isDeletable == 0)
	{
		if (detailIsSmall())
		{
			if (!TimeAdaptivity || (TimeAdaptivity && isEndTimeCycle()))
				for (n = 0; n < ChildNb; n++) Child[n]->combine();

			// Add value 0 to variable isDeletable of the parent
			return 0;

   	}
		else
		{
			if (!TimeAdaptivity || (TimeAdaptivity && isEndTimeCycle()))
				for (n = 0; n < ChildNb; n++) Child[n]->split();

			return 1;
    		}
	}
 	// Add value 1 to variable Deletable of the parent
	return 1;
}
/*
______________________________________________________________________________________________

	predict the cell-average values with a linear interpolation
______________________________________________________________________________________________

*/
Vector Node::predict() const
{
	// --- Local variables ---

	int pi, pj=1, pk=1;	// Parity of Ni, Nj, Nk

	Vector	Result(QuantityNb);

	// --- Init result with the cell-average value of the father ---

	Result = parentCell()->average();

	// --- 1D case ---

	pi = (Ni%2 == 0)?1:-1;
	Result += (pi*-.125) * uncleCell(1,0,0)->average();
	Result -= (pi*-.125) * uncleCell(-1,0,0)->average();

	// --- 2D case ---

 	if (Dimension > 1)
	{
		pj = (Nj%2 == 0)?1:-1;
		Result += (pj*-.125) * uncleCell(0,1,0)->average();
		Result -= (pj*-.125) * uncleCell(0,-1,0)->average();

		Result += (pi*pj*.015625) * uncleCell(1,1,0)->average();
		Result -= (pi*pj*.015625) * uncleCell(1,-1,0)->average();
		Result -= (pi*pj*.015625) * uncleCell(-1,1,0)->average();
		Result += (pi*pj*.015625) * uncleCell(-1,-1,0)->average();
  }

	// --- 3D case ---

 	if (Dimension > 2)
	{
		pk = (Nk%2 == 0)?1:-1;
		Result += (pk*-.125) * uncleCell(0,0,1)->average();
		Result -= (pk*-.125) * uncleCell(0,0,-1)->average();

		Result += (pi*pk*.015625) * uncleCell(1,0,1)->average();
		Result -= (pi*pk*.015625) * uncleCell(1,0,-1)->average();
		Result -= (pi*pk*.015625) * uncleCell(-1,0,1)->average();
		Result += (pi*pk*.015625) * uncleCell(-1,0,-1)->average();

		Result += (pj*pk*.015625) * uncleCell(0,1,1)->average();
		Result -= (pj*pk*.015625) * uncleCell(0,1,-1)->average();
		Result -= (pj*pk*.015625) * uncleCell(0,-1,1)->average();
		Result += (pj*pk*.015625) * uncleCell(0,-1,-1)->average();

		Result += (pi*pj*pk*-.001953125) * uncleCell(1,1,1)->average();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(1,1,-1)->average();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(1,-1,1)->average();
		Result += (pi*pj*pk*-.001953125) * uncleCell(1,-1,-1)->average();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(-1,1,1)->average();
		Result += (pi*pj*pk*-.001953125) * uncleCell(-1,1,-1)->average();
		Result += (pi*pj*pk*-.001953125) * uncleCell(-1,-1,1)->average();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(-1,-1,-1)->average();
  }

	return Result;
}
/*
______________________________________________________________________________________________

	Predict the temporary cell-average values with a linear interpolation
______________________________________________________________________________________________

*/
Vector Node::predictTempAverage() const
{
	// --- Local variables ---

	int pi, pj=1, pk=1;	// Parity of Ni, Nj, Nk

	Vector	Result(QuantityNb);

	// --- Init result with the cell-average value of the father ---

	Result = parentCell()->tempAverage();

	// --- 1D case ---

	pi = (Ni%2 == 0)?1:-1;
	Result += (pi*-.125) * uncleCell(1,0,0)->tempAverage();
	Result -= (pi*-.125) * uncleCell(-1,0,0)->tempAverage();

	// --- 2D case ---

 	if (Dimension > 1)
	{
		pj = (Nj%2 == 0)?1:-1;
		Result += (pj*-.125) * uncleCell(0,1,0)->tempAverage();
		Result -= (pj*-.125) * uncleCell(0,-1,0)->tempAverage();

		Result += (pi*pj*.015625) * uncleCell(1,1,0)->tempAverage();
		Result -= (pi*pj*.015625) * uncleCell(1,-1,0)->tempAverage();
		Result -= (pi*pj*.015625) * uncleCell(-1,1,0)->tempAverage();
		Result += (pi*pj*.015625) * uncleCell(-1,-1,0)->tempAverage();
  }

	// --- 3D case ---

 	if (Dimension > 2)
	{
		pk = (Nk%2 == 0)?1:-1;
		Result += (pk*-.125) * uncleCell(0,0,1)->tempAverage();
		Result -= (pk*-.125) * uncleCell(0,0,-1)->tempAverage();

		Result += (pi*pk*.015625) * uncleCell(1,0,1)->tempAverage();
		Result -= (pi*pk*.015625) * uncleCell(1,0,-1)->tempAverage();
		Result -= (pi*pk*.015625) * uncleCell(-1,0,1)->tempAverage();
		Result += (pi*pk*.015625) * uncleCell(-1,0,-1)->tempAverage();

		Result += (pj*pk*.015625) * uncleCell(0,1,1)->tempAverage();
		Result -= (pj*pk*.015625) * uncleCell(0,1,-1)->tempAverage();
		Result -= (pj*pk*.015625) * uncleCell(0,-1,1)->tempAverage();
		Result += (pj*pk*.015625) * uncleCell(0,-1,-1)->tempAverage();

		Result += (pi*pj*pk*-.001953125) * uncleCell(1,1,1)->tempAverage();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(1,1,-1)->tempAverage();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(1,-1,1)->tempAverage();
		Result += (pi*pj*pk*-.001953125) * uncleCell(1,-1,-1)->tempAverage();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(-1,1,1)->tempAverage();
		Result += (pi*pj*pk*-.001953125) * uncleCell(-1,1,-1)->tempAverage();
		Result += (pi*pj*pk*-.001953125) * uncleCell(-1,-1,1)->tempAverage();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(-1,-1,-1)->tempAverage();
  }

	return Result;
}

/*
______________________________________________________________________________________________

	Predict the values of the gradient with a linear interpolation
______________________________________________________________________________________________

*/
Matrix Node::predictGradient() const
{
	// --- Local variables ---

	int pi, pj=1, pk=1;	// Parity of Ni, Nj, Nk

	Matrix	Result;

	// --- Init result with the cell-average value of the father ---

	Result = parentCell()->gradient();

	// --- 1D case ---

	pi = (Ni%2 == 0)?1:-1;
	Result += (pi*-.125) * uncleCell(1,0,0)->gradient();
	Result -= (pi*-.125) * uncleCell(-1,0,0)->gradient();

	// --- 2D case ---

 	if (Dimension > 1)
	{
		pj = (Nj%2 == 0)?1:-1;
		Result += (pj*-.125) * uncleCell(0,1,0)->gradient();
		Result -= (pj*-.125) * uncleCell(0,-1,0)->gradient();

		Result += (pi*pj*.015625) * uncleCell(1,1,0)->gradient();
		Result -= (pi*pj*.015625) * uncleCell(1,-1,0)->gradient();
		Result -= (pi*pj*.015625) * uncleCell(-1,1,0)->gradient();
		Result += (pi*pj*.015625) * uncleCell(-1,-1,0)->gradient();
  }

	// --- 3D case ---

 	if (Dimension > 2)
	{
		pk = (Nk%2 == 0)?1:-1;
		Result += (pk*-.125) * uncleCell(0,0,1)->gradient();
		Result -= (pk*-.125) * uncleCell(0,0,-1)->gradient();

		Result += (pi*pk*.015625) * uncleCell(1,0,1)->gradient();
		Result -= (pi*pk*.015625) * uncleCell(1,0,-1)->gradient();
		Result -= (pi*pk*.015625) * uncleCell(-1,0,1)->gradient();
		Result += (pi*pk*.015625) * uncleCell(-1,0,-1)->gradient();

		Result += (pj*pk*.015625) * uncleCell(0,1,1)->gradient();
		Result -= (pj*pk*.015625) * uncleCell(0,1,-1)->gradient();
		Result -= (pj*pk*.015625) * uncleCell(0,-1,1)->gradient();
		Result += (pj*pk*.015625) * uncleCell(0,-1,-1)->gradient();

		Result += (pi*pj*pk*-.001953125) * uncleCell(1,1,1)->gradient();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(1,1,-1)->gradient();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(1,-1,1)->gradient();
		Result += (pi*pj*pk*-.001953125) * uncleCell(1,-1,-1)->gradient();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(-1,1,1)->gradient();
		Result += (pi*pj*pk*-.001953125) * uncleCell(-1,1,-1)->gradient();
		Result += (pi*pj*pk*-.001953125) * uncleCell(-1,-1,1)->gradient();
		Result -= (pi*pj*pk*-.001953125) * uncleCell(-1,-1,-1)->gradient();
  }

	return Result;
}

/*
______________________________________________________________________________________________

	Print tree
______________________________________________________________________________________________

*/
void Node::writeTree(const char* FileName) const
{
	// --- Local variables ---

	int n, l;			// Counter

	FILE *output;	// Pointer to output file

	// --- Open file ---

	if ( (Nl == 0) ? (output = fopen(FileName,"w")) : (output = fopen(FileName,"a")) )
	{
		for (l = 1; l <= Nl; l++)
			fprintf(output, "|" );

		fprintf(output, "+ ");
		switch (Dimension)
		{
			case 1:
				fprintf(output, "(%d, %d)", Nl, Ni );
				break;

			case 2:
				fprintf(output, "(%d, %d, %d)", Nl, Ni, Nj);
				break;

			case 3:
				fprintf(output, "(%d, %d, %d, %d)", Nl, Ni, Nj, Nk );
				break;
		};
		switch(Flag)
		{
			case 0:
				fprintf(output," -- node --");
				break;

			case 1:
				fprintf(output," -- leaf --");
				break;

			case 2:
				fprintf(output," -- leaf with virtual children --");
				break;

			case 3:
				fprintf(output," -- virtual leaf --");
				break;

		};
		fprintf(output,"\n");
		fclose(output);
		if (hasChildren())
		{
			for (n = 0; n < ChildNb; n++)
				Child[n]->writeTree(FileName);
		}
	}
	else
	{
		cout << "Node.cpp: In method `void Node::writeText()':\n";
		cout << "Node.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [Node.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
/*
______________________________________________________________________________________________

	write header for Data Explorer visualization
______________________________________________________________________________________________

*/
void Node::writeHeader(const char* FileName) const
{
	// --- Local variables ---

	FILE 	*output;	// Pointer to output file

	// --- Open file ---

	if ((output = fopen(FileName,"w")) != NULL)
	{
		// --- Header ---

		if (Dimension == 1)
		{
			// GNUPLOT
			fprintf(output,"#");
			fprintf(output, TXTFORMAT, " x");
			switch(EquationType)
			{
        			// LINEAR ADVECTION AND BURGERS
				case 1:
				case 2:
					fprintf(output, TXTFORMAT, "u");
					break;

				// FLAME BALL, FLAME FRONT, INTERACTION FLAME-CURL
				case 3:
				case 4:
				case 5:
					fprintf(output, TXTFORMAT, "Temperature");
					fprintf(output, TXTFORMAT, "Concentration");
					fprintf(output, TXTFORMAT, "Reaction rate");
					break;

				// NAVIER-STOKES
				case 6:
					fprintf(output, TXTFORMAT, "Density");
					fprintf(output, TXTFORMAT, "Pressure");
					fprintf(output, TXTFORMAT, "Temperature");
					fprintf(output, TXTFORMAT, "Energy");
					fprintf(output, TXTFORMAT, "Vorticity");
					fprintf(output, TXTFORMAT, "Velocity");
					break;
			};
			fprintf(output, "\n");
		}
		else
    		{
			fprintf(output, "# Data Explorer file\n# generated by Carmen %3.1f\n", CarmenVersion);
			fprintf(output, "points = %d\n",LeafNb);
			fprintf(output, "format = ascii\n");
			fprintf(output, "interleaving = field\n");

			fprintf(output, "field = locations, Q1\n");
			fprintf(output, "structure = %d-vector, scalar\n", Dimension);
			fprintf(output, "type = %s, %s\n", REAL, REAL);
			fprintf(output, "dependency = positions, positions\n");

			fprintf(output, "header = marker \"START_DATA\\n\" \n");
			fprintf(output, "end\n");
			fprintf(output, "START_DATA\n");
    }

		fclose(output);
		return;
	}
	else
	{
		cout << "Node.cpp: In method `void Node::writeHeader()':\n";
		cout << "Node.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [Node.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
/*
______________________________________________________________________________________________

	write cell-average values for graphic visualization
______________________________________________________________________________________________

*/
void Node::writeAverage(const char* FileName)
{
	// --- Local variables ---

	int 		n ;		// Counter on children
	FILE 		*output;	// Pointer to output file
	Vector Qbuf(QuantityNb);  	// Buffer of the vector of conservative quantities
	real a=1., b=1.;		// Weights for the synchronisation of conservative quantities
	int	Nf=1;

	// --- Open file ---

	if ((output = fopen(FileName,"a")) != NULL)
	{
		// If node is not a leaf, recurse to children
		if (isInternalNode())
		{
			fclose(output);
			for (n = 0; n < ChildNb; n++)
				Child[n]->writeAverage(FileName);
		}
		else
		{
			// In case of local time stepping, store value and synchronize it
			if ( TimeAdaptivity && Dimension == 1 && (Nl < (ScaleNb-1)) )
			{
				Qbuf = ThisCell.average();

				Nf = 1<<TimeAdaptivityFactor*(ScaleNb-Nl);
				a = 1.*(IterationNo%Nf - Nf/2);
				b = 1.*(Nf - IterationNo%Nf);

				if ( (a+b) != 0. )
				{
			//		  Qbuf2 = ( b/(a+b) )* ThisCell.average()+ ( a/(a+b) )* ThisCell.lowAverage();

                             //              ThisCell.setAverage(Qbuf2);

//ThisCell.setAverage(( b/(a+b) )* ThisCell.average()+ ( a/(a+b) )* ThisCell.lowAverage()) ;

                                }
			}

			// write cell-average values and details from parent

			fprintf(output, FORMAT, ThisCell.center(1));
			if (Dimension > 1) fprintf(output, FORMAT, ThisCell.center(2));
			if (Dimension > 2) fprintf(output, FORMAT, ThisCell.center(3));

                        if (Dimension > 1)
                          fprintf(output, "%i", Nl);

                        else
                        {
			   switch(EquationType)
			   {
				case 1:
				case 2:
					// LINEAR ADVECTION AND BURGERS
					fprintf(output, FORMAT, ThisCell.average(1));
					break;

				case 3:
				case 4:
					// FLAME FRONT AND FLAME BALL
					fprintf(output, FORMAT, ThisCell.temperature());
                                        fprintf(output, FORMAT, ThisCell.concentration());
				        fprintf(output, FORMAT, ReactionRate(ThisCell.temperature(),ThisCell.concentration()));
					break;

				case 5:
					// INTERACTION FLAME - CURL
					fprintf(output, FORMAT, ThisCell.temperature());
					break;

				case 6:
					// NAVIER-STOKES
					if (ScalarEqNb == 1)
						fprintf(output, FORMAT, ThisCell.average(Dimension+3)/ThisCell.density());
					else
						fprintf(output, FORMAT, ThisCell.density());

					fprintf(output, FORMAT, ThisCell.pressure()*Gamma*Ma*Ma);
					fprintf(output, FORMAT, ThisCell.temperature());
					fprintf(output, FORMAT, ThisCell.energy());
					fprintf(output, FORMAT, 0.);
					fprintf(output, FORMAT, ThisCell.velocity(1));
					break;
                              };
                        }

			fprintf(output,"\n");

			// In case of local time stepping, return to the previous value
			if (TimeAdaptivity && Dimension == 1 && Nl < ScaleNb-1)
				ThisCell.setAverage(Qbuf);

			fclose(output);
		}
	}
	else
	{
		cout << "Node.cpp: In method `void Node::writeAverage()':\n";
		cout << "Node.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [Node.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
/*
______________________________________________________________________________________________

	write cell-average values for graphic visualization on grid of level L
______________________________________________________________________________________________

*/
void Node::writeFineGrid(const char* FileName, const int L) const
{

	// --- Declarations ------------------------------------------------------------------------

	int l=0,i=0,j=0,k=0,n=0; 	// counters

	int ej = (Dimension > 1)? 1:0;
	int ek = (Dimension > 2)? 1:0;

	real x=0., y=0., z=0., t=0.;	// Cell centers and time
	real dx=0., dy=0., dz=0.;	// Cell sizes
	int GridPoints;			// Grid points
	char DependencyType[12];	// positions or connections
	
	PrintGrid FineGrid(L);

	FILE	*output, *outputHeader;
        char	FileNameHeader[25] ;

        sprintf(FileNameHeader,"HeaderAverage.general");
        

	// --- Execution ---------------------------------------------------------------------------

	// --- Compute grid points and set dependency type ---
	
	if (WriteAsPoints)
	{
		GridPoints = (1<<L);
		sprintf(DependencyType,"positions");
	}
	else
	{
		GridPoints = (1<<L)+1;
		sprintf(DependencyType,"connections");
	}
	
	// --- Compute t, dx, dy, dz ---
	
        t = ElapsedTime;

	dx = (XMax[1]-XMin[1])/(1<<L);

	if (Dimension > 1)
		dy = (XMax[2]-XMin[2])/(1<<L);

	if (Dimension > 2)
		dz = (XMax[3]-XMin[3])/(1<<L);


	// --- Compute result on fine mesh ---

  for (l=0; l<=L; l++)
	{
			for (i = 0; i <= ((1<<l)-1); i++)
			for (j = 0; j <= ((1<<l)-1)*ej; j++)
			for (k = 0; k <= ((1<<l)-1)*ek; k++)
      {
				if (node(l,i,j,k)==0)
					FineGrid.predict(l,i,j,k);
				else
					FineGrid.setValue(i,j,k,cell(l,i,j,k)->average());
      }
			FineGrid.refresh();
	}


	// --- Open file ---

	if ((outputHeader = fopen(FileNameHeader,"w")) != NULL)
	{
 	    output = fopen(FileName,"w");
     		// --- Header ---

		switch(PostProcessing)
		{
			// GNUPLOT
			case 1:
				fprintf(outputHeader,"#");
				fprintf(outputHeader, TXTFORMAT, " x");
				switch(EquationType)
				{
        				// LINEAR ADVECTION AND BURGERS
					case 1:
					case 2:
						fprintf(outputHeader, TXTFORMAT, "u");
						break;

					// FLAME BALL, FLAME FRONT, INTERACTION FLAME-CURL
					case 3:
					case 4:
					case 5:
						fprintf(outputHeader, TXTFORMAT, "Temperature");
						fprintf(outputHeader, TXTFORMAT, "Concentration");
						fprintf(outputHeader, TXTFORMAT, "Reaction rate");
						break;

					// NAVIER-STOKES
					case 6:
						fprintf(outputHeader, TXTFORMAT, "Density");
						fprintf(outputHeader, TXTFORMAT, "Pressure");
						fprintf(outputHeader, TXTFORMAT, "Temperature");
						fprintf(outputHeader, TXTFORMAT, "Energy");
						fprintf(outputHeader, TXTFORMAT, "Vorticity");
						fprintf(outputHeader, TXTFORMAT, "Velocity");
						break;
				};
				fprintf(outputHeader, "\n");
				break;

      			// DATA EXPLORTER
			case 2:
				// Header for Data explorer

				fprintf(outputHeader, "# Data Explorer file\n# generated by Carmen %3.1f\n", CarmenVersion);

				switch(Dimension)
				{
					case 2:
						fprintf(outputHeader, "grid = %d x %d\n", GridPoints, GridPoints);
						fprintf(outputHeader, "positions = %f, %f, %f, %f\n#\n", XMin[1], dx, XMin[2], dy);
						break;

					case 3:
						fprintf(outputHeader, "grid = %d x %d x %d\n", GridPoints, GridPoints, GridPoints);
						fprintf(outputHeader, "positions = %f, %f, %f, %f, %f, %f\n#\n", XMin[1], dx, XMin[2], dy, XMin[3], dz );
						break;
    				};
				if (DataIsBinary)
					fprintf(outputHeader, "format = binary\n");
				else
					fprintf(outputHeader, "format = ascii\n");

				fprintf(outputHeader, "interleaving = field\n");

				switch(EquationType)
				{
					case 1:
					case 2:
						// LINEAR ADVECTION AND BURGERS
						fprintf(outputHeader, "field = velocity\n");
						fprintf(outputHeader, "structure = scalar\n");
						fprintf(outputHeader, "type = %s\n", REAL);
						fprintf(outputHeader, "dependency = %s\n", DependencyType);
						break;

					case 3:
					case 4:
						// FLAME FRONT AND FLAME BALL
						fprintf(outputHeader, "field = temperature, concentration, reaction\n");
						fprintf(outputHeader, "structure = scalar, scalar, scalar\n");
						fprintf(outputHeader, "type = %s, %s, %s\n", REAL, REAL, REAL);
						fprintf(outputHeader, "dependency = %s, %s, %s\n", DependencyType, DependencyType, DependencyType);
      						break;

					case 5:
						// INTERACTION FLAME - CURL
						fprintf(outputHeader, "field = temperature, concentration, reaction, velocity\n");
						fprintf(outputHeader, "structure = scalar, scalar, scalar, 2-vector\n");
						fprintf(outputHeader, "type = %s, %s, %s, %s\n", REAL, REAL, REAL, REAL);
						fprintf(outputHeader, "dependency = %s, %s, %s, %s\n", DependencyType, DependencyType, DependencyType, DependencyType);
      						break;

					case 6:
						// NAVIER-STOKES
						fprintf(outputHeader, "field = density, pressure, temperature, energy, velocity\n");
						fprintf(outputHeader, "structure = scalar, scalar, scalar, scalar, %d-vector \n",Dimension);
						fprintf(outputHeader, "type = %s, %s, %s, %s, %s\n", REAL, REAL, REAL, REAL, REAL);
						fprintf(outputHeader, "dependency = %s, %s, %s, %s, %s\n", DependencyType, DependencyType, DependencyType, DependencyType, DependencyType);
      						break;

				};

				fprintf(outputHeader, "header = marker \"START_DATA\\n\" \n");
				fprintf(outputHeader, "end\n");
				fprintf(outputHeader, "START_DATA\n");

				break;

			// TECPLOT
			case 3:
				fprintf(output, "VARIABLES = \"x\"\n");
				if (Dimension > 1)
					fprintf(output,"\"y\"\n");
				if (Dimension > 2)
					fprintf(output,"\"z\"\n");

				switch(EquationType)
				{
					// LINEAR ADVECTION AND BURGERS
					case 1:
					case 2:
						fprintf(output,"\"U\"\n ");
						break;

					// FLAME FRONT AND FLAME BALL
					case 3:
					case 4:
						fprintf(output,"\"T\"\n\"C\"\n\"OMEGA\"\n");
						break;

					// INTERACTION FLAME - CURL
          				case 5:
						fprintf(output,"\"T\"\n\"C\"\n\"OMEGA\"\n\"U\"\n\"V\"");
						break;

					// NAVIER-STOKES
					case 6:
						fprintf(output,"\"RHO\"\n\"P\"\n\"T\"\n\"E\"\n\"U\"\n");
						if (Dimension > 1)
							fprintf(output,"\"V\"\n");
						if (Dimension > 2)
							fprintf(output,"\"W\"\n");
						break;
				};

				fprintf(output,"ZONE T=\"Carmen %3.1f\"\n",CarmenVersion);
				fprintf(output,"I=%i, ",(1<<L));
				if (Dimension > 1)
					fprintf(output,"J=%i, ",(1<<L));
				if (Dimension > 2)
					fprintf(output,"K=%i, ",(1<<L));
				fprintf(output,"F=POINT\n");
				break;

		};

		// --- write values ---

		for (n=0; n < (1<<(Dimension*L)); n++)
   	{

			// -- Compute i, j, k --

			// For Gnuplot and DX, loop order: for i... {for j... {for k...} }
			if (PostProcessing != 3)
			{
				switch(Dimension)
				{
					case 1:
						i = n;
						j = k = 0;
						break;

					case 2:
						j = n%(1<<L);
						i = n/(1<<L);
						k = 0;
						break;

					case 3:
						k = n%(1<<L);
            					j = (n%(1<<(2*L)))/(1<<L);
						i =  n/(1<<(2*L));
						break;
				};
			}
			else
			{
				// For Tecplot, loop order: for k... {for j... {for i...} }
				i = n%(1<<L);
				if (Dimension > 1)
					j = (n%(1<<(2*L)))/(1<<L);
				else
					j = 0;
				if (Dimension > 2)
					k = n/(1<<(2*L));
				else
					k = 0;
			}

			// Compute x, y, z

			x  = XMin[1]+i*dx;
			if (Dimension > 1) y  = XMin[2]+j*dy;
			if (Dimension > 2) z  = XMin[3]+k*dz;

			// For Tecplot, write coordinates
			if (PostProcessing == 3)
			{
					FileWrite(output, FORMAT, x);
					if (Dimension > 1) FileWrite(output, FORMAT, y);
					if (Dimension > 2) FileWrite(output, FORMAT, z);
			}

			switch(EquationType)
			{
				case 1:
				case 2:
				  // ADVECTION-DIFFUSION AND BURGERS
					FileWrite(output, FORMAT, FineGrid.value(i,j,k,1));
					break;

				case 3:
				case 4:
				case 5:
					// FLAME FRONT AND FLAME BALL
					FileWrite(output, FORMAT, FineGrid.temperature(i,j,k));
					FileWrite(output, FORMAT, FineGrid.concentration(i,j,k));
					FileWrite(output, FORMAT, ReactionRate(FineGrid.temperature(i,j,k),FineGrid.concentration(i,j,k)));

					if (EquationType == 5)
					{
						FileWrite(output, FORMAT, CurlVelocity(x, y, t, 1));
						FileWrite(output, FORMAT, CurlVelocity(x, y, t, 2));
					}
					break;

				case 6:
					// NAVIER-STOKES
					if (ScalarEqNb == 1)
						FileWrite(output, FORMAT, FineGrid.value(i,j,k,Dimension+3));
					else
						FileWrite(output, FORMAT, FineGrid.density(i,j,k));

					FileWrite(output, FORMAT, FineGrid.pressure(i,j,k)*Gamma*Ma*Ma);
					FileWrite(output, FORMAT, FineGrid.temperature(i,j,k));
					FileWrite(output, FORMAT, FineGrid.energy(i,j,k));

					if (PostProcessing == 1)
						FileWrite(output, FORMAT, FineGrid.vorticity(i,j,k));

  					for (int AxisNo = 1; AxisNo <= Dimension; AxisNo++)
						FileWrite(output, FORMAT, FineGrid.velocity(i,j,k,AxisNo));

					break;
      };

			// For ASCII data, add a return at the end of the line
			if (!DataIsBinary)
				fprintf(output,"\n");

			// For Gnuplot, add empty lines when j=jmax or k=kmax
			if (PostProcessing == 1)
			{
				if (j==(1<<ScaleNb)-1)
					fprintf(output,"\n");

				if (k==(1<<ScaleNb)-1)
					fprintf(output,"\n");
			}
   	}
		fclose(output);
                fclose(outputHeader);
	}
	else
	{
		cout << "Node.cpp: In method `void writeFineGrid(Node*, char*, int)':\n";
		cout << "Node.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [Node.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}/*
______________________________________________________________________________________________

	write mesh in file for graphic visualization
______________________________________________________________________________________________

*/
void Node::writeMesh(const char* FileName) const
{
	// --- Local variables ---

	int 	n;				// Counter on children
	FILE 	*output;	// Pointer to output file

	// --- Open file ---

	if ( (Nl == 0) ? (output = fopen(FileName,"w")) : (output = fopen(FileName,"a")) )
	{
		if (isInternalNode())
		{
			for (n = 0; n < ChildNb; n++)
				Child[n]->writeMesh(FileName);
		}
		else
		{
			// x-direction
			fprintf(output, FORMAT, ThisCell.center(1)-.5*ThisCell.size(1));
			if (Dimension >1) fprintf(output, FORMAT, ThisCell.center(2)-.5*ThisCell.size(2));
			if (Dimension >2) fprintf(output, FORMAT, ThisCell.center(3)-.5*ThisCell.size(3));
			fprintf(output, "%d", Nl);
			fprintf(output,"\n");

			fprintf(output, FORMAT, ThisCell.center(1)+.5*ThisCell.size(1));
			if (Dimension >1) fprintf(output, FORMAT, ThisCell.center(2)-.5*ThisCell.size(2));
			if (Dimension >2) fprintf(output, FORMAT, ThisCell.center(3)-.5*ThisCell.size(3));
			fprintf(output, "%d", Nl);
			fprintf(output,"\n\n");

			// y-direction
			if (Dimension > 1)
			{
				fprintf(output, FORMAT, ThisCell.center(1)-.5*ThisCell.size(1));
				fprintf(output, FORMAT, ThisCell.center(2)+.5*ThisCell.size(2));
				if (Dimension >2) fprintf(output, FORMAT, ThisCell.center(3)-.5*ThisCell.size(3));
				fprintf(output, "%d", Nl);
				fprintf(output,"\n");

				fprintf(output, FORMAT, ThisCell.center(1)+.5*ThisCell.size(1));
				fprintf(output, FORMAT, ThisCell.center(2)+.5*ThisCell.size(2));
				if (Dimension >2) fprintf(output, FORMAT, ThisCell.center(3)-.5*ThisCell.size(3));
				fprintf(output, "%d", Nl);
				fprintf(output,"\n\n");
			}

			// z-direction
			if (Dimension > 2)
			{
				fprintf(output, FORMAT, ThisCell.center(1)-.5*ThisCell.size(1));
				fprintf(output, FORMAT, ThisCell.center(2)-.5*ThisCell.size(2));
				fprintf(output, FORMAT, ThisCell.center(3)+.5*ThisCell.size(3));
				fprintf(output, "%d", Nl);
				fprintf(output,"\n");

				fprintf(output, FORMAT, ThisCell.center(1)+.5*ThisCell.size(1));
				fprintf(output, FORMAT, ThisCell.center(2)-.5*ThisCell.size(2));
				fprintf(output, FORMAT, ThisCell.center(3)+.5*ThisCell.size(3));
				fprintf(output, "%d", Nl);
				fprintf(output,"\n\n");

				fprintf(output, FORMAT, ThisCell.center(1)-.5*ThisCell.size(1));
				fprintf(output, FORMAT, ThisCell.center(2)+.5*ThisCell.size(2));
				fprintf(output, FORMAT, ThisCell.center(3)+.5*ThisCell.size(3));
				fprintf(output, "%d", Nl);
				fprintf(output,"\n");

				fprintf(output, FORMAT, ThisCell.center(1)+.5*ThisCell.size(1));
				fprintf(output, FORMAT, ThisCell.center(2)+.5*ThisCell.size(2));
				fprintf(output, FORMAT, ThisCell.center(3)+.5*ThisCell.size(3));
				fprintf(output, "%d", Nl);
				fprintf(output,"\n\n\n");
			}
		fprintf(output,"\n");
    }
		fclose(output);
	}
	else
	{
		cout << "Node.cpp: In method `void Node::writeMesh()':\n";
		cout << "Node.cpp: cannot open file " << FileName << '\n';
		cout << "carmen: *** [Node.o] Execution error\n";
		cout << "carmen: abort execution.\n";
		exit(1);
	}
}
/*
______________________________________________________________________________________________

	Store cell-average values into temporary cell-average
______________________________________________________________________________________________

*/
void Node::store()
{
	// --- Local variables ---

	int n=0;	// Counter on children

	// --- Store cell-average value Q into Qs ---

	if ((EquationType==6 && SchemeNb > 5 ) || requiresTimeEvolution() || IterationNo == 1)
	{
		if (UseBoundaryRegions)
		{
			if (IterationNo == 1)
				ThisCell.setOldAverage(ThisCell.average());
			else
				ThisCell.setOldAverage(ThisCell.tempAverage());
		}

		ThisCell.setTempAverage(ThisCell.average());
	}

	// --- Recursion in nodes that have children (real or virtual) ---

	if (hasChildren())
	{
			for (n = 0; n < ChildNb; n++)
				Child[n]->store();
	}
}

/*
______________________________________________________________________________________________

	Store gradient-average values into temporary gradient-average
______________________________________________________________________________________________
*/

void Node::storeGrad()
{
	// --- Local variables ---

	int n=0;	// Counter on children

	// --- Store gradient-average value Grad into Grads ---

	if ((EquationType==6 && SchemeNb > 5) || requiresTimeEvolution() || IterationNo == 1)
		ThisCell.setTempGradient(ThisCell.gradient());

	// --- Recursion in nodes that have children (real or virtual) ---

	if (hasChildren())
	{
		for (n = 0; n < ChildNb; n++)
			Child[n]->storeGrad();
	}
}

/*
______________________________________________________________________________________________

	Compute Divergence
______________________________________________________________________________________________

*/
void Node::computeDivergence()
{
	// --- Local variables ---

	int n=0;									// Counter on children
	real XIn, XOut;						// Borders of the cell
	Vector FluxIn, FluxOut;		// Ingoing and outgoing flux
	Vector InitAverage0;			// Cell-average value of the initial condition

	// --- Computation ---

	if (requiresDivergenceComputation())
	{
		// --- Compute source term -------------------------------------------------------------

		ThisCell.setDivergence(Source(ThisCell));

		// --- Add flux in x-direction ---------------------------------------------------------

		// If the cell is a leaf with virtual children and its left cousin is a node, compute flux on upper level

		if (isLeafWithVirtualChildren() && node(Nl, Ni-1, Nj, Nk) != 0 && node(Nl, Ni-1, Nj, Nk)->isInternalNode() && FluxCorrection)
		{
			FluxIn  = Flux( *childCell(-2,0,0), *childCell(-1,0,0), *childCell(0,0,0) , *childCell(1,0,0), 1 );

			if (Dimension > 1)
				FluxIn += Flux( *childCell(-2,1,0), *childCell(-1,1,0), *childCell(0,1,0), *childCell(1,1,0), 1 );

			if (Dimension > 2)
			{
				FluxIn += Flux( *childCell(-2,0,1), *childCell(-1,0,1), *childCell(0,0,1), *childCell(1,0,1), 1 );
				FluxIn += Flux( *childCell(-2,1,1), *childCell(-1,1,1), *childCell(0,1,1), *childCell(1,1,1), 1 );
        		}

			// Average flux
			FluxIn *= 1./(1<<(Dimension-1));

      		}
		else
			FluxIn  = Flux( *cousinCell(-2,0,0), *cousinCell(-1,0,0), ThisCell, *cousinCell(1,0,0), 1 );


		// If the cell is a leaf with virtual children and its right cousin is a node, compute flux on upper level

		if (isLeafWithVirtualChildren() && node(Nl, Ni+1, Nj, Nk) != 0 && node(Nl, Ni+1, Nj, Nk)->isInternalNode() && FluxCorrection)
		{
			FluxOut  = Flux( *childCell(0,0,0), *childCell(1,0,0), *childCell(2,0,0), *childCell(3,0,0), 1 );

			if (Dimension > 1)
				FluxOut += Flux( *childCell(0,1,0), *childCell(1,1,0), *childCell(2,1,0), *childCell(3,1,0), 1 );

			if (Dimension > 2)
			{
				FluxOut += Flux( *childCell(0,0,1), *childCell(1,0,1), *childCell(2,0,1), *childCell(3,0,1), 1 );
				FluxOut += Flux( *childCell(0,1,1), *childCell(1,1,1), *childCell(2,1,1), *childCell(3,1,1), 1 );
        		}

			// Average flux
			FluxOut *= 1./(1<<(Dimension-1));

     		}
		else
			FluxOut = Flux( *cousinCell(-1,0,0), ThisCell, *cousinCell(1,0,0), *cousinCell(2,0,0), 1 );

		// Add divergence in x-direction

		// Test : use of spherical coordinates (only in 1D)

		if (Dimension == 1 && Coordinate == 2)
      		{
			// Spherical case
			XIn = ThisCell.center(1) - .5*ThisCell.size(1);
			XOut = XIn + ThisCell.size(1);

			ThisCell.setDivergence( ThisCell.divergence() +
                  	3.*(FluxIn*XIn*XIn - FluxOut*XOut*XOut)/(XOut*XOut*XOut - XIn*XIn*XIn) );
		}
		else
			// Cartesian case
			ThisCell.setDivergence( ThisCell.divergence() + (FluxIn - FluxOut)/(ThisCell.size(1)) );


		// --- Add flux in y-direction ---------------------------------------------------------

		if (Dimension > 1)
		{
			// If the cell is a leaf with virtual children and its front cousin is a node, compute flux on upper level

			if (isLeafWithVirtualChildren() && node(Nl, Ni, Nj-1, Nk) != 0 && node(Nl, Ni, Nj-1, Nk)->isInternalNode() && FluxCorrection)
			{
				FluxIn  = Flux( *childCell(0,-2,0), *childCell(0,-1,0), *childCell(0,0,0), *childCell(0,1,0), 2 );
				FluxIn += Flux( *childCell(1,-2,0), *childCell(1,-1,0), *childCell(1,0,0) , *childCell(1,1,0), 2 );

				if (Dimension > 2)
				{
					FluxIn += Flux( *childCell(0,-2,1), *childCell(0,-1,1), *childCell(0,0,1), *childCell(0,1,1), 2 );
					FluxIn += Flux( *childCell(1,-2,1), *childCell(1,-1,1), *childCell(1,0,1), *childCell(1,1,1), 2 );
        			}

				// Average flux
				FluxIn *= 1./(1<<(Dimension-1));

      			}
			else
				FluxIn = Flux( *cousinCell(0,-2,0), *cousinCell(0,-1,0), ThisCell, *cousinCell(0,1,0), 2 );


			// If the cell is a leaf with virtual children and its back cousin is a node, compute flux on upper level

			if (isLeafWithVirtualChildren() && node(Nl, Ni, Nj+1, Nk) != 0 && node(Nl, Ni, Nj+1, Nk)->isInternalNode() && FluxCorrection)
			{
				FluxOut  = Flux( *childCell(0,0,0), *childCell(0,1,0), *childCell(0,2,0), *childCell(0,3,0), 2 );
				FluxOut += Flux( *childCell(1,0,0), *childCell(1,1,0), *childCell(1,2,0), *childCell(1,3,0), 2 );

				if (Dimension > 2)
				{
					FluxOut += Flux( *childCell(0,0,1), *childCell(0,1,1), *childCell(0,2,1), *childCell(0,3,1), 2 );
					FluxOut += Flux( *childCell(1,0,1), *childCell(1,1,1), *childCell(1,2,1), *childCell(1,3,1), 2 );
        			}

				// Average flux
				FluxOut *= 1./(1<<(Dimension-1));

      			}
			else
				FluxOut = Flux( *cousinCell(0,-1,0), ThisCell, *cousinCell(0,1,0), *cousinCell(0,2,0), 2 );

			// Add divergence in y-direction

			ThisCell.setDivergence( ThisCell.divergence() + (FluxIn - FluxOut)/(ThisCell.size(2)) );
		}


		// --- Add flux in z-direction ---------------------------------------------------------

		if (Dimension > 2)
		{
			// If the cell is a leaf with virtual children and its lower cousin is a node, compute flux on upper level

			if (isLeafWithVirtualChildren() && node(Nl, Ni, Nj, Nk-1) != 0 && node(Nl, Ni, Nj, Nk-1)->isInternalNode() && FluxCorrection)
			{
				FluxIn  = Flux( *childCell(0,0,-2), *childCell(0,0,-1), *childCell(0,0,0), *childCell(0,0,1), 3 );
				FluxIn += Flux( *childCell(1,0,-2), *childCell(1,0,-1), *childCell(1,0,0), *childCell(1,0,1), 3 );
				FluxIn += Flux( *childCell(0,1,-2), *childCell(0,1,-1), *childCell(0,1,0), *childCell(0,1,1), 3 );
				FluxIn += Flux( *childCell(1,1,-2), *childCell(1,1,-1), *childCell(1,1,0), *childCell(1,1,1), 3 );

				// Average flux
				FluxIn *= 0.25;

      			}
			else
				FluxIn  = Flux( *cousinCell(0,0,-2), *cousinCell(0,0,-1), ThisCell, *cousinCell(0,0,1), 3 );

			// If the cell is a leaf with virtual children and its upper cousin is a node, compute flux on upper level

			if (isLeafWithVirtualChildren() && node(Nl, Ni, Nj, Nk+1) != 0 && node(Nl, Ni, Nj, Nk+1)->isInternalNode() && FluxCorrection)
			{
				FluxOut  = Flux( *childCell(0,0,0), *childCell(0,0,1), *childCell(0,0,2), *childCell(0,0,3), 3 );
				FluxOut += Flux( *childCell(1,0,0), *childCell(1,0,1), *childCell(1,0,2), *childCell(1,0,3), 3 );
				FluxOut += Flux( *childCell(0,1,0), *childCell(0,1,1), *childCell(0,1,2), *childCell(0,1,3), 3 );
				FluxOut += Flux( *childCell(1,1,0), *childCell(1,1,1), *childCell(1,1,2), *childCell(1,1,3), 3 );

				// Average flux
				FluxOut *= 0.25;

      			}
			else
				FluxOut = Flux( *cousinCell(0,0,-1), ThisCell, *cousinCell(0,0,1), *cousinCell(0,0,2), 3 );

			// Add divergence in z-direction

			ThisCell.setDivergence( ThisCell.divergence() + (FluxIn - FluxOut)/(ThisCell.size(3)) );
		}
  	}

	// --- Recurse on children ---

	if (isInternalNode())
		for (n = 0; n < ChildNb; n++)
			Child[n]->computeDivergence();
}
/*
______________________________________________________________________________________________

	Compute Velocity Gradient (Navier-Stokes only)
______________________________________________________________________________________________

*/
void Node::computeGradient()
{
	// --- Local variables ---

	int n=0;														// Counter on children
	real rho1=0., rho2=0.;          // Densities
	real rhoE1=0., rhoE2=0.;	// Energies
	real V1=0., V2=0.;  		// Velocity
	real dx=0.;			// Cell size
	real dxV=0.;			// Correction of dx for the computation of GradV close to solid walls
	int p=0, q=0;			// Counters
	int ei=0, ej=0, ek=0;		// 1 if this direction is chosen, 0 elsewhere

	// --- Computation ---

//	if (requiresDivergenceComputation() || (CVS && isParentOfLeaf()))

	if (requiresDivergenceComputation())
	{
		if (EquationType != 6)
		{
			cout << "Node.cpp: In method `void Node::computeGradient()':\n";
			cout << "Node.cpp: EquationType not equal to 6 \n";
			cout << "carmen: *** [Node.o] Execution error\n";
			cout << "carmen: abort execution.\n";
			exit(1);
		}

		for (p=1;p <= Dimension; p++)
		{
			ei = (p==1)? 1:0;
			ej = (p==2)? 1:0;
			ek = (p==3)? 1:0;

			dx = ThisCell.size(p);
			dx *= 2.;

			// dxV = correction on dx for the computation of GradV close to solid walls

			if (BoundaryRegion(cell(Nl, Ni+ei,Nj+ej,Nk+ek)->center()) > 3 ||
			    BoundaryRegion(cell(Nl, Ni-ei,Nj-ej,Nk-ek)->center()) > 3 )
			    	dxV = 0.75*dx;
			else
				dxV = dx;

			rho1 = cell(Nl, Ni+ei,Nj+ej,Nk+ek)->density();
			rho2 = cell(Nl, Ni-ei,Nj-ej,Nk-ek)->density();

			ThisCell.setGradient(p, 1, (rho1-rho2)/dx);

			for (q=1; q <= Dimension; q++)
			{
				V1=cell(Nl, Ni+ei, Nj+ej, Nk+ek)->velocity(q);
				V2=cell(Nl, Ni-ei, Nj-ej, Nk-ek)->velocity(q);
				ThisCell.setGradient(p, q+1, (V1-V2)/dxV);
			}

			rhoE1 = cell(Nl, Ni+ei, Nj+ej, Nk+ek)->energy();
			rhoE2 = cell(Nl, Ni-ei, Nj-ej, Nk-ek)->energy();

			ThisCell.setGradient(p, Dimension+2, (rhoE1-rhoE2)/dx);
		}
  }

	// --- Recurse on children ---

	if (isInternalNode())
		for (n = 0; n < ChildNb; n++)
			Child[n]->computeGradient();
}
/*
______________________________________________________________________________________________

	Runge-Kutta step
______________________________________________________________________________________________

*/
void Node::RungeKutta()
{
	// --- Local variables ---

	int n=0;			// Counter on children
	real c1=0., c2=0., c3=0.;	// Runge-Kutta coefficients
	real LocalTimeStep=TimeStep;	// Local time step

	Vector Q(QuantityNb), Qs(QuantityNb), D(QuantityNb);	// Cell-average, temporary cell-average and divergence

	// --- If node is in the tree, recurse on children -----------------------------------------

	if (isInternalNode())
		for (n = 0; n < ChildNb; n++)
			Child[n]->RungeKutta();

	// --- If the leaf is in the boundary, do not perform time evolution -----------------------

	if (UseBoundaryRegions)
	{
		if (isLeaf() && isInsideBoundary())
			return;
	}

	// --- For leaves in the fluid region ------------------------------------------------------

	if (requiresDivergenceComputation())
	{
		// --- Compute local time step in function of the scale ---

		if (TimeAdaptivity)
		{
			// Compute local time step
			LocalTimeStep = TimeStep*(1<<(TimeAdaptivityFactor*(ScaleNb-Nl)));

			// At the end of the time cycle, Q <- Qlow (solution at end of cycle)
			if (isEndTimeCycle() && StepNo == 1 && Nl < ScaleNb)
				ThisCell.setAverage( ThisCell.lowAverage());
		}
		// --- Define Runge-Kutta coefficients ---

		switch(StepNo)
		{
			case 1:
				c1 = 1.; c2 = 0.; c3 = 1.;
				break;
			case 2:
				if (StepNb == 2) {c1 = .5;  c2 = .5;  c3 = .5;  }
				if (StepNb == 3) {c1 = .75; c2 = .25; c3 = .25; }
				break;
			case 3:
				c1 = 1./3.; c2 = 2.*c1; c3 = c2;
				break;
		};

		// --- Runge-Kutta step ---

    		Q  = ThisCell.average();
		Qs = ThisCell.tempAverage();
		D  = ThisCell.divergence();

		// Perform RK step only in fluid region
		ThisCell.setAverage( c1*Qs + c2*Q + (c3 * LocalTimeStep)*D );

		// For the Runge-Kutta-Fehlberg 2(3) method, store second-stage with the RK2 coefficients

		if (!ConstantTimeStep && (StepNo == 2) && (StepNb == 3))   //  MODIFIED 10.12.05
			ThisCell.setLowAverage(0.5*(Qs + Q + LocalTimeStep*D));

		// Correction on concentration: if Y < 0, then Y = 0
    		if ((EquationType >= 3)&&(EquationType <= 5))
			if (ThisCell.average(2) < 0) ThisCell.setAverage(2, 0.);

		// Correction on scalar if Y > 1 => Y = 1, if Y < 0 => Y = 0

		if (EquationType == 6 && ScalarEqNb == 1)
		{
			if (ThisCell.average(Dimension+3) < 0)
				ThisCell.setAverage(Dimension+3, 0);

			if (ThisCell.average(Dimension+3) > ThisCell.average(1))
				ThisCell.setAverage(Dimension+3, ThisCell.average(1));
		}

		// Time adaptivity :
		if (TimeAdaptivity)
		{
			if (isBeginTimeCycle() && StepNo == 1 && Nl < ScaleNb)
			storeTimeEvolution();
		}
	}

}
/*
______________________________________________________________________________________________

	Check stability
______________________________________________________________________________________________

*/
void Node::checkStability()
{
	// --- Local variables ---

	int n, iaux;		// Counter on children
	real x=0., y=0., z=0.;		// Real position

	// --- Recursion ---

	if (isInternalNode())
	{
		for (n = 0; n < ChildNb; n++)
			Child[n]->checkStability();
	}
	else
	{
		// --- Compute x, y, z ---

		x = ThisCell.center(1);
		if (Dimension > 1) y = ThisCell.center(2);
		if (Dimension > 2) z = ThisCell.center(3);

		if (ThisCell.isOverflow())
		{
			iaux=system("echo Unstable computation.>> carmen.prf");
			if (Cluster == 0) iaux=system("echo carmen: unstable computation. >> OUTPUT");
			cout << "carmen: instability detected at iteration no. "<< IterationNo <<"\n";
			cout << "carmen: position ("<< x <<", "<<y<<", "<<z<<")\n";
			cout << "carmen: abort execution.\n";
			exit(1);
		}
	}
}
/*
______________________________________________________________________________________________

	Compute integral values
______________________________________________________________________________________________

*/
void Node::computeIntegral()
{
	// --- Local variables ---

  int 	QuantityNo;			// Quantity number (0 to QuantityNb)
	int	n;				// Counter on children
	int 	AxisNo;				// Counter on dimension
	real 	T, Y;				// Temperature and concentration
	real 	dx, dy=0., dz=0.;		// Cell size
	real 	x, y, z, t;			// position, time
	real  	Omega, Radius;			// local reaction rate
	Vector 	Center(Dimension);  		// local center of the flame ball
	real 	VelocityMax;			// local maximum of the velocity
	real	MemoryCompression = 0.;		// Memory compression

	Vector	GradDensity(Dimension);		// gradient of density
	Vector  GradPressure(Dimension);	// gradient of pressure
	real 	Density=0.;			// density

	real 	X1=0., X2=0.;			// Positions of the center for the computation of GradPressure
	real	P1=0., P2=0.;			// Pressures for the computations of GradPressure

	int 	ei=0, ej=0, ek=0;		// 1 if this direction is chosen, 0 elsewhere


	// --- Init ---

	if (Nl == 0)
	{
		// Only if ExpectedCompression not equal to zero => variable tolerance

		if (ExpectedCompression != 0.)
		{
			MemoryCompression = (1.*CellNb)/(1<<(ScaleNb*Dimension));
			Tolerance = Tolerance*(1.- (ExpectedCompression-MemoryCompression));
			if (Tolerance > 1E+10)
			{
				printf("carmen: ExpectedCompression unreachable\n");
				printf("carmen: maximal compression is %5.2f %%", MemoryCompression*100.);
				printf("carmen: abort execution.\n");
				exit(1);
			}
		}

		// Init integral values

		FlameVelocity	= 0.;
		GlobalMomentum	= 0.;
		GlobalEnergy	= 0.;
		GlobalEnstrophy = 0.;
		ExactMomentum	= 0.;
		ExactEnergy	= 0.;

		GlobalReactionRate 	= 0.;
		AverageRadius		= 0.;
		ReactionRateMax		= 0.;

		for (AxisNo=1; AxisNo <= Dimension; AxisNo++)
			Center.setValue(AxisNo,XCenter[AxisNo]);

		ErrorMax 			= 0.;
		ErrorMid			= 0.;
		ErrorL2				= 0.;
		ErrorNb				= 0;

		RKFError = 0.;

		EigenvalueMax = 0.;
		QuantityMax.setZero();
		QuantityAverage.setZero();

		IntVorticity=0.;
		IntDensity=0.;
		IntMomentum.setZero();
		BaroclinicEffect=0.;

	}

	// --- Recursion ---

	if (isInternalNode())
	{
		for (n = 0; n < ChildNb; n++)
			Child[n]->computeIntegral();
	}
	else if (isLeaf())
	{
		// Whatever the equation, if ConstantTimeStep is false, compute RKFError

		if (!ConstantTimeStep && StepNb == 3)
		{
			for (QuantityNo = 1; QuantityNo <= QuantityNb; QuantityNo++)
			{
				if (Abs(ThisCell.average(QuantityNo)) > RKFAccuracyFactor)
					RKFError = Max(RKFError, Abs(1.-ThisCell.lowAverage(QuantityNo)/ThisCell.average(QuantityNo)));
			}
		}

			switch (EquationType)
			{
				// ADVECTION, BURGERS
				case 1:
				case 2:
					if (Dimension == 1 && (EquationType == 1 || EquationType == 2))
					{
						x  = ThisCell.center(1);
						dx = ThisCell.size(1);
						t  = IterationNo*TimeStep;

						ErrorNb++;
						ErrorGlobalNb++;

						ErrorMax = Max( ErrorMax, fabs(ThisCell.average(1)-AnalyticAverage(x, dx, t) ));
						ErrorGlobalMax = Max( ErrorGlobalMax, ErrorMax );

						ErrorMid       = ( (ErrorNb-1)*ErrorMid + fabs( ThisCell.average(1)-AnalyticAverage(x, dx, t) )) / ErrorNb;
						ErrorGlobalMid = ( (ErrorGlobalNb-1)*ErrorGlobalMid + fabs( ThisCell.average(1)-AnalyticAverage(x, dx, t) )) / ErrorGlobalNb;

						ErrorL2       = sqrt(( (ErrorNb-1)*power2(ErrorL2) + power2( ThisCell.average(1)-AnalyticAverage(x, dx, t) )) / ErrorNb);
						ErrorGlobalL2 = sqrt(( (ErrorGlobalNb-1)*power2(ErrorGlobalL2) + power2( ThisCell.average(1)-AnalyticAverage(x, dx, t) )) / ErrorGlobalNb);

            			// Compute global momentum and energy for MR and exact solutions

						GlobalMomentum += fabs(ThisCell.average(1))*dx;
						ExactMomentum  += fabs(AnalyticAverage(x,dx,t))*dx;

						GlobalEnergy += power2(ThisCell.average(1))*dx;
						ExactEnergy  += power2(AnalyticAverage(x,dx,t))*dx;
					}
					else if (Dimension == 2 && (EquationType == 1 || EquationType == 6))
					{
						x  = ThisCell.center(1);
						dx = ThisCell.size(1);
						y  = ThisCell.center(2);
						dy = ThisCell.size(2);
						t  = IterationNo*TimeStep;

						QuantityNo = (EquationType == 1) ? 1:2;

						ErrorNb++;
						ErrorGlobalNb++;

						ErrorMax = Max( ErrorMax, fabs(ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, t) ));
						ErrorGlobalMax = Max( ErrorGlobalMax, ErrorMax );

						ErrorMid       = ( (ErrorNb-1)*ErrorMid + fabs( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, t) )) / ErrorNb;
						ErrorGlobalMid = ( (ErrorGlobalNb-1)*ErrorGlobalMid + fabs( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, t) )) / ErrorGlobalNb;

						ErrorL2       = sqrt(( (ErrorNb-1)*power2(ErrorL2) + power2( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, t) )) / ErrorNb);
						ErrorGlobalL2 = sqrt(( (ErrorGlobalNb-1)*power2(ErrorGlobalL2) + power2( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, t) )) / ErrorGlobalNb);
					}
					else if (Dimension == 3 && (EquationType == 1 || EquationType == 6))
					{
						x  = ThisCell.center(1);
						dx = ThisCell.size(1);
						y  = ThisCell.center(2);
						dy = ThisCell.size(2);
						z  = ThisCell.center(3);
						dz = ThisCell.size(3);
						t  = IterationNo*TimeStep;

						QuantityNo = (EquationType == 1) ? 1:2;

						ErrorNb++;
						ErrorGlobalNb++;

						ErrorMax = Max( ErrorMax, fabs(ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, z, dz, t) ));
						ErrorGlobalMax = Max( ErrorGlobalMax, ErrorMax );

						ErrorMid       = ( (ErrorNb-1)*ErrorMid + fabs( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, z, dz, t) )) / ErrorNb;
						ErrorGlobalMid = ( (ErrorGlobalNb-1)*ErrorGlobalMid + fabs( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, z, dz, t) )) / ErrorGlobalNb;

						ErrorL2       = sqrt(( (ErrorNb-1)*power2(ErrorL2) + power2( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, z, dz, t) )) / ErrorNb);
						ErrorGlobalL2 = sqrt(( (ErrorGlobalNb-1)*power2(ErrorGlobalL2) + power2( ThisCell.average(QuantityNo)-AnalyticAverage(x, dx, y, dy, z, dz, t) )) / ErrorGlobalNb);
					}
					break;

				// FLAME FRONT AND FLAME-VORTEX INTERACTION
				case 3:
				case 5:
					T  = ThisCell.average(1);
					Y  = ThisCell.average(2);
					dx = ThisCell.size(1);
					if (Dimension > 1) dy = ThisCell.size(2);
					if (Dimension > 2) dz = ThisCell.size(3);

					switch (Dimension)
					{
						case 1:
							FlameVelocity += ReactionRate(T,Y)*dx;
							break;

						case 2:
							if (EquationType==3)
								FlameVelocity += ReactionRate(T,Y)*dx*dy/(XMax[2]-XMin[2]);
							else
								FlameVelocity += ReactionRate(T,Y)*dx*dy;
							break;

						case 3:
							if (EquationType==3)
								FlameVelocity += ReactionRate(T,Y)*dx*dy*dz/((XMax[2]-XMin[2])*(XMax[3]-XMin[3]));
							else
								FlameVelocity += ReactionRate(T,Y)*dx*dy*dz;
							break;
					};
					break;

				// FLAME BALL
				case 4:
					T  = ThisCell.average(1);
					Y  = ThisCell.average(2);

					dx = ThisCell.size(1);
					dy = (Dimension > 1) ? ThisCell.size(2) : 1.;
					dz = (Dimension > 2) ? ThisCell.size(3) : 1.;

					// Compute flame ball average radius

          				Omega = ReactionRate(T,Y);
					ReactionRateMax += Omega;
					if (ReactionRateMax !=0)
					{
						Radius = N2( ThisCell.center() - Center );
						AverageRadius  = ((ReactionRateMax-Omega)*AverageRadius+Omega*Radius)/ReactionRateMax;
          }

					// Compute global reaction rate
					GlobalReactionRate += Omega*dx*dy*dz;   // dy=dz=1 if Dimension=1, dz = 1 if Dimension=2

					break;

        	// NAVIER-STOKES
		case 6:
			dx = ThisCell.size(1);
			dy = (Dimension > 1) ? ThisCell.size(2) : 1.;
			dz = (Dimension > 2) ? ThisCell.size(3) : 1.;

			// --- Compute the global momentum, global energy and global enstrophy ---

			GlobalMomentum 	+= ThisCell.average(2)*dx*dy*dz;
			GlobalEnergy 	+= .5*ThisCell.density()*(ThisCell.velocity()*ThisCell.velocity())*dx*dy*dz;

			if (Dimension > 1)
				GlobalEnstrophy += .5*power2(N2(ThisCell.vorticity()))*dx*dy*dz;

			// --- Compute maximum of the conservative quantities ---

			for (QuantityNo=1; QuantityNo <=QuantityNb; QuantityNo++)
			{
				if ( QuantityMax.value(QuantityNo) < fabs(ThisCell.average(QuantityNo)) )
					QuantityMax.setValue(QuantityNo, fabs(ThisCell.average(QuantityNo)) );
          		}

			// --- Compute the maximal eigenvalue ---

          		VelocityMax = 0.;

			for (AxisNo=1; AxisNo <= Dimension; AxisNo ++)
				VelocityMax = Max( VelocityMax, fabs(ThisCell.velocity(AxisNo)));

			VelocityMax += ThisCell.speedOfSound();

			EigenvalueMax = Max (EigenvalueMax, VelocityMax);

			// --- Compute integral of modulus of vorticity ---

			IntVorticity += N2(ThisCell.vorticity())*dx*dy*dz;
			IntDensity += Abs(ThisCell.density())*dx*dy*dz;
			IntEnergy += Abs(ThisCell.energy())*dx*dy*dz;

			for (AxisNo = 1; AxisNo <= Dimension; AxisNo++)
				IntMomentum.setValue(AxisNo, Abs(ThisCell.average(AxisNo+1))*dx*dy*dz);

			// --- Compute integral of modulus of baroclinic torque ---

			Density = ThisCell.density();

			for (AxisNo = 1; AxisNo <= Dimension; AxisNo ++)
			{
				GradDensity.setValue(AxisNo, ThisCell.gradient(AxisNo,1));

				ei = (AxisNo == 1)? 1:0;
				ej = (AxisNo == 2)? 1:0;
				ek = (AxisNo == 3)? 1:0;

				X1 = cousinCell( ei, ej, ek)->center(AxisNo);
				X2 = cousinCell(-ei,-ej,-ek)->center(AxisNo);

				P1 = cousinCell( ei, ej, ek)->pressure();
				P2 = cousinCell(-ei,-ej,-ek)->pressure();

 				GradPressure.setValue(AxisNo, (P1-P2)/(X1-X2) );
               		}

 			BaroclinicEffect += N2(GradDensity^GradPressure)/(Density*Density)*dx*dy*dz;

		break;
		};
	}
}
/*
______________________________________________________________________________________________

	Check if the tree is graded
______________________________________________________________________________________________

*/
void Node::checkGradedTree()
{
	// --- Local variables ---

	int		n;					// Counter on children
	int   i, j, k;		// Counter in directions
	int   ej, ek;			// 1 if this dimension is existing, 0 else.

	// --- Init ---
	ej = (Dimension > 1)? 1:0;
	ek = (Dimension > 2)? 1:0;


	if (Nl == 0)
	{
		cout << "carmen: testing tree structure ...\n";
		for (n = 0; n < ChildNb; n++)
			Child[n]->checkGradedTree();
		cout << "carmen: tree structure OK. \n";
		return;
	}

	// --- Test if neighbours are existing (eventually virtual) ---

	for (i = -1; i <= 1; i+=1)
	for (j = -1*ej; j <= 1*ej; j+=1)
	for (k = -1*ek; k <= 1*ek; k+=1)
	{
		if (cell(Nl, Ni+i, Nj+j, Nk+k)==0)
		{
			cout << "carmen: Tree not graded':\n";
			cout << "carmen: Node (" << Nl << ", " << Ni << ", "<< Nj << ", "<< Nk << ") \n";
			cout << "carmen: has missing neighbour (" << Nl << ", " << Ni+i << ", "<< Nj+j << ", "<< Nk+k << ") \n";
			cout << "carmen: abort execution.\n";
			exit(1);
		}
	}

  // --- Recurse if it is a node ---

	if (isInternalNode())
	{
		for (n = 0; n < ChildNb; n++)
			Child[n]->checkGradedTree();
	}
}

/*
______________________________________________________________________________________________

	This function return "true" if it is parent of a leaf (with or without virtual children)
______________________________________________________________________________________________

*/
bool Node::isParentOfLeaf() const
{
	// --- Local variables ---

	int		n=0;		// Counter on children

	// --- If the node is a virtual leaf, return false ---

	if (isVirtualLeaf()) return false;

	// --- Test if the node is a leaf

	if (isLeaf()) return true;

	// --- Test if at least one child is a leaf

		for (n = 0; n < ChildNb; n++)
				if (Child[n]->isLeaf()) return true;

	// --- It is not a leaf and all the children are not leaves => return false ---

	return false;
}


/*
______________________________________________________________________________________________

	Back up tree and data
______________________________________________________________________________________________

*/
void Node::backup()
{
	int n;					// Counter on children
	FILE* output;				// Output file
  	int QuantityNo;				// Counter on quantities

	// --- Init ---

	if (Nl==0)
        {
		output = fopen("carmen.bak","w");

                // --- Write header ---

                fprintf(output, "Backup at iteration %i, physical time %e\n", IterationNo, ElapsedTime);
        }
	else
		output = fopen("carmen.bak","a");

	// --- If node is not a leaf, recurse to children ---

	if (isInternalNode())
	{
		fprintf(output,"N\n");
		fclose(output);
		for (n = 0; n < ChildNb; n++)
			Child[n]->backup();
	}
	else
	{
		for (QuantityNo=1; QuantityNo <= QuantityNb; QuantityNo++)
		{
			fprintf(output, FORMAT, ThisCell.average(QuantityNo));
			fprintf(output, "\n");
		}
		fclose(output);
	}
}
/*
______________________________________________________________________________________________

	Restore tree and data
______________________________________________________________________________________________

*/
void Node::restore()
{

	int n=0;		// Counter on children
  	int QuantityNo=0;	// Counter on quantities
	char buf[256];		// Text buffer
        char* caux; 
	// --- Init ---

	if (Nl==0)
        {
  	   GlobalFile = fopen("carmen.bak","r");
           // fgets(buf, 256, GlobalFile);
        }

  	caux=fgets(buf,256,GlobalFile);

	// If the first data is not a 'N', it means that the data has been created using FineMesh

	if (buf[0] != 'N' && Nl==0)
	{
		fclose(GlobalFile);
		restoreFineMesh();
		return;
	}

	// If end of file is reached, close file and return
	if (feof(GlobalFile))
	{
		fclose(GlobalFile);
		return;
	}

	// --- Recurse : if node is not a leaf, split it and restore children ---

	if (buf[0]=='N')
	{
		split();
		for (n = 0; n < ChildNb; n++)
			Child[n]->restore();
	}
	else
	{
		ThisCell.setAverage(1,atof(buf));
		for (QuantityNo=2; QuantityNo <= QuantityNb; QuantityNo++)
		{
			caux=fgets(buf,256,GlobalFile);
			ThisCell.setAverage(QuantityNo, atof(buf));
    		}
		return;
	}
}

/*
______________________________________________________________________________________________

	Restore tree and data from data written using FineMesh
______________________________________________________________________________________________

*/
void Node::restoreFineMesh()
{
	// --- Local variables ---

	int i=0,j=0,k=0;		// Counters in the three directions
	int n=0,iaux;			// Global counter
	int QuantityNo=0;		// Counter on quantities
	FILE* 	input;			// Input file
  	real buf;

	// --- Split the whole tree structure ---

	splitAll();

	// --- Get data from carmen.bak in the FineMesh format


	// -- Open file --

  	input = fopen("carmen.bak","r");

	// -- When there is no back-up file, return --

	if (!input) return;

	// -- Loop on fine-grid cells --

	for (n = 0; n < 1<<(ScaleNb*Dimension); n++)
	{
		// -- Compute i, j, k --

		i = n%(1<<ScaleNb);
		if (Dimension > 1) j = (n%(1<<(2*ScaleNb)))/(1<<ScaleNb);
		if (Dimension > 2) k = n/(1<<(2*ScaleNb));

 		for (QuantityNo=1; QuantityNo <= QuantityNb; QuantityNo++)
		{
        		iaux=fscanf(input,BACKUP_FILE_FORMAT, &buf);
        		cell(ScaleNb,i,j,k)->setAverage(QuantityNo, buf);
		}
    	}

	fclose(input);

}

/*
______________________________________________________________________________________________

	Smoothes the solution by interpolating the value between exact and predicted values.
______________________________________________________________________________________________

*/
void Node::smooth()
{
	int n=0; // Child number

	// --- Recurse on children ---

	if (isInternalNode())
	{
		for (n = 0; n < ChildNb; n++)
			Child[n]->smooth();
	}
	else
	//if (Nl == ScaleNb)
	{
		if (parentCell()->average() == ThisCell.average())
			return;
		else
			ThisCell.setAverage(SmoothCoeff*predict()+(1.-SmoothCoeff)*ThisCell.average());
	}

  return;
}
/*
______________________________________________________________________________________________

Returns true if the detail in the current node is smaller than the required tolerance
______________________________________________________________________________________________

*/
inline bool Node::detailIsSmall() const
{
	bool result=true;

    	if (Nl == 0)
		result = false;
	else
		result = (NormMaxQuantities(ThisCell.average()-predict()) < ComputedTolerance(Nl));

	return result;
}

/*
______________________________________________________________________________________________

Returns true if node is inside a boundary
______________________________________________________________________________________________

*/
bool Node::isInsideBoundary() const
{
	return (ThisCell.isInsideBoundary());
}

/*
______________________________________________________________________________________________

Returns true if node is in the fluid
______________________________________________________________________________________________

*/
bool Node::isInFluid() const
{
	return (!isInsideBoundary());
}

/*
______________________________________________________________________________________________

Returns true if node is on a boundary of the domain
______________________________________________________________________________________________

*/
bool Node::isOnBoundary() const
{
	int i=0,j=0,k=0; 		// Counters in each direction
	int ei = 1;
	int ej = (Dimension > 1)? 1:0;	// 0 in 1D, 1 in 2D and 3D
	int ek = (Dimension > 2)? 1:0;	// 0 in 1D and 2D, 1 in 3D

	Cell *C;

	// Only when boundary regions are used
	if (!UseBoundaryRegions)
		return false;

	// If the node is in the fluid and one of its cousin is inside a boundary, return true
	for (i=-ei; i <= ei; i++)
	for (j=-ej; j <= ej; j++)
	for (k=-ek; k <= ek; k++)
	{
		 C = cousinCell(i,j,k);

		if (ThisCell.isInFluid() && C->isInsideBoundary())
			return true;
		if (ThisCell.isInsideBoundary() && C->isInFluid())
			return true;
		if (BoundaryRegion(ThisCell.center()) != BoundaryRegion(C->center()))
			return true;
	}

	return false;

}
