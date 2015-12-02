/***************************************************************************
                          Node.h  -  description
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

#include "PreProcessor.h"

 /***************************************************************************
 * An object Node is an element of a graded tree structure, used for
 * multiresolution computations. Its contains
 * the following informations:
 *
 * * A pointer to the root node <i>*Root</i> ;<br> * * The corresponding cell <i>ThisCell</i> ;<br>
 * * An array of pointers to the children nodes <i>**Child</i>. Each parent has <i>2**Dimension</i>
 *   children nodes ; <br>
 * * The position of the node <i>Nl, Ni, Nj, Nk</i> into the tree structure (Nl = level) ;<br>
 * * A <i>Flag</i> giving the kind of node :
   0 = not a leaf, 1 = leaf, 2 = leaf with virtual children, 3 = virtual leaf.
 *
 * A leaf is a node without children, a virtual leaf is an artificial leaf created
 * only for the flux computations. No time evolution is made on virtual leaves.
 *
 */


class Node
{
/*
______________________________________________________________________________________________

PUBLIC METHODS

	Constructor and distructor
______________________________________________________________________________________________

*/

 /***************************************************************************
 *  Generates a new node at the position (<i>l, i, j, k</i>) in the tree structure.
 *  A new node is always a leaf. The array of pointers to the children is
 *  allocated, together with the informations on the corresponding cell: cell-center
 *  position and cell size.
 *
 */
public:

 Node(const int l=0, const int i=0, const int j=0, const int k=0);

 /***************************************************************************
 * Removes the node from the tree structure. If the node is not a leaf, all
 * the children are also removed.
 *
 */
	~Node();
/*
______________________________________________________________________________________________

	Get procedures
______________________________________________________________________________________________

*/

 /***************************************************************************
 * Returns the number of cells in the tree.
 *
 */
	inline int		cells() const;

 /***************************************************************************
 * Returns the number of leaves in the tree.
 *
 */
	inline int		leaves() const;

/*
______________________________________________________________________________________________

	Multiresolution procedures
______________________________________________________________________________________________

*/

 /***************************************************************************
 * Computes the details in the leaves and its parent nodes and, in function
 * of the threshold <i>Tolerance</i>, adapt the tree structure.
 *
 */
  	int		adapt();

 /***************************************************************************
 * Checks if the tree is graded. If not, an error is emitted. Only for
 * debugging.
 *
 */
  	void checkGradedTree();


	void		initValue();
	void		addLevel();

 /***************************************************************************
 * Computes the cell-average values of all nodes that are not leaves by projection from
 * the cell-averages values of the leaves. This procedure is required after a time evolution
 * to refresh the internal nodes of the tree.
 *
 */
	Cell* 	project();

 /***************************************************************************
 * Fills the cell-average values of every virtual leaf with values predicted from its parent and
 * uncles. This procedure is required after a time evolution
 * to refresh the virtual leaves of the tree.
 *
 */
	void fillVirtualChildren();
/*
______________________________________________________________________________________________

	Time evolution procedures
______________________________________________________________________________________________

*/

/***************************************************************************
 * Stores cell-average values into temporary cell-average values.
 *
 */
	void store();
               					
/***************************************************************************
 * Stores gradient values into temporary gradient values.
 *
 */
	void storeGrad();
               					
/***************************************************************************
 * Computes the divergence vector with the space discretization scheme.
 *
 */
	void computeDivergence();               

/***************************************************************************
 * Computes one Runge-Kutta step.
 *
 */
	void RungeKutta();           	  				

/***************************************************************************
 * Computes integral values like e.g. flame velocity, global error, etc.
 *
 */
	void computeIntegral();

/***************************************************************************
 * Computes velocity gradient (only for Navier-Stokes).
 *
 */
	void computeGradient();

/***************************************************************************
 * Checks if the computation is numerically unstable, i.e. if one of the
 * cell-averages is overflow. In case of numerical instability, the computation is
 * stopped and a message appears.
 *
 */
	void checkStability();
/*
______________________________________________________________________________________________

	Output procedures
______________________________________________________________________________________________

*/
/***************************************************************************
 * Writes tree structure into file <i>FileName</i>. Only for debugging.
 *
 */
	void writeTree(const char* FileName) const;

/***************************************************************************
 * Writes cell-average values in multiresolution representation and the
 * corresponding mesh into file <i>FileName</i>.
 *
 */
	void writeAverage(const char* FileName);   	

/***************************************************************************
 * Writes mesh data for Gnuplot into file <i>FileName</i>.
 *
 */
	void writeMesh(const char* FileName) const; 

/***************************************************************************
 * Writes header for Data Explorer into file <i>FileName</i>.
 *
 */
	void writeHeader(const char* FileName) const;

/***************************************************************************
 * Writes cell-average values on a regular grid of level <i>L</i> into file <i>FileName</i>.
 *
 */
  	void writeFineGrid(const char* FileName, const int L=ScaleNb) const;

/*
______________________________________________________________________________________________

	Backup-restore procedure (to restart a computation)
______________________________________________________________________________________________

*/
/***************************************************************************
 * Backs up the tree structure and the cell-averages into a file <i>carmen.bak</i>.
 * In further computations, the data can be recovered using <b>Restore()</b>.
 *
 */
	void backup();

/***************************************************************************
 * Restores the tree structure and the cell-averages from the file <i>carmen.bak</i>.
 * This file was created by the method <b>Backup()</b>.
 *
 */
	void restore();


/***************************************************************************
 * Restores the tree structure and the cell-averages from the file <i>carmen.bak</i>
 * in FineMesh format.
 *
 */
	void restoreFineMesh();

/***************************************************************************
 * Deletes the details in the highest level.
 *
 */
	void	smooth();

/*
______________________________________________________________________________________________

PRIVATE METHODS
______________________________________________________________________________________________

	Multiresolution procedures
______________________________________________________________________________________________

*/

/***************************************************************************
 * Splits node with respect to the graded tree structure.
 * If <i>init</i> is true, the new value is computed from the initial condition. Elsewhere
 * it is computed by prediction.
 *
 */
private:
	void split(const bool init=false);
	
/***************************************************************************
 * Splits the complete tree until the smallest scale.
 */
	void splitAll();

/***************************************************************************
 * Combines node, ie remove it if the graded tree structure can be maintained without it.
 *
 */
	void combine();

/***************************************************************************
 * Returns the cell-average value of the current node predicted from the parent node and its nearest neigbours (uncles)
 *
 */
	Vector 	predict() const;

/***************************************************************************
 * Returns the temporary cell-average value of the current node predicted from the parent node and its nearest neigbours (uncles).
 * Only for time adaptivity.
 *
 */
	Vector 	predictTempAverage() const;

/***************************************************************************
 * Returns the gradient of the current node predicted from the parent node and its nearest neigbours (uncles).
 *
 */
	Matrix 	predictGradient() const;

/***************************************************************************
 * Computes the interpolation in time in the current node. Only for time adaptivity.
 *
 */
	inline void computeTimeInterpolation();


/***************************************************************************
 * Computes the extrapolation in time in the current node. Only for time adaptivity.
 *
 */
	inline void computeTimeExtrapolation();

/***************************************************************************
 * Stores the result of the time evolution into temporary cell-average values. Only for time adaptivity.
 *
 */
	inline void storeTimeEvolution();

/***************************************************************************
 * Sets the current node to an internal node.
 *
 */
	inline void setInternalNode();

/***************************************************************************
 * Sets the current node to a simple leaf.
 *
 */
	inline void setSimpleLeaf();

/***************************************************************************
 * Sets the current node to a leaf with virtual children.
 *
 */
	inline void setLeafWithVirtualChildren();

/***************************************************************************
 * Sets the current node to a virtual leaf.
 *
 */
	inline void setVirtualLeaf();

/***************************************************************************
 * Returns true if the current node is an internal node.
 *
 */
	inline bool isInternalNode() const;

/***************************************************************************
 * Returns true if the current node is a leaf (with or without virtual children).
 *
 */
	inline bool isLeaf() const;

/***************************************************************************
 * Returns true if the current node is a simple leaf (no virtual children).
 *
 */
	inline bool isSimpleLeaf() const;

/***************************************************************************
 * Returns true if the current node is a simple leaf with virtual children.
 *
 */
	inline bool isLeafWithVirtualChildren() const;

/***************************************************************************
 * Returns true if the current node is a virtual.
 *
 */
	inline bool isVirtualLeaf() const;

/***************************************************************************
 * Returns true if the current node has children (real or virtual).
 *
 */
	inline bool hasChildren() const;

/***************************************************************************
 * Returns true if the current node or one of its children is a leaf.
 *
 */
	inline bool isParentOfLeaf() const;


/***************************************************************************
 * Returns true if the current node is inside a boundary region.
 *
 */
	bool isInsideBoundary() const;


/***************************************************************************
 * Returns true if the current node is inside the fluid region.
 *
 */
	bool isInFluid() const;


/***************************************************************************
 * Returns true if the current node is on a limit of boundary region.
 *
 */
	bool isOnBoundary() const;


/***************************************************************************
 * Returns true if the current leaf is at the begining of a time cycle.
 * Only useful when TimeAdaptivity = true
 *
 */
	inline bool isBeginTimeCycle() const;


/***************************************************************************
 * Returns true if the current leaf is at the end of a time cycle.
 * Only useful when TimeAdaptivity = true
 *
 */
	inline bool isEndTimeCycle() const;

/***************************************************************************
 * Returns true if a time evolution is required for this node.
 *
 */
	inline bool requiresTimeEvolution() const;

/***************************************************************************
 * Returns true if an interpolation in time is required for this node.
 *
 */
	inline bool requiresTimeInterpolation() const;


/***************************************************************************
 * Returns true if the divergence computation is required for this node.
 *
 */
	inline bool requiresDivergenceComputation() const;

/***************************************************************************
 * Returns true if the detail in the current node is smaller than the required tolerance.
 *
 */
	inline bool detailIsSmall() const;
/*
______________________________________________________________________________________________

	Get procedures
______________________________________________________________________________________________

*/
/***************************************************************************
 * Returns the pointer to the node located at (<i>l, i, j, k</i>). The
 * boundary conditions are not taken into account. If this node does not
 * exist, 0 is returned.
 *
 */
	Node* node(int l, int i, int j = 0, int k = 0) const;

/***************************************************************************
 * Returns the pointer to the cell located at (<i>l, i, j, k</i>), with respect to boundary conditions.
 *
 */
	Cell* cell(int l, int i, int j = 0, int k = 0) const;


/***************************************************************************
 * Returns the pointer to the parent node
 *
 */
	inline Node* parent() const;


/***************************************************************************
 * Returns the pointer to the parent cell
 *
 */
	inline Cell* parentCell() const;

/***************************************************************************
 * Returns the pointer to the uncle node with relative position (<i>i, j, k</i>)
 *
 */
	inline Node* uncle(const int i, const int j=0, const int k=0) const;


/***************************************************************************
 * Returns the pointer to the uncle cell with relative position (<i>i, j, k</i>)
 *
 */
	inline Cell* uncleCell(const int i, const int j=0, const int k=0) const;


/***************************************************************************
 * Returns the pointer to the cousin node with relative position (<i>i, j, k</i>)
 *
 */
	inline Node* cousin(const int i, const int j=0, const int k=0) const;



/***************************************************************************
 * Returns the pointer to the cousin cell with relative position (<i>i, j, k</i>)
 *
 */
	inline Cell* cousinCell(const int i, const int j=0, const int k=0) const;


/***************************************************************************
 * Returns the pointer to the child cell with relative position (<i>i, j, k</i>)
 *
 */
	inline Node* child(const int i, const int j=0, const int k=0) const;


/***************************************************************************
 * Returns the pointer to the child cell with relative position (<i>i, j, k</i>)
 *
 */
	inline Cell* childCell(const int i, const int j=0, const int k=0) const;

/*
______________________________________________________________________________________________

	Procedures on virtual children
______________________________________________________________________________________________

*/

/***************************************************************************
 * Tests and, if possible, deletes the virtual children of the current node.
 *
 */
	void deleteVirtualChildren();

/***************************************************************************
 * Generates virtual children for the current node. If <i>init</i> is true,
 * the cell-average value in the virtual child is computed from the initial
 * condition. Elsewhere, it is computed by prediction from parent and uncles.
 *
 */
	void makeVirtualChildren(bool init=false);

/*
______________________________________________________________________________________________

PRIVATE VARIABLES
______________________________________________________________________________________________

*/
	static Node *Root;	// Pointer on root node
//	static Node Boundary; 	// Boundary cell
	
	int	Nl;		// Scale number
	int Ni, Nj, Nk;		// Position number in x, y, z directions
	
	Node **Child;		// Array of pointers to children (whatever they exist or not)
  	Cell ThisCell;		// Cell corresponding to this node
			
	byte Flag;		// Flag (0 = internal node, 1 = simple leaf, 2 = leaf with virtual children, 3 = virtual leaf
};

/*
______________________________________________________________________________________________

	INLINE FUNCTION
______________________________________________________________________________________________

	Returns the number of cells
______________________________________________________________________________________________

*/
inline int Node::cells() const
{
	return CellNb;
}
/*
______________________________________________________________________________________________

	Get number of leaves
______________________________________________________________________________________________

*/
inline int Node::leaves() const
{
	return LeafNb;
}

/*
______________________________________________________________________________________________

	Returns pointer to child cell using relative position
______________________________________________________________________________________________

*/

inline Cell* Node::childCell(const int i, const int j, const int k) const
{
	return cell(Nl+1, 2*Ni+i, 2*Nj+j, 2*Nk+k);
}

/*
______________________________________________________________________________________________

	Returns pointer to child node using relative position
______________________________________________________________________________________________

*/

inline Node* Node::child(const int i, const int j, const int k) const
{
	return node(Nl+1, 2*Ni+i, 2*Nj+j, 2*Nk+k);
}
/*
______________________________________________________________________________________________

	This function computes the interpolation in time in the node
______________________________________________________________________________________________

*/

inline void Node::computeTimeInterpolation()
{
	ThisCell.setAverage(ThisCell.tempAverage() + ThisCell.average());
}

/*
______________________________________________________________________________________________

	Returns pointer to uncle cell using relative position
______________________________________________________________________________________________

*/

inline Cell* Node::cousinCell(const int i, const int j, const int k) const
{
	return cell(Nl, Ni+i, Nj+j, Nk+k);
}

/*
______________________________________________________________________________________________

	Returns pointer to uncle cell using relative position
______________________________________________________________________________________________

*/

inline Node* Node::cousin(const int i, const int j, const int k) const
{
	return node(Nl, Ni+i, Nj+j, Nk+k);
}

/*
______________________________________________________________________________________________

	Return true if node is has children (real or virtual)
______________________________________________________________________________________________

*/
inline bool Node::hasChildren() const
{
	return (Flag==0 || Flag==2);
}

/*
______________________________________________________________________________________________

	Sets the current node to internal node
______________________________________________________________________________________________

*/
inline void Node::setInternalNode()
{
	Flag = 0;
	return;
}

/*
______________________________________________________________________________________________

	Sets the current node to simple leaf
______________________________________________________________________________________________

*/
inline void Node::setSimpleLeaf()
{
	Flag = 1;
	return;
}
/*
______________________________________________________________________________________________

	Sets the current node to leaf with virtual children
______________________________________________________________________________________________

*/
inline void Node::setLeafWithVirtualChildren()
{
	Flag = 2;
	return;
}
/*
______________________________________________________________________________________________

	Sets the current node to virtual leaf
______________________________________________________________________________________________

*/
inline void Node::setVirtualLeaf()
{
	Flag = 3;
	return;
}
/*
______________________________________________________________________________________________

Returns true if the leaf is at the begining of a time cycle.
Only useful when TimeAdaptivity = true
______________________________________________________________________________________________

*/
inline bool Node::isBeginTimeCycle() const
{
	return ((IterationNo-1)%(1<<(TimeAdaptivityFactor*(ScaleNb-Nl))) == 0);
}

/*
______________________________________________________________________________________________

Returns true if the leaf is at the end of a time cycle.
Only useful when TimeAdaptivity = true
______________________________________________________________________________________________

*/
inline bool Node::isEndTimeCycle() const
{
	return (IterationNo%(1<<(TimeAdaptivityFactor*(ScaleNb-Nl))) == 0);
}

/*
______________________________________________________________________________________________

	Returns pointer to parent node
______________________________________________________________________________________________

*/
inline Node* Node::parent() const
{
	return node(Nl-1, (Ni+4)/2-2, (Nj+4)/2-2, (Nk+4)/2-2);
}

/*
______________________________________________________________________________________________

	Returns pointer to parent cell
______________________________________________________________________________________________

*/
inline Cell* Node::parentCell() const
{
	return cell(Nl-1, (Ni+4)/2-2, (Nj+4)/2-2, (Nk+4)/2-2);
}

/*
______________________________________________________________________________________________

	Returns pointer to uncle node using relative position
______________________________________________________________________________________________

*/
inline Node* Node::uncle(const int i, const int j, const int k) const
{
	return node(Nl-1, (Ni+4)/2-2 + i, (Nj+4)/2-2 + j, (Nk+4)/2-2 + k);
}

/*
______________________________________________________________________________________________

	Returns pointer to uncle cell using relative position
______________________________________________________________________________________________

*/
inline Cell* Node::uncleCell(const int i, const int j, const int k) const
{
	return cell(Nl-1, (Ni+4)/2-2 + i, (Nj+4)/2-2 + j, (Nk+4)/2-2 + k);
}

/*
______________________________________________________________________________________________

	Returns true if node is an internal node (not a leaf)
______________________________________________________________________________________________

*/
inline bool Node::isInternalNode()  const
{
	return (Flag==0);
}

/*
______________________________________________________________________________________________

	Returns true if node is a leaf (with or without virtual children)
______________________________________________________________________________________________

*/
inline bool Node::isLeaf() const
{
	return (Flag==1 || Flag==2);
}

/*
______________________________________________________________________________________________

	Returns true if node is a simple leaf (without virtual children)
______________________________________________________________________________________________

*/
inline bool Node::isSimpleLeaf() const
{
	return (Flag==1);
}

/*
______________________________________________________________________________________________

	Returns true if node is a leaf with virtual children
______________________________________________________________________________________________

*/
inline bool Node::isLeafWithVirtualChildren() const
{
	return (Flag==2);
}

/*
______________________________________________________________________________________________

	Return true if node is a virtual leaf
______________________________________________________________________________________________

*/
inline bool Node::isVirtualLeaf() const
{
	return (Flag==3);
}

/*
______________________________________________________________________________________________

	This function stores the result of the time evolution in Qlow, and interpolates
	the internediary state between [n+1] and [n], which is stored in Q.
______________________________________________________________________________________________

*/

inline void Node::storeTimeEvolution()
{

	// Qlow <- Q
	ThisCell.setLowAverage( ThisCell.average());

	// Q <- 1/2 (Q + Qs)
	ThisCell.setAverage( 0.5*(ThisCell.average() + ThisCell.tempAverage()) );	
}

/*
______________________________________________________________________________________________

	Returns the extrapolation in time. Only for local time stepping
______________________________________________________________________________________________

*/

inline void Node::computeTimeExtrapolation()
{
	Vector A(QuantityNb);

	A = ThisCell.average();
	ThisCell.setAverage(ThisCell.average() + ThisCell.average() - ThisCell.tempAverage());
	ThisCell.setTempAverage(A);
}

/*
______________________________________________________________________________________________

	This function return "true" if this node requires an interpolation in time, instead of a full time evolution
______________________________________________________________________________________________

*/

inline bool Node::requiresDivergenceComputation() const
{
	if (UseBoundaryRegions)
		if (isInsideBoundary()) return false;

	if (TimeAdaptivity)
		// With time adaptivity, compute divergence at the begining and the end of time cycles
		return (isLeaf() && (isBeginTimeCycle()|| isEndTimeCycle()) );
	else
    	// Without time adaptivity, perform TimeEvolution on leaves whatever IterationNo
		return (isLeaf()) ;
}

/*
______________________________________________________________________________________________

	This function return "true" if this node requires a time evolution procedure (no interpolation)
______________________________________________________________________________________________

*/
inline bool Node::requiresTimeEvolution() const
{
	if (UseBoundaryRegions)
		if (isInsideBoundary()) return false;

	if (TimeAdaptivity)
		// With time adaptivity, perform TimeEvolution on leaves every 2^a(L-l) iterations
		return (isLeaf() && isBeginTimeCycle());
	else
    	// Without time adaptivity, perform TimeEvolution on leaves whatever IterationNo
		return (isLeaf()) ;
}

/*
______________________________________________________________________________________________

	This function return "true" if this node requires an interpolation in time, instead of a full time evolution
______________________________________________________________________________________________

*/
inline bool Node::requiresTimeInterpolation() const
{
	if (UseBoundaryRegions)
		if (isInsideBoundary()) return false;

	// When the time adaptivity is not used, never perform interpolation in time
	if (!TimeAdaptivity) return false;

	// With time adaptivity, permform interpolation on leaves where no time evolution is made	
	return (isLeaf() && isEndTimeCycle());
}





