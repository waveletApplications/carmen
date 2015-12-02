/***************************************************************************
                          Parallel.cpp  -  description
                             -------------------
    begin                : Wed Dec 22 2004
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

void CreateMPITopology() {
#if defined PARMPI
 int src;
 int periods[]={1,1,1};
 CartDims[0]=CartDims[1]=CartDims[2]=0;

 MPI_Dims_create(size,Dimension,CartDims);
 MPI_Cart_create(MPI_COMM_WORLD,Dimension,CartDims,periods,1,&comm_cart);
 MPI_Comm_rank(comm_cart, &rank);
 MPI_Cart_coords(comm_cart,rank,Dimension,coords);

 MPI_Cart_shift(comm_cart, 0, -1, &src, &rank_il);
 MPI_Cart_shift(comm_cart, 0, 1, &src, &rank_iu);

 if (Dimension >= 2) {
   MPI_Cart_shift(comm_cart, 1, -1, &src, &rank_jl);
   MPI_Cart_shift(comm_cart, 1, 1, &src, &rank_ju);
 }

 if (Dimension == 3) {
  MPI_Cart_shift(comm_cart, 2, -1, &src, &rank_kl);
  MPI_Cart_shift(comm_cart, 2, 1, &src, &rank_ku);
 }
#endif
}


void FillCellAddr(Cell *Mesh4MPI, int d, int &n) {
#if defined PARMPI

  if ((WhatSend & SendQ) != 0) {
    MPI_Address(Mesh4MPI[d].Q.U, &disp[n]);
   	blocklen[n++]=Mesh4MPI[d].Q.dimension();
  }

  if ((WhatSend & SendQs) != 0) {
  	MPI_Address(Mesh4MPI[d].Qs.U, &disp[n]);
  	blocklen[n++]=Mesh4MPI[d].Qs.dimension();
  }
   
  if ((WhatSend & SendX) != 0) {
    MPI_Address(Mesh4MPI[d].X.U, &disp[n]);
    blocklen[n++]=Mesh4MPI[d].X.dimension();
  }
  
  if ((WhatSend & SenddX) != 0) {
  	MPI_Address(Mesh4MPI[d].dX.U, &disp[n]);
  	blocklen[n++]=Mesh4MPI[d].dX.dimension();
  }
  
  if ((WhatSend & SendD) != 0) {
  	MPI_Address(Mesh4MPI[d].D.U, &disp[n]);
  	blocklen[n++]=Mesh4MPI[d].D.dimension();
  }
   
  if ((WhatSend & SendGrad) != 0) {
   	MPI_Address(Mesh4MPI[d].Grad.U, &disp[n]);
    blocklen[n++]=Mesh4MPI[d].Grad.columns()*Mesh4MPI[d].Grad.lines();
  }  
#endif
}


void FillNbAddr(Cell ***Nb,int l,int i,int j,int &n) {
#if defined PARMPI
  if ((WhatSend & SendQ) != 0) {
    MPI_Address(Nb[l][i][j].Q.U, &disp[n]);
    blocklen[n++]=Nb[l][i][j].Q.dimension();
  }

  if ((WhatSend & SendQs) != 0) {
   MPI_Address(Nb[l][i][j].Qs.U, &disp[n]);
    blocklen[n++]=Nb[l][i][j].Qs.dimension();
  }

  if ((WhatSend & SendX) != 0) {
    MPI_Address(Nb[l][i][j].X.U, &disp[n]);
    blocklen[n++]=Nb[l][i][j].X.dimension();
  }

  if ((WhatSend & SenddX) != 0) {
    MPI_Address(Nb[l][i][j].dX.U, &disp[n]);
    blocklen[n++]=Nb[l][i][j].dX.dimension();
  }

  if ((WhatSend & SendD) != 0) {
	  MPI_Address(Nb[l][i][j].D.U, &disp[n]);
    blocklen[n++]=Nb[l][i][j].D.dimension();
  }


  if ((WhatSend & SendGrad) != 0) {
    MPI_Address(Nb[l][i][j].Grad.U, &disp[n]);
    blocklen[n++]=Nb[l][i][j].Grad.columns()*Nb[l][i][j].Grad.lines();
  }
#endif
}


void CreateMPIType(FineMesh *Root) {
#if defined PARMPI
  int i,j,k;
	int n,d,l;

	Cell *MeshCell;
	MeshCell=Root->MeshCell;

	n=0;
  for (l=0;l<NeighbourNb;l++)
		for (j=0;j<one_D;j++)
			for (k=0;k<two_D;k++)	FillNbAddr(Root->Neighbour_iL,l,j,k,n);
  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeRiL);
	MPI_Type_commit(&MPItypeRiL);

	n=0;
  for (l=0;l<NeighbourNb;l++)
		for (j=0;j<one_D;j++)
			for (k=0;k<two_D;k++)	FillNbAddr(Root->Neighbour_iU,l,j,k,n);
  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeRiU);
	MPI_Type_commit(&MPItypeRiU);

	n=0;
	for (l=0;l<NeighbourNb;l++)
		for (j=0;j<one_D;j++)
			for (k=0;k<two_D;k++)	{
				i=l;
				d=i + (1<<ScaleNb)*(j + (1<<ScaleNb)*k);
				FillCellAddr(MeshCell,d,n);
			}
  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeSiL);
	MPI_Type_commit(&MPItypeSiL);

	n=0;
	for (l=0;l<NeighbourNb;l++)
		for (j=0;j<one_D;j++)
			for (k=0;k<two_D;k++)	{
				i=(1<<ScaleNb)-NeighbourNb+l;
				d=i + (1<<ScaleNb)*(j + (1<<ScaleNb)*k);
				FillCellAddr(MeshCell,d,n);
			}
  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeSiU);
	MPI_Type_commit(&MPItypeSiU);

	if (Dimension >= 2) {
		n=0;
	  for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (k=0;k<two_D;k++)	FillNbAddr(Root->Neighbour_jL,l,i,k,n);
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeRjL);
		MPI_Type_commit(&MPItypeRjL);

		n=0;
	  for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (k=0;k<two_D;k++)	FillNbAddr(Root->Neighbour_jU,l,i,k,n);
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeRjU);
		MPI_Type_commit(&MPItypeRjU);

		n=0;
		for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (k=0;k<two_D;k++)	{
					j=l;
					d=i + (1<<ScaleNb)*(j + (1<<ScaleNb)*k);
					FillCellAddr(MeshCell,d,n);
				}
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeSjL);
		MPI_Type_commit(&MPItypeSjL);

		n=0;
		for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (k=0;k<two_D;k++)	{
					j=(1<<ScaleNb)-NeighbourNb+l;
					d=i + (1<<ScaleNb)*(j + (1<<ScaleNb)*k);
					FillCellAddr(MeshCell,d,n);
				}
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeSjU);
		MPI_Type_commit(&MPItypeSjU);
  }


	if (Dimension == 3) {
		n=0;
	  for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (j=0;j<two_D;j++)	FillNbAddr(Root->Neighbour_kL,l,i,j,n);
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeRkL);
		MPI_Type_commit(&MPItypeRkL);

		n=0;
	  for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (j=0;j<two_D;j++)	FillNbAddr(Root->Neighbour_kU,l,i,j,n);
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeRkU);
		MPI_Type_commit(&MPItypeRkU);

		n=0;
		for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (j=0;j<two_D;j++)	{
					k=l;
					d=i + (1<<ScaleNb)*(j + (1<<ScaleNb)*k);
					FillCellAddr(MeshCell,d,n);
				}
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeSkL);
		MPI_Type_commit(&MPItypeSkL);

		n=0;
		for (l=0;l<NeighbourNb;l++)
			for (i=0;i<one_D;i++)
				for (j=0;j<two_D;j++)	{
					k=(1<<ScaleNb)-NeighbourNb+l;
					d=i + (1<<ScaleNb)*(j + (1<<ScaleNb)*k);
					FillCellAddr(MeshCell,d,n);
				}
	  MPI_Type_hindexed(CellElementsNb*NeighbourNb*one_D*two_D,blocklen,disp,MPI_Type,&MPItypeSkU);
		MPI_Type_commit(&MPItypeSkU);
	}

#endif
}



void FreeMPIType() {
#if defined PARMPI
 	MPI_Type_free(&MPItypeSiL);
 	MPI_Type_free(&MPItypeSiU);
 	MPI_Type_free(&MPItypeRiL);
 	MPI_Type_free(&MPItypeRiU);

	if (Dimension >= 2) {
	 	MPI_Type_free(&MPItypeSjL);
	 	MPI_Type_free(&MPItypeSjU);
	 	MPI_Type_free(&MPItypeRjL);
	 	MPI_Type_free(&MPItypeRjU);
	}

	if (Dimension == 3) {
	 	MPI_Type_free(&MPItypeSkL);
	 	MPI_Type_free(&MPItypeSkU);
	 	MPI_Type_free(&MPItypeRkL);
	 	MPI_Type_free(&MPItypeRkU);
	}
#endif
}


void CreateMPILinks() {
	int exNb;
	exNb=0;
#if defined PARMPI

//Send

	switch (MPISendType) {
	case 0:
	  MPI_Bsend_init(MPI_BOTTOM, 1, MPItypeSiL, rank_il, 100, comm_cart ,&req[exNb++]);
	  MPI_Bsend_init(MPI_BOTTOM, 1, MPItypeSiU, rank_iu, 200, comm_cart, &req[exNb++]);
		break;

	case 10:
	  MPI_Send_init(MPI_BOTTOM, 1, MPItypeSiL, rank_il, 100, comm_cart,&req[exNb++]);
	  MPI_Send_init(MPI_BOTTOM, 1, MPItypeSiU, rank_iu, 200, comm_cart,&req[exNb++]);
		break;

	case 20:
	  MPI_Ssend_init(MPI_BOTTOM, 1, MPItypeSiL, rank_il, 100, comm_cart,&req[exNb++]);
	  MPI_Ssend_init(MPI_BOTTOM, 1, MPItypeSiU, rank_iu, 200, comm_cart,&req[exNb++]);
		break;
	}

	if (Dimension >= 2) {
		switch (MPISendType) {
		case 0:
		  MPI_Bsend_init(MPI_BOTTOM, 1, MPItypeSjL, rank_jl, 300, comm_cart,&req[exNb++]);
	  	MPI_Bsend_init(MPI_BOTTOM, 1, MPItypeSjU, rank_ju, 400, comm_cart,&req[exNb++]);
		  break;

		case 10:
		  MPI_Send_init(MPI_BOTTOM, 1, MPItypeSjL, rank_jl, 300, comm_cart,&req[exNb++]);
		  MPI_Send_init(MPI_BOTTOM, 1, MPItypeSjU, rank_ju, 400, comm_cart,&req[exNb++]);
			break;

		case 20:
		  MPI_Ssend_init(MPI_BOTTOM, 1, MPItypeSjL, rank_jl, 300, comm_cart,&req[exNb++]);
		  MPI_Ssend_init(MPI_BOTTOM, 1, MPItypeSjU, rank_ju, 400, comm_cart,&req[exNb++]);
			break;
		}
	}

	if (Dimension == 3) {
		switch (MPISendType) {
		case 0:
		  MPI_Bsend_init(MPI_BOTTOM, 1, MPItypeSkL, rank_kl, 500, comm_cart,&req[exNb++]);
		  MPI_Bsend_init(MPI_BOTTOM, 1, MPItypeSkU, rank_ku, 600, comm_cart,&req[exNb++]);
			break;

		case 10:
		  MPI_Send_init(MPI_BOTTOM, 1, MPItypeSkL, rank_kl, 500, comm_cart,&req[exNb++]);
		  MPI_Send_init(MPI_BOTTOM, 1, MPItypeSkU, rank_ku, 600, comm_cart,&req[exNb++]);
			break;

		case 20:
  		MPI_Ssend_init(MPI_BOTTOM, 1, MPItypeSkL, rank_kl, 500, comm_cart,&req[exNb++]);
		  MPI_Ssend_init(MPI_BOTTOM, 1, MPItypeSkU, rank_ku, 600, comm_cart,&req[exNb++]);
			break;
		}
	}

//Recv

	MPI_Recv_init(MPI_BOTTOM, 1, MPItypeRiL, rank_il, 200, comm_cart, &req[exNb++]);
	MPI_Recv_init(MPI_BOTTOM, 1, MPItypeRiU, rank_iu, 100, comm_cart, &req[exNb++]);

	if (Dimension >= 2) {
 		MPI_Recv_init(MPI_BOTTOM, 1, MPItypeRjL, rank_jl, 400, comm_cart, &req[exNb++]);
		MPI_Recv_init(MPI_BOTTOM, 1, MPItypeRjU, rank_ju, 300, comm_cart, &req[exNb++]);
	}

	if (Dimension == 3) {
 		MPI_Recv_init(MPI_BOTTOM, 1, MPItypeRkL, rank_kl, 600, comm_cart, &req[exNb++]);
		MPI_Recv_init(MPI_BOTTOM, 1, MPItypeRkU, rank_ku, 500, comm_cart, &req[exNb++]);
	}
#endif
}

void CPUExchange(FineMesh *Root,int WS) {
  CommTimer.start();
#if defined PARMPI
  int i,k;
	int exNb=0;

  WhatSend=WS;
  CellElementsNb=0;

  for (i=0;i<16;i++) {
    k=1<<i;
    if ((WS & k) != 0) CellElementsNb++;
  }
      
  static bool ft=true;
//  if (ft==true) {
    CreateMPIType(Root);
//    CreateMPILinks();
//    ft=false;  
//  }

//  MPI_Startall(4*Dimension,req);

    
//Send
	switch (MPISendType) {
	case 0:
	  MPI_Ibsend(MPI_BOTTOM, 1, MPItypeSiL, rank_il, 100, comm_cart ,&req[exNb++]);
	  MPI_Ibsend(MPI_BOTTOM, 1, MPItypeSiU, rank_iu, 200, comm_cart, &req[exNb++] );
		break;

	case 10:
	  MPI_Isend(MPI_BOTTOM, 1, MPItypeSiL, rank_il, 100, comm_cart,&req[exNb++]);
	  MPI_Isend(MPI_BOTTOM, 1, MPItypeSiU, rank_iu, 200, comm_cart,&req[exNb++]);
		break;

	case 20:
	  MPI_Issend(MPI_BOTTOM, 1, MPItypeSiL, rank_il, 100, comm_cart,&req[exNb++]);
	  MPI_Issend(MPI_BOTTOM, 1, MPItypeSiU, rank_iu, 200, comm_cart,&req[exNb++]);
		break;
	}

	if (Dimension >= 2) {
		switch (MPISendType) {
		case 0:
		  MPI_Ibsend(MPI_BOTTOM, 1, MPItypeSjL, rank_jl, 300, comm_cart,&req[exNb++]);
	  	MPI_Ibsend(MPI_BOTTOM, 1, MPItypeSjU, rank_ju, 400, comm_cart,&req[exNb++]);
		  break;

		case 10:
		  MPI_Isend(MPI_BOTTOM, 1, MPItypeSjL, rank_jl, 300, comm_cart,&req[exNb++]);
		  MPI_Isend(MPI_BOTTOM, 1, MPItypeSjU, rank_ju, 400, comm_cart,&req[exNb++]);
			break;

		case 20:
		  MPI_Issend(MPI_BOTTOM, 1, MPItypeSjL, rank_jl, 300, comm_cart,&req[exNb++]);
		  MPI_Issend(MPI_BOTTOM, 1, MPItypeSjU, rank_ju, 400, comm_cart,&req[exNb++]);
			break;
		}
	}

	if (Dimension == 3) {
		switch (MPISendType) {
		case 0:
		  MPI_Ibsend(MPI_BOTTOM, 1, MPItypeSkL, rank_kl, 500, comm_cart,&req[exNb++]);
		  MPI_Ibsend(MPI_BOTTOM, 1, MPItypeSkU, rank_ku, 600, comm_cart,&req[exNb++]);
			break;

		case 10:
		  MPI_Isend(MPI_BOTTOM, 1, MPItypeSkL, rank_kl, 500, comm_cart,&req[exNb++]);
		  MPI_Isend(MPI_BOTTOM, 1, MPItypeSkU, rank_ku, 600, comm_cart,&req[exNb++]);
			break;

		case 20:
  		MPI_Issend(MPI_BOTTOM, 1, MPItypeSkL, rank_kl, 500, comm_cart,&req[exNb++]);
		  MPI_Issend(MPI_BOTTOM, 1, MPItypeSkU, rank_ku, 600, comm_cart,&req[exNb++]);
			break;
		}
	}

//Recv

	if (MPIRecvType==0) {
		MPI_Recv(MPI_BOTTOM, 1, MPItypeRiL, rank_il, 200, comm_cart, &st[6]);
		MPI_Recv(MPI_BOTTOM, 1, MPItypeRiU, rank_iu, 100, comm_cart, &st[7]);
	} else
	{	
		MPI_Irecv(MPI_BOTTOM, 1, MPItypeRiL, rank_il, 200, comm_cart, &req[exNb++]);
		MPI_Irecv(MPI_BOTTOM, 1, MPItypeRiU, rank_iu, 100, comm_cart, &req[exNb++]);
	}

	if (Dimension >= 2) {
		if (MPIRecvType==0) {
	 		MPI_Recv(MPI_BOTTOM, 1, MPItypeRjL, rank_jl, 400, comm_cart, &st[8]);
			MPI_Recv(MPI_BOTTOM, 1, MPItypeRjU, rank_ju, 300, comm_cart, &st[9]);
		} else
		{
	 		MPI_Irecv(MPI_BOTTOM, 1, MPItypeRjL, rank_jl, 400, comm_cart, &req[exNb++]);
			MPI_Irecv(MPI_BOTTOM, 1, MPItypeRjU, rank_ju, 300, comm_cart, &req[exNb++]);
		}
	}

	if (Dimension == 3) {
		if (MPIRecvType==0) {
	 		MPI_Recv(MPI_BOTTOM, 1, MPItypeRkL, rank_kl, 600, comm_cart, &st[10]);
			MPI_Recv(MPI_BOTTOM, 1, MPItypeRkU, rank_ku, 500, comm_cart, &st[11]);
		} else
		{
	 		MPI_Irecv(MPI_BOTTOM, 1, MPItypeRkL, rank_kl, 600, comm_cart, &req[exNb++]);
			MPI_Irecv(MPI_BOTTOM, 1, MPItypeRkU, rank_ku, 500, comm_cart, &req[exNb++]);
		}
	}     

  FreeMPIType();
#endif
  CommTimer.stop();
}


void ReduceIntegralValues() {
real rb; //Recieve Buffer
rb=0.0;
  CommTimer.start();
#if defined PARMPI
  MPI_Reduce(&ErrorMax,&rb,1,MPI_Type,MPI_MAX,0,MPI_COMM_WORLD);
  ErrorMax=rb;

  MPI_Reduce(&ErrorMid,&rb,1,MPI_Type,MPI_SUM,0,MPI_COMM_WORLD);
  ErrorMid=rb/size;

  MPI_Reduce(&ErrorL2,&rb,1,MPI_Type,MPI_SUM,0,MPI_COMM_WORLD);
  ErrorL2=rb/size;

  MPI_Reduce(&ErrorNb,&rb,1,MPI_Type,MPI_SUM,0,MPI_COMM_WORLD);
  ErrorNb=rb;

  MPI_Allreduce(&FlameVelocity,&rb,1,MPI_Type,MPI_SUM,MPI_COMM_WORLD);
  FlameVelocity=rb;

  MPI_Allreduce(&GlobalMomentum,&rb,1,MPI_Type,MPI_SUM,MPI_COMM_WORLD);
  GlobalMomentum=rb;

  MPI_Allreduce(&GlobalEnergy,&rb,1,MPI_Type,MPI_SUM,MPI_COMM_WORLD);
  GlobalEnergy=rb;

  MPI_Reduce(&ExactMomentum,&rb,1,MPI_Type,MPI_SUM,0,MPI_COMM_WORLD);
  ExactMomentum=rb;

  MPI_Reduce(&ExactEnergy,&rb,1,MPI_Type,MPI_SUM,0,MPI_COMM_WORLD);
  ExactEnergy=rb;
      
  MPI_Allreduce(&GlobalReactionRate,&rb,1,MPI_Type,MPI_SUM,MPI_COMM_WORLD);
  GlobalReactionRate=rb;

  MPI_Allreduce(&EigenvalueMax, &rb,1,MPI_Type,MPI_MAX,MPI_COMM_WORLD);
  EigenvalueMax=rb;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
#endif
  CommTimer.stop();
}
