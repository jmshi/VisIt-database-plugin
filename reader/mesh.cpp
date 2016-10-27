//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in Mesh class
//======================================================================================

// C/C++ headers
#include <cfloat>     // FLT_MAX
#include <cmath>      // std::abs(), pow()
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // sort
#include <iomanip>
#include <stdlib.h>
#include <string.h>  // memcpy

// Athena++ classes headers
#include "athena.hpp"
#include "globals.hpp"
#include "athena_arrays.hpp"
#include "parameter_input.hpp"
#include "io_wrapper.hpp"
//#include "../utils/buffer_utils.hpp"
//#include "../reconstruct/reconstruction.hpp"
//#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"
#include "mesh.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


//--------------------------------------------------------------------------------------
// was Mesh constructor for restarts. Load the restart file
// Now use as a test for suitable reader for VisIt

Mesh::Mesh(ParameterInput *pin, IOWrapper& resfile, int flag)
{
  std::stringstream msg;
  RegionSize block_size;
  enum BoundaryFlag block_bcs[6];
  MeshBlock *pfirst;
  int i, j, dim;
  IOWrapperSize_t *offset, datasize, listsize, headeroffset;


  // --(1) read time and cycle limits from input file; optional
  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  cfl_number = pin->GetReal("time","cfl_number");
  nlim = pin->GetOrAddInteger("time","nlim",-1);
  // other parameter might be needed
  // NGHOST NHYDRO etc. which are imported with defs.hpp for now
  // The future data dump should include those parameters so that
  // we don't have to rely on defs.hpp



  // skip the following mesh_size as it was set later in restart file
  // read physical size of mesh (root level) from input file.
  //mesh_size.x1min = pin->GetReal("mesh","x1min");
  //mesh_size.x2min = pin->GetReal("mesh","x2min");
  //mesh_size.x3min = pin->GetReal("mesh","x3min");
  //mesh_size.x1max = pin->GetReal("mesh","x1max");
  //mesh_size.x2max = pin->GetReal("mesh","x2max");
  //mesh_size.x3max = pin->GetReal("mesh","x3max");
  // read ratios of grid cell size in each direction
  //block_size.x1rat = mesh_size.x1rat = pin->GetOrAddReal("mesh","x1rat",1.0);
  //block_size.x2rat = mesh_size.x2rat = pin->GetOrAddReal("mesh","x2rat",1.0);
  //block_size.x3rat = mesh_size.x3rat = pin->GetOrAddReal("mesh","x3rat",1.0);
  // read MeshBlock parameters
  //block_size.nx1 = pin->GetOrAddInteger("meshblock","nx1",mesh_size.nx1);
  //block_size.nx2 = pin->GetOrAddInteger("meshblock","nx2",mesh_size.nx2);
  //block_size.nx3 = pin->GetOrAddInteger("meshblock","nx3",mesh_size.nx3);

  // user defined mesh input data is not implemented
  nint_user_mesh_data_=0, nreal_user_mesh_data_=0;
  nbnew=0; nbdel=0; int udsize=0;

  // get the end of the header
  headeroffset=resfile.GetPosition();
  // read the restart file
  // the file is already open and the pointer is set to after <par_end>
  IOWrapperSize_t headersize = sizeof(int)*3+sizeof(Real)*2
                             + sizeof(RegionSize)+sizeof(IOWrapperSize_t);
  char *headerdata = new char[headersize];
  if(Globals::my_rank==0) { // the master process reads the header data
    if(resfile.Read(headerdata,1,headersize)!=headersize) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  // --(2) read int             nbtotal, rootlevel, ncycle
  //          Real            time, dt
  //          RegionSize      mesh_size
  //          IOWrapperSize_t datasize
  IOWrapperSize_t hdos = 0;
  memcpy(&nbtotal, &(headerdata[hdos]), sizeof(int));
  hdos+=sizeof(int);
  memcpy(&root_level, &(headerdata[hdos]), sizeof(int));
  hdos+=sizeof(int);
  current_level=root_level;
  memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos+=sizeof(RegionSize);
  memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos+=sizeof(Real);
  memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos+=sizeof(Real);
  memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  hdos+=sizeof(int);
  memcpy(&datasize, &(headerdata[hdos]), sizeof(IOWrapperSize_t));
  hdos+=sizeof(IOWrapperSize_t);

  dim=1;
  if(mesh_size.nx2>1) dim=2;
  if(mesh_size.nx3>1) dim=3;


  // --(3) reconstruct the meshblocks
  // --(3-1) initialize
  loclist=new LogicalLocation[nbtotal];
  offset=new IOWrapperSize_t[nbtotal];
  costlist=new Real[nbtotal];
  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];

  block_size.nx1 = pin->GetOrAddReal("meshblock","nx1",mesh_size.nx1);
  block_size.nx2 = pin->GetOrAddReal("meshblock","nx2",mesh_size.nx2);
  block_size.nx3 = pin->GetOrAddReal("meshblock","nx3",mesh_size.nx3);

  // calculate the number of the blocks
  nrbx1=mesh_size.nx1/block_size.nx1;
  nrbx2=mesh_size.nx2/block_size.nx2;
  nrbx3=mesh_size.nx3/block_size.nx3;

  //initialize user-enrollable functions
  if(mesh_size.x1rat!=1.0)
    user_meshgen_[X1DIR]=true;
  else
    user_meshgen_[X1DIR]=false;
  if(mesh_size.x2rat!=1.0)
    user_meshgen_[X2DIR]=true;
  else
    user_meshgen_[X2DIR]=false;
  if(mesh_size.x3rat!=1.0)
    user_meshgen_[X3DIR]=true;
  else
    user_meshgen_[X3DIR]=false;
  MeshGenerator_[X1DIR]=DefaultMeshGeneratorX1;
  MeshGenerator_[X2DIR]=DefaultMeshGeneratorX2;
  MeshGenerator_[X3DIR]=DefaultMeshGeneratorX3;

  multilevel=false;
  adaptive=false;
  if (pin->GetOrAddString("mesh","refinement","none")=="adaptive")
    adaptive=true, multilevel=true;
  else if (pin->GetOrAddString("mesh","refinement","none")=="static")
    multilevel=true;
  if (adaptive==true) {
    max_level = pin->GetOrAddInteger("mesh","numlevel",1)+root_level-1;
    if(max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63-root_level+1 << "." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  } else {
    max_level = 63;
  }

  // --(3-2) read the ID list: LogicalLocation  loclist[nbtotal]
  //                           Real             costlist[nbtotal]
  listsize=sizeof(LogicalLocation)+sizeof(Real);
  //allocate the idlist buffer
  char *idlist = new char [listsize*nbtotal];
  if(Globals::my_rank==0) { // only the master process reads the ID list
    if(resfile.Read(idlist,listsize,nbtotal)!=nbtotal) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }

  int os=0;
  for(int i=0;i<nbtotal;i++) {
    memcpy(&(loclist[i]), &(idlist[os]), sizeof(LogicalLocation));
    os+=sizeof(LogicalLocation);
    memcpy(&(costlist[i]), &(idlist[os]), sizeof(Real));
    os+=sizeof(Real);
    if(loclist[i].level>current_level) current_level=loclist[i].level;
  }
  delete [] idlist;

  // calculate the header offset and seek
  headeroffset+=headersize+udsize+listsize*nbtotal;
  if(Globals::my_rank!=0)
    resfile.Seek(headeroffset);


  // --(3-3) rebuild the Block Tree
  for(int i=0;i<nbtotal;i++)
    tree.AddMeshBlockWithoutRefine(loclist[i],nrbx1,nrbx2,nrbx3,root_level);
  int nnb;
  // check the tree structure, and assign GID
  tree.GetMeshBlockList(loclist, NULL, nnb);
  if(nnb!=nbtotal) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Tree reconstruction failed. The total numbers of the blocks do not match. ("
        << nbtotal << " != " << nnb << ")" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  LoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // --(4) construct coordinates and read rst data
  // --(4-1) allocate data buffer
  int nb=nblist[Globals::my_rank];
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nb-1;
  char *mbdata = new char [datasize*nb];
  // initialize grid indices
  int is,ie,js,je,ks,ke;
  is = NGHOST;
  ie = is + block_size.nx1 - 1;
  if (block_size.nx2 > 1) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
  } else {
    js = je = 0;
  }
  if (block_size.nx3 > 1) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
  } else {
    ks = ke = 0;
  }
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);
  // allocate the coordinates
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  // allocate the variables
  u.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);
  // determine if MHD enabled
  bool magnetic=false;
  std::cout << u.GetSizeInBytes() << " " << datasize << std::endl;
  if (u.GetSizeInBytes() == datasize) magnetic = false;
  else  magnetic = true;
  if (magnetic) {
    //  Note the extra cell in each longitudinal dirn for interface fields
	b.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
	b.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
	b.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  }

  // --(4-2) load MeshBlocks (parallel)
  if(resfile.Read_at_all(mbdata, datasize, nb, headeroffset+nbs*datasize)!=nb) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  // --(4-3) set x1f x2f x3f
  //             u and b.x1f b.x2f b.x3f for given meshblock
  // we then write out grids and data for each blocks for debug

  // open '*.debug' files
  std::ofstream fp_x1f("x1f.debug");
  std::ofstream fp_x2f("x2f.debug");
  std::ofstream fp_x3f("x3f.debug");
  std::ofstream fp_con("con.debug");
  std::ofstream fp_mag("mag.debug");
  //FILE *fp_x1f, *fp_x2f, *fp_x3f;
  //if ((fp_x1f = fopen("x1f.debug","wb")) == NULL) {
  //  std::cout << "Cannot open x1f.debug" << std::endl;
  //  return;
  //  }
  //if ((fp_x2f = fopen("x2f.debug","wb")) == NULL) {
  //  std::cout << "Cannot open x2f.debug" << std::endl;
  //  return;
  //  }
  //if ((fp_x3f = fopen("x3f.debug","wb")) == NULL) {
  //  std::cout << "Cannot open x3f.debug" << std::endl;
  //  return;
  //  }
  //FILE *fp_cons;
  //if ((fp_x3f = fopen("cons.debug","wb")) == NULL) {
  //  std::cout << "Cannot open cons.debug" << std::endl;
  //  return;
  //  }

  int num;
  for(num=nbs;num<=nbe;num++) {
	char *mbdata_current;
    int buff_os = datasize * (num-nbs);
	mbdata_current = mbdata+buff_os;
	// update block_size.x1min etc. block_bcs etc.
    SetBlockSizeAndBoundaries(loclist[num], block_size, block_bcs);
	// update x1f x2f x3f
    AllocateAndSetBasicCoordinates(is,ie,js,je,ks,ke,loclist[num],block_size,block_bcs);
	// 2) allocate phydro->u  pfield->b
	int os = 0;
	memcpy(u.data(), &mbdata_current[os],u.GetSizeInBytes());
	os += u.GetSizeInBytes();
	if(magnetic) {
	  memcpy(b.x1f.data(), &(mbdata_current[os]), b.x1f.GetSizeInBytes());
	  os += b.x1f.GetSizeInBytes();
	  memcpy(b.x2f.data(), &(mbdata_current[os]), b.x2f.GetSizeInBytes());
	  os += b.x2f.GetSizeInBytes();
	  memcpy(b.x3f.data(), &(mbdata_current[os]), b.x3f.GetSizeInBytes());
	  os += b.x3f.GetSizeInBytes();
	}

	for(int i=0;i<ncells1+1;i++)
	  fp_x1f << x1f(i) << " ";
	fp_x1f << "\n";
	for(int j=0;j<ncells2+1;j++)
	  fp_x2f << x2f(j) << " ";
	fp_x2f << "\n";
	for(int k=0;k<ncells3+1;k++)
	  fp_x3f << x3f(k) << " ";
	fp_x3f << "\n";

	fp_con << "-- meshblock " << (num-nbs) << "--\n\n";
    for(int n=0;n<NHYDRO;n++){
	  fp_con << "u[" << n << "]\n";
	  for(int k=0;k<ncells3;k++){
	  for(int j=0;j<ncells2;j++){
	  for(int i=0;i<ncells1;i++){
        fp_con << u(n,k,j,i) << " ";
	  }
	  fp_con << "\n";
	  }
	  fp_con << "\n\n";
	  }
	}
    if(magnetic) {
      for(int n=0;n<3;n++){
	    fp_mag << "b" << n << "\n";
	    for(int k=0;k<ncells3;k++){
	    for(int j=0;j<ncells2;j++){
	    for(int i=0;i<ncells1;i++){
          if (n==0) {
	        fp_mag << b.x1f(k,j,i) << " ";
	        if (i==ncells1-1) fp_mag << b.x1f(k,j,ncells1) << " ";
	      }
          if (n==1) {
	        fp_mag << b.x2f(k,j,i) << " ";
	        if (j==ncells2-1) fp_mag << b.x2f(k,ncells2,i) << " ";
	      }
          if (n==2) {
	        fp_mag << b.x3f(k,j,i) << " ";
	        if (k==ncells1-1) fp_mag << b.x3f(ncells3,j,i) << " ";
	      }
	    }
	    fp_con << "\n";
	    }
	    fp_con << "\n\n";
	    }
	}}
  }
  delete [] mbdata;
  // clean up
  delete [] offset;
}

//--------------------------------------------------------------------------------------
// destructor

Mesh::~Mesh()
{
  //delete pblock;
  delete [] nslist;
  delete [] nblist;
  delete [] ranklist;
  delete [] costlist;
  delete [] loclist;
  //if(adaptive==true) { // deallocate arrays for AMR
  //  delete [] nref;
  //  delete [] nderef;
  //  delete [] rdisp;
  //  delete [] ddisp;
  //  delete [] bnref;
  //  delete [] bnderef;
  //  delete [] brdisp;
  //  delete [] bddisp;
  //}
  //// delete user Mesh data
  //for(int n=0; n<nreal_user_mesh_data_; n++)
  //  ruser_mesh_data[n].DeleteAthenaArray();
  //if(nreal_user_mesh_data_>0) delete [] ruser_mesh_data;
  //for(int n=0; n<nint_user_mesh_data_; n++)
  //  iuser_mesh_data[n].DeleteAthenaArray();
  //if(nint_user_mesh_data_>0) delete [] iuser_mesh_data;
}


//--------------------------------------------------------------------------------------
// \!fn void Mesh::LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb)
// \brief Calculate distribution of MeshBlocks based on the cost list

void Mesh::LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb)
{
  std::stringstream msg;
  Real totalcost=0, maxcost=0.0, mincost=(FLT_MAX);

  for(int i=0; i<nb; i++) {
    totalcost+=clist[i];
    mincost=std::min(mincost,clist[i]);
    maxcost=std::max(maxcost,clist[i]);
  }
  int j=(Globals::nranks)-1;
  Real targetcost=totalcost/Globals::nranks;
  Real mycost=0.0;
  // create rank list from the end: the master node should have less load
  for(int i=nb-1;i>=0;i--) {
    if(targetcost==0.0) {
      msg << "### FATAL ERROR in LoadBalance" << std::endl
          << "There is at least one process which has no MeshBlock" << std::endl
          << "Decrease the number of processes or use smaller MeshBlocks." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mycost+=clist[i];
    rlist[i]=j;
    if(mycost >= targetcost && j>0) {
      j--;
      totalcost-=mycost;
      mycost=0.0;
      targetcost=totalcost/(j+1);
    }
  }
  slist[0]=0;
  j=0;
  for(int i=1;i<nb;i++) { // make the list of nbstart and nblocks
    if(rlist[i]!=rlist[i-1]) {
      nlist[j]=i-nslist[j];
      slist[++j]=i;
    }
  }
  nlist[j]=nb-slist[j];

#ifdef MPI_PARALLEL
  if(nb % Globals::nranks != 0 && adaptive == false
  && maxcost == mincost && Globals::my_rank==0) {
    std::cout << "### Warning in LoadBalance" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will cause a poor load balance." << std::endl;
  }
#endif
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc,
//                 RegionSize &block_size, enum BundaryFlag *block_bcs)
// \brief Set the physical part of a block_size structure and block boundary conditions

void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                     enum BoundaryFlag *block_bcs)
{
  long int &lx1=loc.lx1;
  long int &lx2=loc.lx2;
  long int &lx3=loc.lx3;
  int &ll=loc.level;
  // calculate physical block size, x1
  if(lx1==0) {
    block_size.x1min=mesh_size.x1min;
    block_bcs[INNER_X1]=mesh_bcs[INNER_X1];
  }
  else {
    Real rx=(Real)lx1/(Real)(nrbx1<<(ll-root_level));
    block_size.x1min=MeshGenerator_[X1DIR](rx,mesh_size);
    block_bcs[INNER_X1]=BLOCK_BNDRY;
  }
  if(lx1==(nrbx1<<(ll-root_level))-1) {
    block_size.x1max=mesh_size.x1max;
    block_bcs[OUTER_X1]=mesh_bcs[OUTER_X1];
  }
  else {
    Real rx=(Real)(lx1+1)/(Real)(nrbx1<<(ll-root_level));
    block_size.x1max=MeshGenerator_[X1DIR](rx,mesh_size);
    block_bcs[OUTER_X1]=BLOCK_BNDRY;
  }

  // calculate physical block size, x2
  if(mesh_size.nx2 == 1) {
    block_size.x2min=mesh_size.x2min;
    block_size.x2max=mesh_size.x2max;
    block_bcs[INNER_X2]=mesh_bcs[INNER_X2];
    block_bcs[OUTER_X2]=mesh_bcs[OUTER_X2];
  }
  else {
    if(lx2==0) {
      block_size.x2min=mesh_size.x2min;
      block_bcs[INNER_X2]=mesh_bcs[INNER_X2];
    }
    else {
      Real rx=(Real)lx2/(Real)(nrbx2<<(ll-root_level));
      block_size.x2min=MeshGenerator_[X2DIR](rx,mesh_size);
      block_bcs[INNER_X2]=BLOCK_BNDRY;
    }
    if(lx2==(nrbx2<<(ll-root_level))-1) {
      block_size.x2max=mesh_size.x2max;
      block_bcs[OUTER_X2]=mesh_bcs[OUTER_X2];
    }
    else {
      Real rx=(Real)(lx2+1)/(Real)(nrbx2<<(ll-root_level));
      block_size.x2max=MeshGenerator_[X2DIR](rx,mesh_size);
      block_bcs[OUTER_X2]=BLOCK_BNDRY;
    }
  }

  // calculate physical block size, x3
  if(mesh_size.nx3 == 1) {
    block_size.x3min=mesh_size.x3min;
    block_size.x3max=mesh_size.x3max;
    block_bcs[INNER_X3]=mesh_bcs[INNER_X3];
    block_bcs[OUTER_X3]=mesh_bcs[OUTER_X3];
  }
  else {
    if(lx3==0) {
      block_size.x3min=mesh_size.x3min;
      block_bcs[INNER_X3]=mesh_bcs[INNER_X3];
    }
    else {
      Real rx=(Real)lx3/(Real)(nrbx3<<(ll-root_level));
      block_size.x3min=MeshGenerator_[X3DIR](rx,mesh_size);
      block_bcs[INNER_X3]=BLOCK_BNDRY;
    }
    if(lx3==(nrbx3<<(ll-root_level))-1) {
      block_size.x3max=mesh_size.x3max;
      block_bcs[OUTER_X3]=mesh_bcs[OUTER_X3];
    }
    else {
      Real rx=(Real)(lx3+1)/(Real)(nrbx3<<(ll-root_level));
      block_size.x3max=MeshGenerator_[X3DIR](rx,mesh_size);
      block_bcs[OUTER_X3]=BLOCK_BNDRY;
    }
  }
  return;
}

inline void Mesh::AllocateAndSetBasicCoordinates(const int is,const int ie,const int js,const int je,const int ks,const int ke, LogicalLocation loc, RegionSize &block_size, enum BoundaryFlag *block_bcs)
{

  // allocate arrays for sizes and positions of cells
  int ng=NGHOST;
  int cflag = 0;
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

  long long nrootmesh, noffset;
  long int &lx1=loc.lx1;
  long int &lx2=loc.lx2;
  long int &lx3=loc.lx3;
  int &ll=loc.level;

  // X1-DIRECTION: initialize sizes and positions of cell FACES (dx1f,x1f)
  nrootmesh=mesh_size.nx1*(1L<<(ll-root_level));
  if(user_meshgen_[X1DIR]==false) { // uniform
    Real dx=(block_size.x1max-block_size.x1min)/(ie-is+1);
    for(int i=is-ng; i<=ie+ng; ++i)
      dx1f(i)=dx;
    x1f(is-ng)=block_size.x1min-ng*dx;
    for(int i=is-ng+1;i<=ie+ng+1;i++)
      x1f(i)=x1f(i-1)+dx;
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
  }
  else {
    for (int i=is-ng; i<=ie+ng+1; ++i) {
      // if there are too many levels, this won't work or be precise enough
      noffset=((i-is)<<cflag)+(long long)lx1*block_size.nx1;
      Real rx=(Real)noffset/(Real)nrootmesh;
      x1f(i)=MeshGenerator_[X1DIR](rx,mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
    for(int i=is-ng; i<=ie+ng; ++i)
      dx1f(i)=x1f(i+1)-x1f(i);
  }

  // correct cell face positions in ghost zones for reflecting boundary condition
  if (block_bcs[INNER_X1] == REFLECTING_BNDRY) {
    for (int i=1; i<=ng; ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }
  if (block_bcs[OUTER_X1] == REFLECTING_BNDRY) {
    for (int i=1; i<=ng; ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

  // X2-DIRECTION: initialize spacing and positions of cell FACES (dx2f,x2f)
  if(ncells2 > 1) {
    nrootmesh=mesh_size.nx2*(1L<<(ll-root_level));
    if(user_meshgen_[X2DIR]==false) { // uniform
      Real dx=(block_size.x2max-block_size.x2min)/(je-js+1);
      for(int j=js-ng; j<=je+ng; ++j)
        dx2f(j)=dx;
      x2f(js-ng)=block_size.x2min-ng*dx;
      for(int j=js-ng+1;j<=je+ng+1;j++)
        x2f(j)=x2f(j-1)+dx;
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
    }
    else {
      for (int j=js-ng; j<=je+ng+1; ++j) {
        // if there are too many levels, this won't work or be precise enough
        noffset=((j-js)<<cflag)+(long long)lx2*block_size.nx2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x2f(j)=MeshGenerator_[X2DIR](rx,mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
      for(int j=js-ng; j<=je+ng; ++j)
        dx2f(j)=x2f(j+1)-x2f(j);
    }

    // correct cell face positions in ghost zones for reflecting boundary condition
    if (block_bcs[INNER_X2] == REFLECTING_BNDRY
     || block_bcs[INNER_X2] == POLAR_BNDRY) { // also polar boundary
      for (int j=1; j<=ng; ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }
    if (block_bcs[OUTER_X2] == REFLECTING_BNDRY
     || block_bcs[OUTER_X2] == POLAR_BNDRY) { // also polar boundary
      for (int j=1; j<=ng; ++j) {
        dx2f(je+j  ) = dx2f(je-j+1);
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    }
  }
  else {
    dx2f(js) = block_size.x2max-block_size.x2min;
    x2f(js  ) = block_size.x2min;
    x2f(je+1) = block_size.x2max;
  }


  // X3-DIRECTION: initialize spacing and positions of cell FACES (dx3f,x3f)
  if(ncells3 > 1) {
    nrootmesh=mesh_size.nx3*(1L<<(ll-root_level));
    if(user_meshgen_[X3DIR]==false) { // uniform
      Real dx=(block_size.x3max-block_size.x3min)/(ke-ks+1);
      for(int k=ks-ng; k<=ke+ng; ++k)
        dx3f(k)=dx;
      x3f(ks-ng)=block_size.x3min-ng*dx;
      for(int k=ks-ng+1;k<=ke+ng+1;k++)
        x3f(k)=x3f(k-1)+dx;
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
    }
    else {
      for (int k=ks-ng; k<=ke+ng+1; ++k) {
        // if there are too many levels, this won't work or be precise enough
        noffset=((k-ks)<<cflag)+(long long)lx3*block_size.nx3;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x3f(k)=MeshGenerator_[X3DIR](rx,mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
      for(int k=ks-ng; k<=ke+ng; ++k)
        dx3f(k)=x3f(k+1)-x3f(k);
    }

    // correct cell face positions in ghost zones for reflecting boundary condition
    if (block_bcs[INNER_X3] == REFLECTING_BNDRY) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }
    if (block_bcs[OUTER_X3] == REFLECTING_BNDRY) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ke+k  ) = dx3f(ke-k+1);
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    }
  }
  else {
    dx3f(ks) = block_size.x3max-block_size.x3min;
    x3f(ks  ) = block_size.x3min;
    x3f(ke+1) = block_size.x3max;
  }

  return;
}
