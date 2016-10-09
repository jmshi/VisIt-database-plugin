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
//! \file meshblocktree.cpp
//  \brief implementation of functions in the MeshBlockTree class
// The MeshBlockTree stores the logical grid structure, and is used for neighbor
// searches, assigning global IDs, etc.  Level is defined as "logical level", where the
// logical root (single block) level is 0.  Note the logical level of the physical root
// grid (user-specified root grid) will be greater than zero if it contains more than
// one MeshBlock
//======================================================================================

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"

// this class header
#include "meshblock_tree.hpp"

//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree()
//  \brief constructor for the logical root

MeshBlockTree::MeshBlockTree()
{
  flag=true;
  gid=-1;
  loc.lx1=0, loc.lx2=0, loc.lx3=0;
  loc.level=0;
  pparent=NULL;
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        pleaf[k][j][i]=NULL;
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz)
//  \brief constructor for a leaf

MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz)
{
  flag=true;
  pparent=parent;
  gid=pparent->gid;
  loc.lx1=(parent->loc.lx1<<1)+ox;
  loc.lx2=(parent->loc.lx2<<1)+oy;
  loc.lx3=(parent->loc.lx3<<1)+oz;
  loc.level=parent->loc.level+1;
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        pleaf[k][j][i]=NULL;
      }
    }
  }
}


//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree::~MeshBlockTree()
//  \brief destructor (for both root and leaves)

MeshBlockTree::~MeshBlockTree()
{
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        if(pleaf[k][j][i]!=NULL)
          delete pleaf[k][j][i];
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CreateRootGrid(long int nx,long int ny, long int nz, int nl)
//  \brief create the root grid; the root grid can be incomplete (less than 8 leaves)

void MeshBlockTree::CreateRootGrid(long int nx, long int ny, long int nz, int nl)
{
  long int mx, my, mz;
  if (loc.level == nl) return;

  for(int k=0; k<=1; k++) {
    if((loc.lx3*2+k)*(1L<<(nl-loc.level-1)) < nz) {
      for(int j=0; j<=1; j++) {
        if((loc.lx2*2+j)*(1L<<(nl-loc.level-1)) < ny) {
          for(int i=0; i<=1; i++) {
            if((loc.lx1*2+i)*(1L<<(nl-loc.level-1)) < nx) {
              flag=false; // if there is a leaf, this is a node
              gid=-1;
              pleaf[k][j][i] = new MeshBlockTree(this, i, j, k);
              pleaf[k][j][i]->CreateRootGrid(nx, ny, nz, nl);
            }
          }
        }
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc,
//   int dim, enum BoundaryFlag* mesh_bcs, long int rbx, long int rby, long int rbz,
//   int rl, int &nnew)
//  \brief add a MeshBlock to the tree, also creates neighboring blocks
/*
void MeshBlockTree::AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc, int dim,
   enum BoundaryFlag* mesh_bcs, long int rbx, long int rby, long int rbz,
   int rl, int &nnew)
{
  int mx, my, mz;
  if(loc.level==rloc.level) return; // done

  if(flag==true) // leaf -> create the finer level
    Refine(root,dim,mesh_bcs,rbx,rby,rbz,rl,nnew);
  // get leaf indexes
  int sh = rloc.level-loc.level-1;
  mx=(int)((rloc.lx1>>sh)&1L);
  my=(int)((rloc.lx2>>sh)&1L);
  mz=(int)((rloc.lx3>>sh)&1L);
  pleaf[mz][my][mx]->AddMeshBlock(root,rloc,dim,mesh_bcs,rbx,rby,rbz,rl,nnew);

  return;
} */

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlockWithoutRefine(LogicalLocation rloc,
//                          long int rbx, long int rby, long int rbz, int rl)
//  \brief add a MeshBlock to the tree without refinement, used in restarting

void MeshBlockTree::AddMeshBlockWithoutRefine(LogicalLocation rloc,
                    long int rbx, long int rby, long int rbz, int rl)
{
  int mx, my, mz;
  if(loc.level==rloc.level) // done
    return;
  if(flag==true) // leaf -> create the finer level
    flag=false;
  // get leaf indexes
  int sh = rloc.level-loc.level-1;
  mx=(int)((rloc.lx1>>sh)&1L);
  my=(int)((rloc.lx2>>sh)&1L);
  mz=(int)((rloc.lx3>>sh)&1L);
  if(pleaf[mz][my][mx]==NULL)
    pleaf[mz][my][mx] = new MeshBlockTree(this, mx, my, mz);
  pleaf[mz][my][mx]->AddMeshBlockWithoutRefine(rloc,rbx,rby,rbz,rl);
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CountMeshBlock(int& count)
//  \brief creates the Location list sorted by Z-ordering

void MeshBlockTree::CountMeshBlock(int& count)
{
  if(pparent==NULL) count=0;
  if(flag==true)
    count++;
  else {
    for(int k=0; k<=1; k++) {
      for(int j=0; j<=1; j++) {
        for(int i=0; i<=1; i++) {
          if(pleaf[k][j][i]!=NULL)
            pleaf[k][j][i]->CountMeshBlock(count);
        }
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::GetMeshBlockList(LogicalLocation *list,
//                                           int *pglist, int& count)
//  \brief creates the Location list sorted by Z-ordering

void MeshBlockTree::GetMeshBlockList(LogicalLocation *list, int *pglist, int& count)
{
  if(pparent==NULL) count=0;
  if(flag==true) {
    list[count]=loc;
    if(pglist!=NULL)
      pglist[count]=gid;
    gid=count;
    count++;
  }
  else {
    for(int k=0; k<=1; k++) {
      for(int j=0; j<=1; j++) {
        for(int i=0; i<=1; i++) {
          if(pleaf[k][j][i]!=NULL)
            pleaf[k][j][i]->GetMeshBlockList(list, pglist, count);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::FindMeshBlock(LogicalLocation tloc)
//  \brief find MeshBlock with LogicalLocation tloc and return a pointer

MeshBlockTree* MeshBlockTree::FindMeshBlock(LogicalLocation tloc)
{
  if(tloc.level==loc.level) return this;
  // get leaf indexes
  int sh = tloc.level-loc.level-1;
  int mx=(int)((tloc.lx1>>sh)&1L);
  int my=(int)((tloc.lx2>>sh)&1L);
  int mz=(int)((tloc.lx3>>sh)&1L);
  if(pleaf[mz][my][mx]==NULL) return NULL;
  else return pleaf[mz][my][mx]->FindMeshBlock(tloc);
}
