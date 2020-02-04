/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: 
------------------------------------------------------------------------- */

#include "pair_mbx.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairMBX::PairMBX(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;
  restartinfo = 0;

  mbx_total_energy = 0;

  // energy terms available to pair compute
  
  nextra = 11;
  pvector = new double[nextra];
}

/* ---------------------------------------------------------------------- */

PairMBX::~PairMBX()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }

  delete [] pvector;
}

/* ---------------------------------------------------------------------- */

void PairMBX::setup()
{

  fix_mbx = NULL;
  int ifix = modify->find_fix_by_style("mbx");
  fix_mbx = (FixMBX *) modify->fix[ifix];

  if(!fix_mbx) error->all(FLERR,"Fix MBX not found");
}

/* ---------------------------------------------------------------------- */

void PairMBX::compute(int eflag, int vflag)
{
  
#if 0

  compute_full();

#else
  
  // compute energy+gradients in parallel

  bblock::System * ptr_mbx = fix_mbx->ptr_mbx;
  
  // compute energy

  fix_mbx->mbxt_start(MBXT_E1B);
  mbx_e1b = ptr_mbx->OneBodyEnergy(true);
  fix_mbx->mbxt_stop(MBXT_E1B);
  accumulate_f();

  fix_mbx->mbxt_start(MBXT_E2B_LOCAL);
  double mbx_e2b_local = ptr_mbx->TwoBodyEnergy(true);
  fix_mbx->mbxt_stop(MBXT_E2B_LOCAL);
  accumulate_f();
  
  fix_mbx->mbxt_start(MBXT_E2B_GHOST);
  double mbx_e2b_ghost = ptr_mbx->TwoBodyEnergy(true, true);
  fix_mbx->mbxt_stop(MBXT_E2B_GHOST);
  accumulate_f();

  mbx_e2b = mbx_e2b_local + mbx_e2b_ghost;

  fix_mbx->mbxt_start(MBXT_E3B_LOCAL);
  double mbx_e3b_local = ptr_mbx->ThreeBodyEnergy(true);
  fix_mbx->mbxt_stop(MBXT_E3B_LOCAL);
  accumulate_f();
  
  fix_mbx->mbxt_start(MBXT_E3B_GHOST);
  double mbx_e3b_ghost = ptr_mbx->ThreeBodyEnergy(true, true);
  fix_mbx->mbxt_stop(MBXT_E3B_GHOST);
  accumulate_f();

  mbx_e3b = mbx_e3b_local + mbx_e3b_ghost;

  // compute energy+gradients in serial on rank 0 for full system

  mbx_disp = 0.0;
  mbx_buck = 0.0;
  mbx_ele  = 0.0;
  
  bblock::System * ptr_mbx_full = fix_mbx->ptr_mbx_full;
  
  fix_mbx->mbxt_start(MBXT_DISP);
  if(comm->me == 0) mbx_disp = ptr_mbx_full->Dispersion(true);
  fix_mbx->mbxt_stop(MBXT_DISP);
  accumulate_f_full();

  
  fix_mbx->mbxt_start(MBXT_BUCK);
  if(comm->me == 0) mbx_buck = ptr_mbx_full->Buckingham(true);
  fix_mbx->mbxt_stop(MBXT_BUCK);
  accumulate_f_full();

  fix_mbx->mbxt_start(MBXT_ELE);
  if(comm->me == 0) mbx_ele = ptr_mbx_full->Electrostatics(true);
  fix_mbx->mbxt_stop(MBXT_ELE);
  accumulate_f_full();

  mbx_total_energy = mbx_e1b + mbx_e2b + mbx_disp + mbx_buck + mbx_e3b + mbx_ele;

  // save total energy from mbx as vdwl
  
  if(eflag) {
    eng_vdwl = mbx_total_energy;

    // generally useful
    
    pvector[0] = mbx_e1b;
    pvector[1] = mbx_e2b;
    pvector[2] = mbx_e3b;
    pvector[3] = mbx_disp;
    pvector[4] = mbx_buck;
    pvector[5] = mbx_ele;
    pvector[6] = mbx_total_energy;
    
    // for debugging
    
    pvector[7]  = mbx_e2b_local;
    pvector[8]  = mbx_e2b_ghost;
    pvector[9]  = mbx_e3b_local;
    pvector[10] = mbx_e3b_ghost;
  }

#if 0
  // debug output
  
  double e1  = 0.0;
  double e2l = 0.0;
  double e2g = 0.0;
  double e3l = 0.0;
  double e3g = 0.0;

  MPI_Reduce(&mbx_e1b,       &e1,  1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&mbx_e2b_local, &e2l, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&mbx_e2b_ghost, &e2g, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&mbx_e3b_local, &e3l, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&mbx_e3b_ghost, &e3g, 1, MPI_DOUBLE, MPI_SUM, 0, world);

  double etot = e1 + e2l + e2g + e3l + e3g + mbx_disp + mbx_buck + mbx_ele;

  if(comm->me == 0) {
    printf("mbx_e1b=   %f  Parallel\n",e1);
    printf("mbx_e2b=   %f (%f, %f)  Parallel\n",e2l+e2g, e2l, e2g);
    printf("mbx_e3b=   %f (%f, %f)  Parallel\n",e3l+e3g, e3l, e3g);
    printf("mbx_disp=  %f  Serial\n",mbx_disp);
    printf("mbx_buck=  %f  Serial\n",mbx_buck);
    printf("mbx_ele=   %f  Serial\n",mbx_ele);
    printf("mbx_total= %f\n",etot);
  }
#endif
  
#endif
  
}

/* ---------------------------------------------------------------------- */
// compute energy of full system on rank 0
/* ---------------------------------------------------------------------- */

void PairMBX::compute_full()
{
  
  // Update coordinates in MBX library

  //  update_mbx_xyz();

  bblock::System * ptr_mbx = fix_mbx->ptr_mbx_full;

#if 0
  
  mbx_total_energy = ptr_mbx->Energy(true);
  accumulate_f_full();

#else
  
  mbx_e1b = ptr_mbx->OneBodyEnergy(true);
  accumulate_f_full();

  mbx_e2b = ptr_mbx->TwoBodyEnergy(true);
  accumulate_f_full();

  mbx_e3b = ptr_mbx->ThreeBodyEnergy(true);
  accumulate_f_full();
  
  mbx_disp = ptr_mbx->Dispersion(true);
  accumulate_f_full();
  
  mbx_buck = ptr_mbx->Buckingham(true);
  accumulate_f_full();

  mbx_ele = ptr_mbx->Electrostatics(true);
  accumulate_f_full();
  
  mbx_total_energy = mbx_e1b + mbx_e2b + mbx_disp + mbx_buck + mbx_e3b + mbx_ele;

  double mbx_e2b_local = mbx_e2b;
  double mbx_e2b_ghost = 0.0;
  double mbx_e3b_local = mbx_e3b;
  double mbx_e3b_ghost = 0.0;

  if(comm->me == 0) {
    printf("mbx_e1b=   %f\n",mbx_e1b);
    printf("mbx_e2b=   %f (%f, %f)\n",mbx_e2b, mbx_e2b_local, mbx_e2b_ghost);
    printf("mbx_e3b=   %f (%f, %f)\n",mbx_e3b, mbx_e3b_local, mbx_e3b_ghost);
    printf("mbx_disp=  %f\n",mbx_disp);
    printf("mbx_buck=  %f\n",mbx_buck);
    printf("mbx_ele=   %f\n",mbx_ele);
    printf("mbx_total= %f\n",mbx_total_energy);
  }
#endif
  
  // save total energy from mbx as vdwl
  
  eng_vdwl = mbx_total_energy;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMBX::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMBX::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMBX::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMBX::init_style()
{
  // request regular or rRESPA neighbor list

  int irequest;
  int respa = 0;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
  }

  // currently request neighbor list, but we don't use it for anything...
  
  irequest = neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMBX::init_one(int i, int j)
{
  if (setflag[i][j] == 0) cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMBX::write_data(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMBX::write_data_all(FILE *fp)
{
}

/* ---------------------------------------------------------------------- */

void *PairMBX::extract(const char *str, int &dim)
{
  dim = 2;
  // if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  // if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}

/* ----------------------------------------------------------------------
   update forces with MBX contribution
------------------------------------------------------------------------- */

void PairMBX::accumulate_f()
{
  fix_mbx->mbxt_start(MBXT_ACCUMULATE_F);
  
  bblock::System * ptr_mbx = fix_mbx->ptr_mbx;

  const int nlocal = atom->nlocal;
  double ** f = atom->f;

  const int * const mol_anchor = fix_mbx->mol_anchor;
  const int * const mol_type = fix_mbx->mol_type;
  char ** mol_names = fix_mbx->mol_names;
  
  std::vector<double> grads = ptr_mbx->GetRealGrads();

  // accumulate forces on local particles
  // -- forces on ghost particles ignored/not needed
  // -- should use a map created from earlier loop loading particles into mbx

  int indx = 0;
  
  //  for(int i=0; i<nall; ++i) {
  for(int i=0; i<nlocal; ++i) {

    if(mol_anchor[i]) {

      const int mtype = mol_type[i];

      // to be replaced with integer comparison
      
      int na = 0;
      if(strcmp("h2o",      mol_names[mtype]) == 0) na = 3;
      else if(strcmp("na",  mol_names[mtype]) == 0) na = 1;
      else if(strcmp("cl",  mol_names[mtype]) == 0) na = 1;
      else if(strcmp("co2", mol_names[mtype]) == 0) na = 3;

      tagint anchor = atom->tag[i];

      for(int j=0; j<na; ++j) {

	const int ii = atom->map(anchor + j);
	f[ii][0] -= grads[indx++];
	f[ii][1] -= grads[indx++];
	f[ii][2] -= grads[indx++];
      }
     
    } // if(anchor)

  }
 
  fix_mbx->mbxt_stop(MBXT_ACCUMULATE_F);
}

/* ----------------------------------------------------------------------
   update forces with MBX contribution from full system
------------------------------------------------------------------------- */

void PairMBX::accumulate_f_full()
{
  fix_mbx->mbxt_start(MBXT_ACCUMULATE_F_FULL);
  
  //  printf("[MBX] Inside accumulate_f_full()\n");
  
  // master rank retrieves forces

  double ** f_full = fix_mbx->f_full;
    
  if(comm->me == 0) {
    
    bblock::System * ptr_mbx = fix_mbx->ptr_mbx_full;

    const int natoms = atom->natoms;
    
    const int * const mol_anchor_full = fix_mbx->mol_anchor_full;
    const int * const mol_type_full = fix_mbx->mol_type_full;
    const tagint * const tag_full = fix_mbx->tag_full;
    const int * const atom_map_full = fix_mbx->atom_map_full;
    char ** mol_names = fix_mbx->mol_names;
    
    std::vector<double> grads = ptr_mbx->GetRealGrads();
    
    // accumulate forces on local particles
    // -- forces on ghost particles ignored/not needed
    // -- should use a map created from earlier loop loading particles into mbx
    
    int indx = 0;
    
    for(int i=0; i<natoms; ++i) {
      
      if(mol_anchor_full[i]) {
	
	const int mtype = mol_type_full[i];
	
	// to be replaced with integer comparison
	
	int na = 0;
	if(strcmp("h2o",      mol_names[mtype]) == 0) na = 3;
	else if(strcmp("na",  mol_names[mtype]) == 0) na = 1;
	else if(strcmp("cl",  mol_names[mtype]) == 0) na = 1;
	else if(strcmp("co2", mol_names[mtype]) == 0) na = 3;
	
	tagint anchor = tag_full[i];
	
	for(int j=0; j<na; ++j) {
	  
	  const int ii = atom_map_full[anchor + j];
	  f_full[ii][0] = -grads[indx++];
	  f_full[ii][1] = -grads[indx++];
	  f_full[ii][2] = -grads[indx++];

	  //	  printf("MASTER:: tag= %i  f= %f %f %f\n",tag_full[ii],f_full[ii][0],f_full[ii][1],f_full[ii][2]);
	}
	
      } // if(anchor)
      
    }
    
  } // if(me == 0)

  // scatter forces to other ranks

  const int nlocal = atom->nlocal;
  double ** f_local = fix_mbx->f_local;
  
  MPI_Scatterv(&(f_full[0][0]), fix_mbx->nlocal_rank3, fix_mbx->nlocal_disp3, MPI_DOUBLE, &(f_local[0][0]), nlocal*3, MPI_DOUBLE, 0, world);
  
  // all ranks accumulate forces into their local arrays
  
  // for(int i=0; i<nlocal; ++i)
  //   printf("(%i):: tag= %i  f= %f %f %f\n",comm->me,atom->tag[i],f_local[i][0],f_local[i][1],f_local[i][2]);

  double ** f = atom->f;
  for(int i=0; i<nlocal; ++i) {
    f[i][0] += f_local[i][0];
    f[i][1] += f_local[i][1];
    f[i][2] += f_local[i][2];
  }
  
  //  printf("[MBX] Leaving accumulate_f_full()\n");
  
  fix_mbx->mbxt_stop(MBXT_ACCUMULATE_F_FULL);
}
