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

#include "angle_cosine_break.h"

#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosineBreak::AngleCosineBreak(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCosineBreak::~AngleCosineBreak()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineBreak::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type,m; //GET
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,a,a11,a12,a22;

  eangle = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    
    if(anglelist[n][3] == 0) continue; //GET
    
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // if either bond breaks, set angle type to 0                               //GET
    //   both in temporary anglelist and permanent angle_type                   //GET
    // if this proc owns all three atoms,                                       //GET
    //   negate angle_type twice if other atoms store it                        //GET
    // if other procs owns 1st and 3rd atom, other procs will also break angle  //GET

    if ((r1 > 1.5) || (r2 > 1.5)) {
       anglelist[n][3]=0;
       for (m = 0; m < atom->num_angle[i2]; m++)
        if ( (atom->angle_atom1[i2][m] == atom->tag[i1]) && (atom->angle_atom3[i2][m] == atom->tag[i3]) )
          atom->angle_type[i2][m] = 0;
       if (i1 < atom->nlocal)
        for (m = 0; m < atom->num_angle[i1]; m++)
          if( (atom->angle_atom2[i1][m] == atom->tag[i2]) && (atom->angle_atom3[i1][m] == atom->tag[i3]))
            atom->angle_type[i1][m] = 0;
       if (i3 < atom->nlocal)
        for (m = 0; m < atom->num_angle[i3]; m++)
          if( (atom->angle_atom1[i3][m] == atom->tag[i1]) && (atom->angle_atom2[i3][m] == atom->tag[i2]))
            atom->angle_type[i3][m] = 0;
       continue;
    }

    // c = cosine of angle

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // force & energy

    if (eflag) eangle = k[type]*(1.0+c);

    a = k[type];
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;
        
    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
			 delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineBreak::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void AngleCosineBreak::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  double k_one = utils::numeric(FLERR,arg[1],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleCosineBreak::equilibrium_angle(int /*i*/)
{
  return MY_PI;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCosineBreak::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCosineBreak::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) utils::sfread(FLERR,&k[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleCosineBreak::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g\n",i,k[i]);
}

/* ---------------------------------------------------------------------- */

double AngleCosineBreak::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  return k[type]*(1.0+c);
}
