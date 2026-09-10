/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Repulsive (WCA-like) Gay-Berne potential based on pair_gayberne.

   Contributing author for original GB: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "pair_gayberne_wca.h"

#include <mpi.h>
#include <cmath>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

// precompute 2^(1/6) once
namespace {
constexpr double TWO_1_6 = 1.1224620483093729814;
}

static const char cite_pair_gayberne_wca[] =
  "pair gayberne/wca command:\n\n"
  "@Article{Brown09,\n"
  " author =  {W. M. Brown, M. K. Petersen, S. J. Plimpton, and G. S. Grest},\n"
  " title =   {Liquid crystal nanodroplets in solution},\n"
  " journal = {J.~Chem.~Phys.},\n"
  " year =    2009,\n"
  " volume =  130,\n"
  " pages =   {044901}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairGayBerneWCA::PairGayBerneWCA(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_gayberne_wca);

  single_enable = 0;
  writedata = 1;
  orientation = nullptr;
  orientation_nmax = 0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairGayBerneWCA::~PairGayBerneWCA()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(form);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(shape1);
    memory->destroy(shape2);
    memory->destroy(well);
    memory->destroy(cut);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(chi_iso);
    delete [] lshape;
    delete [] setwell;
  }
  memory->destroy(orientation);
}

/* ---------------------------------------------------------------------- */

void PairGayBerneWCA::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  double a[3][3],temp[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // An atom can occur as j in many neighbor pairs.  Cache its orientation-
  // dependent matrices once per timestep instead of rebuilding them for
  // every occurrence in the inner loop.

  if (atom->nmax > orientation_nmax) {
    orientation_nmax = atom->nmax;
    memory->grow(orientation,orientation_nmax,"pair:orientation");
  }
  const int nall = nlocal + atom->nghost;
  for (i = 0; i < nall; i++) {
    const int ktype = type[i];
    if (form[ktype][ktype] != ELLIPSE_ELLIPSE) {
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++) {
          orientation[i].g[k][l] = k == l ? shape2[ktype][0] : 0.0;
          orientation[i].b[k][l] = k == l ? well[ktype][0] : 0.0;
        }
      orientation[i].planar = true;
      continue;
    }
    double *quat = bonus[ellipsoid[i]].quat;
    orientation[i].planar = quat[1] == 0.0 && quat[2] == 0.0;
    if (orientation[i].planar) {
      // Exact planar rotation, without constructing/multiplying 3x3 matrices.
      const double w2 = quat[0]*quat[0], z2 = quat[3]*quat[3];
      const double c = w2-z2, s = 2.0*quat[0]*quat[3];
      const double cc = c*c, ss = s*s, cs = c*s;
      const double zz = (w2+z2)*(w2+z2);
      double (*g)[3] = orientation[i].g;
      double (*b)[3] = orientation[i].b;
      g[0][0] = shape2[ktype][0]*cc+shape2[ktype][1]*ss;
      g[1][1] = shape2[ktype][0]*ss+shape2[ktype][1]*cc;
      g[0][1] = g[1][0] = (shape2[ktype][0]-shape2[ktype][1])*cs;
      g[2][2] = shape2[ktype][2]*zz;
      b[0][0] = well[ktype][0]*cc+well[ktype][1]*ss;
      b[1][1] = well[ktype][0]*ss+well[ktype][1]*cc;
      b[0][1] = b[1][0] = (well[ktype][0]-well[ktype][1])*cs;
      b[2][2] = well[ktype][2]*zz;
      g[0][2] = g[2][0] = g[1][2] = g[2][1] = 0.0;
      b[0][2] = b[2][0] = b[1][2] = b[2][1] = 0.0;
      continue;
    }
    MathExtra::quat_to_mat_trans(quat,a);
    MathExtra::diag_times3(well[ktype],a,temp);
    MathExtra::transpose_times3(a,temp,orientation[i].b);
    MathExtra::diag_times3(shape2[ktype],a,temp);
    MathExtra::transpose_times3(a,temp,orientation[i].g);
  }

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    OrientationMatrices &oi = orientation[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      if (factor_lj == 0.0) continue;

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {
        if (epsilon[itype][jtype] == 0.0) continue;

        switch (form[itype][jtype]) {
        case SPHERE_SPHERE:
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= -r2inv;
          if (eflag) one_eng =
                       r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) +
                       epsilon[itype][jtype];
          fforce[0] = r12[0]*forcelj;
          fforce[1] = r12[1]*forcelj;
          fforce[2] = r12[2]*forcelj;
          ttor[0] = ttor[1] = ttor[2] = 0.0;
          rtor[0] = rtor[1] = rtor[2] = 0.0;
          break;

        default:
          // The contact surface is the ellipsoid r^T G^-1 r = 2.
          // Its axis-aligned bounds are sqrt(2*G_kk).  Adding a positive
          // radial WCA offset expands each bound by at most that offset.
          // These tests are conservative for arbitrary particle rotations.
          {
            const double pad = MAX(0.0,(TWO_1_6-gamma)*sigma[itype][jtype]);
            bool outside = false;
            for (int k = 0; k < 3; k++) {
              const double distance = std::abs(r12[k])-pad;
              if (distance > 0.0 && distance*distance >
                  2.0*(oi.g[k][k]+orientation[j].g[k][k])) {
                outside = true;
                break;
              }
            }
            if (outside) continue;
          }
          if (oi.planar && orientation[j].planar && r12[2] == 0.0)
            one_eng = gayberne_fast<true>(itype,jtype,oi,orientation[j],
                                          r12,rsq,fforce,ttor,rtor);
          else
            one_eng = gayberne_fast<false>(itype,jtype,oi,orientation[j],
                                           r12,rsq,fforce,ttor,rtor);
          break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
        tor[i][1] += ttor[1];
        tor[i][2] += ttor[2];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
          tor[j][1] += rtor[1];
          tor[j][2] += rtor[2];
        }

        if (eflag) evdwl = factor_lj*one_eng;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                 -r12[0],-r12[1],-r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGayBerneWCA::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(form,n+1,n+1,"pair:form");
  // Diagonal entries may be unused in pair_style hybrid.
  for (int i = 1; i <= n; i++) form[i][i] = SPHERE_SPHERE;
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(shape1,n+1,3,"pair:shape1");
  memory->create(shape2,n+1,3,"pair:shape2");
  memory->create(well,n+1,3,"pair:well");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(chi_iso,n+1,n+1,"pair:chi_iso");
  for (int i = 1; i <= n; i++)
    well[i][0] = well[i][1] = well[i][2] = 1.0;
  lshape = new double[n+1];
  setwell = new int[n+1];
  for (int i = 1; i <= n; i++) setwell[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGayBerneWCA::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  gamma      = utils::numeric(FLERR,arg[0],false,lmp);
  upsilon    = utils::numeric(FLERR,arg[1],false,lmp)/2.0;
  mu         = utils::numeric(FLERR,arg[2],false,lmp);
  cut_global = utils::numeric(FLERR,arg[3],false,lmp);

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

void PairGayBerneWCA::coeff(int narg, char **arg)
{
  if (narg < 10 || narg > 11)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one   = utils::numeric(FLERR,arg[3],false,lmp);
  double eia_one     = utils::numeric(FLERR,arg[4],false,lmp);
  double eib_one     = utils::numeric(FLERR,arg[5],false,lmp);
  double eic_one     = utils::numeric(FLERR,arg[6],false,lmp);
  double eja_one     = utils::numeric(FLERR,arg[7],false,lmp);
  double ejb_one     = utils::numeric(FLERR,arg[8],false,lmp);
  double ejc_one     = utils::numeric(FLERR,arg[9],false,lmp);

  double cut_one = cut_global;
  if (narg == 11) cut_one = utils::numeric(FLERR,arg[10],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      if (eia_one != 0.0 || eib_one != 0.0 || eic_one != 0.0) {
        well[i][0] = std::pow(eia_one,-1.0/mu);
        well[i][1] = std::pow(eib_one,-1.0/mu);
        well[i][2] = std::pow(eic_one,-1.0/mu);
        if (eia_one == eib_one && eib_one == eic_one) setwell[i] = 2;
        else setwell[i] = 1;
      }
      if (eja_one != 0.0 || ejb_one != 0.0 || ejc_one != 0.0) {
        well[j][0] = std::pow(eja_one,-1.0/mu);
        well[j][1] = std::pow(ejb_one,-1.0/mu);
        well[j][2] = std::pow(ejc_one,-1.0/mu);
        if (eja_one == ejb_one && ejb_one == ejc_one) setwell[j] = 2;
        else setwell[j] = 1;
      }
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGayBerneWCA::init_style()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Pair gayberne/wca requires atom style ellipsoid");

  neighbor->request(this,instance_me);

  // per-type shape precalculations
  // require that atom shapes are identical within each type
  // if shape = 0 for point particle, set shape = 1 as required by Gay-Berne

  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->shape_consistency(i,shape1[i][0],shape1[i][1],shape1[i][2]))
      error->all(FLERR,
                 "Pair gayberne/wca requires atoms with same type have same shape");
    if (shape1[i][0] == 0.0)
      shape1[i][0] = shape1[i][1] = shape1[i][2] = 1.0;
    shape2[i][0] = shape1[i][0]*shape1[i][0];
    shape2[i][1] = shape1[i][1]*shape1[i][1];
    shape2[i][2] = shape1[i][2]*shape1[i][2];
    lshape[i] = (shape1[i][0]*shape1[i][1]+shape1[i][2]*shape1[i][2]) *
      std::sqrt(shape1[i][0]*shape1[i][1]);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
------------------------------------------------------------------------- */

double PairGayBerneWCA::init_one(int i, int j)
{
  if (setwell[i] == 0 || setwell[j] == 0)
    error->all(FLERR,"Pair gayberne/wca epsilon a,b,c coeffs are not all set");

  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  const double sij   = sigma[i][j];
  const double epsij = epsilon[i][j];

  const double sij2  = sij*sij;
  const double sij4  = sij2*sij2;
  const double sij6  = sij4*sij2;
  const double sij12 = sij6*sij6;

  lj1[i][j] = 48.0 * epsij * sij12;
  lj2[i][j] = 24.0 * epsij * sij6;
  lj3[i][j] =  4.0 * epsij * sij12;
  lj4[i][j] =  4.0 * epsij * sij6;

  int ishape = 0;
  if (shape1[i][0] != shape1[i][1] ||
      shape1[i][0] != shape1[i][2] ||
      shape1[i][1] != shape1[i][2]) ishape = 1;
  if (setwell[i] == 1) ishape = 1;
  int jshape = 0;
  if (shape1[j][0] != shape1[j][1] ||
      shape1[j][0] != shape1[j][2] ||
      shape1[j][1] != shape1[j][2]) jshape = 1;
  if (setwell[j] == 1) jshape = 1;

  if (ishape == 0 && jshape == 0)
    form[i][i] = form[j][j] = form[i][j] = form[j][i] = SPHERE_SPHERE;
  else if (ishape == 0) {
    form[i][i] = SPHERE_SPHERE; form[j][j] = ELLIPSE_ELLIPSE;
    form[i][j] = SPHERE_ELLIPSE;  form[j][i] = ELLIPSE_SPHERE;
  } else if (jshape == 0) {
    form[j][j] = SPHERE_SPHERE; form[i][i] = ELLIPSE_ELLIPSE;
    form[j][i] = SPHERE_ELLIPSE;  form[i][j] = ELLIPSE_SPHERE;
  } else
    form[i][i] = form[j][j] = form[i][j] = form[j][i] = ELLIPSE_ELLIPSE;

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i]   = sigma[i][j];
  lj1[j][i]     = lj1[i][j];
  lj2[j][i]     = lj2[i][j];
  lj3[j][i]     = lj3[i][j];
  lj4[j][i]     = lj4[i][j];
  chi_iso[i][j] = chi_iso[j][i] = (setwell[i] == 2 && setwell[j] == 2) ?
    std::pow(2.0/(well[i][0]+well[j][0]),mu) : -1.0;
  // Do not put pairs in the neighbor list beyond the largest possible WCA
  // range.  The user cutoff can still impose a shorter range explicitly.

  double cut_wca;
  if (form[i][j] == SPHERE_SPHERE) {
    cut_wca = TWO_1_6*sigma[i][j];
  } else {
    const double imax = MAX(shape1[i][0],MAX(shape1[i][1],shape1[i][2]));
    const double jmax = MAX(shape1[j][0],MAX(shape1[j][1],shape1[j][2]));
    cut_wca = std::sqrt(2.0*(imax*imax+jmax*jmax)) +
      (TWO_1_6-gamma)*sigma[i][j];
  }

  return MIN(cut[i][j],MAX(0.0,cut_wca));
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGayBerneWCA::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    fwrite(&setwell[i],sizeof(int),1,fp);
    if (setwell[i]) fwrite(&well[i][0],sizeof(double),3,fp);
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGayBerneWCA::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) utils::sfread(FLERR,&setwell[i],sizeof(int),1,fp,NULL,error);
    MPI_Bcast(&setwell[i],1,MPI_INT,0,world);
    if (setwell[i]) {
      if (me == 0) utils::sfread(FLERR,&well[i][0],sizeof(double),3,fp,NULL,error);
      MPI_Bcast(&well[i][0],3,MPI_DOUBLE,0,world);
    }
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGayBerneWCA::write_restart_settings(FILE *fp)
{
  fwrite(&gamma,sizeof(double),1,fp);
  fwrite(&upsilon,sizeof(double),1,fp);
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGayBerneWCA::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&gamma,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&upsilon,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&mu,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&upsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairGayBerneWCA::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g\n",i,
            epsilon[i][i],sigma[i][i],
            std::pow(well[i][0],-mu),std::pow(well[i][1],-mu),std::pow(well[i][2],-mu),
            std::pow(well[i][0],-mu),std::pow(well[i][1],-mu),std::pow(well[i][2],-mu));
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairGayBerneWCA::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g\n",i,j,
              epsilon[i][j],sigma[i][j],
              std::pow(well[i][0],-mu),std::pow(well[i][1],-mu),std::pow(well[i][2],-mu),
              std::pow(well[j][0],-mu),std::pow(well[j][1],-mu),std::pow(well[j][2],-mu),
              cut[i][j]);
}

/* ----------------------------------------------------------------------
   Symmetric positive definite matrix operations.  LDL^T avoids the general
   pivoting solver and supplies the inverse and determinant together.
------------------------------------------------------------------------- */

namespace {
struct Symmetric3 {
  double xx, yy, zz, xy, xz, yz;

  static Symmetric3 sum(const double a[3][3], const double b[3][3])
  {
    return {a[0][0]+b[0][0],a[1][1]+b[1][1],a[2][2]+b[2][2],
            a[0][1]+b[0][1],a[0][2]+b[0][2],a[1][2]+b[1][2]};
  }

  template<bool PLANAR>
  double invert(Symmetric3 &v) const
  {
    const double ix = 1.0/xx;
    const double l10 = xy*ix;
    const double d1 = yy-xy*l10;
    const double iy = 1.0/d1;
    if (PLANAR) {
      v = {ix+l10*l10*iy,iy,1.0/zz,-l10*iy,0.0,0.0};
      return xx*d1*zz;
    }
    const double l20 = xz*ix;
    const double l21 = (yz-xz*l10)*iy;
    const double d2 = zz-xz*l20-d1*l21*l21;
    const double iz = 1.0/d2;
    const double t = l10*l21-l20;
    v = {ix+l10*l10*iy+t*t*iz,iy+l21*l21*iz,iz,
         -l10*iy-t*l21*iz,t*iz,-l21*iz};
    return xx*d1*d2;
  }

  template<bool PLANAR>
  void multiply(const double *x, double *y) const
  {
    y[0] = xx*x[0]+xy*x[1];
    y[1] = xy*x[0]+yy*x[1];
    y[2] = 0.0;
    if (!PLANAR) {
      y[0] += xz*x[2];
      y[1] += yz*x[2];
      y[2] = xz*x[0]+yz*x[1]+zz*x[2];
    }
  }
};

// Exact reductions for the common exponents, no approximate math/fast-math.
inline double gb_power(double x, double exponent)
{
  if (exponent == 0.0) return 1.0;
  if (exponent == 1.0) return x;
  if (exponent == 2.0) return x*x;
  if (exponent == 0.5) return std::sqrt(x);
  if (exponent == 1.5) return x*std::sqrt(x);
  return std::pow(x,exponent);
}
}

/* ----------------------------------------------------------------------
   Same GB/WCA energy and derivatives as the original analytic kernel.
   PLANAR uses the exact xy block when both orientations and r lie in xy.
   It retains the zz determinant factor (this is not a different 2-D model).
------------------------------------------------------------------------- */

template<bool PLANAR>
double PairGayBerneWCA::gayberne_fast(int itype, int jtype,
                                    const OrientationMatrices &oi,
                                    const OrientationMatrices &oj,
                                    const double *r12, double rsq,
                                    double *fforce, double *ttor, double *rtor)
{
  const double r = std::sqrt(rsq);
  const double rinv = 1.0/r;
  const double n[3] = {r12[0]*rinv,r12[1]*rinv,r12[2]*rinv};
  const Symmetric3 g = Symmetric3::sum(oi.g,oj.g);
  Symmetric3 invg;
  const double detg = g.invert<PLANAR>(invg);
  if (!(detg > 0.0))
    error->one(FLERR,"Bad shape matrix in pair gayberne/wca");

  double k[3];
  invg.multiply<PLANAR>(n,k);
  const double nk = MathExtra::dot3(n,k);
  const double contact = std::sqrt(2.0/nk);
  const double sigma0 = sigma[itype][jtype];
  if (r >= contact+(TWO_1_6-gamma)*sigma0) {
    fforce[0] = fforce[1] = fforce[2] = 0.0;
    ttor[0] = ttor[1] = ttor[2] = 0.0;
    rtor[0] = rtor[1] = rtor[2] = 0.0;
    return 0.0;
  }

  const double rho = sigma0/(r-contact+gamma*sigma0);
  const double rho2 = rho*rho;
  const double rho6 = rho2*rho2*rho2;
  const double eps = epsilon[itype][jtype];
  // Algebraically 4*eps*(rho^12-rho^6)+eps, stable near the minimum.
  const double ur = 4.0*eps*(rho6-0.5)*(rho6-0.5);
  const double derivative = 24.0*eps*rho*rho6*(2.0*rho6-1.0)/sigma0;
  const double eta = gb_power(2.0*lshape[itype]*lshape[jtype]/detg,upsilon);
  double chi = chi_iso[itype][jtype];
  double v[3] = {0.0,0.0,0.0};
  double nv = 0.0, chi_derivative = 0.0;
  if (chi < 0.0) {
    const Symmetric3 b = Symmetric3::sum(oi.b,oj.b);
    Symmetric3 invb;
    if (!(b.invert<PLANAR>(invb) > 0.0))
      error->one(FLERR,"Bad well matrix in pair gayberne/wca");
    invb.multiply<PLANAR>(n,v);
    nv = MathExtra::dot3(n,v);
    chi = gb_power(2.0*nv,mu);
    chi_derivative = 2.0*mu/nv;
  }
  const double energy = eta*chi*ur;
  const double radial = eta*chi*derivative;
  const double shape_derivative = radial*contact*contact*contact*0.5;
  const double angular = energy*chi_derivative;
  constexpr int ndim = PLANAR ? 2 : 3;
  for (int m = 0; m < ndim; m++)
    fforce[m] = -radial*n[m]-shape_derivative*rinv*(k[m]-nk*n[m])+
      angular*rinv*(v[m]-nv*n[m]);

  // d(log det G)/d(theta_i) = 2*axial(G_i G^-1).
  // This replaces the repeated symbolic 3x3 eta-torque derivatives.
  const double eta_derivative = 2.0*energy*upsilon;
  if (PLANAR) {
    const double kgx = oi.g[0][0]*k[0]+oi.g[0][1]*k[1];
    const double kgy = oi.g[0][1]*k[0]+oi.g[1][1]*k[1];
    const double vb_x = oi.b[0][0]*v[0]+oi.b[0][1]*v[1];
    const double vb_y = oi.b[0][1]*v[0]+oi.b[1][1]*v[1];
    const double axial = (oi.g[0][0]-oi.g[1][1])*invg.xy+
      oi.g[0][1]*(invg.yy-invg.xx);
    ttor[2] = shape_derivative*(k[0]*kgy-k[1]*kgx)+
      angular*(vb_x*v[1]-vb_y*v[0])+eta_derivative*axial;
    fforce[2] = ttor[0] = ttor[1] = rtor[0] = rtor[1] = 0.0;
    rtor[2] = r12[0]*fforce[1]-r12[1]*fforce[0]-ttor[2];
  } else {
    double gk[3], bv[3], cross_g[3], cross_b[3];
    MathExtra::matvec(oi.g,k,gk);
    MathExtra::matvec(oi.b,v,bv);
    MathExtra::cross3(k,gk,cross_g);
    MathExtra::cross3(bv,v,cross_b);
    const double axial[3] = {
      oi.g[0][1]*invg.xz-oi.g[0][2]*invg.xy+
        (oi.g[1][1]-oi.g[2][2])*invg.yz+oi.g[1][2]*(invg.zz-invg.yy),
      oi.g[0][2]*(invg.xx-invg.zz)+(oi.g[2][2]-oi.g[0][0])*invg.xz+
        oi.g[1][2]*invg.xy-oi.g[0][1]*invg.yz,
      (oi.g[0][0]-oi.g[1][1])*invg.xy+oi.g[0][1]*(invg.yy-invg.xx)+
        oi.g[0][2]*invg.yz-oi.g[1][2]*invg.xz};
    for (int m = 0; m < 3; m++)
      ttor[m] = shape_derivative*cross_g[m]+angular*cross_b[m]+eta_derivative*axial[m];
    // Rotational invariance: tau_i + tau_j = r_ij cross F_i.
    // This also applies with Newton off; unused ghost torques are not tallied.
    MathExtra::cross3(r12,fforce,rtor);
    for (int m = 0; m < 3; m++) rtor[m] -= ttor[m];
  }
  // An isotropic sphere has no torque, including floating-point residue.
  if (form[itype][jtype] == SPHERE_ELLIPSE)
    ttor[0] = ttor[1] = ttor[2] = 0.0;
  if (form[itype][jtype] == ELLIPSE_SPHERE)
    rtor[0] = rtor[1] = rtor[2] = 0.0;
  return energy;
}

/* ---------------------------------------------------------------------- */

double PairGayBerneWCA::memory_usage()
{
  return (double) orientation_nmax*sizeof(OrientationMatrices);
}
