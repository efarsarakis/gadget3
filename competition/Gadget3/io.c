#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "utilities.h"
#include "proto.h"

#ifdef LT_METAL_COOLING_WAL
#include "lt_sfr/lt_wal_cooling.h"
#endif

#ifdef VORONOI
#include "voronoi.h"
#endif

#ifdef JD_DPP
#include "cr_electrons.h"
#endif

/*! \file io.c
 *  \brief Output of a snapshot file to disk.
 */

static int n_type[6];
static long long ntot_type_all[6];

#ifdef LT_SEvDbg
double metals[2 * LT_NMetP + 3], tot_metals[2 * LT_NMetP + 3];
#endif

#ifdef LT_SEv_INFO_DETAILS
double DetailsWSum[LT_NMetP], DetailsWoSum[LT_NMetP];
double my_weight_sum_details, My_weight_sum_details;
unsigned int *my_weight_sum_details_N;
#endif


#if (defined(LT_STELLAREVOLUTION) || defined(LT_SMOOTH_Z) || defined(LT_METAL_COOLING) || defined(LT_METAL_COOLING_WAL))
double XH, yhelium;

#if defined(LT_LOG_SMOOTH)
double metalmass;
void Spy_Z(int);
#endif
#endif

#ifdef LT_METAL_COOLING_WAL
double Redshift, DZ;
#endif


static int n_info;

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savepositions(int num)
{
  size_t bytes;
  char buf[500];
  int n, filenr, gr, ngroups, masterTask, lastTask;

  CPU_Step[CPU_MISC] += measure_time();

#if defined(SFR) || defined(BLACK_HOLES)
  rearrange_particle_sequence();
  /* ensures that new tree will be constructed */
  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TreeDomainUpdateFrequency * All.TotNumPart);
#endif

#if defined(JD_DPP) && defined(JD_DPPONSNAPSHOTONLY)
  compute_Dpp(0);		/* Particle loop already inside, just add water */
#endif

  if(DumpFlag == 1)
    {
      if(ThisTask == 0)
	printf("\nwriting snapshot file #%d... \n", num);

#ifdef ORDER_SNAPSHOTS_BY_ID
      double t0, t1;

      if(All.TotN_gas > 0)
	{
	  if(ThisTask == 0)
	    printf
	      ("\nThe option ORDER_SNAPSHOTS_BY_ID does not work yet with gas particles, only simulations with collisionless particles are allowed.\n\n");
	  endrun(0);
	}

      t0 = second();
      for(n = 0; n < NumPart; n++)
	{
	  P[n].GrNr = ThisTask;
	  P[n].SubNr = n;
	}
      parallel_sort(P, NumPart, sizeof(struct particle_data), io_compare_P_ID);
      t1 = second();
      if(ThisTask == 0)
	printf("Reordering of particle-data in ID-sequence took = %g sec\n", timediff(t0, t1));
#endif

#ifdef VORONOI_MESHOUTPUT
#ifdef ORDER_SNAPSHOTS_BY_ID
      endrun(112);
#endif

      domain_Decomposition(0, 0);

      force_treebuild(NumPart, NULL);
      voronoi_mesh();
      voronoi_setup_exchange();
#endif

#ifdef WRITE_KEY_FILES
      for(n = 0; n < 6; n++)
	NKeys[n] = NKeysTotal[n] = 0;

      domain_Decomposition(0, 1);
#endif

      if(!(CommBuffer = mymalloc("CommBuffer", bytes = All.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}


      if(NTask < All.NumFilesPerSnapshot)
	{
	  if(ThisTask == 0)
	    printf
	      ("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
	  endrun(0);
	}
      if(All.SnapFormat < 1 || All.SnapFormat > 3)
	{
	  if(ThisTask == 0)
	    printf("Unsupported File-Format\n");
	  endrun(0);
	}
#ifndef  HAVE_HDF5
      if(All.SnapFormat == 3)
	{
	  if(ThisTask == 0)
	    printf("Code wasn't compiled with HDF5 support enabled!\n");
	  endrun(0);
	}
#endif


      /* determine global and local particle numbers */
      for(n = 0; n < 6; n++)
	n_type[n] = 0;

      for(n = 0; n < NumPart; n++)
	n_type[P[n].Type]++;

      sumup_large_ints(6, n_type, ntot_type_all);

#ifdef WRITE_KEY_FILES
      sumup_large_ints(6, NKeys, NKeysTotal);
#endif

#ifdef LT_SEvDbg
      int i;

      for(i = 0; i < 2 * LT_NMetP + 3; i++)
	metals[i] = tot_metals[i] = 0;
#endif

#ifdef LT_SEv_INFO_DETAILS
      MPI_Reduce(DetailsW, DetailsWSum, LT_NMetP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(DetailsWo, DetailsWoSum, LT_NMetP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      int i;
      if(ThisTask == 0)
	for(i = 0; i < LT_NMetP; i++)
	  printf("%g DETAILS: %4s :: %12.8g  *  %12.8g\n", All.Time, MetNames[i], DetailsWSum[i],
		 DetailsWoSum[i]);
#endif

#ifdef LT_METAL_COOLING_WAL
      Redshift = 1.0 / All.Time - 1;	/* note that ComovingIntegration is alway on with wal cooling */
      get_cool_redshift(Redshift, &DZ);
#endif

      /* assign processors to output files */
      distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(All.NumFilesPerSnapshot > 1)
	{
	  if(ThisTask == 0)
	    {
	      sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
	      mkdir(buf, 02755);
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      if(All.NumFilesPerSnapshot > 1)
	sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
      else
	sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);

#if (defined(LT_STELLAREVOLUTION) || defined(LT_SMOOTH_Z) || defined(LT_METAL_COOLING) || defined(LT_METAL_COOLING_WAL))
      XH = HYDROGEN_MASSFRAC;
      yhelium = (1 - XH) / (4 * XH);

#if defined(LT_LOG_SMOOTH)
      if(ThisTask == 0)
	printf("logging z-smooth effect..\n");
      fflush(stdout);
      Spy_Z(num);
#endif
#endif

      ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
      if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	    {
	      write_file(buf, masterTask, lastTask);

#ifdef VORONOI_MESHOUTPUT
	      if(All.NumFilesPerSnapshot > 1)
		sprintf(buf, "%s/snapdir_%03d/voronoi_mesh_%03d.%d", All.OutputDir, num, num, filenr);
	      else
		sprintf(buf, "%s/voronoi_mesh_%03d", All.OutputDir, num);
	      write_voronoi_mesh(buf, masterTask, lastTask);
#endif
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      myfree(CommBuffer);

#ifdef VORONOI_MESHOUTPUT
      myfree(List_P);
      myfree(ListExports);
      myfree(DT);
      myfree(DP - 5);
      myfree(VF);		/* free the list of faces */
#endif

#ifdef ORDER_SNAPSHOTS_BY_ID
      t0 = second();
      parallel_sort(P, NumPart, sizeof(struct particle_data), io_compare_P_GrNr_SubNr);
      t1 = second();
      if(ThisTask == 0)
	printf("Restoring order of particle-data took = %g sec\n", timediff(t0, t1));
#endif

#ifdef JD_RELAX_CLUSTERS
      if(ThisTask == 0)
	printf("Setting vel=0 \n");
      size_t ipart;
      for(ipart = 0; ipart < N_gas; ipart++)
	P[0].Vel[0] = P[0].Vel[1] = P[0].Vel[2] = 0;
#endif

      if(ThisTask == 0)
	printf("done with snapshot.\n");

      All.Ti_lastoutput = All.Ti_Current;

      CPU_Step[CPU_SNAPSHOT] += measure_time();

#ifdef SUBFIND_RESHUFFLE_CATALOGUE
      endrun(0);
#endif
    }

#ifdef WRITE_KEY_FILES
  /*
     domain_Decomposition(0, 0);
   */
#endif

#ifdef FOF
  if(RestartFlag != 4)
    {
      if(ThisTask == 0)
	printf("\ncomputing group catalogue...\n");

      fof_fof(num);

      if(ThisTask == 0)
	printf("done with group catalogue.\n");

      CPU_Step[CPU_FOF] += measure_time();
    }
#endif

#ifdef POWERSPEC_ON_OUTPUT
  if(RestartFlag != 4)
    {
      if(ThisTask == 0)
	printf("\ncomputing power spectra...\n");

      calculate_power_spectra(num, &ntot_type_all[0]);

      if(ThisTask == 0)
	printf("done with power spectra.\n");

      CPU_Step[CPU_MISC] += measure_time();
    }
#endif
#ifdef LT_SEvDbg
  MPI_Reduce(&metals[0], &tot_metals[0], 2 * LT_NMetP + 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      fprintf(FdMetSumCheck, "@.@ %g ", All.Time);
      int i;

      for(i = 0; i < LT_NMetP * 2 + 3; i++)
	fprintf(FdMetSumCheck, " %g ", tot_metals[i] / All.HubbleParam);
      fprintf(FdMetSumCheck, " \n ");
    }
#endif

}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex;
  MyOutputFloat *fp;
  MyIDType *ip;
  int *ip_int;
  float *fp_single;
  integertime dt_step;

#ifdef MAGNETIC
  double a2_inv, gadget2gauss = 1.0;
#endif

#ifdef PERIODIC
  MyFloat boxSize;
#endif
#ifdef PMGRID
  double dt_gravkick_pm = 0;
#endif
  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

#ifdef OUTPUTCOOLRATE
  double tcool, u;
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#ifdef BP_REAL_CRs
  int Nbin;
#endif

#if defined(DISTORTIONTENSORPS) && defined(GDE_LOGOUTPUT)
  MyBigFloat tmpval1, tmpval2, zero;
  zero = 0.0;
#endif
#if (defined(OUTPUT_DISTORTIONTENSORPS) || defined(OUTPUT_TIDALTENSORPS))
  MyBigFloat half_kick_add[6][6];
  int l;
#endif

#ifdef LT_STELLAREVOLUTION
  int control[LT_NMetP];
  float Temperature, xclouds;

#ifdef LT_TRACK_CONTRIBUTES
  Contrib *contrib;
#endif
#endif

#ifdef LT_MV_CHEMICALDIFFUSION
  int el, cc;
  double maxgrad, myD, tscale, myAvgA;
#endif


  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
#ifdef MAGNETIC
#ifndef MU0_UNITY
      gadget2gauss /= sqrt(All.UnitTime_in_s * All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g /
			   (All.HubbleParam * All.HubbleParam));
#endif
      a2_inv = 1.0 / (All.Time * All.Time);
#endif
    }
  else
    {
      a3inv = fac1 = fac2 = 1;
#ifdef MAGNETIC
#ifndef MU0_UNITY
      gadget2gauss /= sqrt(All.UnitTime_in_s * All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g);
#endif
      a2_inv = 1;
#endif
    }

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkick_pm =
      get_gravkick_factor(All.PM_Ti_begstep,
			  All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif

#ifdef DISTORTIONTENSORPS
  MyBigFloat flde, psde;
#endif

  fp = (MyOutputFloat *) CommBuffer;
  fp_single = (float *) CommBuffer;
  ip = (MyIDType *) CommBuffer;
  ip_int = (int *) CommBuffer;
#ifdef LT_TRACK_CONTRIBUTES
  contrib = (Contrib *) CommBuffer;
#endif

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
#ifdef PERIODIC
		boxSize = All.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = All.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#endif
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    dt_step = (P[pindex].TimeBin ? (((integertime) 1) << P[pindex].TimeBin) : 0);

	    if(All.ComovingIntegrationOn)
	      {
		dt_gravkick =
		  get_gravkick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
		  get_gravkick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
		dt_hydrokick =
		  get_hydrokick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
		  get_hydrokick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
	      }
	    else
	      dt_gravkick = dt_hydrokick =
		(All.Ti_Current - (P[pindex].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k] + P[pindex].g.GravAccel[k] * dt_gravkick;
		if(P[pindex].Type == 0)
		  fp[k] += SphP[pindex].a.HydroAccel[k] * dt_hydrokick;
	      }
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += P[pindex].GravPM[k] * dt_gravkick_pm;
#endif
	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv);

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Mass;
	    n++;
	  }
      break;

    case IO_U:			/* internal energy */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifdef VORONOI_MESHRELAX
#ifdef VORONOI_MESHRELAX_KEEPRESSURE
	    SphP[pindex].Entropy = SphP[pindex].Pressure / (GAMMA_MINUS1 * SphP[pindex].d.Density);
#endif
	    *fp++ = SphP[pindex].Entropy;
#else

#if !defined(EOS_DEGENERATE) && !defined(TRADITIONAL_SPH_FORMULATION) && !defined(ISOTHERM_EQS)
	    *fp++ =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1));
#else
	    *fp++ = SphP[pindex].Entropy;
#endif

#endif
	    n++;
	  }
      break;

    case IO_RHO:		/* density */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
      break;

    case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(RADTRANSFER)
	    *fp++ = SphP[pindex].elec;
#else
	    *fp++ = SphP[pindex].Ne;
#endif
	    n++;
	  }
#endif
      break;

    case IO_NH:		/* neutral hydrogen fraction */
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(RADTRANSFER)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(RADTRANSFER)
	    *fp++ = SphP[pindex].HI;
#else
	    double u, ne, nh0 = 0, nHeII = 0;
	    double fac_entr_to_u;

	    fac_entr_to_u = pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
	    u =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy * fac_entr_to_u) * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

#ifndef LT_METAL_COOLING_WAL

	    ne = SphP[pindex].Ne;

	    AbundanceRatios(u, SphP[pindex].d.Density * a3inv, &ne, &nh0, &nHeII);
#else
	    double rho, temp;
	    double Metallicities[LT_NMet];

	    rho = SphP[pindex].d.Density * a3inv;

	    set_metallicities(pindex, &Metallicities[0], a3inv);

	    temp = convert_u_to_temp(u, Redshift, DZ, &Metallicities[0]);
	    nh0 = AbundanceRatios(temp, rho, &ne, &nh0, &nHeII);
#endif

	    *fp++ = (MyOutputFloat) nh0;
#endif
	    n++;
	  }
#endif
      break;

    case IO_HII:		/* ionized hydrogen abundance */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)  || defined(RADTRANSFER)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HII;
	    n++;
	  }
#endif
      break;

    case IO_HeI:		/* neutral Helium */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || (defined(RADTRANSFER) && defined(RT_INCLUDE_HE))
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeI;
	    n++;
	  }
#endif
      break;

    case IO_HeII:		/* ionized Heluum */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || (defined(RADTRANSFER) && defined(RT_INCLUDE_HE))
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeII;
	    n++;
	  }
#endif
      break;

    case IO_HeIII:		/* double ionised Helium */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeIII;
	    n++;
	  }
#endif
      break;

    case IO_H2I:		/* H2 molecule */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].H2I;
	    n++;
	  }
#endif
      break;

    case IO_H2II:		/* ionised H2 molecule */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].H2II;
	    n++;
	  }
#endif
      break;

    case IO_HM:		/* H minus */
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HM;
	    n++;
	  }
#endif
      break;

    case IO_HD:		/* HD */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HD;
	    n++;
	  }
#endif
      break;

    case IO_DI:		/* deuterium */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DI;
	    n++;
	  }
#endif
      break;

    case IO_DII:		/* deuteriumII */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DII;
	    n++;
	  }
#endif
      break;

    case IO_HeHII:		/* HeH+ */
#ifdef UM_CHEMISTRY
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeHII;
	    n++;
	  }
#endif
      break;

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = PPP[pindex].Hsml;
	    n++;
	  }
      break;


    case IO_VALPHA:		/* artificial viscosity of particle  */
#ifdef VORONOI_TIME_DEP_ART_VISC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].alpha;
	    n++;
	  }
#endif
      break;

    case IO_SFR:		/* star formation rate */
#ifdef SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifndef LT_STELLAREVOLUTION
	    *fp++ = get_starformation_rate(pindex);
#else
	    *fp++ = get_starformation_rate(pindex, &Temperature, &xclouds);
#endif
	    n++;
	  }
#endif
      break;

    case IO_AGE:		/* stellar formation time */
#ifdef STELLARAGE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].StellarAge;
	    n++;
	  }
#endif
      break;

    case IO_VSTURB_DISS:
#if defined (VS_TURB) || defined (AB_TURB)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DuDt_diss;
	    n++;
	  }
#endif
      break;

    case IO_VSTURB_DRIVE:
#if defined (VS_TURB) || defined (AB_TURB)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DuDt_drive;
	    n++;
	  }
#endif
      break;

    case IO_HSMS:		/* smoothing length for star particles */
#ifdef SUBFIND
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].DM_Hsml;
	    n++;
	  }
#endif
      break;

    case IO_ACRS:		/* accreation length for star particles */
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Hsml;
	    n++;
	  }
#endif
      break;

    case IO_Z:			/* gas and star metallicity */
#ifdef METALS
#ifndef CS_MODEL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Metallicity;
	    n++;
	  }
#else
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(P[pindex].Zm[0] <= 0 || P[pindex].Zm[6] <= 0)
	      {
		printf("NEGATIVE METALLICITY: H=%7.3e, He=%7.3e\n", P[pindex].Zm[6], P[pindex].Zm[0]);
		endrun(3758698);
	      }
	    for(k = 0; k < 12; k++)
	      *fp++ = P[pindex].Zm[k];
	    n++;
	  }
#endif
#endif
      break;

    case IO_POT:		/* gravitational potential */
#if defined(OUTPUTPOTENTIAL)  || defined(SUBFIND_RESHUFFLE_AND_POTENTIAL)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].p.Potential;
	    n++;
	  }
#endif
      break;

    case IO_ACCEL:		/* acceleration */
#ifdef OUTPUTACCELERATION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].g.GravAccel[k];
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += fac1 * P[pindex].GravPM[k];
#endif
	    if(P[pindex].Type == 0)
	      for(k = 0; k < 3; k++)
		fp[k] += fac2 * SphP[pindex].a.HydroAccel[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_DTENTR:		/* rate of change of entropy */
#ifdef OUTPUTCHANGEOFENTROPY
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].e.DtEntropy;
	    n++;
	  }
#endif
      break;

    case IO_STRESSDIAG:	/* Diagonal components of viscous shear tensor */
#ifdef OUTPUTSTRESS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = SphP[pindex].u.s.StressDiag[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_STRESSOFFDIAG:	/* Offdiagonal components of viscous shear tensor */
#ifdef OUTPUTSTRESS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = SphP[pindex].u.s.StressOffDiag[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_STRESSBULK:	/* Viscous bulk tensor */
#ifdef OUTPUTBULKSTRESS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].u.s.StressBulk;
	    n++;
	  }
#endif
      break;

    case IO_DELAYTIME:
#ifdef WINDS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DelayTime;
	    n++;
	  }
#endif
      break;

    case IO_SHEARCOEFF:	/* Shear viscosity coefficient */
#ifdef OUTPUTSHEARCOEFF
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = get_shear_viscosity(pindex) *
	      pow((SphP[pindex].Entropy * pow(SphP[pindex].d.Density * a3inv,
					      GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);
	    n++;
	  }
#endif
      break;

    case IO_TSTP:		/* timestep  */
#ifdef OUTPUTTIMESTEP

      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ =
	      (P[pindex].TimeBin ? (((integertime) 1) << P[pindex].TimeBin) : 0) * All.Timebase_interval;
	    n++;
	  }
#endif
      break;

    case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].b2.BPred[k] * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_BSMTH:		/* smoothed magnetic field */
#ifdef OUTPUTBSMOOTH
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].BSmooth[k] * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_DBDT:		/* rate of change of magnetic field  */
#ifdef DBOUTPUT
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].DtB[k] * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_VRMS:		/* turbulent velocity around mean */
#ifdef JD_VTURB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Vrms * sqrt(a3inv);
	    n++;
	  }
#endif
      break;

    case IO_VBULK:		/* mean velocity in kernel */
#if defined(JD_VTURB) && !defined(FS_TURB_ESTIM)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = SphP[pindex].Vbulk[k] * sqrt(a3inv);
	    n++;
	  }
#endif
      break;

    case IO_VTAN:		/* turbulent velocity around mean, tangential part */
#ifdef JD_DECOMPOSE_VTURB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Vtan * sqrt(a3inv);
	    n++;
	  }
#endif
      break;

    case IO_VRAD:		/* turbulent velocity around mean, radial part */
#ifdef JD_DECOMPOSE_VTURB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Vrad * sqrt(a3inv);
	    n++;
	  }
#endif
      break;


    case IO_TRUENGB:		/* True Number of Neighbours  */
#ifdef JD_VTURB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip_int++ = P[pindex].TrueNGB;
	    n++;
	  }
#endif
      break;

    case IO_DPP:		/* Reacceleration Coefficient */
#ifdef JD_DPP
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Dpp;
	    n++;
	  }
#endif
      break;

    case IO_VDIV:		/* Divergence of Vel */
#if defined(OUTPUT_DIV_CURL)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].v.DivVel;
	    n++;
	  }
#endif
      break;

    case IO_VROT:		/* Velocity Curl */
#if defined(OUTPUT_DIV_CURL)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].r.CurlVel;
	    n++;
	  }
#endif
      break;

    case IO_VORT:		/* Vorticity */
#if defined(OUTPUT_VORTICITY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Vorticity[0];
	    *fp++ = SphP[pindex].Vorticity[1];
	    *fp++ = SphP[pindex].Vorticity[2];
	    n++;
	  }
#endif
      break;

    case IO_DIVB:		/* divergence of magnetic field  */
#ifdef TRACEDIVB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (SphP[pindex].divB * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_ABVC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_ART_VISC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].alpha;
	    n++;
	  }
#endif
      break;


    case IO_AMDC:		/* artificial magnetic dissipation of particle  */
#ifdef TIME_DEP_MAGN_DISP
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Balpha;
	    n++;
	  }
#endif
      break;

    case IO_VTURB:		/* Local Turbulence estimation from SF  */
#ifdef FS_TURB_ESTIM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Vturb;
	    n++;
	  }
#endif
      break;

    case IO_LTURB:		/* Local Coherence lenght estimation from SF  */
#ifdef FS_TURB_ESTIM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Lturb;
	    n++;
	  }
#endif
      break;

    case IO_ALFA2_DYN:		/* alfa^2 Dynamo  */
#ifdef FS_ALFA2_DYN
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].alfa2;
	    n++;
	  }
#endif
      break;

    case IO_ETA2_DYN:		/* Eta in Dynamo  */
#ifdef FS_ETA2_DYN
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].eta2;
	    n++;
	  }
#endif
      break;

    case IO_PHI:		/* divBcleaning fuction of particle  */
#ifdef OUTPUTDEDNER
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (SphP[pindex].PhiPred * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_XPHI:		/* Cold fraction in SF  */
#if (defined(OUTPUT_XPHI) || defined(OUTPUTDEDNER))
#ifdef SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].XColdCloud;
	    n++;
	  }
#endif
#endif
      break;


    case IO_GRADPHI:		/* divBcleaning fuction of particle  */
#ifdef OUTPUTDEDNER
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].GradPhi[k] * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_ROTB:		/* rot of magnetic field  */
#ifdef OUTPUT_ROTB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].RotB[k] * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_SROTB:		/* smoothed rot of magnetic field  */
#ifdef OUTPUT_SROTB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].SmoothedRotB[k] * a2_inv * gadget2gauss);
	    n++;
	  }
#endif
      break;

    case IO_EULERA:		/* magnetic field  */
#ifdef EULERPOTENTIALS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (SphP[pindex].EulerA * sqrt(a2_inv * gadget2gauss));
	    n++;
	  }
#endif
      break;

    case IO_EULERB:		/* magnetic field  */
#ifdef EULERPOTENTIALS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (SphP[pindex].EulerB * sqrt(a2_inv * gadget2gauss));
	    n++;
	  }
#endif
      break;


    case IO_COOLRATE:		/* current cooling rate of particle  */
#ifdef OUTPUTCOOLRATE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifndef UM_CHEMISTRY
	    double ne = SphP[pindex].Ne;
#else
	    double ne = SphP[pindex].elec;
#endif
	    /* get cooling time */
	    u = SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1);

#ifndef LT_METAL_COOLING
	    tcool = GetCoolingTime(u, SphP[pindex].d.Density * a3inv, &ne);
#else

#ifndef LT_METAL_COOLING_on_SMOOTH_Z
	    Zcool = get_metallicity_solarunits(get_metallicity(pindex, Iron));
#else
	    if((metalmass = get_metalmass(SphP[pindex].Metals)) > 0)
	      Zcool =
		get_metallicity_solarunits(SphP[pindex].Zsmooth * SphP[pindex].Metals[Iron] / metalmass);
	    else
	      Zcool = NO_METAL;
#endif // LT_METAL_COOLING_on_SMOOTH_Z

#ifndef UM_CHEMISTRY
#if defined (UM_MET_IN_LT_COOLING)
/* :: This consider the possibility to include low T metal cooling without using non-eq calculations (only in Luca's cooling):: */
	    um_ZsPoint = SphP[pindex].Metals;
	    um_mass = P[pindex].Mass;
	    //      printf("::: set um_ZsPoint\n\n");
	    /*
	       printf(" um_ZsPoint[Iron] %g \n", um_ZsPoint[Iron]);
	       printf(" um_ZsPoint[Oxygen]  %g \n", um_ZsPoint[Oxygen]);
	       printf(" um_ZsPoint[Silicon] %g \n", um_ZsPoint[Silicon]);
	       printf(" um_ZsPoint[Carbon]  %g \n", um_ZsPoint[Carbon]);
	     */
#endif // UM_MET_IN_LT_COOLING
	    tcool = GetCoolingTime(u, SphP[pindex].a2.Density * a3inv, &ne, Zcool);
#else
	    tcool = Um_GetCoolingTime(u, SphP[pindex].a2.Density * a3inv, &ne, Zcool, pindex);
#endif // UM_CHEMISTRY
	    tcool = GetCoolingTime(u, SphP[pindex].a2.Density * a3inv, &ne);
#endif // LT_METAL_COOLING

	    /* convert cooling time with current thermal energy to du/dt */
	    if(tcool != 0)
	      *fp++ = u / tcool;
	    else
	      *fp++ = 0;
	    n++;
	  }
#endif // OUTPUTCOOLRATE
      break;

    case IO_CONDRATE:		/* current heating/cooling due to thermal conduction  */
      break;

    case IO_DENN:		/* density normalization factor */
#ifdef OUTPUTDENSNORM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DensityNorm;
	    n++;
	  }
#endif
      break;

    case IO_EGYPROM:
#ifdef CS_FEEDBACK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].EnergySN;
	    n++;
	  }
#endif
      break;

    case IO_EGYCOLD:
#ifdef CS_FEEDBACK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].EnergySNCold;
	    n++;
	  }
#endif
      break;

    case IO_CR_C0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = SphP[pindex].CR_C0[CRpop];
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_CR_Q0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = SphP[pindex].CR_q0[CRpop];
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_CR_P0:
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_THERMO_VARIABLES)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = CR_Physical_Pressure(&SphP[pindex], CRpop);
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_CR_E0:
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_THERMO_VARIABLES)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = SphP[pindex].CR_E0[CRpop];
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_CR_n0:
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_THERMO_VARIABLES)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = SphP[pindex].CR_n0[CRpop];
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_CR_ThermalizationTime:
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_TIMESCALES)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = CR_Tab_GetThermalizationTimescale(SphP[pindex].CR_q0[CRpop] *
							    pow(SphP[pindex].d.Density * a3inv, 0.333333),
							    SphP[pindex].d.Density * a3inv, CRpop);
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_CR_DissipationTime:
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_TIMESCALES)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      fp[CRpop] = CR_Tab_GetDissipationTimescale(SphP[pindex].CR_q0[CRpop] *
							 pow(SphP[pindex].d.Density * a3inv, 0.333333),
							 SphP[pindex].d.Density * a3inv, CRpop);
	    n++;
	    fp += NUMCRPOP;
	  }
#endif
      break;

    case IO_BPCR_pNORM:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
	      fp[Nbin] = SphP[pindex].CRpNorm[Nbin];
	    n++;
	    fp += BP_REAL_CRs;
	  }
#endif
      break;

    case IO_BPCR_eNORM:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
	      fp[Nbin] = SphP[pindex].CReNorm[Nbin];
	    n++;
	    fp += BP_REAL_CRs;
	  }
#endif
      break;

    case IO_BPCR_pSLOPE:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
	      fp[Nbin] = SphP[pindex].CRpSlope[Nbin];
	    n++;
	    fp += BP_REAL_CRs;
	  }
#endif
      break;

    case IO_BPCR_eSLOPE:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
	      fp[Nbin] = SphP[pindex].CReSlope[Nbin];
	    n++;
	    fp += BP_REAL_CRs;
	  }
#endif
      break;

    case IO_BPCR_pE:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
	      fp[Nbin] = SphP[pindex].CRpE[Nbin];
	    n++;
	    fp += BP_REAL_CRs;
	  }
#endif
      break;

    case IO_BPCR_pN:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
	      fp[Nbin] = SphP[pindex].CRpN[Nbin];
	    n++;
	    fp += BP_REAL_CRs;
	  }
#endif

    case IO_BPCR_pPRESSURE:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CRpPressure;
	    n++;
	  }
#endif
      break;

    case IO_BPCR_ePRESSURE:
#if defined BP_REAL_CRs
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CRePressure;
	    n++;
	  }
#endif
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = BPP(pindex).BH_Mass;
	    n++;
	  }
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = BPP(pindex).BH_Mdot;
	    n++;
	  }
#endif
      break;

    case IO_BHPROGS:
#ifdef BH_COUNTPROGS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip_int++ = BPP(pindex).BH_CountProgs;
	    n++;
	  }
#endif
      break;

    case IO_BHMBUB:
#ifdef BH_BUBBLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mass_bubbles;
	    n++;
	  }
#endif
      break;

    case IO_BHMINI:
#ifdef BH_BUBBLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mass_ini;
	    n++;
	  }
#endif
      break;

    case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mass_radio;
	    n++;
	  }
#endif
      break;

    case IO_ACRB:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Hsml;
	    n++;
	  }
#endif
      break;

    case IO_MACH:
#if defined(MACHNUM) && !defined(AB_SHOCK)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_MachNumber;
	    n++;
	  }
#endif

#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Mach;
	    n++;
	  }
#endif
      break;

    case IO_MACH1:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Mach1;
	    n++;
	  }
#endif
      break;

    case IO_MACH2:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Mach2;
	    n++;
	  }
#endif
      break;

    case IO_MACH3:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Mach3;
	    n++;
	  }
#endif
      break;

    case IO_RUP:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Up_Rho;
	    n++;
	  }
#endif
      break;

    case IO_RDOWN:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Down_Rho;
	    n++;
	  }
#endif
      break;

    case IO_PUP:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Up_Pressure;
	    n++;
	  }
#endif
      break;

    case IO_PDOWN:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Down_Pressure;
	    n++;
	  }
#endif
      break;

    case IO_SHSPEED:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Speed;
	    n++;
	  }
#endif
      break;

    case IO_SHCOMPRESS:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_Compress;
	    n++;
	  }
#endif
      break;

    case IO_SHNORMAL:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].Shock_N1[k]);
	    n++;
	  }
#endif
      break;

    case IO_VUP:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].Shock_Up_V1[k]);
	    n++;
	  }
#endif
      break;

    case IO_VDOWN:
#ifdef AB_SHOCK
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = (SphP[pindex].Shock_Down_V1[k]);
	    n++;
	  }
#endif
      break;

    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_DtEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_CSND:
#ifdef OUTPUT_PRESHOCK_CSND
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalSoundSpeed;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_DENSITY:
#if defined(CR_OUTPUT_JUMP_CONDITIONS) || defined(OUTPUT_PRESHOCK_CSND)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalDensity;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_ENERGY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_XCR:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_XCR;
	    n++;
	  }
#endif
      break;

    case IO_DENSITY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_DensityJump;
	    n++;
	  }
#endif
      break;

    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_EnergyJump;
	    n++;
	  }
#endif
      break;

    case IO_CRINJECT:
#if defined( COSMIC_RAYS ) && defined( CR_OUTPUT_INJECTION )
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CR_Specific_SupernovaHeatingRate;
	    n++;
	  }
#endif
      break;

    case IO_TIDALTENSORPS:
      /* 3x3 configuration-space tidal tensor that is driving the GDE */
#ifdef OUTPUT_TIDALTENSORPS
      for(n = 0; n < pc; pindex++)

	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		for(l = 0; l < 3; l++)
		  {
		    fp[k * 3 + l] = (MyOutputFloat) P[pindex].tidal_tensorps[k][l];
#if defined(PMGRID)
		    fp[k * 3 + l] += (MyOutputFloat) P[pindex].tidal_tensorpsPM[k][l];
#endif

		  }
	      }

	    fflush(stderr);
	    n++;
	    fp += 9;
	  }
#endif
      break;


    case IO_DISTORTIONTENSORPS:
      /* full 6D phase-space distortion tensor from GDE integration */
#ifdef OUTPUT_DISTORTIONTENSORPS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    get_half_kick_distortion(pindex, half_kick_add);
	    for(k = 0; k < 6; k++)
	      {
		for(l = 0; l < 6; l++)
		  {
#ifdef GDE_LOGOUTPUT
		    tmpval1 = P[pindex].distortion_tensorps[k][l] + half_kick_add[k][l];
		    tmpval2 = GDE_ABS(P[pindex].distortion_tensorps[k][l] + half_kick_add[k][l]) + 1.0;
		    if(tmpval1 == zero)
		      fp[k * 6 + l] = (MyOutputFloat) 0;
		    if(tmpval1 > zero)
		      fp[k * 6 + l] = +(MyOutputFloat) (GDE_LOG(tmpval2) / log(10));
		    if(tmpval1 < zero)
		      fp[k * 6 + l] = -(MyOutputFloat) (GDE_LOG(tmpval2) / log(10));
#else
		    fp[k * 6 + l] =
		      (MyOutputFloat) (P[pindex].distortion_tensorps[k][l] + half_kick_add[k][l]);
#endif

		  }
	      }
	    n++;
	    fp += 36;

	  }
#endif
      break;

    case IO_CAUSTIC_COUNTER:
      /* caustic counter */
#ifdef DISTORTIONTENSORPS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (MyOutputFloat) P[pindex].caustic_counter;
	    n++;
	  }
#endif
      break;

    case IO_FLOW_DETERMINANT:
      /* physical NON-CUTOFF corrected stream determinant = 1.0/normed stream density * 1.0/initial stream density */
#if defined(DISTORTIONTENSORPS) && !defined(GDE_LEAN)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    get_current_ps_info(pindex, &flde, &psde);
	    *fp++ = (MyOutputFloat) flde;
	    n++;
	  }
#endif
      break;

    case IO_STREAM_DENSITY:
      /* physical CUTOFF corrected stream density = normed stream density * initial stream density */
#ifdef DISTORTIONTENSORPS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifdef GDE_LOGOUTPUT
	    tmpval1 = GDE_ABS(P[pindex].stream_density);
	    if(tmpval1 > zero)
	      *fp++ = (MyOutputFloat) (GDE_LOG(tmpval1) / log(10));
	    else
	      *fp++ = (MyOutputFloat) 0;
#else
	    *fp++ = (MyOutputFloat) (P[pindex].stream_density);
#endif
	    n++;
	  }
#endif
      break;

    case IO_PHASE_SPACE_DETERMINANT:
      /* determinant of phase-space distortion tensor -> should be 1 due to Liouville theorem */
#ifdef DISTORTIONTENSORPS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    get_current_ps_info(pindex, &flde, &psde);
	    *fp++ = (MyOutputFloat) psde;
	    n++;
	  }
#endif
      break;

    case IO_ANNIHILATION_RADIATION:
      /* time integrated stream density in physical units */
#if defined(DISTORTIONTENSORPS) && !defined(GDE_LEAN)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (MyOutputFloat) (P[pindex].annihilation * GDE_INITDENSITY(pindex));
	    *fp++ = (MyOutputFloat) (P[pindex].analytic_caustics);
	    *fp++ = (MyOutputFloat) (P[pindex].analytic_annihilation * GDE_INITDENSITY(pindex));
	    n++;
	  }
#endif
      break;

    case IO_LAST_CAUSTIC:
      /* extensive information on the last caustic the particle has passed */
#ifdef OUTPUT_LAST_CAUSTIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (MyOutputFloat) P[pindex].lc_Time;
	    *fp++ = (MyOutputFloat) P[pindex].lc_Pos[0];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Pos[1];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Pos[2];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Vel[0];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Vel[1];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Vel[2];
	    *fp++ = (MyOutputFloat) P[pindex].lc_rho_normed_cutoff;

	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_x[0];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_x[1];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_x[2];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_y[0];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_y[1];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_y[2];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_z[0];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_z[1];
	    *fp++ = (MyOutputFloat) P[pindex].lc_Dir_z[2];

	    *fp++ = (MyOutputFloat) P[pindex].lc_smear_x;
	    *fp++ = (MyOutputFloat) P[pindex].lc_smear_y;
	    *fp++ = (MyOutputFloat) P[pindex].lc_smear_z;
	    n++;
	  }
#endif
      break;

    case IO_SHEET_ORIENTATION:
      /* initial orientation of the CDM sheet where the particle started */
#if defined(DISTORTIONTENSORPS) && (!defined(GDE_LEAN) || defined(GDE_READIC))
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 0, 0);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 0, 1);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 0, 2);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 1, 0);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 1, 1);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 1, 2);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 2, 0);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 2, 1);
	    *fp++ = (MyOutputFloat) GDE_VMATRIX(pindex, 2, 2);
	    n++;
	  }
#endif
      break;

    case IO_INIT_DENSITY:
      /* initial stream density in physical units  */
#if defined(DISTORTIONTENSORPS) && (!defined(GDE_LEAN) || defined(GDE_READIC))
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(All.ComovingIntegrationOn)
	      *fp++ =
		GDE_INITDENSITY(pindex) / (GDE_TIMEBEGIN(pindex) * GDE_TIMEBEGIN(pindex) *
					   GDE_TIMEBEGIN(pindex));
	    else
	      *fp++ = GDE_INITDENSITY(pindex);
	    n++;
	  }
#endif
      break;

    case IO_SECONDORDERMASS:
      break;

    case IO_EOSTEMP:
#ifdef EOS_DEGENERATE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].temp;
	    n++;
	  }
#endif
      break;

    case IO_EOSXNUC:
#ifdef EOS_DEGENERATE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < EOS_NSPECIES; k++)
	      {
		*fp++ = SphP[pindex].xnuc[k];
	      }
	    n++;
	  }
#endif
      break;

    case IO_PRESSURE:
#if defined(EOS_DEGENERATE)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Pressure;
	    n++;
	  }
#endif
      break;

    case IO_RADGAMMA:
#ifdef RADTRANSFER
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < N_BINS; k++)
	      fp[k] = SphP[pindex].n_gamma[k] / 1e53;

	    n++;
	    fp += N_BINS;
	  }
#endif
      break;

    case IO_RAD_ACCEL:
#if defined(RADTRANSFER) && defined(RT_OUTPUT_RAD_ACCEL)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = SphP[pindex].RadAccel[k];

	    n++;
	    fp += 3;
	  }
#endif
      break;

    case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
#ifdef RT_OUTPUT_ET
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 6; k++)
	      {
		fp[k] = SphP[pindex].ET[k];
	      }
	    n++;
	    fp += 6;
	  }
#endif
#endif
      break;

    case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].DM_Hsml;
	    n++;
	  }
#endif
      break;

    case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].u.DM_Density;
	    n++;
	  }
#endif
      break;

    case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].v.DM_VelDisp;
	    n++;
	  }
#endif
      break;

    case IO_DMHSML_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].DM_Hsml_V;
	    n++;
	  }
#endif
      break;

    case IO_DMDENSITY_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].DM_Density_V;
	    n++;
	  }
#endif
      break;

    case IO_Zs:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < LT_NMetP; k++)
	      control[k] = 0;
	    for(k = 0; k < LT_NMetP; k++)
	      {
		if(type == 4)
		  *fp_single++ = MetP[P[pindex].pt.MetID].Metals[k];
		else
		  *fp_single++ = SphP[pindex].Metals[k];
		if(*(fp_single - 1) > 1)
		  control[k]++;
#ifdef LT_SEvDbg
		metals[k + LT_NMetP * (type == 4)] += *(fp - 1);
		metals[2 * LT_NMetP + 2] += *(fp - 1);
		metals[2 * LT_NMetP + (type == 4)] += *(fp - 1);
#endif
	      }
	    if(*(fp_single - 1) > 1)
	      {
		printf(" ------- [%d][%d] ", ThisTask, pindex);
		for(k = 0; k < LT_NMetP; k++)
		  printf(" %d ", control[k]);
		printf("\n");
	      }
	    n++;
	  }
#endif
      break;
    case IO_iMass:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = MetP[P[pindex].pt.MetID].iMass;
	    n++;
	  }
#endif
      break;
    case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].XColdCloud;
	    n++;
	  }
#endif
      break;
    case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    double u, Ne, mu, rho;
#ifdef LT_METAL_COOLING_WAL
	    double Metallicities[LT_NMet];
#endif


	    u =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv,
							     GAMMA_MINUS1) * All.UnitPressure_in_cgs /
		   All.UnitDensity_in_cgs);
#ifndef LT_METAL_COOLING_WAL

#ifndef UM_CHEMISTRY
	    Ne = (double) SphP[pindex].Ne;
	    mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[pindex].Ne);
#else
	    Ne = (double) SphP[pindex].elec;
	    mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[pindex].elec);
#endif

#endif

	    rho = SphP[pindex].d.Density * a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
	    if(SphP[pindex].XColdCloud > 0)
	      {
		u -= SphP[pindex].XColdCloud * All.EgySpecCold;	/* obtain the energy density of the hot phase */
		u /= (1 - SphP[pindex].XColdCloud);
		rho *= (1 - SphP[pindex].XColdCloud);
	      }

#ifndef LT_METAL_COOLING_WAL
	    *fp++ = (MyFloat) convert_u_to_temp(u, rho, &Ne);
#else
	    set_metallicities(pindex, &Metallicities[0], a3inv);
	    Metallicities[HydPos] *= (1 - SphP[pindex].XColdCloud);

	    *fp++ = (MyFloat) convert_u_to_temp(u, Redshift, DZ, &Metallicities[0]);
#endif
	    n++;
	  }
#endif
      break;

    case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    double u, Ne, mu, rho;
#ifdef LT_METAL_COOLING_WAL
	    double Metallicities[LT_NMet];
#endif

	    u =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv,
							     GAMMA_MINUS1) * All.UnitPressure_in_cgs /
		   All.UnitDensity_in_cgs);
#ifndef LT_METAL_COOLING_WAL

	    Ne = (double) SphP[pindex].Ne;
	    mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[pindex].Ne);

	    rho =
	      SphP[pindex].a2.Density * a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
	    *fp++ = (MyFloat) convert_u_to_temp(u, rho, &Ne);
#else
	    set_metallicities(pindex, &Metallicities[0], a3inv);
	    *fp++ = (MyFloat) convert_u_to_temp(u, Redshift, DZ, &Metallicities[0]);
#endif

	    n++;
	  }
#endif
      break;

    case IO_ZAGE:
#ifdef LT_ZAGE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(type == 0)
	      {
		if(SphP[pindex].ZAgeW > 0)
#ifndef LT_LOGZAGE
		  *fp++ = SphP[pindex].ZAge / SphP[pindex].ZAgeW;
#else
		  *fp++ = exp(10, SphP[pindex].ZAge) / SphP[pindex].ZAgeW;
#endif
		else
		  *fp++ = 0;
	      }
	    else
	      *fp++ = MetP[P[pindex].pt.MetID].ZAge;
	    n++;
	  }
#endif
      break;

    case IO_ZAGE_LLV:
#ifdef LT_ZAGE_LLV
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(type == 0)
	      {
		if(SphP[pindex].ZAgeW > 0)
#ifndef LT_LOGZAGE
		  *fp++ = SphP[pindex].ZAge_llv / SphP[pindex].ZAgeW_llv;
#else
		  *fp++ = exp(10, SphP[pindex].ZAge_llv) / SphP[pindex].ZAgeW_llv;
#endif
		else
		  *fp++ = 0;
	      }
	    else
	      *fp++ = MetP[P[pindex].pt.MetID].ZAge_llv;
	    n++;
	  }
#endif
      break;
    case IO_CONTRIB:
#ifdef LT_TRACK_CONTRIBUTES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(type == 0)
	      *contrib++ = SphP[pindex].contrib;
	    else if(type == 4)
	      *contrib++ = MetP[P[pindex].pt.MetID].contrib;
	    n++;
	  }
#endif
      break;
    case IO_ZSMOOTH:
#if defined(LT_SMOOTH_Z) && !defined(LT_SMOOTH_ALLMETALS)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = SphP[pindex].Zsmooth;
	    n++;
	  }
#endif
      break;
    case IO_allZSMOOTH:
#if defined(LT_SMOOTH_ALLMETALS)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < LT_NMet; k++)
	      *fp_single++ = SphP[pindex].Zsmooth[k];
	    n++;
	  }
#endif

    case IO_ABUNDANCE:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < LT_NMetP; k++)
	      *fp++ = SphP[pindex].Abundance[k];
	    n++;
	  }
#endif
      break;

    case IO_ABND_GRAD:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < LT_NMetP; k++)
	      for(cc = 0; cc < 3; cc++)
		*fp++ = SphP[pindex].GradAbnd[k][cc];
	    n++;
	  }
#endif
      break;

    case IO_ABND_SPH:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < LT_NMetP; k++)
	      *fp++ = SphP[pindex].AvgAbnd[k];
	    n++;
	  }
#endif
      break;

    case IO_DIFFUSING_CB:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < LT_NMetP; k++)
	      *fp++ = SphP[pindex].B[k];
	    n++;
	  }
#endif
      break;

    case IO_DIFFUSING_CA:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].A;
	    n++;
	  }
#endif
      break;

    case IO_DIFF_TSTEP:
#ifdef LT_MV_CHEMICALDIFFUSION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifdef LT_MV_CHEMICALDIFFUSION_CONSTANT_D
	    myD = All.ConstantD;
#else
	    myD = SphP[pindex].D;
#endif
	    if(All.StandardCourantDiffusion)
	      tscale = PPP[pindex].Hsml * PPP[pindex].Hsml * All.DiffusionCourantFac / myD;
	    else
	      {
		maxgrad = 0;
		double AbndGradN;
		for(el = 0; el < LT_NMetP; el++)
		  {
		    AbndGradN = VECT_NORM(3, SphP[pindex].GradAbnd[el]);
		    if(maxgrad < AbndGradN)
		      {
			maxgrad = AbndGradN;
			//myAvgA = SphP[pindex].AvgAbnd[el];
			myAvgA = SphP[pindex].Abundance[el];
		      }
		  }

		if(maxgrad != 0)
		  {
		    if((tscale = myAvgA / (myD * fabs(maxgrad))) > All.DiffusionCourantFac)
		      tscale = All.DiffusionCourantFac;
		    if(tscale < All.DiffusionCourantMinFac)
		      tscale = All.DiffusionCourantMinFac;
		    tscale *= 2 * PPP[pindex].Hsml;
		  }
		else
		  tscale = PPP[pindex].Hsml * PPP[pindex].Hsml * 0.05 / myD;

	      }

	    *fp++ = tscale;
	    n++;
	  }
#endif
      break;

    case IO_DIFFUSING_CD:
#if defined(LT_MV_CHEMICALDIFFUSION) && !defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#if !defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
#if !defined(LT_MV_CHEMICALDIFFUSION_NO_DENSITY_IN_D)
	    *fp++ = 2 * sqrt(SphP[pindex].v_rms.vel_rms) * SphP[pindex].d.Density * PPP[pindex].Hsml;
#else
	    *fp++ = 2 * sqrt(SphP[pindex].v_rms.vel_rms);
#endif
#else
	    *fp++ = All.ConstantD;
#endif
	    n++;
	  }
#endif
      break;


      break;
    case IO_CHEM:
#ifdef CHEMCOOL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < TRAC_NUM; k++)
	      *fp++ = SphP[pindex].TracAbund[k];

	    n++;
	  }
#endif
      break;
    case IO_AGS_SOFT:		/* Adaptive Gravitational Softening: softening */
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTGRAVSOFT)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AGS_Hsml;
	    n++;
	  }
#endif
      break;
    case IO_AGS_DENS:		/* Adaptive Gravitational Softening: (number) density */
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTGRAVNUMDENS)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AGS_Density;
	    n++;
	  }
#endif
      break;
    case IO_AGS_ZETA:		/* Adaptive Gravitational Softening: zeta */
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTZETA)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AGS_zeta;
	    n++;
	  }
#endif
      break;
    case IO_AGS_OMEGA:		/* Adaptive Gravitational Softening: omega */
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTOMEGA)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AGS_omega;
	    n++;
	  }
#endif
      break;
    case IO_AGS_NGBS:		/* Adaptive Gravitational Softening: neighbours */
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTNGBS)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AGS_REALNumNgb;
	    n++;
	  }
#endif
      break;
    case IO_AGS_CORR:		/* Adaptive Gravitational Softening: correction */
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AGS_corr;
	    n++;
	  }
#endif
      break;
    case IO_PSUM:
#ifdef SIDM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].sidm_PSum;
	    n++;
	  }
#endif
      break;

    case IO_SIDMNUMNGB:
#ifdef SIDM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].sidm_NumNgb;
	    n++;
	  }
#endif
      break;

    case IO_NUMTOTALSCATTER:
#ifdef SIDM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].sidm_NumTotalScatter;
	    n++;
	  }
#endif
      break;

    case IO_SIDMHSML:
#ifdef SIDM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].sidm_Hsml;
	    n++;
	  }
#endif
      break;

    case IO_SIDMDENSITY:
#ifdef SIDM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].sidm_Density;
	    n++;
	  }
#endif
      break;

    case IO_SIDMVELDISP:
#ifdef SIDM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp_single++ = P[pindex].sidm_VelDisp;
	    n++;
	  }
#endif
      break;

    case IO_MG_PHI:		/* scalar field fluctuations */
#ifdef MODGRAV
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].mg_phi;
	    n++;
	  }
#endif
      break;

    case IO_MG_ACCEL:		/* modified gravity acceleration */
#if defined(OUTPUTACCELERATION) && defined(MODGRAV)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].ModGravAccel[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_LASTENTRY:
      endrun(213);
      break;
    }

  *startindex = pindex;
}




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_VBULK:
    case IO_BFLD:
    case IO_GRADPHI:
    case IO_BSMTH:
    case IO_DBDT:
    case IO_ROTB:
    case IO_SROTB:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_RAD_ACCEL:
    case IO_VORT:
    case IO_SHNORMAL:
    case IO_VUP:
    case IO_VDOWN:
    case IO_MG_ACCEL:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_ID:
      bytes_per_blockelement = sizeof(MyIDType);
      break;

    case IO_BHPROGS:
    case IO_TRUENGB:
    case IO_AGS_NGBS:
      bytes_per_blockelement = sizeof(int);
      break;

    case IO_MASS:
    case IO_SECONDORDERMASS:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
    case IO_HSML:
    case IO_VALPHA:
    case IO_SFR:
    case IO_AGE:
    case IO_DELAYTIME:
    case IO_HSMS:
    case IO_ACRS:
    case IO_POT:
    case IO_DTENTR:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DIVB:
    case IO_VRMS:
    case IO_VRAD:
    case IO_VTAN:
    case IO_VDIV:
    case IO_VROT:
    case IO_DPP:
    case IO_ABVC:
    case IO_AMDC:
    case IO_VTURB:
    case IO_LTURB:
    case IO_ALFA2_DYN:
    case IO_ETA2_DYN:
    case IO_PHI:
    case IO_XPHI:
    case IO_EULERA:
    case IO_EULERB:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_EGYPROM:
    case IO_EGYCOLD:
    case IO_BHMASS:
    case IO_ACRB:
    case IO_BHMDOT:
    case IO_BHMBUB:
    case IO_BHMINI:
    case IO_BHMRAD:
    case IO_MACH:
    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_DTENERGY:
    case IO_PRESHOCK_CSND:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_CAUSTIC_COUNTER:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_EOSTEMP:
    case IO_PRESSURE:
    case IO_INIT_DENSITY:
    case IO_iMass:
    case IO_CLDX:
    case IO_HTEMP:
    case IO_TEMP:
    case IO_ZSMOOTH:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
    case IO_MG_PHI:
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      if(mode)
	bytes_per_blockelement = NUMCRPOP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = NUMCRPOP * sizeof(MyOutputFloat);
      break;

    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
#ifndef BP_REAL_CRs
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;
#else
      if(mode)
	bytes_per_blockelement = BP_REAL_CRs * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = BP_REAL_CRs * sizeof(MyOutputFloat);
      break;
#endif

    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
      bytes_per_blockelement = sizeof(float);
      break;

    case IO_RADGAMMA:
#ifdef RADTRANSFER
      if(mode)
	bytes_per_blockelement = N_BINS * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = N_BINS * sizeof(MyOutputFloat);
#endif
      break;

    case IO_EDDINGTON_TENSOR:
      if(mode)
	bytes_per_blockelement = 6 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 6 * sizeof(MyOutputFloat);


    case IO_Z:
#ifndef CS_MODEL
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;
#else
      if(mode)
	bytes_per_blockelement = 12 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 12 * sizeof(MyOutputFloat);
      break;
#endif

    case IO_TIDALTENSORPS:
      if(mode)
	bytes_per_blockelement = 9 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
      break;

    case IO_DISTORTIONTENSORPS:
      if(mode)
	bytes_per_blockelement = 36 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 36 * sizeof(MyOutputFloat);
      break;

    case IO_ANNIHILATION_RADIATION:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_LAST_CAUSTIC:
      if(mode)
	bytes_per_blockelement = 20 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 20 * sizeof(MyOutputFloat);
      break;

    case IO_SHEET_ORIENTATION:
      if(mode)
	bytes_per_blockelement = 9 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
      break;

    case IO_EOSXNUC:
#ifdef EOS_DEGENERATE
      if(mode)
	bytes_per_blockelement = EOS_NSPECIES * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = EOS_NSPECIES * sizeof(MyOutputFloat);
      break;
#else
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;
#endif

    case IO_Zs:
#ifdef LT_STELLAREVOLUTION
      if(mode)
	bytes_per_blockelement = LT_NMetP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_allZSMOOTH:
#ifdef LT_SMOOTH_ALLMETALS
      if(mode)
	bytes_per_blockelement = LT_NMetP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_ZAGE:
    case IO_ZAGE_LLV:
#if defined(LT_ZAGE) || defined(LT_ZAGE_LLV)
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_CONTRIB:
#ifdef LT_TRACK_CONTRIBUTES
      bytes_per_blockelement = sizeof(Contrib);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_DIFFUSING_CB:
    case IO_ABUNDANCE:
    case IO_ABND_SPH:
#ifdef LT_MV_CHEMICALDIFFUSION
      if(mode)
	bytes_per_blockelement = LT_NMetP * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CD:
    case IO_DIFF_TSTEP:
#ifdef LT_MV_CHEMICALDIFFUSION
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;

    case IO_ABND_GRAD:
#ifdef LT_MV_CHEMICALDIFFUSION
      if(mode)
	bytes_per_blockelement = LT_NMetP * 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = LT_NMetP * 3 * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;



    case IO_CHEM:
#ifdef CHEMCOOL
      if(mode)
	bytes_per_blockelement = TRAC_NUM * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = TRAC_NUM * sizeof(MyOutputFloat);
#else
      bytes_per_blockelement = 0;
#endif
      break;


    case IO_LASTENTRY:
      endrun(214);
      break;
    }

  return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    case IO_TRUENGB:
    case IO_BHPROGS:
      typekey = 0;		/* native int */
      break;

    default:
      typekey = 1;		/* native MyOutputFloat */
      break;
    }

  return typekey;
}



int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_BFLD:
    case IO_BSMTH:
    case IO_GRADPHI:
    case IO_DBDT:
    case IO_ROTB:
    case IO_SROTB:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_VBULK:
    case IO_RAD_ACCEL:
    case IO_VORT:
    case IO_MG_ACCEL:
    case IO_SHNORMAL:
    case IO_VUP:
    case IO_VDOWN:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_SECONDORDERMASS:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
    case IO_HSML:
    case IO_VALPHA:
    case IO_SFR:
    case IO_AGE:
    case IO_DELAYTIME:
    case IO_HSMS:
    case IO_ACRS:
    case IO_POT:
    case IO_DTENTR:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_VRMS:
    case IO_TRUENGB:
    case IO_VTAN:
    case IO_VRAD:
    case IO_VDIV:
    case IO_VROT:
    case IO_DPP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_VTURB:
    case IO_LTURB:
    case IO_ALFA2_DYN:
    case IO_ETA2_DYN:
    case IO_PHI:
    case IO_XPHI:
    case IO_EULERA:
    case IO_EULERB:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_EGYPROM:
    case IO_EGYCOLD:
    case IO_BHMASS:
    case IO_ACRB:
    case IO_BHMDOT:
    case IO_BHPROGS:
    case IO_BHMBUB:
    case IO_BHMINI:
    case IO_BHMRAD:
    case IO_MACH:
    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_DTENERGY:
    case IO_PRESHOCK_CSND:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_CAUSTIC_COUNTER:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_EOSTEMP:
    case IO_PRESSURE:
    case IO_INIT_DENSITY:
    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
    case IO_iMass:
    case IO_CLDX:
    case IO_HTEMP:
    case IO_TEMP:
    case IO_ZAGE:
    case IO_ZAGE_LLV:
    case IO_ZSMOOTH:
    case IO_CONTRIB:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_AGS_NGBS:
    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
    case IO_MG_PHI:
      values = 1;
      break;

    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      values = NUMCRPOP;
      break;

    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
#ifdef BP_REAL_CRs
      values = BP_REAL_CRs;
#endif
      break;

    case IO_EDDINGTON_TENSOR:
      values = 6;
      break;

    case IO_RADGAMMA:
#ifdef RADRANSFER
      values = N_BINS;
#else
      values = 0;
#endif
      break;

    case IO_Z:
#ifndef CS_MODEL
      values = 1;
#else
      values = 12;
#endif
      break;

    case IO_Zs:
#ifdef LT_STELLAREVOLUTION
      values = LT_NMetP;
#else
      values = 0;
#endif
      break;

    case IO_allZSMOOTH:
#ifdef LT_SMOOTH_ALLMETALS
      values = LT_NMetP;
#else
      values = 0;
#endif
      break;

    case IO_ABUNDANCE:
    case IO_DIFFUSING_CB:
    case IO_ABND_SPH:
#ifdef LT_MV_CHEMICALDIFFUSION
      values = LT_NMetP;
#else
      values = 0;
#endif
      break;

    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CD:
    case IO_DIFF_TSTEP:
#ifdef LT_MV_CHEMICALDIFFUSION
      values = 1;
#else
      values = 0;
#endif
      break;

    case IO_ABND_GRAD:
#ifdef LT_MV_CHEMICALDIFFUSION
      values = LT_NMetP * 3;
#else
      values = 0;
#endif
      break;
    case IO_TIDALTENSORPS:
      values = 9;
      break;
    case IO_DISTORTIONTENSORPS:
      values = 36;
      break;
    case IO_ANNIHILATION_RADIATION:
      values = 3;
      break;
    case IO_LAST_CAUSTIC:
      values = 20;
      break;
    case IO_SHEET_ORIENTATION:
      values = 9;
      break;


    case IO_EOSXNUC:
#ifndef EOS_DEGENERATE
      values = 1;
#else
      values = EOS_NSPECIES;
#endif
      break;


    case IO_CHEM:
#ifdef CHEMCOOL
      values = TRAC_NUM;
#else
      values = 0;
#endif
      break;


    case IO_LASTENTRY:
      endrun(215);
      break;
    }
  return values;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, nsel, ntot_withmasses, ngas, ngasAlpha, nstars, nngb;

  nall = 0;
  nsel = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  ngasAlpha = NUMCRPOP * header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
    case IO_SECONDORDERMASS:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_AGS_NGBS:
    case IO_MG_PHI:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_RAD_ACCEL:
    case IO_RADGAMMA:
    case IO_EDDINGTON_TENSOR:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
    case IO_HSML:
    case IO_DELAYTIME:
    case IO_VALPHA:
    case IO_SFR:
    case IO_DTENTR:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_BSMTH:
    case IO_BFLD:
    case IO_DBDT:
    case IO_VRMS:
    case IO_VBULK:
    case IO_VTAN:
    case IO_VRAD:
    case IO_VDIV:
    case IO_VROT:
    case IO_VORT:
    case IO_DPP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_VTURB:
    case IO_LTURB:
    case IO_ALFA2_DYN:
    case IO_ETA2_DYN:
    case IO_PHI:
    case IO_XPHI:
    case IO_GRADPHI:
    case IO_ROTB:
    case IO_SROTB:
    case IO_EULERA:
    case IO_EULERB:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_MACH:
    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_SHNORMAL:
    case IO_VUP:
    case IO_VDOWN:
    case IO_DTENERGY:
    case IO_PRESHOCK_CSND:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_EOSTEMP:
    case IO_EOSXNUC:
    case IO_PRESSURE:
    case IO_CHEM:
    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
    case IO_ABUNDANCE:
    case IO_DIFFUSING_CD:
    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CB:
    case IO_ABND_GRAD:
    case IO_ABND_SPH:
    case IO_DIFF_TSTEP:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;

    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngasAlpha;
      break;

    case IO_AGE:
      for(i = 0; i < 6; i++)
#ifdef BLACK_HOLES
	if(i != 4 && i != 5)
	  typelist[i] = 0;
      return nstars + header.npart[5];
#else
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
#endif
      break;

    case IO_TRUENGB:
      nngb = ngas;
      for(i = 1; i < 4; i++)
	typelist[i] = 0;
#ifdef LT_STELLAREVOLUTION
      nngb += header.npart[4];
#else
      typelist[4] = 0;
#endif
#ifdef BLACK_HOLES
      nngb += header.npart[5];
#else
      typelist[5] = 0;
#endif
      return nngb;
      break;

    case IO_HSMS:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;

    case IO_ACRS:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;

    case IO_Z:
    case IO_EGYPROM:
    case IO_EGYCOLD:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;

    case IO_BHMASS:
    case IO_ACRB:
    case IO_BHMDOT:
    case IO_BHMBUB:
    case IO_BHMINI:
    case IO_BHMRAD:
    case IO_BHPROGS:
      for(i = 0; i < 6; i++)
	if(i != 5)
	  typelist[i] = 0;
      return header.npart[5];
      break;

    case IO_TIDALTENSORPS:
    case IO_DISTORTIONTENSORPS:
    case IO_CAUSTIC_COUNTER:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_ANNIHILATION_RADIATION:
    case IO_LAST_CAUSTIC:
    case IO_SHEET_ORIENTATION:
    case IO_INIT_DENSITY:
      for(i = 0; i < 6; i++)
	if(((1 << i) & (GDE_TYPES)))
	  nsel += header.npart[i];
	else
	  typelist[i] = 0;
      return nsel;
      break;


    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
      for(i = 0; i < 6; i++)
	if(((1 << i) & (FOF_PRIMARY_LINK_TYPES)))
	  nsel += header.npart[i];
	else
	  typelist[i] = 0;
      return nsel;
      break;

    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
      for(i = 0; i < 6; i++)
#ifdef SIDM
	if(((1 << i) & (SIDM)))
	  nsel += header.npart[i];
	else
	  typelist[i] = 0;
#else
	typelist[i] = 0;
#endif
      return nsel;
      break;


    case IO_Zs:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;
    case IO_ZAGE:
    case IO_ZAGE_LLV:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;
    case IO_iMass:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;
      break;
    case IO_CLDX:
    case IO_HTEMP:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;
    case IO_TEMP:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;
    case IO_ZSMOOTH:
      for(i = 0; i < 6; i++)
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
	if(i != 0)
#else
	if(i != 0 && i != 4)
#endif
	  typelist[i] = 0;
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
      return ngas;
#else
      return ngas + nstars;
#endif
      break;
    case IO_allZSMOOTH:
      for(i = 0; i < 6; i++)
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
	if(i != 0)
#else
	if(i != 0 && i != 4)
#endif
	  typelist[i] = 0;
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
      return ngas;
#else
      return ngas + nstars;
#endif
      break;
    case IO_CONTRIB:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;

    case IO_MG_ACCEL:
      for(i = 2; i < 6; i++)
	typelist[i] = 0;
      return ngas + header.npart[1];
      break;

    case IO_LASTENTRY:
      endrun(216);
      break;
    }

  endrun(212);
  return 0;
}



/*! This function tells whether a block in the output file is present or not.
 */
int blockpresent(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
      return 1;			/* always present */

    case IO_NE:
    case IO_NH:
      if(All.CoolingOn == 0)
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY) || defined(RADTRANSFER)
	return 1;
#else
	return 0;
#endif
      else
	return 1;
      break;

    case IO_RADGAMMA:
#ifdef RADTRANSFER
      return N_BINS;
#else
      return 0;
#endif
      break;

    case IO_RAD_ACCEL:
#ifdef RT_OUTPUT_RAD_ACCEL
      return 3;
#else
      return 0;
#endif

    case IO_VBULK:		/* mean velocity in kernel */
#if defined(JD_VTURB) && !defined(FS_TURB_ESTIM)
      return 1;
#else
      return 0;
#endif

    case IO_HSMS:
#if defined(SUBFIND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ACRS:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SFR:
    case IO_AGE:
    case IO_Z:
      if(All.StarformationOn == 0)
	return 0;
      else
	{
#ifdef SFR
	  if(blocknr == IO_SFR)
	    return 1;
#endif
#ifdef STELLARAGE
	  if(blocknr == IO_AGE)
	    return 1;
#endif
#ifdef METALS
	  if(blocknr == IO_Z)
	    return 1;
#endif
	}
      return 0;
      break;

    case IO_DELAYTIME:
#ifdef WINDS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VALPHA:
#ifdef VORONOI_TIME_DEP_ART_VISC
      return 1;
#else
      return 0;
#endif
      break;

    case IO_EGYPROM:
    case IO_EGYCOLD:
#ifdef CS_FEEDBACK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HeI:
    case IO_HeII:
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY) ||defined(RADTRANSFER)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HII:
    case IO_HeIII:
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
      return 1;
#else
      return 0;
#endif
      break;



    case IO_H2I:
    case IO_H2II:
    case IO_HM:
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HD:
    case IO_DI:
    case IO_DII:
    case IO_HeHII:
#if defined (UM_CHEMISTRY)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_POT:
#if defined(OUTPUTPOTENTIAL) || defined(SUBFIND_RESHUFFLE_AND_POTENTIAL)
      return 1;
#else
      return 0;
#endif


    case IO_VSTURB_DISS:
    case IO_VSTURB_DRIVE:
#if defined (VS_TURB) || defined (AB_TURB)
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ACCEL:
#ifdef OUTPUTACCELERATION
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DTENTR:
#ifdef OUTPUTCHANGEOFENTROPY
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSDIAG:
#ifdef OUTPUTSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSOFFDIAG:
#ifdef OUTPUTSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSBULK:
#ifdef OUTPUTBULKSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SHEARCOEFF:
#ifdef OUTPUTSHEARCOEFF
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TSTP:
#ifdef OUTPUTTIMESTEP
      return 1;
#else
      return 0;
#endif


    case IO_BFLD:
#ifdef MAGNETIC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_BSMTH:
#ifdef OUTPUTBSMOOTH
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DBDT:
#ifdef DBOUTPUT
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VRMS:
#if defined(JD_VTURB) && !defined(FS_TURB_ESTIM)
    case IO_VBULK:
#endif
    case IO_TRUENGB:
#if defined(JD_VTURB)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VRAD:
    case IO_VTAN:
#if defined(JD_DECOMPOSE_VTURB)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VDIV:
    case IO_VROT:
#if defined(OUTPUT_DIV_CURL)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VORT:
#ifdef OUTPUT_VORTICITY
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DPP:
#ifdef JD_DPP
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DIVB:
#ifdef TRACEDIVB
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ABVC:
#ifdef TIME_DEP_ART_VISC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_AMDC:
#ifdef TIME_DEP_MAGN_DISP
      return 1;
#else
      return 0;
#endif
      break;

    case IO_LTURB:
#ifdef FS_TURB_ESTIM
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VTURB:
#ifdef FS_TURB_ESTIM
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ALFA2_DYN:
#ifdef FS_ALFA2_DYN
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ETA2_DYN:
#ifdef FS_ETA2_DYN
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PHI:
#ifdef OUTPUTDEDNER
      return 1;
#else
      return 0;
#endif
      break;

    case IO_XPHI:
#if (defined(SFR) && (defined(OUTPUT_XPHI) || defined(OUTPUTDEDNER)))
      return 1;
#else
      return 0;
#endif
      break;

    case IO_GRADPHI:
#ifdef OUTPUTDEDNER
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ROTB:
#ifdef OUTPUT_ROTB
      return 1;
#else
      return 0;
#endif
      break;


    case IO_SROTB:
#ifdef OUTPUT_SROTB
      return 1;
#else
      return 0;
#endif
      break;

    case IO_EULERA:
#ifdef EULERPOTENTIALS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_EULERB:
#ifdef EULERPOTENTIALS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_COOLRATE:
#ifdef OUTPUTCOOLRATE
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CONDRATE:
#ifdef OUTPUTCONDRATE
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DENN:
#ifdef OUTPUTDENSNORM
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ACRB:
    case IO_BHMASS:
    case IO_BHMDOT:
#ifdef BLACK_HOLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BHPROGS:
#ifdef BH_COUNTPROGS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BHMBUB:
    case IO_BHMINI:
#ifdef BH_BUBBLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MACH:
#if defined(MACHNUM) || defined(AB_SHOCK)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SHSPEED:
    case IO_SHCOMPRESS:
    case IO_SHNORMAL:
#ifdef AB_SHOCK
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MACH1:
    case IO_MACH2:
    case IO_MACH3:
    case IO_RUP:
    case IO_RDOWN:
    case IO_PUP:
    case IO_PDOWN:
    case IO_VUP:
    case IO_VDOWN:
#if defined(AB_SHOCK_EXTRAOUTPUT)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PRESHOCK_CSND:
#ifdef OUTPUT_PRESHOCK_CSND
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_C0:
    case IO_CR_Q0:
#ifdef COSMIC_RAYS
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
#ifdef CR_OUTPUT_THERMO_VARIABLES
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
#ifdef CR_OUTPUT_TIMESCALES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PRESHOCK_DENSITY:
#if defined(CR_OUTPUT_JUMP_CONDITIONS) || defined(OUTPUT_PRESHOCK_CSND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_BPCR_pNORM:
    case IO_BPCR_eNORM:
    case IO_BPCR_pSLOPE:
    case IO_BPCR_eSLOPE:
    case IO_BPCR_pE:
    case IO_BPCR_pN:
    case IO_BPCR_pPRESSURE:
    case IO_BPCR_ePRESSURE:
#ifdef BP_REAL_CRs
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CRINJECT:
#ifdef CR_OUTPUT_INJECTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TIDALTENSORPS:
#ifdef OUTPUT_TIDALTENSORPS
      return 1;
#else
      return 0;
#endif
    case IO_DISTORTIONTENSORPS:
#ifdef OUTPUT_DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_CAUSTIC_COUNTER:
#ifdef DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_FLOW_DETERMINANT:
#if defined(DISTORTIONTENSORPS) && !defined(GDE_LEAN)
      return 1;
#else
      return 0;
#endif

    case IO_STREAM_DENSITY:
#ifdef DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_PHASE_SPACE_DETERMINANT:
#ifdef DISTORTIONTENSORPS
      return 1;
#else
      return 0;
#endif

    case IO_ANNIHILATION_RADIATION:
#if defined(DISTORTIONTENSORPS) && !defined(GDE_LEAN)
      return 1;
#else
      return 0;
#endif

    case IO_LAST_CAUSTIC:
#ifdef OUTPUT_LAST_CAUSTIC
      return 1;
#else
      return 0;
#endif

    case IO_SHEET_ORIENTATION:
#if defined(DISTORTIONTENSORPS) && (!defined(GDE_LEAN) || defined(GDE_READIC))
      return 1;
#else
      return 0;
#endif

    case IO_INIT_DENSITY:
#if defined(DISTORTIONTENSORPS) && (!defined(GDE_LEAN) || defined(GDE_READIC))
      return 1;
#else
      return 0;
#endif


    case IO_SECONDORDERMASS:
      if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
	return 1;
      else
	return 0;

    case IO_EOSTEMP:
    case IO_EOSXNUC:
    case IO_PRESSURE:
#ifdef EOS_DEGENERATE
      return 1;
#else
      return 0;
#endif

    case IO_EDDINGTON_TENSOR:
#if defined(RADTRANSFER) && defined(RT_OUTPUT_ET)
      return 1;
#else
      return 0;
#endif

    case IO_DMHSML:
    case IO_DMDENSITY:
    case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
#ifdef SIDM
      return 1;
#else
      return 0;
#endif
      break;





    case IO_DMHSML_V:
    case IO_DMDENSITY_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_Zs:
    case IO_iMass:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ZAGE:
#if defined(LT_ZAGE)
      return 1;
#else
      return 0;
#endif
      break;
    case IO_ZAGE_LLV:
#if defined(LT_ZAGE_LLV)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CONTRIB:
#ifdef LT_TRACK_CONTRIBUTES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ZSMOOTH:
#if defined(LT_SMOOTH_Z) && !defined(LT_SMOOTH_ALLMETALS)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_allZSMOOTH:
#ifdef LT_SMOOTH_ALLMETALS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ABUNDANCE:
    case IO_ABND_GRAD:
    case IO_ABND_SPH:
    case IO_DIFFUSING_CB:
    case IO_DIFFUSING_CA:
    case IO_DIFFUSING_CD:
    case IO_DIFF_TSTEP:
#ifdef LT_MV_CHEMICALDIFFUSION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_CHEM:
#ifdef CHEMCOOL
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_SOFT:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTGRAVSOFT)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_DENS:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTGRAVNUMDENS)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_ZETA:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTZETA)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_OMEGA:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTOMEGA)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_CORR:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_AGS_NGBS:
#if defined (ADAPTGRAVSOFT) && defined(AGS_OUTPUTNGBS)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MG_PHI:
#ifdef MODGRAV
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MG_ACCEL:
#if defined(OUTPUTACCELERATION) && defined(MODGRAV)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_LASTENTRY:
      return 0;			/* will not occur */
    }


  return 0;			/* default: not present */
}




/*! This function associates a short 4-character block name with each block number.
 *  This is stored in front of each block for snapshot FileFormat=2.
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
  switch (blocknr)
    {
    case IO_POS:
      strncpy(label, "POS ", 4);
      break;
    case IO_VEL:
      strncpy(label, "VEL ", 4);
      break;
    case IO_ID:
      strncpy(label, "ID  ", 4);
      break;
    case IO_MASS:
      strncpy(label, "MASS", 4);
      break;
    case IO_U:
      strncpy(label, "U   ", 4);
      break;
    case IO_RHO:
      strncpy(label, "RHO ", 4);
      break;
    case IO_NE:
      strncpy(label, "NE  ", 4);
      break;
    case IO_NH:
      strncpy(label, "NH  ", 4);
      break;
    case IO_HII:
      strncpy(label, "HII ", 4);
      break;
    case IO_HeI:
      strncpy(label, "HeI ", 4);
      break;
    case IO_HeII:
      strncpy(label, "HeII", 4);
      break;
    case IO_HeIII:
      strncpy(label, "He3 ", 4);
      break;
    case IO_H2I:
      strncpy(label, "H2I ", 4);
      break;
    case IO_H2II:
      strncpy(label, "H2II", 4);
      break;
    case IO_HM:
      strncpy(label, "HM  ", 4);
      break;
    case IO_HD:
      strncpy(label, "HD  ", 4);
      break;
    case IO_DI:
      strncpy(label, "DI  ", 4);
      break;
    case IO_DII:
      strncpy(label, "DII ", 4);
      break;
    case IO_HeHII:
      strncpy(label, "HeHp", 4);
      break;
    case IO_HSML:
      strncpy(label, "HSML", 4);
      break;
    case IO_VALPHA:
      strncpy(label, "VALP", 4);
      break;
    case IO_SFR:
      strncpy(label, "SFR ", 4);
      break;
    case IO_AGE:
      strncpy(label, "AGE ", 4);
      break;
    case IO_DELAYTIME:
      strncpy(label, "DETI", 4);
      break;
    case IO_HSMS:
      strncpy(label, "HSMS", 4);
      break;
    case IO_ACRS:
      strncpy(label, "ACRS", 4);
      break;
    case IO_Z:
      strncpy(label, "Z   ", 4);
      break;
    case IO_POT:
      strncpy(label, "POT ", 4);
      break;
    case IO_ACCEL:
      strncpy(label, "ACCE", 4);
      break;
    case IO_DTENTR:
      strncpy(label, "ENDT", 4);
      break;
    case IO_STRESSDIAG:
      strncpy(label, "STRD", 4);
      break;
    case IO_STRESSOFFDIAG:
      strncpy(label, "STRO", 4);
      break;
    case IO_STRESSBULK:
      strncpy(label, "STRB", 4);
      break;
    case IO_SHEARCOEFF:
      strncpy(label, "SHCO", 4);
      break;
    case IO_TSTP:
      strncpy(label, "TSTP", 4);
      break;
    case IO_BFLD:
      strncpy(label, "BFLD", 4);
      break;
    case IO_BSMTH:
      strncpy(label, "BFSM", 4);
      break;
    case IO_DBDT:
      strncpy(label, "DBDT", 4);
      break;
    case IO_VBULK:
      strncpy(label, "VBLK", 4);
      break;
    case IO_VRMS:
      strncpy(label, "VRMS", 4);
      break;
    case IO_VTAN:
      strncpy(label, "VTAN", 4);
      break;
    case IO_VRAD:
      strncpy(label, "VRAD", 4);
      break;
    case IO_TRUENGB:
      strncpy(label, "TNGB", 4);
      break;
    case IO_DPP:
      strncpy(label, "DPP ", 4);
      break;
    case IO_VDIV:
      strncpy(label, "VDIV", 4);
      break;
    case IO_VROT:
      strncpy(label, "VROT", 4);
      break;
    case IO_VORT:
      strncpy(label, "VORT", 4);
      break;
    case IO_DIVB:
      strncpy(label, "DIVB", 4);
      break;
    case IO_ABVC:
      strncpy(label, "ABVC", 4);
      break;
    case IO_AMDC:
      strncpy(label, "AMDC", 4);
      break;
    case IO_VTURB:
      strncpy(label, "VTRB", 4);
      break;
    case IO_LTURB:
      strncpy(label, "LTRB", 4);
      break;
    case IO_ALFA2_DYN:
      strncpy(label, "ADYN", 4);
      break;
    case IO_ETA2_DYN:
      strncpy(label, "EDYN", 4);
      break;
    case IO_PHI:
      strncpy(label, "PHI ", 4);
      break;
    case IO_XPHI:
      strncpy(label, "XPHI", 4);
      break;
    case IO_GRADPHI:
      strncpy(label, "GPHI", 4);
      break;
    case IO_ROTB:
      strncpy(label, "ROTB", 4);
      break;
    case IO_SROTB:
      strncpy(label, "SRTB", 4);
      break;
    case IO_EULERA:
      strncpy(label, "EULA", 4);
      break;
    case IO_EULERB:
      strncpy(label, "EULB", 4);
      break;
    case IO_COOLRATE:
      strncpy(label, "COOR", 4);
      break;
    case IO_CONDRATE:
      strncpy(label, "CONR", 4);
      break;
    case IO_DENN:
      strncpy(label, "DENN", 4);
      break;
    case IO_EGYPROM:
      strncpy(label, "EGYP", 4);
      break;
    case IO_EGYCOLD:
      strncpy(label, "EGYC", 4);
      break;
    case IO_CR_C0:
      strncpy(label, "CRC0", 4);
      break;
    case IO_CR_Q0:
      strncpy(label, "CRQ0", 4);
      break;
    case IO_CR_P0:
      strncpy(label, "CRP0", 4);
      break;
    case IO_CR_E0:
      strncpy(label, "CRE0", 4);
      break;
    case IO_CR_n0:
      strncpy(label, "CRn0", 4);
      break;
    case IO_CR_ThermalizationTime:
      strncpy(label, "CRco", 4);
      break;
    case IO_CR_DissipationTime:
      strncpy(label, "CRdi", 4);
      break;
    case IO_BHMASS:
      strncpy(label, "BHMA", 4);
      break;
    case IO_ACRB:
      strncpy(label, "ACRB", 4);
      break;
    case IO_BHMDOT:
      strncpy(label, "BHMD", 4);
      break;
    case IO_BHPROGS:
      strncpy(label, "BHPC", 4);
      break;
    case IO_BHMBUB:
      strncpy(label, "BHMB", 4);
      break;
    case IO_BHMINI:
      strncpy(label, "BHMI", 4);
      break;
    case IO_BHMRAD:
      strncpy(label, "BHMR", 4);
      break;
    case IO_MACH:
      strncpy(label, "MACH", 4);
      break;
    case IO_MACH1:
      strncpy(label, "MVEL", 4);
      break;
    case IO_MACH2:
      strncpy(label, "MPRS", 4);
      break;
    case IO_MACH3:
      strncpy(label, "MRHO", 4);
      break;
    case IO_RUP:
      strncpy(label, "RUP ", 4);
      break;
    case IO_RDOWN:
      strncpy(label, "RDO ", 4);
      break;
    case IO_PUP:
      strncpy(label, "PUP ", 4);
      break;
    case IO_PDOWN:
      strncpy(label, "PDO ", 4);
      break;
    case IO_SHSPEED:
      strncpy(label, "SHSP", 4);
      break;
    case IO_SHCOMPRESS:
      strncpy(label, "SHCP", 4);
      break;
    case IO_SHNORMAL:
      strncpy(label, "SHNR", 4);
      break;
    case IO_VUP:
      strncpy(label, "VUP ", 4);
      break;
    case IO_VDOWN:
      strncpy(label, "VDO ", 4);
      break;
    case IO_DTENERGY:
      strncpy(label, "DTEG", 4);
      break;
    case IO_PRESHOCK_CSND:
      strncpy(label, "PSCS", 4);
      break;
    case IO_PRESHOCK_DENSITY:
      strncpy(label, "PSDE", 4);
      break;
    case IO_PRESHOCK_ENERGY:
      strncpy(label, "PSEN", 4);
      break;
    case IO_PRESHOCK_XCR:
      strncpy(label, "PSXC", 4);
      break;
    case IO_DENSITY_JUMP:
      strncpy(label, "DJMP", 4);
      break;
    case IO_ENERGY_JUMP:
      strncpy(label, "EJMP", 4);
      break;
    case IO_CRINJECT:
      strncpy(label, "CRDE", 4);
      break;
    case IO_TIDALTENSORPS:
      strncpy(label, "TIPS", 4);
      break;
    case IO_DISTORTIONTENSORPS:
      strncpy(label, "DIPS", 4);
      break;
    case IO_CAUSTIC_COUNTER:
      strncpy(label, "CACO", 4);
      break;
    case IO_FLOW_DETERMINANT:
      strncpy(label, "FLDE", 4);
      break;
    case IO_STREAM_DENSITY:
      strncpy(label, "STDE", 4);
      break;
    case IO_SECONDORDERMASS:
      strncpy(label, "SOMA", 4);
      break;
    case IO_PHASE_SPACE_DETERMINANT:
      strncpy(label, "PSDE", 4);
      break;
    case IO_ANNIHILATION_RADIATION:
      strncpy(label, "ANRA", 4);
      break;
    case IO_LAST_CAUSTIC:
      strncpy(label, "LACA", 4);
      break;
    case IO_SHEET_ORIENTATION:
      strncpy(label, "SHOR", 4);
      break;
    case IO_INIT_DENSITY:
      strncpy(label, "INDE", 4);
      break;
    case IO_EOSTEMP:
      strncpy(label, "TEMP", 4);
      break;
    case IO_EOSXNUC:
      strncpy(label, "XNUC", 4);
      break;
    case IO_PRESSURE:
      strncpy(label, "P   ", 4);
      break;
    case IO_RADGAMMA:
      strncpy(label, "RADG", 4);
      break;
    case IO_RAD_ACCEL:
      strncpy(label, "RADA", 4);
      break;
    case IO_EDDINGTON_TENSOR:
      strncpy(label, "ET", 4);
      break;
    case IO_PSUM:
      strncpy(label, "PTSU", 4);
      break;
    case IO_SIDMNUMNGB:
      strncpy(label, "DMNB", 4);
      break;
    case IO_NUMTOTALSCATTER:
      strncpy(label, "NTSC", 4);
      break;
    case IO_SIDMHSML:
      strncpy(label, "SHSM", 4);
      break;
    case IO_SIDMDENSITY:
      strncpy(label, "SRHO", 4);
      break;
    case IO_SIDMVELDISP:
      strncpy(label, "SVEL", 4);
      break;
    case IO_DMHSML:
      strncpy(label, "DMHS", 4);
      break;
    case IO_DMDENSITY:
      strncpy(label, "DMDE", 4);
      break;
    case IO_DMVELDISP:
      strncpy(label, "DMVD", 4);
      break;
    case IO_DMHSML_V:
      strncpy(label, "VDMH", 4);
      break;
    case IO_DMDENSITY_V:
      strncpy(label, "VDMD", 4);
      break;
    case IO_Zs:
      strncpy(label, "Zs  ", 4);
      break;
    case IO_ZAGE:
      strncpy(label, "ZAge", 4);
      break;
    case IO_ZAGE_LLV:
      strncpy(label, "ZAlv", 4);
      break;
    case IO_iMass:
      strncpy(label, "iM  ", 4);
      break;
    case IO_CLDX:
      strncpy(label, "CLDX", 4);
      break;
    case IO_HTEMP:
      strncpy(label, "HOTT", 4);
      break;
    case IO_TEMP:
      strncpy(label, "TEMP", 4);
      break;
    case IO_CONTRIB:
      strncpy(label, "TRCK", 4);
      break;
    case IO_ZSMOOTH:
      strncpy(label, "ZsMT", 4);
      break;
    case IO_allZSMOOTH:
      strncpy(label, "ZSMT", 4);
      break;
    case IO_ABUNDANCE:
      strncpy(label, "DABN", 4);
      break;
    case IO_ABND_GRAD:
      strncpy(label, "DABG", 4);
      break;
    case IO_ABND_SPH:
      strncpy(label, "DABS", 4);
      break;
    case IO_DIFFUSING_CB:
      strncpy(label, "DCFB", 4);
      break;
    case IO_DIFFUSING_CA:
      strncpy(label, "DCFA", 4);
      break;
    case IO_DIFFUSING_CD:
      strncpy(label, "DCFD", 4);
      break;
    case IO_DIFF_TSTEP:
      strncpy(label, "DTST", 4);
      break;
    case IO_CHEM:
      strncpy(label, "CHEM", 4);
      break;
    case IO_BPCR_pNORM:
      strncpy(label, "CRpN", 4);
      break;
    case IO_BPCR_eNORM:
      strncpy(label, "CReN", 4);
      break;
    case IO_BPCR_pSLOPE:
      strncpy(label, "CRpS", 4);
      break;
    case IO_BPCR_eSLOPE:
      strncpy(label, "CReS", 4);
      break;
    case IO_BPCR_pE:
      strncpy(label, "CRpE", 4);
      break;
    case IO_BPCR_pN:
      strncpy(label, "CRpD", 4);
      break;
    case IO_BPCR_pPRESSURE:
      strncpy(label, "CRpP", 4);
      break;
    case IO_BPCR_ePRESSURE:
      strncpy(label, "CReP", 4);
      break;
    case IO_AGS_SOFT:
      strncpy(label, "AGSH", 4);
      break;
    case IO_AGS_DENS:
      strncpy(label, "AGSD", 4);
      break;
    case IO_AGS_ZETA:
      strncpy(label, "AGSZ", 4);
      break;
    case IO_AGS_OMEGA:
      strncpy(label, "AGSO", 4);
      break;
    case IO_AGS_CORR:
      strncpy(label, "AGSC", 4);
      break;
    case IO_AGS_NGBS:
      strncpy(label, "AGSN", 4);
      break;
    case IO_VSTURB_DISS:
      strncpy(label, "VSDI", 4);
      break;
    case IO_VSTURB_DRIVE:
      strncpy(label, "VSDR", 4);
      break;
    case IO_MG_PHI:
      strncpy(label, "MGPH", 4);
      break;
    case IO_MG_ACCEL:
      strncpy(label, "MGAC", 4);
      break;

    case IO_LASTENTRY:
      endrun(217);
      break;
    }
}


void get_dataset_name(enum iofields blocknr, char *buf)
{
  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_NE:
      strcpy(buf, "ElectronAbundance");
      break;
    case IO_NH:
      strcpy(buf, "NeutralHydrogenAbundance");
      break;
    case IO_RADGAMMA:
      strcpy(buf, "photon number density");
      break;
    case IO_RAD_ACCEL:
      strcpy(buf, "rad acceleration");
      break;
    case IO_HII:
      strcpy(buf, "HII");
      break;
    case IO_HeI:
      strcpy(buf, "HeI");
      break;
    case IO_HeII:
      strcpy(buf, "HeII");
      break;
    case IO_HeIII:
      strcpy(buf, "HeIII");
      break;
    case IO_H2I:
      strcpy(buf, "H2I");
      break;
    case IO_H2II:
      strcpy(buf, "H2II");
      break;
    case IO_HM:
      strcpy(buf, "HM");
      break;
    case IO_HD:
      strcpy(buf, "HD  ");
      break;
    case IO_DI:
      strcpy(buf, "DI  ");
      break;
    case IO_DII:
      strcpy(buf, "DII ");
      break;
    case IO_HeHII:
      strcpy(buf, "HeHp");
      break;
    case IO_DELAYTIME:
      strcpy(buf, "DelayTime");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_VALPHA:
      strcpy(buf, "ArtificialViscosityV");
      break;
    case IO_SFR:
      strcpy(buf, "StarFormationRate");
      break;
    case IO_AGE:
      strcpy(buf, "StellarFormationTime");
      break;
    case IO_HSMS:
      strcpy(buf, "StellarSmoothingLength");
      break;
    case IO_ACRS:
      strcpy(buf, "StellarSpreadingLength");
      break;
    case IO_Z:
      strcpy(buf, "Metallicity");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_STRESSDIAG:
      strcpy(buf, "DiagonalStressTensor");
      break;
    case IO_STRESSOFFDIAG:
      strcpy(buf, "OffDiagonalStressTensor");
      break;
    case IO_STRESSBULK:
      strcpy(buf, "BulkStressTensor");
      break;
    case IO_SHEARCOEFF:
      strcpy(buf, "ShearCoefficient");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
    case IO_BFLD:
      strcpy(buf, "MagneticField");
      break;
    case IO_BSMTH:
      strcpy(buf, "SmoothedMagneticField");
      break;
    case IO_DBDT:
      strcpy(buf, "RateOfChangeOfMagneticField");
      break;
    case IO_VRMS:
      strcpy(buf, "RMSVelocity");
      break;
    case IO_VBULK:
      strcpy(buf, "BulkVelocity");
      break;
    case IO_VRAD:
      strcpy(buf, "RMSRadialVelocity");
      break;
    case IO_VTAN:
      strcpy(buf, "RMSTangentialVelocity");
      break;
    case IO_TRUENGB:
      strcpy(buf, "TrueNumberOfNeighbours");
      break;
    case IO_DPP:
      strcpy(buf, "MagnetosReaccCoefficient");
      break;
    case IO_VDIV:
      strcpy(buf, "VelocityDivergence");
      break;
    case IO_VROT:
      strcpy(buf, "VelocityCurl");
      break;
    case IO_VORT:
      strcpy(buf, "Vorticity");
      break;
    case IO_DIVB:
      strcpy(buf, "DivergenceOfMagneticField");
      break;
    case IO_ABVC:
      strcpy(buf, "ArtificialViscosity");
      break;
    case IO_AMDC:
      strcpy(buf, "ArtificialMagneticDissipatio");
      break;
    case IO_VTURB:
      strcpy(buf, "TurbVelociy");
      break;
    case IO_LTURB:
      strcpy(buf, "TurbLengh");
      break;
    case IO_ALFA2_DYN:
      strcpy(buf, "AlfaDynamo");
      break;
    case IO_ETA2_DYN:
      strcpy(buf, "EtaDynamo");
      break;
    case IO_PHI:
      strcpy(buf, "DivBcleaningFunctionPhi");
      break;
    case IO_XPHI:
      strcpy(buf, "ColdGasFraction_PHI");
      break;
    case IO_GRADPHI:
      strcpy(buf, "DivBcleaningFunctionGadPhi");
      break;
    case IO_ROTB:
      strcpy(buf, "RotationB");
      break;
    case IO_SROTB:
      strcpy(buf, "SmoothedRotationB");
      break;
    case IO_EULERA:
      strcpy(buf, "EulerPotentialA");
      break;
    case IO_EULERB:
      strcpy(buf, "EulerPotentialB");
      break;
    case IO_COOLRATE:
      strcpy(buf, "CoolingRate");
      break;
    case IO_CONDRATE:
      strcpy(buf, "ConductionRate");
      break;
    case IO_DENN:
      strcpy(buf, "Denn");
      break;
    case IO_EGYPROM:
      strcpy(buf, "EnergyReservoirForFeeback");
      break;
    case IO_EGYCOLD:
      strcpy(buf, "EnergyReservoirForColdPhase");
      break;
    case IO_CR_C0:
      strcpy(buf, "CR_C0");
      break;
    case IO_CR_Q0:
      strcpy(buf, "CR_q0");
      break;
    case IO_CR_P0:
      strcpy(buf, "CR_P0");
      break;
    case IO_CR_E0:
      strcpy(buf, "CR_E0");
      break;
    case IO_CR_n0:
      strcpy(buf, "CR_n0");
      break;
    case IO_CR_ThermalizationTime:
      strcpy(buf, "CR_ThermalizationTime");
      break;
    case IO_CR_DissipationTime:
      strcpy(buf, "CR_DissipationTime");
      break;
    case IO_BHMASS:
      strcpy(buf, "BH_Mass");
      break;
    case IO_ACRB:
      strcpy(buf, "BH_AccreationLength");
      break;
    case IO_BHMDOT:
      strcpy(buf, "BH_Mdot");
      break;
    case IO_BHPROGS:
      strcpy(buf, "BH_NProgs");
      break;
    case IO_BHMBUB:
      strcpy(buf, "BH_Mass_bubbles");
      break;
    case IO_BHMINI:
      strcpy(buf, "BH_Mass_ini");
      break;
    case IO_BHMRAD:
      strcpy(buf, "BH_Mass_radio");
      break;
    case IO_MACH:
      strcpy(buf, "MachNumber");
      break;
    case IO_MACH1:
      strcpy(buf, "MachNumber Vel");
      break;
    case IO_MACH2:
      strcpy(buf, "MachNumber Press");
      break;
    case IO_MACH3:
      strcpy(buf, "MachNumber Rho");
      break;
    case IO_RUP:
      strcpy(buf, "Rho Upstream");
      break;
    case IO_RDOWN:
      strcpy(buf, "Rho Downstream");
      break;
    case IO_PUP:
      strcpy(buf, "Pressure Upstream");
      break;
    case IO_PDOWN:
      strcpy(buf, "Pressure Downstream");
      break;
    case IO_SHSPEED:
      strcpy(buf, "ShockSpeed");
      break;
    case IO_SHCOMPRESS:
      strcpy(buf, "ShockCompressionRatio");
      break;
    case IO_SHNORMAL:
      strcpy(buf, "ShockNormal");
      break;
    case IO_VUP:
      strcpy(buf, "Velocity Upstream");
      break;
    case IO_VDOWN:
      strcpy(buf, "Velocity Downstream");
      break;
    case IO_DTENERGY:
      strcpy(buf, "DtEnergy");
      break;
    case IO_PRESHOCK_CSND:
      strcpy(buf, "Preshock_SoundSpeed");
      break;
    case IO_PRESHOCK_DENSITY:
      strcpy(buf, "Preshock_Density");
      break;
    case IO_PRESHOCK_ENERGY:
      strcpy(buf, "Preshock_Energy");
      break;
    case IO_PRESHOCK_XCR:
      strcpy(buf, "Preshock_XCR");
      break;
    case IO_DENSITY_JUMP:
      strcpy(buf, "DensityJump");
      break;
    case IO_ENERGY_JUMP:
      strcpy(buf, "EnergyJump");
      break;
    case IO_CRINJECT:
      strcpy(buf, "CR_DtE");
      break;
    case IO_TIDALTENSORPS:
      strcpy(buf, "TidalTensorPS");
      break;
    case IO_DISTORTIONTENSORPS:
      strcpy(buf, "DistortionTensorPS");
      break;
    case IO_CAUSTIC_COUNTER:
      strcpy(buf, "CausticCounter");
      break;
    case IO_FLOW_DETERMINANT:
      strcpy(buf, "FlowDeterminant");
      break;
    case IO_STREAM_DENSITY:
      strcpy(buf, "StreamDensity");
      break;
    case IO_SECONDORDERMASS:
      strcpy(buf, "2lpt-mass");
      break;
    case IO_PHASE_SPACE_DETERMINANT:
      strcpy(buf, "PhaseSpaceDensity");
      break;
    case IO_ANNIHILATION_RADIATION:
      strcpy(buf, "AnnihilationRadiation");
      break;
    case IO_LAST_CAUSTIC:
      strcpy(buf, "LastCaustic");
      break;
    case IO_SHEET_ORIENTATION:
      strcpy(buf, "SheetOrientation");
      break;
    case IO_INIT_DENSITY:
      strcpy(buf, "InitDensity");
      break;
    case IO_EOSTEMP:
      strcpy(buf, "Temperature");
      break;
    case IO_EOSXNUC:
      strcpy(buf, "Nuclear mass fractions");
      break;
    case IO_PRESSURE:
      strcpy(buf, "Pressure");
      break;
    case IO_EDDINGTON_TENSOR:
      strcpy(buf, "EddingtonTensor");
      break;
    case IO_PSUM:
      strcpy(buf, "PSum");
      break;
    case IO_SIDMNUMNGB:
      strcpy(buf, "DMNumNgb");
      break;
    case IO_NUMTOTALSCATTER:
      strcpy(buf, "NumTotalScatter");
      break;
    case IO_SIDMHSML:
      strcpy(buf, "SIDMHsml");
      break;
    case IO_SIDMDENSITY:
      strcpy(buf, "SIDMRho");
      break;
    case IO_SIDMVELDISP:
      strcpy(buf, "SVelDisp");
      break;
    case IO_DMHSML:
      strcpy(buf, "DM Hsml");
      break;
    case IO_DMDENSITY:
      strcpy(buf, "DM Density");
      break;
    case IO_DMVELDISP:
      strcpy(buf, "DM Velocity Dispersion");
      break;
    case IO_DMHSML_V:
      strcpy(buf, "DM Hsml Voronoi");
      break;
    case IO_DMDENSITY_V:
      strcpy(buf, "DM Density Voronoi");
      break;
    case IO_Zs:
      strcpy(buf, "Mass of Metals");
      break;
    case IO_ZAGE:
      strcpy(buf, "Metallicity-averaged time");
      break;
    case IO_ZAGE_LLV:
      strcpy(buf, "long-living-Metallicity-averaged time");
      break;
    case IO_iMass:
      strcpy(buf, "SSPInitialMass");
      break;
    case IO_CLDX:
      strcpy(buf, "CloudFraction");
      break;
    case IO_HTEMP:
      strcpy(buf, "HotPhaseTemperature");
      break;
    case IO_TEMP:
      strcpy(buf, "Temperature");
      break;
    case IO_CONTRIB:
      strcpy(buf, "TrackContributes");
      break;
    case IO_ZSMOOTH:
      strcpy(buf, "smoothed metallicity");
      break;
    case IO_allZSMOOTH:
      strcpy(buf, "smoothed metallicities (all elements)");
      break;
    case IO_ABUNDANCE:
      strcpy(buf, "abundance by diffusion");
      break;
    case IO_ABND_GRAD:
      strcpy(buf, "abundance gradient");
      break;
    case IO_ABND_SPH:
      strcpy(buf, "SPH-smoothed abundance");
      break;
    case IO_DIFFUSING_CB:
      strcpy(buf, "Diffusion B-coefficient");
      break;
    case IO_DIFFUSING_CA:
      strcpy(buf, "Diffusion A-coeffient");
      break;
    case IO_DIFFUSING_CD:
      strcpy(buf, "Diffusion D-coefficient");
      break;
    case IO_DIFF_TSTEP:
      strcpy(buf, "Diffusion Time-Step");
      break;
    case IO_CHEM:
      strcpy(buf, "ChemicalAbundances");
      break;
    case IO_BPCR_pNORM:
      strcpy(buf, "BPCRpNormalization");
      break;
    case IO_BPCR_eNORM:
      strcpy(buf, "BPCReNormalization");
      break;
    case IO_BPCR_pSLOPE:
      strcpy(buf, "BPCRpSlope");
      break;
    case IO_BPCR_eSLOPE:
      strcpy(buf, "BPCReSlope");
      break;
    case IO_BPCR_pE:
      strcpy(buf, "BPCRpEnergy");
      break;
    case IO_BPCR_pN:
      strcpy(buf, "BPCRpNumber");
      break;
    case IO_BPCR_pPRESSURE:
      strcpy(buf, "BPCRpPressure");
      break;
    case IO_BPCR_ePRESSURE:
      strcpy(buf, "BPCRePressure");
      break;
    case IO_AGS_SOFT:
      strcpy(buf, "AGS-Softening");
      break;
    case IO_AGS_DENS:
      strcpy(buf, "AGS-Density");
      break;
    case IO_AGS_ZETA:
      strcpy(buf, "AGS-Zeta");
      break;
    case IO_AGS_OMEGA:
      strcpy(buf, "AGS-Omega");
      break;
    case IO_AGS_CORR:
      strcpy(buf, "AGS-Correction");
      break;
    case IO_AGS_NGBS:
      strcpy(buf, "AGS-Neighbours");
      break;
    case IO_VSTURB_DISS:
      strcpy(buf, "TurbulenceDissipation");
      break;
    case IO_VSTURB_DRIVE:
      strcpy(buf, "TurbulenceDriving");
      break;
    case IO_MG_PHI:
      strcpy(buf, "ModifiedGravityPhi");
      break;
    case IO_MG_ACCEL:
      strcpy(buf, "ModifiedGravityAcceleration");
      break;

    case IO_LASTENTRY:
      endrun(218);
      break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the SPH smoothing
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 */
void write_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];
  int n_for_this_task, n, p, pc, offset = 0, task, i;
  size_t blockmaxlen;
  int ntot_type[6], nn[6];
  enum iofields blocknr;
  char label[8];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;

#ifdef WRITE_KEY_FILES
  int NKeys_file[6];
  unsigned int TotNumKeysFile;
  struct io_header key_header;
  FILE *FdKey;
  char buf_key[300];
  struct info_block *KeyInfoBlock;
  int *comm_point_int, *data_point_int, size_of_element, local_offset, write_offset[6];
  peanokey *comm_point_key, *data_point_key;
#define SKIP_KEY  {my_fwrite(&blksize,sizeof(int),1,FdKey);}
#endif

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0, pcsum = 0;
  char buf[500];
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#if defined(LT_SEvDbg) || defined(LT_TRACK_CONTRIBUTES) || defined(LT_POPIII_FLAGS) || defined(LT_SMOOTH_Z)
#ifndef HAVE_HDF5
  char buf[300];
#endif

#ifdef LT_SEvDbg
  FILE *fd_dbg = 0x0;
#endif

#if defined(LT_TRACK_CONTRIBUTES)
  FILE *FdTrck;
  FILE *fd_trck_back;
#define SKIP_TRCK  {my_fwrite(&blksize,sizeof(int),1,FdTrck);}
#endif
#if defined(LT_SMOOTH_Z)
  FILE *FdZsmooth;
  FILE *fd_zsmooth_back;
#define SKIP_ZSMOOTH  {my_fwrite(&blksize,sizeof(int),1,FdZsmooth);}
#endif
#endif

#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
  FILE *fd_pot = 0;

#define SKIP_POT  {my_fwrite(&blksize,sizeof(int),1,fd_pot);}
#endif

  /* determine particle numbers of each type in file */

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	ntot_type[n] = n_type[n];
#ifdef WRITE_KEY_FILES
      for(n = 0; n < 6; n++)
	write_offset[n] = 0;
#endif

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
#ifdef WRITE_KEY_FILES
	  MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_WRT_OFF, MPI_COMM_WORLD);
#endif
	  for(n = 0; n < 6; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
#ifdef WRITE_KEY_FILES
      MPI_Recv(&write_offset[0], 6, MPI_INT, writeTask, TAG_WRT_OFF, MPI_COMM_WORLD, &status);
#endif
      MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

#ifdef WRITE_KEY_FILES

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	NKeys_file[n] = NKeys[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  for(n = 0; n < 6; n++)
	    NKeys_file[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&NKeys_file[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&NKeys[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&NKeys_file[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

  TotNumKeysFile = 0;
  for(n = 0; n < 6; n++)
    TotNumKeysFile += NKeys_file[n];

#endif

  /* fill file header */

  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

#if defined(NEUTRINOS) && defined(OMIT_NEUTRINOS_IN_SNAPS)
  if(DumpFlag != 2)
    {
      header.npart[2] = 0;
      header.npartTotal[2] = 0;
      header.npartTotalHighWord[2] = 0;
    }
#endif

  if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    header.flag_ic_info = FLAG_EVOLVED_2LPT;

  if(header.flag_ic_info == FLAG_ZELDOVICH_ICS)
    header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;

  if(header.flag_ic_info == FLAG_NORMALICS_2LPT)
    header.flag_ic_info = FLAG_EVOLVED_2LPT;

  if(header.flag_ic_info == 0 && All.ComovingIntegrationOn != 0)
    header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];

#ifdef COSMIC_RAYS
  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
    header.SpectralIndex_CR_Pop[CRpop] = All.CR_Alpha[CRpop];
#endif


  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.flag_entropy_instead_u = 0;

#ifdef COOLING
  header.flag_cooling = 1;
#endif

#ifdef SFR
  header.flag_sfr = 1;
  header.flag_feedback = 1;
#ifdef STELLARAGE
  header.flag_stellarage = 1;
#endif
#ifdef METALS
  header.flag_metals = 1;
#endif
#ifdef LT_STELLAREVOLUTION
  //header.flag_stellarevolution = 2;
#endif
#endif

  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  header.flag_doubleprecision = 1;
#else
  header.flag_doubleprecision = 0;
#endif

#ifdef LT_STELLAREVOLUTION
  for(n = 0; n < LT_NMet; n++)
    {
      header.names[n][0] = MetNames[n][0];
      header.names[n][1] = MetNames[n][1];
    }
  for(n = LT_NMet; n < 15; n++)
    {
      header.names[n][0] = ' ';
      header.names[n][1] = ' ';
    }
#endif

  /* open file and write header */

#ifdef WRITE_KEY_FILES
  key_header = header;
  for(n = 0; n < 6; n++)
    {
      key_header.npart[n] = NKeys_file[n];
      key_header.npartTotal[n] = (unsigned int) NKeysTotal[n];
      key_header.npartTotalHighWord[n] = (unsigned int) (NKeysTotal[n] >> 32);
    }
  for(n = 0; n < 3; n++)
    key_header.mass[n] = DomainCorner[n];
  key_header.mass[3] = DomainFac;
  for(n = 4; n < 6; n++)
    key_header.mass[n] = 0;
  key_header.flag_feedback = BITS_PER_DIMENSION_SAVE_KEYS;
#endif

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  sprintf(buf, "%s.hdf5", fname);
	  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
		}
	    }

	  write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      endrun(123);
	    }

	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
#ifdef WRITE_KEY_FILES
	  sprintf(buf_key, "%s.key", fname);
	  FdKey = fopen(buf_key, "w");
	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP_KEY;
	      my_fwrite((void *) "HEAD", sizeof(char), 4, FdKey);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, FdKey);
	      SKIP_KEY;
	    }

	  blksize = sizeof(key_header);
	  SKIP_KEY;
	  my_fwrite(&key_header, sizeof(header), 1, FdKey);
	  SKIP_KEY;
#endif
	}
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
      char fname_pot[1000], *s_split_1, *s_split_2;

      s_split_1 = strtok(fname, ".");
      s_split_2 = strtok(NULL, " ,");

      sprintf(fname_pot, "%s_pot.%s", s_split_1, s_split_2);
      if(!(fd_pot = fopen(fname_pot, "w")))
	{
	  printf("can't open file `%s' for writing pot_snapshot.\n", fname_pot);
	  endrun(1234);
	}
#endif
#ifdef LT_SEvDbg
      sprintf(buf, "%s.dbg", fname);
      fd_dbg = fopen(buf, "w");
#endif
#ifdef LT_TRACK_CONTRIBUTES
      sprintf(buf, "%s.trck", fname);
      if((FdTrck = fopen(buf, "w")) == 0x0)
	{
	  printf("Task %d : it was impossible to open file %s\n", ThisTask, buf);
	  endrun(125);
	}
      blksize = 6 * sizeof(int);
      SKIP_TRCK;

      my_fwrite(&header.npart[0], sizeof(int), 1, FdTrck);
      my_fwrite(&header.npart[4], sizeof(int), 1, FdTrck);
      my_fwrite(&PowerBase, sizeof(int), 1, FdTrck);
      blksize = LT_Nbits;
      my_fwrite(&blksize, sizeof(int), 1, FdTrck);
      blksize = LT_power10_Nbits;
      my_fwrite(&blksize, sizeof(int), 1, FdTrck);
      blksize = SFs_dim;
      my_fwrite(&blksize, sizeof(int), 1, FdTrck);

      blksize = 6 * sizeof(int);
      SKIP_TRCK;
      fd_trck_back = fd;
#endif
#ifdef LT_SMOOTH_Z
      sprintf(buf, "%s.zsmooth", fname);
      if((FdZsmooth = fopen(buf, "w")) == 0x0)
	{
	  printf("Task %d : it was impossible to open file %s\n", ThisTask, buf);
	  endrun(125);
	}
#ifndef LT_SMOOTHZ_IN_IMF_SWITCH
      blksize = 4;
      SKIP_ZSMOOTH;
      my_fwrite(&header.npart[0], sizeof(int), 1, FdZsmooth);
#else
      blksize = 8;
      SKIP_ZSMOOTH;
      my_fwrite(&header.npart[0], sizeof(int), 1, FdZsmooth);
      my_fwrite(&header.npart[4], sizeof(int), 1, FdZsmooth);
      blksize = 8;
#endif
      SKIP_ZSMOOTH;
      fd_zsmooth_back = fd;
#endif
    }

  if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
      n_info = 0;
      InfoBlock = (struct info_block *) mymalloc("InfoBlock", sizeof(struct info_block) * 1000);

      for(bnr = 0; bnr < 1000; bnr++)
	{
	  blocknr = (enum iofields) bnr;

	  if(blocknr == IO_SECONDORDERMASS)
	    continue;

	  if(blocknr == IO_LASTENTRY)
	    break;

	  if(blockpresent(blocknr))
	    {
	      bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
	      npart = get_particles_in_block(blocknr, &typelist[0]);

	      if(npart > 0)
		{
		  for(type = 0; type < 6; type++)
		    InfoBlock[n_info].is_present[type] = typelist[type];
		  InfoBlock[n_info].ndim = get_values_per_blockelement(blocknr);
		  get_Tab_IO_Label(blocknr, label);
		  for(i = 0; i < 4; i++)
		    InfoBlock[n_info].label[i] = label[i];
		  switch (get_datatype_in_block(blocknr))
		    {
		    case 0:
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "LONG    ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "LONGN   ", 8);
		      break;
		    case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
		      break;
#else
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
		      break;
#endif
		    case 2:
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "LLONG   ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "LLONGN  ", 8);
		      break;
		    }
		  n_info++;
		}
	    }
	}
#ifdef WRITE_INFO_BLOCK
      if(All.SnapFormat == 2)
	{
	  blksize = sizeof(int) + 4 * sizeof(char);
	  SKIP;
	  my_fwrite((void *) "INFO", sizeof(char), 4, fd);
	  nextblock = n_info * sizeof(struct info_block) + 2 * sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  SKIP;
	}
      blksize = n_info * sizeof(struct info_block);
      SKIP;
      my_fwrite(InfoBlock, sizeof(struct info_block), n_info, fd);
      SKIP;
#endif
    }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_SECONDORDERMASS)
	continue;

      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);

	  blockmaxlen = (size_t) ((All.BufferSize * 1024 * 1024) / bytes_per_blockelement);

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == 0)
		{
		  char buf[1000];

		  get_dataset_name(blocknr, buf);
		  printf("writing block %d (%s)...\n", bnr, buf);
		}

	      if(ThisTask == writeTask)
		{

		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
#if defined(LT_TRACK_CONTRIBUTES)
		      if(blocknr == IO_CONTRIB)
			fd = FdTrck;
#endif
#if defined(LT_SMOOTH_Z)
#if !defined(LT_SMOOTH_ALLMETALS)
		      if(blocknr == IO_ZSMOOTH)
#else
		      if(blocknr == IO_allZSMOOTH)
#endif
			fd = FdZsmooth;
#endif

		      if(All.SnapFormat == 2)
			{
			  blksize = sizeof(int) + 4 * sizeof(char);
			  SKIP;
			  get_Tab_IO_Label(blocknr, label);
			  my_fwrite(label, sizeof(char), 4, fd);
			  nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
			  my_fwrite(&nextblock, sizeof(int), 1, fd);
			  SKIP;
			}
		      blksize = npart * bytes_per_blockelement;
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
		      if(blocknr == IO_POT)
			SKIP_POT;
		      if(blocknr != IO_POT)
			SKIP;
#else
		      SKIP;
#endif
		    }
		}

	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {
#ifdef HAVE_HDF5
#if defined(LT_TRACK_CONTRIBUTES) || defined(LT_SMOOTH_Z)
		      if(blocknr != IO_CONTRIB && blocknr != IO_ZSMOOTH && blocknr != IO_allZSMOOTH)
#endif
			if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			  {
			    switch (get_datatype_in_block(blocknr))
			      {
			      case 0:
				hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
				break;
			      case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
				hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
				hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
				break;
			      case 2:
				hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
				break;
			      }

			    dims[0] = header.npart[type];
			    dims[1] = get_values_per_blockelement(blocknr);
			    if(dims[1] == 1)
			      rank = 1;
			    else
			      rank = 2;

			    get_dataset_name(blocknr, buf);

			    hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
			    hdf5_dataset =
			      H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file,
					H5P_DEFAULT);

			    pcsum = 0;
			  }
#endif

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > (int) blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_write_buffer(blocknr, &offset, pc, type);

			      if(ThisTask == writeTask && task != writeTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask != writeTask && task == ThisTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					  TAG_PDATA, MPI_COMM_WORLD);

			      if(ThisTask == writeTask)
				{
				  if(All.SnapFormat == 3)
				    {
#ifdef HAVE_HDF5
#if defined(LT_TRACK_CONTRIBUTES) || defined(LT_SMOOTH_Z)
				      if(blocknr != IO_CONTRIB && blocknr != IO_ZSMOOTH
					 && blocknr != IO_allZSMOOTH)
					{
#endif
					  start[0] = pcsum;
					  start[1] = 0;

					  count[0] = pc;
					  count[1] = get_values_per_blockelement(blocknr);
					  pcsum += pc;

					  H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							      start, NULL, count, NULL);

					  dims[0] = pc;
					  dims[1] = get_values_per_blockelement(blocknr);
					  hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

					  hdf5_status =
					    H5Dwrite(hdf5_dataset, hdf5_datatype,
						     hdf5_dataspace_memory,
						     hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

					  H5Sclose(hdf5_dataspace_memory);
#if defined(LT_TRACK_CONTRIBUTES) || defined(LT_SMOOTH_Z)
					}
#endif
#endif
				    }
				  else
				    {
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
				      if(blocknr == IO_POT)
					my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd_pot);
				      else
					my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
#else
				      my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
#endif
				    }
				}

			      n_for_this_task -= pc;
			    }
			}

#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  if(All.SnapFormat == 3)
			    {
			      H5Dclose(hdf5_dataset);
			      H5Sclose(hdf5_dataspace_in_file);
			      H5Tclose(hdf5_datatype);
			    }
			}
#endif
		    }
		}

	      if(ThisTask == writeTask)
		{
		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
		      if(blocknr == IO_POT)
			SKIP_POT;
		      if(blocknr != IO_POT)
			SKIP;
#else
		      SKIP;
#endif
		    }
#if defined(LT_TRACK_CONTRIBUTES)
		  if(blocknr == IO_CONTRIB)
		    fd = fd_trck_back;
#endif
#if defined(LT_SMOOTH_Z)
#if !defined(LT_SMOOTH_ALLMETALS)
		  if(blocknr == IO_ZSMOOTH)
#else
		  if(blocknr == IO_allZSMOOTH)
#endif
		    fd = fd_zsmooth_back;
#endif
		}
	    }
	}
    }

  if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
#ifndef WRITE_INFO_BLOCK
      if(All.SnapFormat == 2)
	{
	  blksize = sizeof(int) + 4 * sizeof(char);
	  SKIP;
	  my_fwrite((void *) "INFO", sizeof(char), 4, fd);
	  nextblock = n_info * sizeof(struct info_block) + 2 * sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  SKIP;
	}
      blksize = n_info * sizeof(struct info_block);
      SKIP;
      my_fwrite(InfoBlock, sizeof(struct info_block), n_info, fd);
      SKIP;
#endif
      myfree(InfoBlock);
    }

#ifdef WRITE_KEY_FILES
  if(ThisTask == writeTask)
    {
      KeyInfoBlock = (struct info_block *) mymalloc("KeyInfoBlock", sizeof(struct info_block) * 3);

      for(bnr = 0; bnr < 3; bnr++)
	{
	  switch (bnr)
	    {
	    case 0:
	      strncpy(label, "KEY ", 4);
	      break;
	    case 1:
	      strncpy(label, "NKEY", 4);
	      break;
	    case 2:
	      strncpy(label, "OKEY", 4);
	      break;
	    default:
	      printf("Undefined block to write to KEY file ...\n");
	      endrun(76112541);
	    }
	  strncpy(KeyInfoBlock[bnr].label, label, 4);

	  if(bnr == 0)
	    if(BITS_PER_DIMENSION_SAVE_KEYS > 10)
	      strncpy(KeyInfoBlock[bnr].type, "LLONG   ", 8);
	    else
	      strncpy(KeyInfoBlock[bnr].type, "LONG    ", 8);
	  else
	    strncpy(KeyInfoBlock[bnr].type, "LONG    ", 8);

	  for(type = 0; type < 6; type++)
	    if(NKeys_file[type] > 0)
	      KeyInfoBlock[bnr].is_present[type] = 1;
	    else
	      KeyInfoBlock[bnr].is_present[type] = 0;

	  KeyInfoBlock[bnr].ndim = 1;
	}

      if(All.SnapFormat == 2)
	{
	  blksize = sizeof(int) + 4 * sizeof(char);
	  SKIP_KEY;
	  my_fwrite((void *) "INFO", sizeof(char), 4, FdKey);
	  nextblock = 3 * sizeof(struct info_block) + 2 * sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, FdKey);
	  SKIP_KEY;
	}
      blksize = 3 * sizeof(struct info_block);
      SKIP_KEY;
      my_fwrite(KeyInfoBlock, sizeof(struct info_block), 3, FdKey);
      SKIP_KEY;

    }

  if(npart > 0)
    for(bnr = 0; bnr < 3; bnr++)
      {
	switch (bnr)
	  {
	  case 0:
	    data_point_key = KeyIndex;
	    if(BITS_PER_DIMENSION_SAVE_KEYS > 10)
	      size_of_element = sizeof(peanokey);
	    else
	      size_of_element = sizeof(int);
	    break;
	  case 1:
	    data_point_int = NPartPerKey;
	    size_of_element = sizeof(int);
	    break;
	  case 2:
	    data_point_int = PartKeyOffset;
	    size_of_element = sizeof(int);
	    break;
	  default:
	    printf("Undefined block to write to KEY file ...\n");
	    endrun(76112541);
	  }

	if(ThisTask == 0)
	  printf("writing block (%c%c%c%c) ...\n", label[0], label[1], label[2], label[3]);

	if(ThisTask == writeTask)
	  {
	    strncpy(label, KeyInfoBlock[bnr].label, 4);

	    if(All.SnapFormat == 2)
	      {
		blksize = sizeof(int) + 4 * sizeof(char);
		SKIP_KEY;
		my_fwrite(label, sizeof(char), 4, FdKey);
		nextblock = TotNumKeysFile * size_of_element + 2 * sizeof(int);
		my_fwrite(&nextblock, sizeof(int), 1, FdKey);
		SKIP_KEY;
	      }
	    blksize = TotNumKeysFile * size_of_element;
	    SKIP_KEY;
	  }

	if(ThisTask == 0)
	  printf(" ... %d bytes = %u x %d (%g MB) ...\n", blksize, TotNumKeysFile, size_of_element,
		 blksize / 1024. / 1024.);

	for(type = 0; type < 6; type++)
	  {
	    for(task = writeTask; task <= lastTask; task++)
	      {
		for(n = 0, local_offset = 0; n < type; n++)
		  local_offset += NKeys[n];

		if(task == ThisTask)
		  {
		    n_for_this_task = NKeys[type];

		    for(p = writeTask; p <= lastTask; p++)
		      if(p != ThisTask)
			MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
		  }
		else
		  MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

		while(n_for_this_task > 0)
		  {
		    pc = n_for_this_task;

		    if(pc > blockmaxlen)
		      pc = blockmaxlen;

		    if(ThisTask == task)
		      {

			comm_point_int = (int *) CommBuffer;
			comm_point_key = (peanokey *) CommBuffer;

			switch (bnr)
			  {
			  case 0:
			    for(n = 0; n < pc; n++)
			      if(BITS_PER_DIMENSION_SAVE_KEYS > 10)
				*comm_point_key++ = data_point_key[local_offset++];
			      else
				*comm_point_int++ = (int) data_point_key[local_offset++];
			    break;
			  case 1:
			    for(n = 0; n < pc; n++)
			      *comm_point_int++ = data_point_int[local_offset++];
			    break;
			  case 2:
			    for(n = 0; n < pc; n++)
			      *comm_point_int++ = data_point_int[local_offset++] + write_offset[type];
			    break;
			  }
		      }

		    if(ThisTask == writeTask && task != writeTask)
		      MPI_Recv(CommBuffer, size_of_element * pc, MPI_BYTE, task,
			       TAG_PDATA, MPI_COMM_WORLD, &status);

		    if(ThisTask != writeTask && task == ThisTask)
		      MPI_Ssend(CommBuffer, size_of_element * pc, MPI_BYTE, writeTask,
				TAG_PDATA, MPI_COMM_WORLD);

		    if(ThisTask == writeTask)
		      my_fwrite(CommBuffer, size_of_element, pc, FdKey);

		    n_for_this_task -= pc;
		  }
	      }
	  }
	if(ThisTask == writeTask)
	  SKIP_KEY;
      }

  if(ThisTask == writeTask)
    {
      myfree(KeyInfoBlock);
      fclose(FdKey);
    }
#endif


  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Gclose(hdf5_headergrp);
	  H5Fclose(hdf5_file);
#endif
#ifdef LT_SEvDbg
	  fclose(fd_dbg);
#endif
#if defined(LT_TRACK_CONTRIBUTES)
	  if(FdTrck != 0x0)
	    fclose(FdTrck);
#endif
#if defined(LT_SMOOTH_Z)
	  if(FdZsmooth != 0x0)
	    fclose(FdZsmooth);
#endif
	}
      else
	{
	  fclose(fd);
	}
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
      fclose(fd_pot);
#endif
    }
}




#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };

#ifdef COSMIC_RAYS
  hsize_t adim_alpha[1] = { NUMCRPOP };
#endif
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

#ifdef COSMIC_RAYS
  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim_alpha, NULL);
  hdf5_attribute = H5Acreate(handle, "SpectralIndex_CR_Pop", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.SpectralIndex_CR_Pop);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
#endif


  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}
#endif





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
	  printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
	  fflush(stdout);
	  endrun(777);
	}
    }
  else
    nwritten = 0;

  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(size * nmemb == 0)
    return 0;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
      else
	printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(778);
    }
  return nread;
}

/** This is so we don't have to preface every printf call with the if
    statement to only make it print once. */
void mpi_printf(const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vprintf(fmt, l);
      fflush(stdout);
      va_end(l);
    }
}

#if defined(ORDER_SNAPSHOTS_BY_ID) || defined(SUBFIND_READ_FOF) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
int io_compare_P_ID(const void *a, const void *b)
{
  if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID))
    return -1;

  if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID))
    return +1;

  return 0;
}

int io_compare_P_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
    return -1;

  if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
    return +1;

  if(((struct particle_data *) a)->SubNr < (((struct particle_data *) b)->SubNr))
    return -1;

  if(((struct particle_data *) a)->SubNr > (((struct particle_data *) b)->SubNr))
    return +1;

  return 0;
}

int io_compare_P_GrNr_ID(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
    return -1;

  if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
    return +1;

  if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID))
    return -1;

  if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID))
    return +1;

  return 0;
}

#endif

#ifdef LT_SEv_INFO_DETAILS
void write_IDZ(void)
{
  char buf[500];

  FILE *myFile;

  fpos_t IDZpos;

  float myZ;

  int j, totn_written = 0;

  unsigned long long int i;

  for(j = 0; j < NTask; j++)
    {
      if(ThisTask == j && n_type[4] > 0)
	{
	  sprintf(buf, "%s%s%03d", All.OutputDir, "IDZ.dat", ThisTask);
	  if((myFile = fopen(buf, "w")) == 0x0)
	    {
	      printf("error in opening file '%s'\n", buf);
	      endrun(1);
	    }
	  if(totn_written == 0)
	    {
	      fwrite(&ntot_type_all[4], sizeof(long long), 1, myFile);
	      totn_written = 1;
	    }

	  for(i = 0; i < NumPart; i++)
	    if(P[i].Type == 4)
	      {
		myZ = (float) get_metallicity(i, -1);
#ifndef LONGIDS
		fwrite(&P[i].ID, sizeof(int), 1, myFile);
#else
		fwrite(&P[i].ID, sizeof(long long), 1, myFile);
#endif
		fwrite(&MetP[P[i].pt.MetID].iMass, sizeof(float), 1, myFile);
		fwrite(&myZ, sizeof(float), 1, myFile);
	      }
	  fclose(myFile);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }

  return;
}

#endif
#if ( (defined(LT_SMOOTH_Z) && !defined(LT_SMOOTH_ALLMETALS)) || defined(LT_METAL_COOLING)) && defined(LT_LOG_SMOOTH)

void Spy_Z(int num)
{
  int nwrite;
  int TaskStart, TaskStop;
  int i, k, j, b, ndata;
  unsigned int Ndone, *Ndone_all;
  unsigned long long int NdoneTot;
  double a3inv;
  double T, Z, rho, u, mu, Ne;
  double Ls, Ls_a, Ls_b, Lt, L0, Zs, Zs_a, Zs_b, Zt, metal_mass;
  float *buffer;

#define BUNCH 500000

  FILE *file;
  char fname[500];
  char head[1000];

  ignore_failure_in_convert_u = 1;

  a3inv = 1 / (All.Time * All.Time * All.Time);

  XH = HYDROGEN_MASSFRACTION;
  yhelium = (1 - XH) / (4 * XH);

  nwrite = (NTask / All.NumFilesWrittenInParallel) + (NTask % All.NumFilesWrittenInParallel);


  if((buffer = (float *) malloc(BUNCH * sizeof(float))) != 0x0)
    {


      for(k = 0; k < nwrite; k++)
	{

	  TaskStart = k * All.NumFilesWrittenInParallel;
	  TaskStop = (k + 1) * All.NumFilesWrittenInParallel;
	  if(TaskStop > NTask)
	    TaskStop = NTask;

	  if(ThisTask >= TaskStart && ThisTask < TaskStop)
	    {

	      Ndone = 0;
	      for(i = 0; i < NumPart; i++)
		if(P[i].Type == 0 && SphP[i].a2.Density > 0)
		  {
		    metal_mass = get_metalmass(SphP[i].Metals);
#ifndef LT_SMOOTH_Z
		    if(metal_mass > 0)
#else
		    if((metal_mass > 0) || (SphP[i].Zsmooth > 0))
#endif
		      Ndone++;
		  }

	      if(Ndone > 0)
		continue;

	      sprintf(fname, "zsmooth.%d.%d.txt", num, ThisTask);
	      if((file = fopen(fname, "w")) == 0x0)
		{
		  free(buffer);
		  ignore_failure_in_convert_u = 0;
		  return;
		}

	      sprintf(head, "temp dens Z zeroZ_cool ");
	      ndata = 4;
#ifdef LT_METAL_COOLING
	      sprintf(&head[strlen(head)], "Zcool ");
	      ndata++;
#endif

#ifdef LT_SMOOTH_Z
	      sprintf(&head[strlen(head)], "smoothZ ");
	      ndata++;
#ifdef LT_SMOOTH_SIZE
	      sprintf(&head[strlen(head)], "sZ_Nngb ");
	      ndata++;
#endif
#ifdef LT_METAL_COOLING
	      sprintf(&head[strlen(head)], "sZ_cool ");
	      ndata++;
#endif
#endif
	      i = strlen(head);
	      fwrite(&i, sizeof(int), 1, file);
	      fprintf(file, "%s", head);

	      b = 0;

	      for(i = 0; i < NumPart; i++)
		if(P[i].Type == 0 && SphP[i].a2.Density > 0)
		  {
		    metal_mass = get_metalmass(SphP[i].Metals);
#ifndef LT_SMOOTH_Z
		    if(metal_mass == 0)
		      continue;
#else
		    if((metal_mass == 0) && (SphP[i].Zsmooth == 0))
		      continue;
#endif

		    Ne = (double) SphP[i].Ne;

		    mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[i].Ne);

		    u = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1) *
		      All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
		    rho =
		      SphP[i].a2.Density * a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

		    T = convert_u_to_temp(u, rho, &Ne);

		    if((Zt = get_metallicity(i, -1)) / 0.02 < 1e-8)
		      Zt = 1.77e-3 * 1e-8;
#ifdef LT_SMOOTH_Z
		    Zs = SphP[i].Zsmooth;
#endif

#ifdef LT_METAL_COOLING
		    Lt = CoolingRateFromU(u, rho, &Ne, get_metallicity_solarunits(get_metallicity(i, Iron)));
#endif
#ifdef LT_SMOOTH_Z
		    if(metal_mass > 0)
		      {
			Ls =
			  CoolingRateFromU(u, rho, &Ne,
					   get_metallicity_solarunits(Zs * SphP[i].Metals[Iron] /
								      metal_mass));
		      }
		    else
		      Ls = Ls_a = Ls_b = 0;
#endif
		    L0 = CoolingRateFromU(u, rho, &Ne, NO_METAL);

		    /*
		       Ls = CoolingRateGetMetalLambda(T, get_metallicity_solarunits(SphP[i].Zsmooth));
		       Lt = GetMetalLambda(T, get_metallicity(i, -1));
		     */

		    buffer[b++] = (float) T;
		    buffer[b++] = (float) SphP[i].a2.Density * a3inv / SFs[0].PhysDensThresh[0];
		    buffer[b++] = (float) Zt;
		    buffer[b++] = (float) L0;
#ifdef LT_METAL_COOLING
		    buffer[b++] = (float) Lt;
#endif
#ifdef LT_SMOOTH_Z
		    buffer[b++] = (float) Zs;
#ifdef LT_SMOOTH_SIZE
		    buffer[b++] = (float) SphP[i].SmoothNgb;
#endif
#ifdef LT_METAL_COOLING
		    buffer[b++] = (float) Ls;
#endif
#endif

		    if(b > BUNCH - ndata)
		      {
			fwrite(buffer, sizeof(float), b, file);
			b = 0;
		      }
		  }
	      if(b > 0)
		{
		  fwrite(buffer, sizeof(float), b, file);
		  b = 0;
		}
	      fclose(file);
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }


  free(buffer);
  if(ThisTask == 0)
    Ndone_all = (unsigned int *) calloc(NTask, sizeof(unsigned));
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Gather(&Ndone, 1, MPI_INT, Ndone_all, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      NdoneTot = Ndone;
      for(i = 1; i < NTask; i++)
	NdoneTot += Ndone_all[i];
      free(Ndone_all);
      printf("done for %llu particles\n", NdoneTot);
    }

  ignore_failure_in_convert_u = 0;
  return;
}

#endif
