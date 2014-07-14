#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */


void compute_grav_accelerations(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  if(ThisTask == 0)
    {
      printf("Start gravity force computation...\n");
      fflush(stdout);
    }


#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      long_range_force();

      CPU_Step[CPU_MESH] += measure_time();
    }
#endif


#ifndef ONLY_PM

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
#endif


  if(All.Ti_Current == 0 && RestartFlag == 0 && header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    second_order_ics();		/* produces the actual ICs from the special second order IC file */


#ifdef FORCETEST
  gravity_forcetest();
#endif

  if(ThisTask == 0)
    {
      printf("gravity force computation done.\n");
      fflush(stdout);
    }
}



void compute_densities(void)
{
  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

#ifdef CS_MODEL
      CPU_Step[CPU_MISC] += measure_time();

#if defined(CS_SNI) || defined(CS_SNII)
      cs_flag_SN_starparticles();	/* mark SNI star particles */
#endif
      cs_copy_densities();
      cs_find_low_density_tail();

      CPU_Step[CPU_CSMISC] += measure_time();
#endif


#ifndef VORONOI
#ifdef FS_ALFA2_QUENCH
      compute_global_quantities_of_system();
#endif
      density();		/* computes density, and pressure */
#else
      voronoi_mesh();
      voronoi_setup_exchange();

      voronoi_density();
#endif


#if (defined(SMOOTH_DIVB)  || defined(SMOOTH_ROTB) || defined(BSMOOTH) || (defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)) || defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_XCLD) || defined(LT_TRACK_WINDS) || defined(BSMOOTH_TIME)) || defined(LT_MV_CHEMICALDIFFUSION) || defined(VSMOOTH)
      smoothed_values();
#endif


#if defined(SNIA_HEATING)
      snIa_heating();
#endif


#if defined(CS_MODEL) && defined(CS_FEEDBACK)

      CPU_Step[CPU_CSMISC] += measure_time();

      cs_find_hot_neighbours();

      cs_promotion();
      cs_copy_densities();
      CPU_Step[CPU_CSMISC] += measure_time();
      density();		/* recalculate densities again */
      CPU_Step[CPU_CSMISC] += measure_time();
#endif

    }
}

void compute_hydro_accelerations(void)
{
  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

#ifndef VORONOI
      force_update_hmax();	/* update smoothing lengths in tree */
      hydro_force();		/* adds hydrodynamical accelerations  and computes du/dt  */
#else
      voronoi_hydro_force();
#endif

#ifdef WINDTUNNEL
      set_injection_accel();
#endif

      if(ThisTask == 0)
	{
	  printf("hydro force computation done.\n");
	  fflush(stdout);
	}
    }
}
