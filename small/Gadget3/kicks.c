#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

void apply_long_range_kick(integertime, integertime);

void do_first_halfstep_kick(void)
{
  int i;
  integertime ti_step, tstart, tend;

#if defined (VS_TURB) || defined (AB_TURB) || defined (TURB_DRIVING)
  do_turb_driving_step_first_half();
#endif

#ifdef PMGRID
  if(All.PM_Ti_begstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
      tstart = All.PM_Ti_begstep;
      tend = tstart + ti_step / 2;

      apply_long_range_kick(tstart, tend);
    }
#endif

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  int il;
#pragma omp parallel for private(i)
  for(il = 0; il < NActivePart; il++)
    {
      i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
      ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;

      tstart = P[i].Ti_begstep;	/* beginning of step */
      tend = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */

#ifdef BLACK_HOLES
      if(P[i].Mass > 0)
#endif
	do_the_kick(i, tstart, tend, P[i].Ti_current);
    }
}

void do_second_halfstep_kick(void)
{
  int i, j;
  integertime ti_step, tstart, tend;

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
      tstart = All.PM_Ti_begstep + ti_step / 2;
      tend = tstart + ti_step / 2;

      apply_long_range_kick(tstart, tend);
    }
#endif
#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  int il;
#pragma omp parallel for private(i)
  for(il = 0; il < NActivePart; il++)
    {
      i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
      ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;

      tstart = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
      tend = P[i].Ti_begstep + ti_step;	/* end of step */

#ifdef BLACK_HOLES
      if(P[i].Mass > 0)
#endif
	do_the_kick(i, tstart, tend, P[i].Ti_current);

      /* after completion of a full step, set the predication values of SPH quantities
       * to the current values. They will then predicted along in drift operations
       */
#ifdef BLACK_HOLES
      if(P[i].Type == 0 && P[i].Mass > 0)
#else
      if(P[i].Type == 0)
#endif
	{
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] = P[i].Vel[j];

	  SphP[i].EntropyPred = SphP[i].Entropy;

#ifdef EOS_DEGENERATE
	  for(j = 0; j < 3; j++)
	    SphP[i].xnucPred[j] = SphP[i].xnuc[j];
#endif

	  SphP[i].Pressure = get_pressure(i);

	  set_predicted_sph_quantities_for_extra_physics(i);
	}
    }

#if defined (VS_TURB) || defined (AB_TURB) || defined (TURB_DRIVING)
  do_turb_driving_step_second_half();
#endif
}


#ifdef PMGRID
void apply_long_range_kick(integertime tstart, integertime tend)
{
  int i, j;
  double dt_gravkick, dvel[3];

  if(All.ComovingIntegrationOn)
    dt_gravkick = get_gravkick_factor(tstart, tend);
  else
    dt_gravkick = (tend - tstart) * All.Timebase_interval;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < NumPart; i++)
    {
#ifdef BLACK_HOLES
      if(P[i].Mass > 0)
#endif
	for(j = 0; j < 3; j++)	/* do the kick, only collisionless particles */
	  {
	    dvel[j] = P[i].GravPM[j] * dt_gravkick;
#ifndef SIDM_FREEZE
	    P[i].Vel[j] += dvel[j];
	    P[i].dp[j] += P[i].Mass * dvel[j];
#endif
	  }
#ifdef STATICBRANDT_KICK
      P[i].Vel[0] = -OmegaR(i, BRANDT_MODE) * (P[i].Pos[1] - LONG_Y / 2.);
      P[i].Vel[1] = OmegaR(i, BRANDT_MODE) * (P[i].Pos[0] - LONG_X / 2.);
      P[i].Vel[2] = 0.0;
      P[i].dp[0] = P[i].Vel[0] * P[i].Mass;
      P[i].dp[1] = P[i].Vel[1] * P[i].Mass;
      P[i].dp[2] = P[i].Vel[2] * P[i].Mass;
#endif

#ifdef DISTORTIONTENSORPS
      do_long_range_phase_space_kick(i, dt_gravkick);
#endif
    }
}
#endif


void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent)
{
  int j;
  double dv[3];
  double dt_entr, dt_gravkick, dt_hydrokick;

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (tend - tstart) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(tstart, tend);
      dt_hydrokick = get_hydrokick_factor(tstart, tend);
    }
  else
    {
      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
    }

  /* do the kick */
  for(j = 0; j < 3; j++)
    {
      dv[j] = P[i].g.GravAccel[j] * dt_gravkick;
#ifdef RELAXOBJECT
      dv[j] -= P[i].Vel[j] * All.RelaxFac * dt_gravkick;
#endif
#ifndef SIDM_FREEZE
      P[i].Vel[j] += dv[j];
      P[i].dp[j] += P[i].Mass * dv[j];
#endif
    }
#ifdef STATICBRANDT_KICK
  P[i].Vel[0] = -OmegaR(i, BRANDT_MODE) * (P[i].Pos[1] - LONG_Y / 2.);
  P[i].Vel[1] = OmegaR(i, BRANDT_MODE) * (P[i].Pos[0] - LONG_X / 2.);
  P[i].Vel[2] = 0.0;
  P[i].dp[0] = P[i].Vel[0] * P[i].Mass;
  P[i].dp[1] = P[i].Vel[1] * P[i].Mass;
  P[i].dp[2] = P[i].Vel[2] * P[i].Mass;
#endif

#ifdef DISTORTIONTENSORPS
  do_the_phase_space_kick(i, dt_gravkick);
#endif


  if(P[i].Type == 0)		/* kick for SPH quantities */
    {
      for(j = 0; j < 3; j++)
	{
	  dv[j] = SphP[i].a.HydroAccel[j] * dt_hydrokick;
#ifdef RT_RAD_PRESSURE
	  dv[j] += SphP[i].RadAccel[j] * dt_hydrokick;
#endif
#ifndef STATICBRANDT_KICK
	  P[i].Vel[j] += dv[j];
	  P[i].dp[j] += P[i].Mass * dv[j];
#endif
	}
#ifdef STATICBRANDT_KICK
      P[i].Vel[0] = -OmegaR(i, BRANDT_MODE) * (P[i].Pos[1] - LONG_Y / 2.);
      P[i].Vel[1] = OmegaR(i, BRANDT_MODE) * (P[i].Pos[0] - LONG_X / 2.);
      P[i].Vel[2] = 0.0;
      P[i].dp[0] = P[i].Vel[0] * P[i].Mass;
      P[i].dp[1] = P[i].Vel[1] * P[i].Mass;
      P[i].dp[2] = P[i].Vel[2] * P[i].Mass;
#endif

      double dEntr = SphP[i].e.DtEntropy * dt_entr;

#if defined(EOS_DEGENERATE)
      dEntr *= All.UnitTime_in_s;
#endif

      SphP[i].Entropy = DMAX(SphP[i].Entropy + dEntr, 0.5 * SphP[i].Entropy);

#ifdef JD_RELAX_CLUSTERS
      SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].U / pow(SphP[i].d.Density, GAMMA_MINUS1);
#endif

      check_particle_for_temperature_minimum(i);

      do_sph_kick_for_extra_physics(i, tstart, tend, dt_entr);
    }
}

void set_predicted_sph_quantities_for_extra_physics(int i)
{
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
  int j1;
  for(j1 = 0; j1 < 3; j1++)
    SphP[i].b2.BPred[j1] = SphP[i].b1.B[j1];
#endif
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
  SphP[i].PhiPred = SphP[i].Phi;
#endif
}


void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
  int j;

#ifdef MAGNETIC
  double dt_mag;
  if(All.ComovingIntegrationOn)
    dt_mag = get_magkick_factor(tstart, tend);
  else
    dt_mag = (tend - tstart) * All.Timebase_interval;
#endif

  for(j = 0; j < 3; j++)
    {
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
      SphP[i].b1.B[j] += SphP[i].DtB[j] * dt_mag;
#endif
    }

#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
  SphP[i].Phi += SphP[i].DtPhi * dt_mag;
#endif
#ifdef TIME_DEP_ART_VISC
  SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
  SphP[i].alpha = DMIN(SphP[i].alpha, All.ArtBulkViscConst);
  if(SphP[i].alpha < All.AlphaMin)
    SphP[i].alpha = All.AlphaMin;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
  SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
#ifdef VORONOI_RELAX_VIA_VISC
  if(SphP[i].alpha < All.ArtBulkViscConst / 128.0 / 128.0)
    SphP[i].alpha = All.ArtBulkViscConst / 128.0 / 128.0;
#else
  if(SphP[i].alpha < All.ArtBulkViscConst / 128.0)
    SphP[i].alpha = All.ArtBulkViscConst / 128.0;
#endif
#endif
#ifdef TIME_DEP_MAGN_DISP
  SphP[i].Balpha += SphP[i].DtBalpha * dt_entr;
  SphP[i].Balpha = DMIN(SphP[i].Balpha, All.ArtMagDispConst);
  if(SphP[i].Balpha < All.ArtMagDispMin)
    SphP[i].Balpha = All.ArtMagDispMin;
#endif

#ifdef NUCLEAR_NETWORK
  for(j = 0; j < EOS_NSPECIES; j++)
    SphP[i].xnuc[j] += SphP[i].dxnuc[j] * dt_entr * All.UnitTime_in_s;

  network_normalize(SphP[i].xnuc, &SphP[i].Entropy, &All.nd, &All.nw);
#endif

#ifdef CHEMISTRY
  /* update the chemical abundances for the new density and temperature */
  double a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
  double a_end = All.TimeBegin * exp(tend * All.Timebase_interval);
  int mode;
  /* time in cosmic expansion parameter */
  compute_abundances(mode = 1, i, a_start, a_end);
#endif
}
