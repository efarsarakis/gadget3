#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"


void reconstruct_timebins(void)
{
  int i, bin;
  long long glob_sum;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinCount[bin] = 0;
      TimeBinCountSph[bin] = 0;
      FirstInTimeBin[bin] = -1;
      LastInTimeBin[bin] = -1;
#ifdef SFR
      TimeBinSfr[bin] = 0;
#endif
#ifdef BLACK_HOLES
      TimeBin_BH_mass[bin] = 0;
      TimeBin_BH_dynamicalmass[bin] = 0;
      TimeBin_BH_Mdot[bin] = 0;
      TimeBin_BH_Medd[bin] = 0;
#endif
#ifdef LT_STELLAREVOLUTION
      TimeBinCountStars[bin] = 0;
#endif
    }

  for(i = 0; i < NumPart; i++)
    {
      bin = P[i].TimeBin;

      if(TimeBinCount[bin] > 0)
	{
	  PrevInTimeBin[i] = LastInTimeBin[bin];
	  NextInTimeBin[i] = -1;
	  NextInTimeBin[LastInTimeBin[bin]] = i;
	  LastInTimeBin[bin] = i;
	}
      else
	{
	  FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
	  PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	}
      TimeBinCount[bin]++;
      if(P[i].Type == 0)
	TimeBinCountSph[bin]++;

#ifdef SFR
      if(P[i].Type == 0)
	TimeBinSfr[bin] += SphP[i].Sfr;
#endif
#ifdef BLACK_HOLES
      if(P[i].Type == 5)
	{
	  TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
	  TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	  TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
	  TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
	}
#endif
#ifdef LT_STELLAREVOLUTION
      if(P[i].Type == 4)
	TimeBinCountStars[MetP[P[i].pt.MetID].ChemTimeBin]++;
#endif
    }

  make_list_of_active_particles();

  NumForceUpdate = 0;
#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  int il;
#pragma omp parallel for reduction(+:NumForceUpdate) private(i)
  for(il = 0; il < NActivePart; il++)
    {
      i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
      NumForceUpdate++;
      if(i >= NumPart)
	{
	  printf("Bummer i=%d\n", i);
	  terminate("inconsistent list");
	}
    }

  sumup_large_ints(1, &NumForceUpdate, &glob_sum);

  GlobNumForceUpdate = glob_sum;
}





void drift_particle(int i, integertime time1)
{
  int j;
  double dt_drift;

  integertime time0 = P[i].Ti_current;

  if(time1 < time0)
    {
      printf("i=%d time0=%d time1=%d\n", i, (int) time0, (int) time1);
      terminate("no prediction into past allowed");
    }

  if(time1 == time0)
    return;

  if(All.ComovingIntegrationOn)
    dt_drift = get_drift_factor(time0, time1);
  else
    dt_drift = (time1 - time0) * All.Timebase_interval;

#ifndef SIDM_FREEZE
  for(j = 0; j < 3; j++)
    P[i].Pos[j] += P[i].Vel[j] * dt_drift;
#endif

#ifdef DISTORTIONTENSORPS
  do_phase_space_drift(i, dt_drift);
#endif

#if defined(BLACK_HOLES) || defined(CA_BH_ACCRETION)
  if(P[i].Type == 0 && P[i].Mass > 0)
#else
  if(P[i].Type == 0)
#endif
    {
      double dt_gravkick, dt_hydrokick, dt_entr;

      if(All.ComovingIntegrationOn)
	{
	  dt_entr = (time1 - time0) * All.Timebase_interval;
	  dt_gravkick = get_gravkick_factor(time0, time1);
	  dt_hydrokick = get_hydrokick_factor(time0, time1);
	}
      else
	{
	  dt_entr = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
	}

#ifdef PMGRID
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] +=
	  (P[i].g.GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#else
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] += P[i].g.GravAccel[j] * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] += SphP[i].TurbAccel[j] * dt_gravkick;
#endif

      SphP[i].EntropyPred += SphP[i].e.DtEntropy * dt_entr;

#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
      bp_cr_diff_update(i, dt_entr);
#endif

      SphP[i].d.Density *= exp(-SphP[i].v.DivVel * dt_drift);

#ifdef EOS_DEGENERATE
      {
	double maxfac = dt_entr;
	double xnuc, tmpfac;
	for(j = 0; j < EOS_NSPECIES; j++)
	  {
	    xnuc = SphP[i].xnucPred[j] + SphP[i].dxnuc[j] * dt_entr;
	    if(xnuc > 1.0)
	      {
		tmpfac = (1.0 - xnuc) / SphP[i].dxnuc[j];
		if(tmpfac < maxfac)
		  maxfac = tmpfac;
	      }
	    if(xnuc < 0.0)
	      {
		tmpfac = (0.0 - xnuc) / SphP[i].dxnuc[j];
		if(tmpfac < maxfac)
		  maxfac = tmpfac;
	      }
	  }

	if(maxfac > 0)
	  {
	    for(j = 0; j < EOS_NSPECIES; j++)
	      {
		SphP[i].xnucPred[j] += SphP[i].dxnuc[j] * maxfac;
	      }
	  }
      }
#endif

      SphP[i].Pressure = get_pressure(i);


#ifdef KD_RESTRICT_NEIGHBOURS	/* prevent particles from getting too many true neighbours */
      double fak = exp(0.333333333333 * SphP[i].v.DivVel * dt_drift);
      if(fak > 1 && P[i].TrueNGB > All.DesNumNgb * 3)
	fak = 1;
      if(fak < 1 && P[i].TrueNGB < All.DesNumNgb * 0.4)
	fak = 1;
      PPP[i].Hsml *= fak;
#else
      PPP[i].Hsml *= exp(0.333333333333 * SphP[i].v.DivVel * dt_drift);
#endif

      if(PPP[i].Hsml < All.MinGasHsml)
	PPP[i].Hsml = All.MinGasHsml;

#ifdef MAXHSML
      if(PPP[i].Hsml > All.MaxHsml)
	PPP[i].Hsml = All.MaxHsml;
#endif

      drift_sph_extra_physics(i, time0, time1, dt_entr);
    }

  P[i].Ti_current = time1;
}

void check_particle_for_temperature_minimum(int i)
{
  double minentropy;
  if(All.MinEgySpec)
    {
#ifndef TRADITIONAL_SPH_FORMULATION
      minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].d.Density * All.cf_a3inv, GAMMA_MINUS1);
#else
      minentropy = All.MinEgySpec;
#endif
      if(SphP[i].Entropy < minentropy)
	{
	  SphP[i].Entropy = minentropy;
	  SphP[i].e.DtEntropy = 0;
	}
    }
}



void move_particles(integertime time1)
{
  int i;

//  if(ThisTask == 0)
//    printf("MOVE\n");
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < NumPart; i++)
    drift_particle(i, time1);
}


/* return the pressure of particle i */
double get_pressure(int i)
{
  MyFloat press = 0;

#ifndef EOS_DEGENERATE
#ifndef MHM
#ifndef SOFTEREQS
#ifndef VORONOI_MESHRELAX

#ifndef TRADITIONAL_SPH_FORMULATION
  press = SphP[i].EntropyPred * pow(SphP[i].d.Density, GAMMA);
#else
  press = GAMMA_MINUS1 * SphP[i].EntropyPred * SphP[i].d.Density;
#endif

#endif
#else
  if(SphP[i].d.Density * All.cf_a3inv >= All.PhysDensThresh)
    press =
      All.FactorForSofterEQS * SphP[i].EntropyPred * pow(SphP[i].d.Density,
							 GAMMA) + (1 -
								   All.FactorForSofterEQS) *
      All.cf_afac1 * GAMMA_MINUS1 * SphP[i].d.Density * All.InitGasU;
  else
    press = SphP[i].EntropyPred * pow(SphP[i].d.Density, GAMMA);
#endif
#else
  /* Here we use an isothermal equation of state */
  press = All.cf_afac1 * GAMMA_MINUS1 * SphP[i].d.Density * All.InitGasU;

#endif
#else
  /* call tabulated eos with physical units */
  struct eos_result res;
  eos_calc_egiven(SphP[i].d.Density * All.UnitDensity_in_cgs, SphP[i].xnucPred, SphP[i].EntropyPred,
		  &SphP[i].temp, &res);
  press = res.p.v / All.UnitPressure_in_cgs;
  SphP[i].dpdr =
    (res.p.drho +
     res.temp * gsl_pow_2(res.p.dtemp / (SphP[i].d.Density * All.UnitDensity_in_cgs)) / res.e.dtemp);
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#if defined( CR_UPDATE_PARANOIA )
  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
    CR_Particle_Update(SphP + i, CRpop);
#endif
#ifndef CR_NOPRESSURE
  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
    press += CR_Comoving_Pressure(SphP + i, CRpop);
#endif
#endif

#ifdef BP_REAL_CRs
  press += SphP[i].CRpPressure;
#endif


  return press;
}



void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
#ifdef MAGNETIC
  double dt_mag;
  if(All.ComovingIntegrationOn)
    dt_mag = get_magkick_factor(tstart, tend);
  else
    dt_mag = (tend - tstart) * All.Timebase_interval;
#endif

#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
#ifdef DIVBCLEANING_DEDNER
  SphP[i].PhiPred += SphP[i].DtPhi * dt_mag;
#endif
  int j;
  for(j = 0; j < 3; j++)
    SphP[i].b2.BPred[j] += SphP[i].DtB[j] * dt_mag;
#endif
#ifdef EULER_DISSIPATION
  double dea = SphP[i].DtEulerA * dt_entr;
  double deb = SphP[i].DtEulerB * dt_entr;
#ifdef WINDS
  if(SphP[i].DelayTime <= 0)
#endif
#ifdef BLACK_HOLES
    if(P[i].Mass > 0)
#endif
      if(dea < 10.0 && deb < 10.0)
	{
	  SphP[i].EulerA += dea;
	  SphP[i].EulerB += deb;
	}
      else
	{
	  printf("Task %d: Large DtEuler DtEulA/B=%e/%e deulA/B=%e/%e Dt_Entr=%e\n",
		 ThisTask, SphP[i].DtEulerA, SphP[i].DtEulerB, dea, deb, dt_entr);
	}
#endif
}




/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < NumPart; i++)
#ifdef WINDTUNNEL
    if(P[i].Type == 0)
#endif
      for(j = 0; j < 3; j++)
	{
	  while(P[i].Pos[j] < 0)
	    P[i].Pos[j] += boxsize[j];

	  while(P[i].Pos[j] >= boxsize[j])
	    P[i].Pos[j] -= boxsize[j];
	}
}
#endif
