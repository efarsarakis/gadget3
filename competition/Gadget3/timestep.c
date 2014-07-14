#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

/*! \file timestep.c
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps
 */

static double dt_displacement = 0;

#ifdef LT_STELLAREVOLUTION
static double time_convert_factor;
#endif

void set_cosmo_factors_for_current_time(void)
{

  if(All.ComovingIntegrationOn)
    {
      All.cf_atime = All.Time;
      All.cf_a2inv = 1 / (All.Time * All.Time);
      All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
      All.cf_afac1 = pow(All.Time, 3 * GAMMA_MINUS1);
      All.cf_afac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      All.cf_afac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      All.cf_hubble_a = hubble_function(All.Time);
    }
  else
    {
      All.cf_atime = 1;
      All.cf_a2inv = 1;
      All.cf_a3inv = 1;
      All.cf_afac1 = 1;
      All.cf_afac2 = 1;
      All.cf_afac3 = 1;
      All.cf_hubble_a = 1;
    }
}


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void find_timesteps(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  int i, bin, binold, prev, next;
  integertime ti_step, ti_step_old, ti_min;
  double aphys;

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin || dt_displacement == 0)
    find_dt_displacement_constraint(All.cf_hubble_a * All.cf_atime * All.cf_atime);

#ifdef LT_STELLAREVOLUTION
  time_convert_factor = (1000 * SEC_PER_MEGAYEAR) /	/* chemical time is in Gyr */
    All.UnitTime_in_s *		/* convert in code units   */
    All.HubbleParam;
#endif

#ifdef MAKEGLASS
  do_glass_making_step();
#endif


#ifdef FORCE_EQUAL_TIMESTEPS
  for(i = FirstActiveParticle, ti_min = TIMEBASE; i >= 0; i = NextActiveParticle[i])
    {
      ti_step = get_timestep(i, &aphys, 0);

      if(ti_step < ti_min)
	ti_min = ti_step;
    }

  if(ti_min > (dt_displacement / All.Timebase_interval))
    ti_min = (dt_displacement / All.Timebase_interval);

  ti_step = TIMEBASE;
  while(ti_step > ti_min)
    ti_step >>= 1;

  integertime ti_min_glob;

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &ti_step, &ti_min_glob);
#else
  MPI_Allreduce(&ti_step, &ti_min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef RELAXOBJECT
  determine_relaxfac();
#endif


  /* Now assign new timesteps  */

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef FORCE_EQUAL_TIMESTEPS
      ti_step = ti_min_glob;
#else
      ti_step = get_timestep(i, &aphys, 0);
#endif

      /* make it a power 2 subdivision */
      ti_min = TIMEBASE;
      while(ti_min > ti_step)
	ti_min >>= 1;
      ti_step = ti_min;

      bin = get_timestep_bin(ti_step);
      binold = P[i].TimeBin;

      if(bin > binold)		/* timestep wants to increase */
	{
	  while(TimeBinActive[bin] == 0 && bin > binold)	/* make sure the new step is synchronized */
	    bin--;

	  ti_step = bin ? (((integertime) 1) << bin) : 0;
	}

      if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
	{
	  ti_step = 0;
	  bin = 0;
	}

      if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	{
	  terminate("we are beyond the end of the timeline");	/* should not happen */
	  ti_step = TIMEBASE - All.Ti_Current;
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;
	}

      if(bin != binold)
	{
	  TimeBinCount[binold]--;
	  if(P[i].Type == 0)
	    {
	      TimeBinCountSph[binold]--;
#ifdef SFR
	      TimeBinSfr[binold] -= SphP[i].Sfr;
	      TimeBinSfr[bin] += SphP[i].Sfr;
#endif
	    }

#ifdef BLACK_HOLES
	  if(P[i].Type == 5)
	    {
	      TimeBin_BH_mass[binold] -= BPP(i).BH_Mass;
	      TimeBin_BH_dynamicalmass[binold] -= P[i].Mass;
	      TimeBin_BH_Mdot[binold] -= BPP(i).BH_Mdot;
	      if(BPP(i).BH_Mass > 0)
		TimeBin_BH_Medd[binold] -= BPP(i).BH_Mdot / BPP(i).BH_Mass;
	      TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
	      TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	      TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
	      if(BPP(i).BH_Mass > 0)
		TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
	    }
#endif

	  prev = PrevInTimeBin[i];
	  next = NextInTimeBin[i];

	  if(FirstInTimeBin[binold] == i)
	    FirstInTimeBin[binold] = next;
	  if(LastInTimeBin[binold] == i)
	    LastInTimeBin[binold] = prev;
	  if(prev >= 0)
	    NextInTimeBin[prev] = next;
	  if(next >= 0)
	    PrevInTimeBin[next] = prev;

	  if(TimeBinCount[bin] > 0)
	    {
	      PrevInTimeBin[i] = LastInTimeBin[bin];
	      NextInTimeBin[LastInTimeBin[bin]] = i;
	      NextInTimeBin[i] = -1;
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

	  P[i].TimeBin = bin;
	}

#ifndef WAKEUP
      ti_step_old = binold ? (((integertime) 1) << binold) : 0;
#else
      ti_step_old = P[i].dt_step;
#endif

      P[i].Ti_begstep += ti_step_old;

#ifdef WAKEUP
      P[i].dt_step = ti_step;
#endif
    }


#ifdef CONDUCTION
  if(All.Conduction_Ti_endstep == All.Ti_Current)
    {
      ti_step = TIMEBASE;
      while(ti_step > (All.MaxSizeConductionStep / All.Timebase_interval))
	ti_step >>= 1;
      while(ti_step > (All.MaxSizeTimestep / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.Conduction_Ti_endstep - All.Conduction_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.Conduction_Ti_endstep) % ti_step) > 0)
	    ti_step = All.Conduction_Ti_endstep - All.Conduction_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      All.Conduction_Ti_begstep = All.Conduction_Ti_endstep;
      All.Conduction_Ti_endstep = All.Conduction_Ti_begstep + ti_step;
    }
#endif

#ifdef CR_DIFFUSION
  if(All.CR_Diffusion_Ti_endstep == All.Ti_Current)
    {
      if(All.CR_Diffusion_Ti_endstep < All.Ti_Current)
	endrun(1231);

      ti_step = TIMEBASE;
      while(ti_step > (All.CR_DiffusionMaxSizeTimestep / All.Timebase_interval))
	ti_step >>= 1;
      while(ti_step > (All.MaxSizeTimestep / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.CR_Diffusion_Ti_endstep - All.CR_Diffusion_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.CR_Diffusion_Ti_endstep) % ti_step) > 0)
	    ti_step = All.CR_Diffusion_Ti_endstep - All.CR_Diffusion_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      All.CR_Diffusion_Ti_begstep = All.CR_Diffusion_Ti_endstep;
      All.CR_Diffusion_Ti_endstep = All.CR_Diffusion_Ti_begstep + ti_step;
    }
#endif


#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  bin = get_timestep_bin(ti_step);
	  binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);

	  while(TimeBinActive[bin] == 0 && bin > binold)	/* make sure the new step is synchronized */
	    bin--;

	  ti_step = bin ? (((integertime) 1) << bin) : 0;
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;
    }
#endif


#ifdef WAKEUP
  process_wake_ups();
#endif

  CPU_Step[CPU_TIMELINE] += measure_time();
}



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
integertime get_timestep(int p,	/*!< particle index */
			 double *aphys,	/*!< acceleration (physical units) */
			 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
					   aphys */ )
{
  double ax, ay, az, ac;
  double csnd = 0, dt = 0, dt_courant = 0, dt_divv = 0;
  integertime ti_step;
#ifdef CHEMCOOL
  double hubble_param;

  if(All.ComovingIntegrationOn)
    hubble_param = All.HubbleParam;
  else
    hubble_param = 1.0;
#endif

#ifdef BLACK_HOLES
  double dt_accr;

#ifdef UNIFIED_FEEDBACK
  double meddington = 0;
#endif
#endif

#ifdef NS_TIMESTEP
  double dt_NS = 0;
#endif

#ifdef NONEQUILIBRIUM
  double dt_cool, dt_elec;
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#ifdef NUCLEAR_NETWORK
  double dt_network, dt_species;
  int k;
#endif

#ifdef SIDM
  double dt_sidm;
#endif
#ifdef DEDNER_TIMESTEP
  double dt_ded;
#endif
#ifdef FS_ETA2_DYN
  double dt_eta;
#endif

  if(flag == 0)
    {
      ax = All.cf_a2inv * P[p].g.GravAccel[0];
      ay = All.cf_a2inv * P[p].g.GravAccel[1];
      az = All.cf_a2inv * P[p].g.GravAccel[2];

#ifdef PMGRID
      ax += All.cf_a2inv * P[p].GravPM[0];
      ay += All.cf_a2inv * P[p].GravPM[1];
      az += All.cf_a2inv * P[p].GravPM[2];
#endif

#if defined(AB_TURB) || defined(TURB_DRIVING)
      if(P[p].Type == 0)
	{
	  ax += SphP[p].TurbAccel[0];
	  ay += SphP[p].TurbAccel[1];
	  az += SphP[p].TurbAccel[2];
	}
#endif

      if(P[p].Type == 0)
	{
	  ax += All.cf_afac2 * SphP[p].a.HydroAccel[0];
	  ay += All.cf_afac2 * SphP[p].a.HydroAccel[1];
	  az += All.cf_afac2 * SphP[p].a.HydroAccel[2];
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;


  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;

	  dt /= All.cf_hubble_a;	/* convert dloga to physical timestep  */
#ifdef ADAPTGRAVSOFT
	  ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * P[p].AGS_Hsml / 2.8 / (dt * dt);
#else
	  ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / (dt * dt);
#endif
	  *aphys = ac;
	  return flag;
	}
#ifdef ADAPTGRAVSOFT
      dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * P[p].AGS_Hsml / 2.8 / ac);
#else
      dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * PPP[p].Hsml / 2.8 / ac);
#else
      if(P[p].Type == 0)
	dt =
	  sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] *
	       pow(P[p].Mass / All.ReferenceGasMass, 1.0 / 3) / ac);
#endif
#endif
      break;

    default:
      fprintf(stderr, "\n !!!2@@@!!! \n");
      endrun(888);
      fprintf(stderr, "\n !!!2@@@!!! \n");
      break;
    }

#ifdef BLACK_HOLES
  if(P[p].Type == 0 && P[p].Mass > 0)
#else
  if(P[p].Type == 0)
#endif
    {
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].d.Density);

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      double dmax1, dmax2;

      if(All.ComovingIntegrationOn)
	dt_courant =
	  All.CourantFac * All.Time * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / (All.cf_afac3 * csnd);
      else
	dt_courant = All.CourantFac * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / csnd;

      if(dt_courant > 2 * All.CourantFac * SphP[p].MinViscousDt)
	dt_courant = 2 * All.CourantFac * SphP[p].MinViscousDt;
#else
      if(All.ComovingIntegrationOn)
	dt_courant = 2 * All.CourantFac * All.Time * PPP[p].Hsml / (All.cf_afac3 * SphP[p].MaxSignalVel);
      else
	dt_courant = 2 * All.CourantFac * PPP[p].Hsml / SphP[p].MaxSignalVel;
#endif

#ifdef VORONOI
#ifdef TWODIMS
      dt_courant = All.CourantFac * sqrt(SphP[p].Volume / M_PI) / SphP[p].MaxSignalVel;
#else
      dt_courant = All.CourantFac * pow(SphP[p].Volume / (4.0 / 3 * M_PI), 1.0 / 3) / SphP[p].MaxSignalVel;
#endif
#endif

      if(dt_courant < dt)
	dt = dt_courant;

#ifdef FS_ETA2_DYN

      if(SphP[p].eta2 != 0)
	{
	  dt_eta =
	    All.CourantFac * All.CourantFac * PPP[p].Hsml * PPP[p].Hsml / (All.cf_afac2 * SphP[p].eta2);

	  if(dt_eta < dt)
	    dt = dt_eta;
	}
#endif
#ifdef DEDNER_TIMESTEP
      if(SphP[p].divB != 0)
	{
	  dt_ded = sqrt(SphP[p].d.Density / (MU0_1 * SphP[p].divB * SphP[p].divB));
	  if(P[p].ID != 0)
	    if(P[p].TimeBin > 10 && TimeBinCountSph[P[p].TimeBin] > All.TotN_gas / 100.)
	      dt = dt_ded < dt ? dt_ded : dt;
	}
#endif

      /* make sure that the velocity divergence does not imply a too large change of density or smoothing length in the step */
      if(SphP[p].v.DivVel != 0)
	{
	  dt_divv = 1.5 / fabs(All.cf_a2inv * SphP[p].v.DivVel);
	  if(dt_divv < dt)
	    dt = dt_divv;
	}

#ifdef MYFALSE
      dt_viscous = All.CourantFac * SphP[p].MaxViscStep / All.cf_hubble_a;	/* to convert dloga to physical dt */

      if(dt_viscous < dt)
	dt = dt_viscous;
#endif

#ifdef NS_TIMESTEP
      if(fabs(SphP[p].ViscEntropyChange))
	{
	  dt_NS = VISC_TIMESTEP_PARAMETER * SphP[p].Entropy / SphP[p].ViscEntropyChange / All.cf_hubble_a;

	  if(dt_NS < dt)
	    dt = dt_NS;
	}
#endif


#ifdef NUCLEAR_NETWORK
      if(SphP[p].temp > 1e7)
	{
	  /* check if the new timestep blows up our abundances */
	  dt_network = dt * All.UnitTime_in_s;
	  for(k = 0; k < EOS_NSPECIES; k++)
	    {
	      if(SphP[p].dxnuc[k] > 0)
		{
		  dt_species = (1.0 - SphP[p].xnuc[k]) / SphP[p].dxnuc[k];
		  if(dt_species < dt_network)
		    dt_network = dt_species;
		}
	      else if(SphP[p].dxnuc[k] < 0)
		{
		  dt_species = (0.0 - SphP[p].xnuc[k]) / SphP[p].dxnuc[k];
		  if(dt_species < dt_network)
		    dt_network = dt_species;
		}

	    }

	  dt_network /= All.UnitTime_in_s;
	  if(dt_network < dt)
	    dt = dt_network;
	}
#endif

#ifdef CHEMCOOL
      dt_courant = do_chemcool(p, dt_courant * All.UnitTime_in_s / hubble_param);

      if(dt_courant < dt)
	dt = dt_courant;
#endif
    }


#ifdef SIDM
  /* constraint time step to avoid too many multiple scatter events (prefactor 0.1 can be adopted) */
  if((1 << P[p].Type) & (SIDM))
    {
      if(P[p].sidm_Density > 0.0 && P[p].sidm_VelDisp > 0.0)
	{
	  dt_sidm =
	    SIDM_DTIME_FAC / ((All.HubbleParam * All.HubbleParam * All.cf_a3inv * P[p].sidm_Density) *
			      sidm_cross_sigma(P[p].sidm_VelDisp / All.cf_atime) * (P[p].sidm_VelDisp /
										    All.cf_atime));
	  if(dt_sidm < dt)
	    dt = dt_sidm;
	}
    }
#endif

#ifdef BLACK_HOLES
  if(P[p].Type == 5)
    {
      if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
	{
	  dt_accr = 0.25 * BPP(p).BH_Mass / BPP(p).BH_Mdot;
	  if(dt_accr < dt)
	    dt = dt_accr;
	}
    }
#endif

#ifdef  BLACK_HOLES
  if(P[p].Type == 5)
    {
      double dt_ngbs =
	(BPP(p).BH_TimeBinGasNeighbor ? (1 << BPP(p).BH_TimeBinGasNeighbor) : 0) * All.Timebase_interval /
	All.cf_hubble_a;

      if(dt > dt_ngbs && dt_ngbs > 0)
	dt = 1.01 * dt_ngbs;
    }
#endif

#ifdef BH_BUBBLES
  if(P[p].Type == 5)
    {
      if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
	{
#ifdef UNIFIED_FEEDBACK
	  meddington = (4 * M_PI * GRAVITY * C * PROTONMASS /
			(All.BlackHoleRadiativeEfficiency * C * C * THOMPSON)) * BPP(p).BH_Mass *
	    All.UnitTime_in_s;
	  if(BPP(p).BH_Mdot < All.RadioThreshold * meddington)
#endif
	    dt_accr = (All.BlackHoleRadioTriggeringFactor - 1) * BPP(p).BH_Mass / BPP(p).BH_Mdot;
	  if(dt_accr < dt)
	    dt = dt_accr;
	}
    }
#endif

#ifdef NONEQUILIBRIUM
  /* another criterion given by the local cooling time */

  if(P[p].Type == 0)
    {
      dt_cool = fabs(SphP[p].t_cool);	/* still in yrs */
      dt_cool *= SEC_PER_YEAR;	/* in seconds */
      dt_cool /= All.UnitTime_in_s;
      dt_cool *= All.HubbleParam;	/* internal units */

      dt_cool = All.Epsilon * dt_cool;


#ifndef UM_CONTINUE
      if(dt_cool > 0 && dt_cool < dt)
	dt = dt_cool;
#else
#ifdef WINDS
      if(dt_cool > 0 && dt_cool < dt && SphP[p].DelayTime < 0)
	dt = dt_cool;
#else
      if(dt_cool > 0 && dt_cool < dt)
	dt = dt_cool;
#endif
#endif


      /* yet another criterion given by the electron number density change */

      dt_elec = fabs(SphP[p].t_elec);	/* still in yrs */
      dt_elec *= SEC_PER_YEAR;	/* in seconds */
      dt_elec /= All.UnitTime_in_s;
      dt_elec *= All.HubbleParam;	/* internal units */


      dt_elec = All.Epsilon * dt_elec;

#ifndef UM_CONTINUE
      if(dt_elec > 0 && dt_elec < dt)
	dt = dt_elec;
#else
      if(dt_elec > 0 && dt_elec < dt && SphP[p].DelayTime < 0)
	dt = dt_elec;
#endif
    }
#endif

#if defined(LT_STELLAREVOLUTION) && defined(LT_CONSTRAIN_DYNAMICAL_TIMESTEP)

  double dt_chem;

  if(P[i].Type == 4)
    {
      dt_chem = MetP[P[i].MetID].NextChemTime - All.Time_Age;	/* this is the next chemical time in Gyr */
      dt_chem *= time_convert_factor;	/* convert in code physical units        */

      if(dt_chem > 0 && dt_chem < dt)
	dt = dt_chem;
    }
#endif


  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     All.cf_hubble_a=1.
   */
  dt *= All.cf_hubble_a;

#ifdef ONLY_PM
  dt = All.MaxSizeTimestep;
#endif



  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;


  if(dt >= dt_displacement)
    dt = dt_displacement;


#ifdef CONDUCTION
  if(P[p].Type == 0)
    if(dt >= All.MaxSizeConductionStep)
      dt = All.MaxSizeConductionStep;
#endif
#ifdef CR_DIFFUSION
  if(P[p].Type == 0)
    if(dt >= All.CR_DiffusionMaxSizeTimestep)
      dt = All.CR_DiffusionMaxSizeTimestep;
#endif

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
#ifndef LONGIDS
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
	     PPP[p].Hsml, csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) *
	     All.cf_hubble_a, All.SofteningTable[P[p].Type]);
#else
	  printf
	    ("Part-ID=%llu  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (MyIDType) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
	     PPP[p].Hsml, csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) *
	     All.cf_hubble_a, All.SofteningTable[P[p].Type]);
#endif

#ifdef NS_TIMESTEP
#ifndef LONGIDS
	  printf
	    ("Part-ID=%d  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     (int) P[p].ID, dt_NS * All.cf_hubble_a, SphP[p].Entropy, SphP[p].d.Density,
	     SphP[p].ViscEntropyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval,
	     All.IonMeanFreePath *
	     pow((SphP[p].Entropy * pow(SphP[p].d.Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1),
		 2.0) / SphP[p].d.Density);
#else
	  printf
	    ("Part-ID=%llu  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     (MyIDType) P[p].ID, dt_NS * All.cf_hubble_a, SphP[p].Entropy, SphP[p].d.Density,
	     SphP[p].ViscEntropyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval,
	     All.IonMeanFreePath *
	     pow((SphP[p].Entropy * pow(SphP[p].d.Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1),
		 2.0) / SphP[p].d.Density);
#endif

	  printf("Stressd=(%g|%g|%g) \n", SphP[p].u.s.StressDiag[0], SphP[p].u.s.StressDiag[1],
		 SphP[p].u.s.StressDiag[2]);
	  printf("Stressoffd=(%g|%g|%g) \n", SphP[p].u.s.StressOffDiag[0], SphP[p].u.s.StressOffDiag[1],
		 SphP[p].u.s.StressOffDiag[2]);
#endif


	}
      else
	{
#ifndef LONGIDS
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
#else
	  printf("Part-ID=%llu  dt=%g ac=%g xyz=(%g|%g|%g)\n", (MyIDType) P[p].ID, dt, ac, P[p].Pos[0],
		 P[p].Pos[1], P[p].Pos[2]);
#endif
	}
      fflush(stdout);
      fprintf(stderr, "\n @ fflush \n");
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = (integertime) (dt / All.Timebase_interval);

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  if(ti_step == 0)
    {
#ifndef LONGIDS
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dt_elec=%g dt_cool=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g)\n\n",
	     ThisTask, P[p].ID, dt, SphP[p].t_elec, SphP[p].t_cool, All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
#else
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%llu dt=%g dt_elec=%g dt_cool=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g)\n\n",
	     ThisTask, P[p].ID, dt, SphP[p].t_elec, SphP[p].t_cool, All.Timebase_interval, (int) ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
#endif
      fflush(stdout);
      endrun(818);
    }
#endif

  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%llu dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
	     ThisTask, (unsigned long long) P[p].ID, dt, dt_courant, dt_divv, dt_displacement,
	     All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].g.GravAccel[0], P[p].g.GravAccel[1],
	     P[p].g.GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g) dens=%g hsml=%g  dTdE=%g\n", SphP[p].a.HydroAccel[0],
	       SphP[p].a.HydroAccel[1], SphP[p].a.HydroAccel[2], SphP[p].d.Density, PPP[p].Hsml,
	       SphP[p].e.dDtEntropy);

#ifdef DIVBFORCE3
      if(P[p].Type == 0)
	printf("mhd-frc=(%g|%g|%g) mhd-corr=(%g|%g|%g)\n",
	       SphP[p].magacc[0], SphP[p].magacc[1], SphP[p].magacc[2],
	       SphP[p].magcorr[0], SphP[p].magcorr[1], SphP[p].magcorr[2]);
#endif
#if defined(FS_ETA2_DYN) && defined(DEDNER_TIMESTEP)
      if(P[p].Type == 0)
	printf("dt_eta= %g mhd-dt_ded= %g \n", dt_eta, dt_ded);
#endif

#ifdef COSMIC_RAYS
      if(P[p].Type == 0)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  printf("Cosmic Ray Properties: C0: %g -- q0  : %g -- P  : %g\n"
		 "                       Rho: %g\n\n",
		 SphP[p].CR_C0[CRpop], SphP[p].CR_q0[CRpop], CR_Particle_Pressure(SphP + p, CRpop),
		 SphP[p].d.Density);
#endif

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, type;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(P[i].Mass > 0)
	    {
	      if(mim[P[i].Type] > P[i].Mass)
		mim[P[i].Type] = P[i].Mass;
	    }
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      sumup_large_ints(6, count, count_sum);

#ifdef SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      min_mass[5] = min_mass[0];
#endif
#endif

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0 || (type == 4 && All.StarformationOn))
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

#ifdef BLACK_HOLES
	      if(type == 5)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
#endif
	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);


#ifdef NEUTRINOS
	      if(type != 2)	/* don't constrain the step to the neutrinos */
#endif
		if(dt < dt_displacement)
		  dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}

int get_timestep_bin(integertime ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}




#ifdef MAKEGLASS
void do_glass_making_step(void)
{
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].g.GravAccel[j] *= -1;
#ifdef PMGRID
	  P[i].GravPM[j] *= -1;
	  P[i].g.GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
#endif
	}

      disp = sqrt(P[i].g.GravAccel[0] * P[i].g.GravAccel[0] + P[i].g.GravAccel[1] * P[i].g.GravAccel[1]
		  + P[i].g.GravAccel[2] * P[i].g.GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n", dmean, globmax,
	     sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].g.GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].g.GravAccel[j] = 0;
	}
    }
}
#endif

#ifdef RELAXOBJECT
void determine_relaxfac(void)
{
  if(All.Time < 0.2 * All.TimeMax)
    {
      All.RelaxFac = 1. / All.RelaxBaseFac;
    }
  else if(All.Time > 0.8 * All.TimeMax)
    {
      All.RelaxFac = 0.;
    }
  else
    {
      All.RelaxFac =
	1. / (All.RelaxBaseFac * pow(10., (All.Time - 0.2 * All.TimeMax) / (0.6 * All.TimeMax) * 3.));
    }
}
#endif



#ifdef WAKEUP
void process_wake_ups(void)
{
  int i, k, n, dt_bin, dt_step;
  int ti_step_old, ti_next_for_bin, ti_next_kick, ti_next_kick_global, max_time_bin_active;
  int bin, binold, prev, next;
  int time0, time1_old, time1_new;
  double dt_entr, dt_gravkick, dt_hydrokick;
  long long ntot;

  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      dt_bin = (((integertime) 1) << n);
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_glob);
#else
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    printf("predicting next timestep: %g\n", (ti_next_kick_global - All.Ti_Current) * All.Timebase_interval);

  max_time_bin_active = 0;
  /* get the highest bin, that is active next time */
  for(n = 0; n < TIMEBINS; n++)
    {
      dt_bin = (((integertime) 1) << n);

      if((ti_next_kick_global % dt_bin) == 0)
	max_time_bin_active = n;
    }

  /* move the particle on the highest bin, that is active in the next timestep and that is lower than its last timebin */
  bin = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n] > 0)
	{
	  bin = n;
	  break;
	}
    }
  n = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type != 0)
	continue;

      if(!SphP[i].wakeup)
	continue;

      binold = P[i].TimeBin;
      if(TimeBinActive[binold])
	continue;

      bin = max_time_bin_active < binold ? max_time_bin_active : binold;

      if(bin != binold)
	{
	  TimeBinCount[binold]--;
	  if(P[i].Type == 0)
	    TimeBinCountSph[binold]--;

	  prev = PrevInTimeBin[i];
	  next = NextInTimeBin[i];

	  if(FirstInTimeBin[binold] == i)
	    FirstInTimeBin[binold] = next;
	  if(LastInTimeBin[binold] == i)
	    LastInTimeBin[binold] = prev;
	  if(prev >= 0)
	    NextInTimeBin[prev] = next;
	  if(next >= 0)
	    PrevInTimeBin[next] = prev;

	  if(TimeBinCount[bin] > 0)
	    {
	      PrevInTimeBin[i] = LastInTimeBin[bin];
	      NextInTimeBin[LastInTimeBin[bin]] = i;
	      NextInTimeBin[i] = -1;
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

	  P[i].TimeBin = bin;

	  if(TimeBinActive[bin])
	    NumForceUpdate++;

	  /* correct quantities predicted for a longer timestep */
	  ti_step_old = P[i].dt_step;
	  dt_step = ti_next_kick_global - P[i].Ti_begstep;
	  P[i].dt_step = dt_step;
	  /*
	     dt_entr = (-ti_step_old / 2 + dt_step / 2) * All.Timebase_interval;
	   */
	  time0 = P[i].Ti_begstep;
	  time1_old = P[i].Ti_begstep + ti_step_old;
	  time1_new = P[i].Ti_begstep + dt_step;

	  /* This part has still to be adapted ...
	     #ifdef PMGRID
	     if(All.ComovingIntegrationOn)
	     dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
	     get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
	     else
	     dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
	     #endif
	   */
	  if(All.ComovingIntegrationOn)
	    {
	      dt_entr = dt_gravkick = dt_hydrokick = (-(time1_old - time0) / 2
						      + (time1_new - time0) / 2) * All.Timebase_interval;
	      dt_gravkick = -get_gravkick_factor(time0, time1_old) / 2
		+ get_gravkick_factor(time0, time1_new) / 2;
	      dt_hydrokick = -get_hydrokick_factor(time0, time1_old) / 2
		+ get_hydrokick_factor(time0, time1_new) / 2;
	    }
	  else
	    {
	      dt_entr = dt_gravkick = dt_hydrokick = (-(time1_old - time0) / 2
						      + (time1_new - time0) / 2) * All.Timebase_interval;
	    }

	  /* This may now work in comoving runs */
	  /* WARNING: this velocity correction is inconsistent,
	   * as the position of the particle was calculated with a "wrong" velocity before  */
	  for(k = 0; k < 3; k++)
	    {
	      P[i].Vel[k] += P[i].g.GravAccel[k] * dt_gravkick;
	    }

	  for(k = 0; k < 3; k++)
	    {
	      P[i].Vel[k] += SphP[i].a.HydroAccel[k] * dt_hydrokick;
	    }
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
	  printf("Not yet addapted to new integration scheme for B !!\n");
	  endrun(87654210);
	  for(k = 0; k < 3; k++)
	    {
	      SphP[i].B[k] += SphP[i].DtB[k] * dt_entr;
	      SphP[i].BPred[k] = SphP[i].B[k] - SphP[i].DtB[k] * dt_entr;
	    }
#endif
#if !defined(EOS_DEGENERATE)
	  SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
#else
	  SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr * All.UnitTime_in_s;
#endif

#ifdef NUCLEAR_NETWORK
	  for(k = 0; k < EOS_NSPECIES; k++)
	    {
	      SphP[i].xnuc[k] += SphP[i].dxnuc[k] * dt_entr * All.UnitTime_in_s;
	    }
	  network_normalize(SphP[i].xnuc, &SphP[i].Entropy, &All.nd, &All.nw);
#endif

	  n++;
	}
    }

  sumup_large_ints(1, &n, &ntot);
  if(ThisTask == 0)
    printf("%d%09d particles woken up.\n", (int) (ntot / 1000000000), (int) (ntot % 1000000000));
}
#endif
