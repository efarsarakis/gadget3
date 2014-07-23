#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"


/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */
void run(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  if(RestartFlag != 1)		/* need to compute forces at initial synchronization time, unless we restarted from restart files */
    {
      output_log_messages();

      domain_Decomposition(0, 0);

      set_non_standard_physics_for_current_time();

      compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

      compute_densities();	/* densities for synchronous particles */

      compute_hydro_accelerations();	/* hydro-accels for synchronous particles */

      calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */
    }

  CPU_Step[CPU_MISC] += measure_time();
  int i;
  for(i = 0; i < CPU_PARTS; i++)
    CPU_Step[i] = 0;

  while(1)			/* main timestep iteration loop */
    {
#ifdef GDE_BIGFLOAT
      xErrNo = 0;
#endif
      compute_statistics();	/* regular statistics outputs (like total energy) */

      //create_snapshot_if_desired();

      write_cpu_log();		/* output some CPU usage log-info (accounts for everything needed up to the current sync-point) */

      if(All.NumCurrentTiStep == 3)
	endrun(0);
      
      if(All.Ti_Current >= TIMEBASE)	/* check whether we reached the final time */
	{
	  if(ThisTask == 0)
	    printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);

	  restart(0);		/* write a restart file to allow continuation of the run for a larger value of TimeMax */

	  if(All.Ti_lastoutput != All.Ti_Current)	/* make a snapshot at the final time in case none has produced at this time */
	    savepositions(All.SnapshotFileCount++);	/* this will be overwritten if All.TimeMax is increased and the run is continued */

	  break;
	}

      find_timesteps();		/* find-timesteps */

      do_first_halfstep_kick();	/* half-step kick at beginning of timestep for synchronous particles */

      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time.
					 */

      output_log_messages();	/* write some info to log-files */

      set_non_standard_physics_for_current_time();	/* update auxiliary physics for current time */


      if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)	/* check whether we have a big step */
	{
	  domain_Decomposition(0, 0);	/* do domain decomposition if step is big enough, and set new list of active particles  */
	}
      else
	{
	  force_update_tree();	/* update tree dynamically with kicks of last step so that it can be reused */

	  make_list_of_active_particles();	/* now we can set the new chain list of active particles */
	}

      compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

      compute_densities();	/* densities for synchronous particles */

      compute_hydro_accelerations();	/* hydro-accels for synchronous particles */

      do_second_halfstep_kick();	/* this does the half-step kick at the end of the timestep */

      calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */

      /* Check whether we need to interrupt the run */
      int stopflag = 0;
      if(ThisTask == 0)
	{
	  FILE *fd;
	  char stopfname[1000];
	  sprintf(stopfname, "%sstop", All.OutputDir);
	  if((fd = fopen(stopfname, "r")))	/* Is the stop-file present? If yes, interrupt the run. */
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  if(CPUThisRun > 0.85 * All.TimeLimitCPU)	/* are we running out of CPU-time ? If yes, interrupt run. */
	    {
	      printf("reaching time-limit. stopping.\n");
	      stopflag = 2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
	  restart(0);		/* write restart file */
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag == 2 && ThisTask == 0)
	    {
	      FILE *fd;
	      char contfname[1000];
	      sprintf(contfname, "%scont", All.OutputDir);
	      if((fd = fopen(contfname, "w")))
		fclose(fd);

	      if(All.ResubmitOn)
		execute_resubmit_command();
	    }
	  return;
	}

      if(ThisTask == 0)
	{
	  /* is it time to write one of the regularly space restart-files? */
	  if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	    {
	      All.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag == 3)
	{
	  restart(0);		/* write an occasional restart file */
	  stopflag = 0;
	  All.TimeLastRestartFile += report_time();
	}

      set_random_numbers();	/* draw a new list of random numbers */

      report_memory_usage(&HighMark_run, "RUN");
#ifdef GDE_BIGFLOAT
      if(xErrNo != 0)
	terminate("GDE_BIGFLOAT error\n");
#endif
    }

}



void set_non_standard_physics_for_current_time(void)
{
#if defined(RADIATIVE_RATES) || defined(RADIATION)
  init_rad(All.Time);
#endif

#ifdef COOLING
  IonizeParams();		/* set UV background for the current time */
#endif
}



void calculate_non_standard_physics(void)
{

#ifdef SIDM
  if(ThisTask == 0)
    printf("Starting sidm...\n");
  sidm_DoScatter();
  if(ThisTask == 0)
    printf("done.\n");
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  reset_turb_temp();
#if defined(POWERSPEC_GRID)
  if(All.Time >= All.TimeNextTurbSpectrum)
    {
      powerspec_turb(All.FileNumberTurbSpectrum++);
      All.TimeNextTurbSpectrum += All.TimeBetTurbSpectrum;
    }
#endif
#endif

#ifdef REIONIZATION
  heating();
#endif

#ifdef CONDUCTION
  if(All.Conduction_Ti_endstep == All.Ti_Current)
    conduction();
#endif

#ifdef CR_DIFFUSION
  if(All.CR_Diffusion_Ti_endstep == All.Ti_Current)
    cosmic_ray_diffusion();
#endif

#ifdef RADTRANSFER
  double timeeach = 0, timeall = 0, tstart = 0, tend = 0;

  if(Flag_FullStep)		/* only do it for full timesteps */
    {
      All.Radiation_Ti_endstep = All.Ti_Current;

      if(ThisTask == 0)
	{
	  printf("Start Eddington tensor computation...\n");
	  fflush(stdout);
	}

      eddington();

#ifdef RT_RAD_PRESSURE
      n();
#endif

      if(ThisTask == 0)
	{
	  printf("done Eddington tensor! \n");
	  fflush(stdout);
	}

#ifdef EDDINGTON_TENSOR_SFR
      density_sfr();
      sfr_lum();
#endif

#ifdef EDDINGTON_TENSOR_STARS
      star_lum();
#endif

#ifdef EDDINGTON_TENSOR_GAS
      gas_lum();
#endif

#ifdef EDDINGTON_TENSOR_BH
      bh_lum();
#endif

      /***** evolve the transport of radiation *****/
      if(ThisTask == 0)
	{
	  printf("start radtransfer...\n");
	  fflush(stdout);
	}

      tstart = second();

      radtransfer();

      tend = second();
      timeeach = timediff(tstart, tend);
      MPI_Allreduce(&timeeach, &timeall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(ThisTask == 0)
	{
#ifndef RT_SILENT
	  printf("time consumed is %g \n", timeall);
#endif
	  printf("done with radtransfer! \n");
	  fflush(stdout);
	}

      All.Radiation_Ti_begstep = All.Radiation_Ti_endstep;
    }

  /* chemistry updated at sub-stepping as well */
#ifndef UM_CHEMISTRY
  radtransfer_update_chemistry();
#endif

  if(Flag_FullStep)		/* only do it for full timesteps */
    rt_write_stats();
#endif

#ifdef MHM
  /***** kinetic feedback *****/
  kinetic_feedback_mhm();
#endif

#ifdef BLACK_HOLES
  /***** black hole accretion and feedback *****/
  blackhole_accretion();
#endif

#if defined(BLACK_HOLES) || defined(VARIABLE_WINDS)
#ifdef FOF
  /* this will find new black hole seed halos and/or assign host halo masses for the variable wind model */
  if(All.Time >= All.TimeNextOnTheFlyFoF)
    {
#ifdef KD_SEED_STAR_MASS_FRACTION
      fof_fof(-2);
#else
      fof_fof(-1);
#endif
      if(All.ComovingIntegrationOn)
	All.TimeNextOnTheFlyFoF *= All.TimeBetOnTheFlyFoF;
      else
	All.TimeNextOnTheFlyFoF += All.TimeBetOnTheFlyFoF;
    }
#endif
#endif

#ifdef COOLING	/**** radiative cooling and star formation *****/

#ifdef CS_MODEL
  cs_cooling_and_starformation();
#else

#ifdef SFR
  cooling_and_starformation();

#ifdef LT_STELLAREVOLUTION
  if(ThisTask == 0)
    {
      printf("Start supernovae computation...\n");
      fflush(stdout);
    }

  //tstartSn = second();

  evolve_SN();

  //tendSn = second();
  //      CPU_sev = timediff(tstartSn, tendSn);

  //MPI_Reduce(&CPU_sev, &sumCPU_sev, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //if(ThisTask == 0)
  //All.CPU_SEv += sumCPU_sev / NTask;
#endif

#else
  cooling_only();
#endif

#endif
  CPU_Step[CPU_COOLINGSFR] += measure_time();
#endif /*ends COOLING */

#ifdef CHEMCOOL
  do_chemcool(-1, 0);
#endif

#if defined(CS_MODEL) && defined(CS_ENRICH)
#ifndef CS_FEEDBACK
  Flag_phase = 0;		/* no destinction between phases */

  cs_update_weights();
  CPU_Step[CPU_WEIGHTS_HOT] += measure_time();
  cs_enrichment();
  CPU_Step[CPU_ENRICH_HOT] += measure_time();
#else

  Flag_phase = 1;		/* COLD phase Flag_phase  = 1 */

  cs_update_weights();
  CPU_Step[CPU_WEIGHTS_COLD] += measure_time();
  cs_enrichment();
  CPU_Step[CPU_ENRICH_COLD] += measure_time();

  Flag_phase = 2;		/* HOT phase Flag_phase = 2 */

  cs_update_weights();
  CPU_Step[CPU_WEIGHTS_HOT] += measure_time();
  cs_enrichment();
  CPU_Step[CPU_ENRICH_HOT] += measure_time();

  Flag_phase = 0;
#endif
#endif

#ifdef CS_TESTS
  cs_energy_test();
#endif

#ifndef BH_BUBBLES
#ifdef BUBBLES
  double hubble_a;

  /**** bubble feedback *****/
  if(All.Time >= All.TimeOfNextBubble)
    {
#ifdef FOF
      fof_fof(-1);
      bubble();
#else
      bubble();
#endif
      if(All.ComovingIntegrationOn)
	{
	  hubble_a = hubble_function(All.Time);
	  All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	}
      else
	All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

      if(ThisTask == 0)
	printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
    }
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
  if(All.Time >= All.TimeOfNextBubble)
    {
      fof_fof(-1);

      if(All.ComovingIntegrationOn)
	{
	  hubble_a = hubble_func(All.Time);
	  All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	}
      else
	All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

      if(ThisTask == 0)
	printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
    }
#endif

#ifdef SINKS
  do_sinks();
#endif

#ifdef INVARIANCETEST
  compare_partitions();
#endif

#ifdef MOL_CLOUDS
  if(ThisTask == 0)
    {
      printf("Starting MOL_CLOUDS...\n");
      fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();

  do_mol_clouds();

  if(ThisTask == 0)
    {
      printf("done.\n");
      fflush(stdout);
    }
#endif

#ifdef BP_REAL_CRs
  bp_cr_evol();
#endif

#ifdef SCF_HYBRID
  SCF_do_center_of_mass_correction(0.75, 10.0 * SCF_HQ_A, 0.01, 1000);
#endif
}


void compute_statistics(void)
{
  if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
    {
#ifdef COMPUTE_POTENTIAL_ENERGY
      compute_potential();
#endif
      energy_statistics();	/* compute and output energy statistics */

#ifdef SCFPOTENTIAL
      SCF_write(0);
#endif

      All.TimeLastStatistics += All.TimeBetStatistics;
    }
}


void execute_resubmit_command(void)
{
  char buf[1000];
  sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
  system(buf);
#endif
}


void create_snapshot_if_desired(void)
{
#ifdef SYNCRONIZ_OUTPUT
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      {
#ifdef OUTPUTPOTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(RECOMPUTE_POTENTIAL_ON_OUTPUT))
	domain_Decomposition(0, 0);

	compute_potential();
#endif
#endif
	savepositions(All.SnapshotFileCount++);
	All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
      }
#endif
}


/*! This function finds the next synchronization point of the system
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.  If the
 * system dirfts over the desired time of a snapshot file, the
 * function will drift to this moment, generate an output, and then
 * resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, i, prev;
  integertime dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;
  int highest_active_bin, highest_occupied_bin;
  double timeold;

  timeold = All.Time;

  All.NumCurrentTiStep++;	/* we are now moving to the next sync point */

  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE, highest_occupied_bin = 0; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      highest_occupied_bin = n;
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
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_global);
#else
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif


#ifndef SYNCRONIZ_OUTPUT
  while(ti_next_kick_global >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      All.Ti_Current = All.Ti_nextoutput;

      if(All.ComovingIntegrationOn)
	All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
	All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

      set_cosmo_factors_for_current_time();

#ifdef TIMEDEPGRAV
      All.G = All.Gini * dGfak(All.Time);
#endif

      move_particles(All.Ti_nextoutput);

      CPU_Step[CPU_DRIFT] += measure_time();

#ifdef OUTPUTPOTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(RECOMPUTE_POTENTIAL_ON_OUTPUT))
      domain_Decomposition(0, 0);

      compute_potential();
#endif
#endif

#ifdef LT_SEvDbg
      get_metals_sumcheck(9);
#endif

      savepositions(All.SnapshotFileCount++);	/* write snapshot file */

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }
#endif



  All.Previous_Ti_Current = All.Ti_Current;
  All.Ti_Current = ti_next_kick_global;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  set_cosmo_factors_for_current_time();

#ifdef TIMEDEPGRAV
  All.G = All.Gini * dGfak(All.Time);
#endif

  All.TimeStep = All.Time - timeold;

#ifdef SUB_TURB_DRIVING
  sub_turb_move_perturbers(timeold, All.Time);
#endif


#ifdef WINDTUNNEL
  All.WindCurrentX += All.WindVel * All.TimeStep;
#endif

#ifdef LT_STELLAREVOLUTION
  All.Time_Age = get_age(All.Time);
#endif


  /* mark the bins that will be active */
  for(n = 1, TimeBinActive[0] = 1, NumForceUpdate = TimeBinCount[0], highest_active_bin = 0; n < TIMEBINS;
      n++)
    {
      dt_bin = (((integertime) 1) << n);
      if((ti_next_kick_global % dt_bin) == 0)
	{
	  TimeBinActive[n] = 1;
	  NumForceUpdate += TimeBinCount[n];
	  if(TimeBinCount[n])
	    highest_active_bin = n;
	}
      else
	TimeBinActive[n] = 0;
    }

  sumup_large_ints(1, &NumForceUpdate, &GlobNumForceUpdate);
  MPI_Allreduce(&highest_active_bin, &All.HighestActiveTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&highest_occupied_bin, &All.HighestOccupiedTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(GlobNumForceUpdate == All.TotNumPart)
    {
      Flag_FullStep = 1;
      if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
	terminate("Something is wrong with the time bins.\n");
    }
  else
    Flag_FullStep = 0;


  /* move the new set of active/synchronized particles */
  /* Note: We do not yet call make_list_of_active_particles(), since we
   * may still need to old list in the dynamic tree update
   */
  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	    drift_particle(i, All.Ti_Current);
	}
    }
}


void make_list_of_active_particles(void)
{
  int i, n, prev;
  /* make a link list with the particles in the active time bins */
  FirstActiveParticle = -1;

  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	    {
#ifdef BLACK_HOLES
	      if(P[i].Mass == 0)
		continue;
#endif
	      if(prev == -1)
		FirstActiveParticle = i;

	      if(prev >= 0)
		NextActiveParticle[prev] = i;

	      prev = i;
	    }
	}
    }

  if(prev >= 0)
    NextActiveParticle[prev] = -1;

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  build_active_particle_list();
#endif
}

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
void build_active_particle_list()
{
  int i;
  NActivePart = 0;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      ActiveParticleList[NActivePart] = i;
      NActivePart++;
    }
}
#endif

/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
integertime find_next_outputtime(integertime ti_curr)
{
  int i, iter = 0;
  integertime ti, ti_next;
  double next, time;

  DumpFlag = 1;
  ti_next = -1;


  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
	      else
		ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

#ifdef PROCESS_TIMES_OF_OUTPUTLIST
	      /* first, determine maximum output interval based on All.MaxSizeTimestep */
	      integertime timax = (integertime) (All.MaxSizeTimestep / All.Timebase_interval);

	      /* make it a power 2 subdivision */
	      integertime ti_min = TIMEBASE;
	      while(ti_min > timax)
		ti_min >>= 1;
	      timax = ti_min;

	      double multiplier = ti / ((double) timax);

	      /* now round this to the nearest multiple of timax */
	      ti = ((integertime) (multiplier + 0.5)) * timax;
#endif

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount)
			All.SnapshotFileCount = i;
		    }

		  if(ti_next > ti)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount)
			All.SnapshotFileCount = i;
		    }
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}
      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
	  else
	    ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }


  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n\n", next, DumpFlag);

    }

  return ti_next;
}




/*! This routine writes for every synchronisation point in the timeline information to two log-files:
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdTimebins we inform about the distribution of particles over the timebins, and which timebins are active on this step.
 * code is stored.
 */
void output_log_messages(void)
{
  double z;
  int i, j;
  long long tot, tot_sph;
  long long tot_count[TIMEBINS];
  long long tot_count_sph[TIMEBINS];
  long long tot_cumulative[TIMEBINS];
  int weight, corr_weight;
  double sum, avg_CPU_TimeBin[TIMEBINS], frac_CPU_TimeBin[TIMEBINS];

  sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);
  sumup_large_ints(TIMEBINS, TimeBinCountSph, tot_count_sph);

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nSync-Point %d, Time: %g, Redshift: %g, Nf = %d%09d, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z,
		  (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  printf("\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z, All.TimeStep,
		  log(All.Time) - log(All.Time - All.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nSync-Point %d, Time: %g, Nf = %d%09d, Systemstep: %g\n", All.NumCurrentTiStep,
		  All.Time, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep);
	  printf("\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
	  fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
		  All.TimeStep);
	  fflush(FdInfo);
	}

      for(i = 1, tot_cumulative[0] = tot_count[0]; i < TIMEBINS; i++)
	tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];


      for(i = 0; i < TIMEBINS; i++)
	{
	  for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
	    sum += All.CPU_TimeBinMeasurements[i][j];
	  if(All.CPU_TimeBinCountMeasurements[i])
	    avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
	  else
	    avg_CPU_TimeBin[i] = 0;
	}

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count[i] > 0; i--, weight *= 2)
	{
	  if(weight > 1)
	    corr_weight = weight / 2;
	  else
	    corr_weight = weight;

	  frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
	  sum += frac_CPU_TimeBin[i];
	}

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count[i] > 0; i--)
	{
	  if(sum)
	    frac_CPU_TimeBin[i] /= sum;
	}


      printf
	("Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
      fprintf(FdTimebin,
	      "Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
      for(i = TIMEBINS - 1, tot = tot_sph = 0; i >= 0; i--)
	if(tot_count_sph[i] > 0 || tot_count[i] > 0)
	  {
	    printf(" %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n",
		   TimeBinActive[i] ? 'X' : ' ',
		   i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
		   i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
		   (i == All.HighestActiveTimeBin) ? '<' : ' ',
		   (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
		   avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

	    fprintf(FdTimebin,
		    " %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n",
		    TimeBinActive[i] ? 'X' : ' ', i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
		    i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
		    (i == All.HighestActiveTimeBin) ? '<' : ' ',
		    (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
		    avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

	    if(TimeBinActive[i])
	      {
		tot += tot_count[i];
		tot_sph += tot_count_sph[i];
	      }
	  }
      printf("               ------------------------\n");
      fprintf(FdTimebin, "               ------------------------\n");

#ifdef PMGRID
      if(All.PM_Ti_endstep == All.Ti_Current)
	{
	  printf("PM-Step. Total: %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	  fprintf(FdTimebin, "PM-Step. Total: %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	}
      else
#endif
	{
	  printf("Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	  fprintf(FdTimebin, "Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);

	}
      fprintf(FdTimebin, "\n");
      fflush(FdTimebin);

    }

  output_extra_log_messages();
}




void write_cpu_log(void)
{
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum;
  int i;

  CPU_Step[CPU_MISC] += measure_time();

  for(i = 1, CPU_Step[0] = 0; i < CPU_PARTS; i++)
    CPU_Step[0] += CPU_Step[i];

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	avg_CPU_Step[i] /= NTask;

      put_symbol(0.0, 1.0, '#');

      for(i = 1, tsum = 0.0; i < CPU_PARTS; i++)
	{
	  if(max_CPU_Step[i] > 0)
	    {
	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_Symbol[i]);
	      tsum += t1 - t0;

	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_SymbolImbalance[i]);
	      tsum += t1 - t0;
	    }
	}

      put_symbol(tsum / max_CPU_Step[0], 1.0, '-');

      fprintf(FdBalance, "Step=%7d  sec=%10.3f  Nf=%2d%09d  %s\n", All.NumCurrentTiStep, max_CPU_Step[0],
	      (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), CPU_String);
      fflush(FdBalance);

      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
	{
	  All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
	  memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0],
		  &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1],
		  (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
	}

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements
							    [All.HighestActiveTimeBin]++] = max_CPU_Step[0];
    }

  CPUThisRun += CPU_Step[0];

  for(i = 0; i < CPU_PARTS; i++)
    CPU_Step[i] = 0;

  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	All.CPU_Sum[i] += avg_CPU_Step[i];

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU,
	      "total         %10.2f  %5.1f%%\n"
	      "treegrav      %10.2f  %5.1f%%\n"
	      "   treebuild  %10.2f  %5.1f%%\n"
	      "   treeupdate %10.2f  %5.1f%%\n"
	      "   treewalk   %10.2f  %5.1f%%\n"
	      "   treecomm   %10.2f  %5.1f%%\n"
	      "   treeimbal  %10.2f  %5.1f%%\n"
	      "pmgrav        %10.2f  %5.1f%%\n"
	      "sph           %10.2f  %5.1f%%\n"
	      "   density    %10.2f  %5.1f%%\n"
	      "   denscomm   %10.2f  %5.1f%%\n"
	      "   densimbal  %10.2f  %5.1f%%\n"
	      "   hydrofrc   %10.2f  %5.1f%%\n"
	      "   hydcomm    %10.2f  %5.1f%%\n"
	      "   hydmisc    %10.2f  %5.1f%%\n"
	      "   hydnetwork %10.2f  %5.1f%%\n"
	      "   hydimbal   %10.2f  %5.1f%%\n"
	      "   hmaxupdate %10.2f  %5.1f%%\n"
	      "domain        %10.2f  %5.1f%%\n"
	      "potential     %10.2f  %5.1f%%\n"
	      "predict       %10.2f  %5.1f%%\n"
	      "kicks         %10.2f  %5.1f%%\n"
	      "i/o           %10.2f  %5.1f%%\n"
	      "peano         %10.2f  %5.1f%%\n"
	      "sfrcool       %10.2f  %5.1f%%\n"
	      "blackholes    %10.2f  %5.1f%%\n"
	      "fof/subfind   %10.2f  %5.1f%%\n" "smoothing     %10.2f  %5.1f%%\n"
#ifdef CS_MODEL
	      "hotngbs       %10.2f  %5.1f%%\n"
	      "weights_hot   %10.2f  %5.1f%%\n"
	      "enrich_hot    %10.2f  %5.1f%%\n"
	      "weights_cold  %10.2f  %5.1f%%\n"
	      "enrich_cold   %10.2f  %5.1f%%\n" "cs_misc       %10.2f  %5.1f%%\n"
#endif
	      "misc          %10.2f  %5.1f%%\n",
	      All.CPU_Sum[CPU_ALL], 100.0,
	      All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]
	      + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
	      + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]
	      + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
	      + All.CPU_Sum[CPU_TREEMISC],
	      (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]
	       + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
	       + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]
	       + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
	       + All.CPU_Sum[CPU_TREEMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEBUILD],
	      (All.CPU_Sum[CPU_TREEBUILD]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEUPDATE],
	      (All.CPU_Sum[CPU_TREEUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2],
	      (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV],
	      (All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2],
	      (All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MESH],
	      (All.CPU_Sum[CPU_MESH]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT]
	      + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
	      + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
	      + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC] + All.CPU_Sum[CPU_HYDNETWORK],
	      (All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT]
	       + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
	       + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
	       + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC] +
	       All.CPU_Sum[CPU_HYDNETWORK]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DENSCOMPUTE],
	      (All.CPU_Sum[CPU_DENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DENSCOMM],
	      (All.CPU_Sum[CPU_DENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DENSWAIT],
	      (All.CPU_Sum[CPU_DENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDCOMPUTE],
	      (All.CPU_Sum[CPU_HYDCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDCOMM],
	      (All.CPU_Sum[CPU_HYDCOMM]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDMISC],
	      (All.CPU_Sum[CPU_HYDMISC]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDNETWORK],
	      (All.CPU_Sum[CPU_HYDNETWORK]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDWAIT],
	      (All.CPU_Sum[CPU_HYDWAIT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_TREEHMAXUPDATE],
	      (All.CPU_Sum[CPU_TREEHMAXUPDATE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DOMAIN],
	      (All.CPU_Sum[CPU_DOMAIN]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_POTENTIAL],
	      (All.CPU_Sum[CPU_POTENTIAL]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DRIFT],
	      (All.CPU_Sum[CPU_DRIFT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_TIMELINE],
	      (All.CPU_Sum[CPU_TIMELINE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_SNAPSHOT],
	      (All.CPU_Sum[CPU_SNAPSHOT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_PEANO],
	      (All.CPU_Sum[CPU_PEANO]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_COOLINGSFR],
	      (All.CPU_Sum[CPU_COOLINGSFR]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_BLACKHOLES],
	      (All.CPU_Sum[CPU_BLACKHOLES]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_FOF],
	      (All.CPU_Sum[CPU_FOF]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SMTHCOMPUTE] + All.CPU_Sum[CPU_SMTHWAIT] + All.CPU_Sum[CPU_SMTHCOMM] +
	      All.CPU_Sum[CPU_SMTHMISC],
	      (All.CPU_Sum[CPU_SMTHCOMPUTE] + All.CPU_Sum[CPU_SMTHWAIT] + All.CPU_Sum[CPU_SMTHCOMM] +
	       All.CPU_Sum[CPU_SMTHMISC]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef CS_MODEL
	      All.CPU_Sum[CPU_HOTNGBS],
	      (All.CPU_Sum[CPU_HOTNGBS]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_WEIGHTS_HOT],
	      (All.CPU_Sum[CPU_WEIGHTS_HOT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_ENRICH_HOT],
	      (All.CPU_Sum[CPU_ENRICH_HOT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_WEIGHTS_COLD],
	      (All.CPU_Sum[CPU_WEIGHTS_COLD]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_ENRICH_COLD],
	      (All.CPU_Sum[CPU_ENRICH_COLD]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_CSMISC],
	      (All.CPU_Sum[CPU_CSMISC]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
	      All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);
#ifdef MODGRAV
      fprintf(FdCPU,
	      "modgrav       %10.2f  %5.1f%%\n"
	      "   cic        %10.2f  %5.1f%%\n"
	      "   fieldsolve %10.2f  %5.1f%%\n"
	      "   eff_mass   %10.2f  %5.1f%%\n",
	      All.CPU_Sum[CPU_MG_CIC] + All.CPU_Sum[CPU_MG_FIELDSOLVE] + All.CPU_Sum[CPU_MG_EFF_MASS],
	      (All.CPU_Sum[CPU_MG_CIC] + All.CPU_Sum[CPU_MG_FIELDSOLVE] +
	       All.CPU_Sum[CPU_MG_EFF_MASS]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_MG_CIC],
	      All.CPU_Sum[CPU_MG_CIC] / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_MG_FIELDSOLVE],
	      All.CPU_Sum[CPU_MG_FIELDSOLVE] / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_MG_EFF_MASS],
	      All.CPU_Sum[CPU_MG_EFF_MASS] / All.CPU_Sum[CPU_ALL] * 100);
#endif
      fprintf(FdCPU, "\n");
      fflush(FdCPU);
    }
}


void put_symbol(double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j < 0)
    j = 0;
  if(i >= CPU_STRING_LEN)
    i = CPU_STRING_LEN;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    CPU_String[i++] = c;

  CPU_String[CPU_STRING_LEN] = 0;
}


/*! This routine first calls a computation of various global
 * quantities of the particle distribution, and then writes some
 * statistics about the energies in the various particle components to
 * the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

#if defined(MAGNETIC_STATISTICS)
      fprintf(FdEnergy, " %g", SysState.EnergyMag);
#if defined(TRACEDIVB)
      fprintf(FdEnergy, " %g", SysState.DivBerr);
#endif
#endif

      fprintf(FdEnergy, " \n");
      fflush(FdEnergy);
    }
}




void output_extra_log_messages(void)
{

#if defined (VS_TURB) || defined (AB_TURB)
  log_turb_temp();
#endif

  if(ThisTask == 0)
    {
#ifdef LT_STELLAREVOLUTION_NOT_YET_FIXED
      long long tot, tot_sph, tot_stars;
      long long tot_count[TIMEBINS], tot_count_sph[TIMEBINS], tot_count_stars[TIMEBINS];
      int i;

      sumup_large_ints(TIMEBINS, TimeBinCountStars, tot_count_stars);

      printf("Occupied timebins: non-sph         sph       stars    dt\n");
      for(i = TIMEBINS - 1, tot = tot_sph = tot_stars = 0; i >= 0; i--)
	if(tot_count_sph[i] > 0 || tot_count[i] > 0 || tot_count_stars[i] > 0)
	  {
	    printf(" %c  bin=%2d     %2d%09d %2d%09d %2d%09d   %6g\n",
		   TimeBinActive[i] ? 'X' : ' ',
		   i,
		   (int) ((tot_count[i] - tot_count_sph[i]) / 1000000000),
		   (int) ((tot_count[i] - tot_count_sph[i]) % 1000000000),
		   (int) (tot_count_sph[i] / 1000000000),
		   (int) (tot_count_sph[i] % 1000000000),
		   (int) (tot_count_stars[i] / 1000000000),
		   (int) (tot_count_stars[i] % 1000000000),
		   i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0);
	    if(TimeBinActive[i])
	      {
		tot += tot_count[i];
		tot_sph += tot_count_sph[i];
		tot_stars += tot_count_stars[i];
	      }
	  }
#endif


#ifdef CHEMISTRY
      printf("Abundances elec: %g, HM: %g, H2I: %g, H2II: %g\n",
	     SphP[1].elec, SphP[1].HM, SphP[1].H2I, SphP[1].H2II);
#endif

#ifdef XXLINFO
      if(Flag_FullStep == 1)
	{
	  fprintf(FdXXL, "%d %g ", All.NumCurrentTiStep, All.Time);
#ifdef MAGNETIC
	  fprintf(FdXXL, "%e ", MeanB);
#ifdef TRACEDIVB
	  fprintf(FdXXL, "%e ", MaxDivB);
#endif
#endif
#ifdef TIME_DEP_ART_VISC
	  fprintf(FdXXL, "%f ", MeanAlpha);
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
#endif

#ifdef MODGRAV
      if(All.ComovingIntegrationOn == 1)
	{
	  double hubble_a;

	  hubble_a = hubble_function(All.Time);
	  fprintf(FdMOG, "%d %g %e ", All.NumCurrentTiStep, All.Time, hubble_a);
	  fprintf(FdMOG, "\n");
	  fflush(FdMOG);
	}
#endif

#ifdef DARKENERGY
      if(All.ComovingIntegrationOn == 1)
	{
	  double hubble_a;

	  hubble_a = hubble_function(All.Time);
	  fprintf(FdDE, "%d %g %e ", All.NumCurrentTiStep, All.Time, hubble_a);
#ifndef TIMEDEPDE
	  fprintf(FdDE, "%e ", All.DarkEnergyParam);
#else
	  fprintf(FdDE, "%e %e ", get_wa(All.Time), DarkEnergy_a(All.Time));
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "%e %e", dHfak(All.Time), dGfak(All.Time));
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
#endif
    }
}


void check_particles_info(const char *func, const char *file, int linenr)
{
  int i, k, vok = 0, pok = 0, vsph = 0;
  double vv;

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Checking particle data (function %s in file %s at line %d) ...\n", func, file, linenr);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < NumPart; i++)
    {

      for(k = 0; k < 3; k++)
	{
	  if(P[i].Vel[k] > -1e8 && P[i].Vel[k] < 1e8)
	    {
	      vv = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
	      if(vv > 15000)
		{
		  printf
		    ("task=%d: WARNING: Large velocity for particle %d ID %llu v[%d]=%g, renormalizing it !!\n",
		     ThisTask, i, (unsigned long long) P[i].ID, k, vv);
		  fflush(stdout);
		  P[i].Vel[0] = P[i].Vel[0] / vv * 10000;
		  P[i].Vel[1] = P[i].Vel[1] / vv * 10000;
		  P[i].Vel[2] = P[i].Vel[2] / vv * 10000;
		}
	      vok++;
	    }
	  else
	    {
	      printf
		("task=%d:  strange value in velocity in for particle %d ID %llu , type=%d, mass=%g, v[%d]=%g\n",
		 ThisTask, i, (unsigned long long) P[i].ID, P[i].Type, P[i].Mass, k, P[i].Vel[k]);
	      fflush(stdout);
	      endrun(712401);
	    }

	  if(P[i].Pos[k] > -10000 && P[i].Pos[k] < All.BoxSize + 10000)
	    pok++;
	  else
	    {
	      printf("task=%d:  strange value in position in for particle %d ID %llu x[%d]=%g\n", ThisTask, i,
		     (unsigned long long) P[i].ID, k, P[i].Pos[k]);
	      fflush(stdout);
	      endrun(712402);
	    }
	}

      if(P[i].Type == 0)
	{
	  if((SphP[i].EntropyPred > -1e20 && SphP[i].EntropyPred < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange Entropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].EntropyPred);

	  if((SphP[i].e.DtEntropy > -1e20 && SphP[i].e.DtEntropy < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange DtEntropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].e.DtEntropy);

	  if((SphP[i].Pressure > -1e20 && SphP[i].Pressure < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange Pressure value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].Pressure);
	}

      if(P[i].Type > 5 || P[i].Type < 0)
	{
	  printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
	  endrun(712411);
	}

#ifdef LT_STELLAREVOLUTION
      if(P[i].Type == 4)
	if(MetP[P[i].pt.MetID].PID != i)
	  {
	    printf("task=%d:  error in cross-indexes for star-particle %d ID %llu\n", ThisTask, i,
		   (unsigned long long) P[i].ID);
	    fflush(stdout);
	    endrun(712412);
	  }
#endif

#if defined(BLACK_HOLES) && defined(DETACH_BLACK_HOLES)
      if(P[i].Type == 5)
	if(BHP[P[i].pt.BHID].PID != i)
	  {
	    printf("task=%d:  error in cross-indexes for bh-particle %d ID %llu\n", ThisTask, i,
		   (unsigned long long) P[i].ID);
	    fflush(stdout);
	    endrun(712413);
	  }
#endif
    }

  if(ThisTask == 0)
    printf("Positions and Velocities fine for (%d,%d,%d) of %d cases on task 0...\n", pok, vok, vsph,
	   NumPart * 3);

  MPI_Barrier(MPI_COMM_WORLD);
}
