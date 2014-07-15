#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


void force_update_tree(void)
{
  int i, j;

  if(ThisTask == 0)
    printf("kicks will prepare for dynamic update of tree\n");

  GlobFlag++;
  DomainNumChanged = 0;
  DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));


  /* note: the current list of active particles still refers to that
   * synchronized at the previous time.
   */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      force_kick_node(i, P[i].dp);	/* kick the parent nodes with this momentum
					   difference, also updated maximum velocity, softening and soundspeed, if needed */
      for(j = 0; j < 3; j++)
	P[i].dp[j] = 0;
    }

  force_finish_kick_nodes();
  myfree(DomainList);

  if(ThisTask == 0)
    printf("Tree has been updated dynamically.\n");
}









void force_kick_node(int i, MyDouble * dp)
{
  int j, no;
  MyFloat v, vmax;

#ifdef SCALARFIELD
  MyFloat dp_dm[3];
#endif

#ifdef NEUTRINOS
  if(P[i].Type == 2)
    return;
#endif

  for(j = 0; j < 3; j++)
    {
#ifdef SCALARFIELD
      if(P[i].Type != 0)
	dp_dm[j] = dp[j];
      else
	dp_dm[j] = 0;
#endif
    }

  for(j = 0, vmax = 0; j < 3; j++)
    if((v = fabs(P[i].Vel[j])) > vmax)
      vmax = v;

  no = Father[i];

  while(no >= 0)
    {
      force_drift_node(no, All.Ti_Current);

      for(j = 0; j < 3; j++)
	{
	  Extnodes[no].dp[j] += dp[j];
#ifdef SCALARFIELD
	  Extnodes[no].dp_dm[j] += dp_dm[j];
#endif
	}

      if(Extnodes[no].vmax < vmax)
	Extnodes[no].vmax = vmax;

      Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);

      Extnodes[no].Ti_lastkicked = All.Ti_Current;

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* top-level tree-node reached */
	{
	  if(Extnodes[no].Flag != GlobFlag)
	    {
	      Extnodes[no].Flag = GlobFlag;
	      DomainList[DomainNumChanged++] = no;
	    }
	  break;
	}

      no = Nodes[no].u.d.father;
    }
}






void force_finish_kick_nodes(void)
{
  int i, j, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *counts_dp, *offset_list, *offset_dp, *offset_vmax;
  MyLongDouble *domainDp_loc, *domainDp_all;

#ifdef SCALARFIELD
  MyLongDouble *domainDp_dm_loc, *domainDp_dm_all;
#endif
  MyFloat *domainVmax_loc, *domainVmax_all;

  /* share the momentum-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  counts_dp = (int *) mymalloc("counts_dp", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_dp = (int *) mymalloc("offset_dp", sizeof(int) * NTask);
  offset_vmax = (int *) mymalloc("offset_vmax", sizeof(int) * NTask);

  domainDp_loc = (MyLongDouble *) mymalloc("domainDp_loc", DomainNumChanged * 3 * sizeof(MyLongDouble));
#ifdef SCALARFIELD
  domainDp_dm_loc = (MyLongDouble *) mymalloc("domainDp_dm_loc", DomainNumChanged * 3 * sizeof(MyLongDouble));
#endif
  domainVmax_loc = (MyFloat *) mymalloc("domainVmax_loc", DomainNumChanged * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  domainDp_loc[i * 3 + j] = Extnodes[DomainList[i]].dp[j];
#ifdef SCALARFIELD
	  domainDp_dm_loc[i * 3 + j] = Extnodes[DomainList[i]].dp_dm[j];
#endif
	}
      domainVmax_loc[i] = Extnodes[DomainList[i]].vmax;
    }

  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_dp[0] = 0, offset_vmax[0] = 0; ta < NTask;
      ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_dp[ta] = offset_dp[ta - 1] + counts[ta - 1] * 3 * sizeof(MyLongDouble);
	  offset_vmax[ta] = offset_vmax[ta - 1] + counts[ta - 1] * sizeof(MyFloat);
	}
    }

  if(ThisTask == 0)
    {
      printf("I exchange kick momenta for %d top-level nodes out of %d\n", totDomainNumChanged, NTopleaves);
    }

  domainDp_all = (MyLongDouble *) mymalloc("domainDp_all", totDomainNumChanged * 3 * sizeof(MyLongDouble));
#ifdef SCALARFIELD
  domainDp_dm_all =
    (MyLongDouble *) mymalloc("domainDp_dm_all", totDomainNumChanged * 3 * sizeof(MyLongDouble));
#endif
  domainVmax_all = (MyFloat *) mymalloc("domainVmax_all", totDomainNumChanged * sizeof(MyFloat));

  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    {
      counts_dp[ta] = counts[ta] * 3 * sizeof(MyLongDouble);
      counts[ta] *= sizeof(MyFloat);
    }


  MPI_Allgatherv(domainDp_loc, DomainNumChanged * 3 * sizeof(MyLongDouble), MPI_BYTE,
		 domainDp_all, counts_dp, offset_dp, MPI_BYTE, MPI_COMM_WORLD);

#ifdef SCALARFIELD
  MPI_Allgatherv(domainDp_dm_loc, DomainNumChanged * 3 * sizeof(MyLongDouble), MPI_BYTE,
		 domainDp_dm_all, counts_dp, offset_dp, MPI_BYTE, MPI_COMM_WORLD);
#endif

  MPI_Allgatherv(domainVmax_loc, DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainVmax_all, counts, offset_vmax, MPI_BYTE, MPI_COMM_WORLD);


  /* construct momentum kicks in top-level tree */
  for(i = 0; i < totDomainNumChanged; i++)
    {
      no = domainList_all[i];

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the local one is kicked twice */
	no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);

	  for(j = 0; j < 3; j++)
	    {
	      Extnodes[no].dp[j] += domainDp_all[3 * i + j];
#ifdef SCALARFIELD
	      Extnodes[no].dp_dm[j] += domainDp_dm_all[3 * i + j];
#endif
	    }

	  if(Extnodes[no].vmax < domainVmax_all[i])
	    Extnodes[no].vmax = domainVmax_all[i];

	  Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);
	  Extnodes[no].Ti_lastkicked = All.Ti_Current;

	  no = Nodes[no].u.d.father;
	}
    }

  myfree(domainList_all);
  myfree(domainVmax_all);
#ifdef SCALARFIELD
  myfree(domainDp_dm_all);
#endif
  myfree(domainDp_all);
  myfree(domainVmax_loc);
#ifdef SCALARFIELD
  myfree(domainDp_dm_loc);
#endif
  myfree(domainDp_loc);
  myfree(offset_vmax);
  myfree(offset_dp);
  myfree(offset_list);
  myfree(counts_dp);
  myfree(counts);
}



void force_drift_node(int no, int time1)
{
  int j, time0;
  double dt_drift, dt_drift_hmax, fac;

  if(time1 == Nodes[no].Ti_current)
    return;

  time0 = Extnodes[no].Ti_lastkicked;

  if(Nodes[no].u.d.bitflags & (1 << BITFLAG_NODEHASBEENKICKED))
    {
      if(Extnodes[no].Ti_lastkicked != Nodes[no].Ti_current)
	{
	  printf("Task=%d Extnodes[no].Ti_lastkicked=%d  Nodes[no].Ti_current=%d\n",
		 ThisTask, (int) Extnodes[no].Ti_lastkicked, (int) Nodes[no].Ti_current);
	  terminate("inconsistency in drift node");
	}

      if(Nodes[no].u.d.mass)
	fac = 1 / Nodes[no].u.d.mass;
      else
	fac = 0;

#ifdef SCALARFIELD
      double fac_dm;

      if(Nodes[no].mass_dm)
	fac_dm = 1 / Nodes[no].mass_dm;
      else
	fac_dm = 0;
#endif

      for(j = 0; j < 3; j++)
	{
	  Extnodes[no].vs[j] += fac * FLT(Extnodes[no].dp[j]);
	  Extnodes[no].dp[j] = 0;
#ifdef SCALARFIELD
	  Extnodes[no].vs_dm[j] += fac_dm * FLT(Extnodes[no].dp_dm[j]);
	  Extnodes[no].dp_dm[j] = 0;
#endif

#ifdef FLTROUNDOFFREDUCTION
	  Extnodes[no].s_base[j] = Nodes[no].u.d.s[j];
#ifdef SCALARFIELD
	  Extnodes[no].s_dm_base[j] = Nodes[no].s_dm[j];
#endif
#endif
	}
#ifdef FLTROUNDOFFREDUCTION
      Extnodes[no].len_base = Nodes[no].len;
#endif
      Nodes[no].u.d.bitflags &= (~(1 << BITFLAG_NODEHASBEENKICKED));
    }

  if(All.ComovingIntegrationOn)
    {
      dt_drift_hmax = get_drift_factor(Nodes[no].Ti_current, time1);
#ifdef FLTROUNDOFFREDUCTION
      dt_drift = get_drift_factor(time0, time1);
#else
      dt_drift = dt_drift_hmax;
#endif
    }
  else
    {

      dt_drift_hmax = (time1 - Nodes[no].Ti_current) * All.Timebase_interval;
#ifdef FLTROUNDOFFREDUCTION
      dt_drift = (time1 - time0) * All.Timebase_interval;
#else
      dt_drift = dt_drift_hmax;
#endif
    }

#ifdef FLTROUNDOFFREDUCTION
  for(j = 0; j < 3; j++)
    Nodes[no].u.d.s[j] = Extnodes[no].s_base[j] + Extnodes[no].vs[j] * dt_drift;
  Nodes[no].len = Extnodes[no].len_base + 2 * Extnodes[no].vmax * dt_drift;

#ifdef SCALARFIELD
  for(j = 0; j < 3; j++)
    Nodes[no].s_dm[j] = Extnodes[no].s_dm_base[j] + Extnodes[no].vs_dm[j] * dt_drift;
#endif

#else
  for(j = 0; j < 3; j++)
    Nodes[no].u.d.s[j] += Extnodes[no].vs[j] * dt_drift;
  Nodes[no].len += 2 * Extnodes[no].vmax * dt_drift;

#ifdef SCALARFIELD
  for(j = 0; j < 3; j++)
    Nodes[no].s_dm[j] += Extnodes[no].vs_dm[j] * dt_drift;
#endif

#endif

#ifndef KD_HMAX_ESTIMATE
  Extnodes[no].hmax *= exp(0.333333333333 * Extnodes[no].divVmax * dt_drift_hmax);
#endif

  Nodes[no].Ti_current = time1;
}





/*! This function updates the hmax-values in tree nodes that hold SPH
 *  particles. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the Hsml-values are potentially changed
 *  in the SPH-denity computation, force_update_hmax() should be carried
 *  out just before the hydrodynamical SPH forces are computed, i.e. after
 *  density().
 */
void force_update_hmax(void)
{
  int i, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *offset_list, *offset_hmax;
  MyFloat *domainHmax_loc, *domainHmax_all;
#ifdef KD_HMAX_ESTIMATE
  int OffsetSIZE = 1;
  integertime ti_step, time0, time1;
  double dt_drift_hmax;
  MyFloat hmax_est;
#else
  int OffsetSIZE = 2;
#endif

  GlobFlag++;

  DomainNumChanged = 0;
  DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
#ifdef BLACK_HOLES
    if(P[i].Type == 0 && P[i].Mass > 0)
#else
    if(P[i].Type == 0)
#endif
      {
	no = Father[i];

	while(no >= 0)
	  {
	    force_drift_node(no, All.Ti_Current);

#ifdef KD_HMAX_ESTIMATE
	    ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
	    time0 = All.Ti_Current;
	    time1 = P[i].Ti_begstep + ti_step;
	    dt_drift_hmax = get_drift_factor(time0, time1);
	    hmax_est = PPP[i].Hsml * exp(0.333333333333 * SphP[i].v.DivVel * dt_drift_hmax);
	    if(PPP[i].Hsml > Extnodes[no].hmax || hmax_est > Extnodes[no].hmax)
	      {
		if(PPP[i].Hsml > Extnodes[no].hmax)
		  Extnodes[no].hmax = PPP[i].Hsml;

		if(hmax_est > Extnodes[no].hmax)
		  Extnodes[no].hmax = hmax_est;
#else
	    if(PPP[i].Hsml > Extnodes[no].hmax || SphP[i].v.DivVel > Extnodes[no].divVmax)
	      {
		if(PPP[i].Hsml > Extnodes[no].hmax)
		  Extnodes[no].hmax = PPP[i].Hsml;

		if(SphP[i].v.DivVel > Extnodes[no].divVmax)
		  Extnodes[no].divVmax = SphP[i].v.DivVel;
#endif
		if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node */
		  {
		    if(Extnodes[no].Flag != GlobFlag)
		      {
			Extnodes[no].Flag = GlobFlag;
			DomainList[DomainNumChanged++] = no;
		      }
		    break;
		  }
	      }
	    else
	      break;

	    no = Nodes[no].u.d.father;
	  }
      }

  /* share the hmax-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_hmax = (int *) mymalloc("offset_hmax", sizeof(int) * NTask);

  domainHmax_loc = (MyFloat *) mymalloc("domainHmax_loc", DomainNumChanged * OffsetSIZE * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      domainHmax_loc[OffsetSIZE * i] = Extnodes[DomainList[i]].hmax;
#ifndef KD_HMAX_ESTIMATE
      domainHmax_loc[OffsetSIZE * i + 1] = Extnodes[DomainList[i]].divVmax;
#endif
    }


  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_hmax[0] = 0; ta < NTask; ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_hmax[ta] = offset_hmax[ta - 1] + counts[ta - 1] * OffsetSIZE * sizeof(MyFloat);
	}
    }

  if(ThisTask == 0)
    printf("Hmax exchange: %d topleaves out of %d\n", totDomainNumChanged, NTopleaves);

  domainHmax_all = (MyFloat *) mymalloc("domainHmax_all", totDomainNumChanged * OffsetSIZE * sizeof(MyFloat));
  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    counts[ta] *= OffsetSIZE * sizeof(MyFloat);

  MPI_Allgatherv(domainHmax_loc, OffsetSIZE * DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainHmax_all, counts, offset_hmax, MPI_BYTE, MPI_COMM_WORLD);


  for(i = 0; i < totDomainNumChanged; i++)
    {
      no = domainList_all[i];

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the hmax is updated twice */
	no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);


#ifdef KD_HMAX_ESTIMATE
	  if(domainHmax_all[OffsetSIZE * i] > Extnodes[no].hmax)
	    Extnodes[no].hmax = domainHmax_all[OffsetSIZE * i];
#else
	  if(domainHmax_all[OffsetSIZE * i] > Extnodes[no].hmax
	     || domainHmax_all[OffsetSIZE * i + 1] > Extnodes[no].divVmax)
	    {
	      if(domainHmax_all[OffsetSIZE * i] > Extnodes[no].hmax)
		Extnodes[no].hmax = domainHmax_all[OffsetSIZE * i];

	      if(domainHmax_all[OffsetSIZE * i + 1] > Extnodes[no].divVmax)
		Extnodes[no].divVmax = domainHmax_all[OffsetSIZE * i + 1];
	    }
#endif
	  else
	    break;

	  no = Nodes[no].u.d.father;
	}
    }


  myfree(domainList_all);
  myfree(domainHmax_all);
  myfree(domainHmax_loc);
  myfree(offset_hmax);
  myfree(offset_list);
  myfree(counts);
  myfree(DomainList);

  CPU_Step[CPU_TREEHMAXUPDATE] += measure_time();
}
