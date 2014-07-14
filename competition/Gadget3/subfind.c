#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#if defined(HAVE_HDF5)&& defined(WRITE_SUB_IN_SNAP_FORMAT)
#include <hdf5.h>
#endif

#include "fof.h"

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "subfind.h"

#ifdef SUBFIND

static struct subfind_id_list
{
  MyIDType ID;
  int GrNr;
  int SubNr;
  float BindingEgy;

#ifdef SUBFIND_SAVE_PARTICLELISTS
  float Pos[3];
  float Vel[3];
  int Type;
#ifdef STELLARAGE
  float Mass;
  float StellarAge;
#endif
#endif
}
 *ID_list;

static int Nids;


void subfind(int num)
{
  double t0, t1, tstart, tend;
  int i, gr, nlocid, offset, limit, ncount, ntotingrouplocal, nminingrouplocal, nmaxingrouplocal;

#ifdef DENSITY_SPLIT_BY_TYPE
  struct unbind_data *d;
  int j, n, count[6], countall[6];
  double a3inv;
#endif

  if(ThisTask == 0)
    printf("\nWe now execute a parallel version of SUBFIND.\n");

  tstart = second();

#ifdef LT_ADD_GAL_TO_SUB
  if(ThisTask == 0)
    printf("  loading CB07 tables ...\n");

  int IMFt = 1;

  Filters_Effective_Wavelenght = mymalloc("Filters_Effective_Wavelenght", sizeof(float) * LT_ADD_GAL_TO_SUB);
  tempiAS = mymalloc("tempiAS", sizeof(float) * 219);
  CB07 = mymalloc("CB07", sizeof(float) * 7 * 219 * LT_ADD_GAL_TO_SUB);
#ifdef OBSERVER_FRAME
  CB07obs = mymalloc("CB07obs", sizeof(float) * 7 * 219 * LT_ADD_GAL_TO_SUB);
#endif
  if(ThisTask == 0)
    load_CB07_table(IMFt, num);

  MPI_Bcast(&Filters_Effective_Wavelenght[0], sizeof(float) * LT_ADD_GAL_TO_SUB, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tempiAS[0], sizeof(float) * 219, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&CB07[0], sizeof(float) * 7 * 219 * LT_ADD_GAL_TO_SUB, MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef OBSERVER_FRAME
  MPI_Bcast(&CB07obs[0], sizeof(float) * 7 * 219 * LT_ADD_GAL_TO_SUB, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
#endif

#ifdef DENSITY_SPLIT_BY_TYPE
  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  for(j = 0; j < 6; j++)
    count[j] = 0;

  /* let's count number of particles of selected species */
  for(i = 0; i < NumPart; i++)
    count[P[i].Type]++;

  MPI_Allreduce(count, countall, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* do first loop: basically just defining the hsml for different species */
  for(j = 0; j < 6; j++)
    {
      if((1 << j) & (DENSITY_SPLIT_BY_TYPE))
	{
#ifdef BLACK_HOLES
	  if(j == 5)
	    countall[j] = 0;	/* this will prevent that the black holes are treated separately */
#endif

	  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

	  if(countall[j] > All.DesNumNgb)
	    {
	      /* build index list of particles of selectes species */
	      d = (struct unbind_data *) mymalloc("	      d", count[j] * sizeof(struct unbind_data));
	      for(i = 0, n = 0; i < NumPart; i++)
		if(P[i].Type == j)
		  d[n++].index = i;

	      t0 = second();
	      if(ThisTask == 0)
		printf("Tree construction for species %d (%d).\n", j, countall[j]);

	      CPU_Step[CPU_FOF] += measure_time();

	      force_treebuild(count[j], d);

	      myfree(d);

	      t1 = second();
	      if(ThisTask == 0)
		printf("tree build for species %d took %g sec\n", j, timediff(t0, t1));
	    }
	  else
	    {
	      t0 = second();
	      if(ThisTask == 0)
		printf("Tree construction.\n");

	      CPU_Step[CPU_FOF] += measure_time();

	      force_treebuild(NumPart, NULL);

	      t1 = second();
	      if(ThisTask == 0)
		printf("tree build took %g sec\n", timediff(t0, t1));
	    }


	  /* let's determine the local densities */
	  t0 = second();
	  subfind_setup_smoothinglengths(j);
	  subfind_density(j);
	  t1 = second();
	  if(ThisTask == 0)
	    printf("density and smoothing length for species %d took %g sec\n", j, timediff(t0, t1));

	  force_treefree();

	  /* let's save density contribution of own species */
	  for(i = 0; i < NumPart; i++)
	    if(P[i].Type == j)
	      P[i].w.density_sum = P[i].u.DM_Density;

	}
    }

  /* do second loop: now calculate all density contributions */
  for(j = 0; j < 6; j++)
    {
      if((1 << j) & (DENSITY_SPLIT_BY_TYPE))
	{
	  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

	  /* build index list of particles of selectes species */
	  d = (struct unbind_data *) mymalloc("	  d", count[j] * sizeof(struct unbind_data));
	  for(i = 0, n = 0; i < NumPart; i++)
	    if(P[i].Type == j)
	      d[n++].index = i;

	  t0 = second();
	  if(ThisTask == 0)
	    printf("Tree construction for species %d (%d).\n", j, countall[j]);

	  CPU_Step[CPU_FOF] += measure_time();

	  force_treebuild(count[j], d);

	  myfree(d);

	  t1 = second();
	  if(ThisTask == 0)
	    printf("tree build for species %d took %g sec\n", j, timediff(t0, t1));

	  /* let's determine the local densities */
	  t0 = second();
	  for(i = 0; i < 6; i++)
	    if((1 << i) & (DENSITY_SPLIT_BY_TYPE))
	      if(j != i)
		{
		  if(countall[i] > All.DesNumNgb)
		    {
		      if(ThisTask == 0)
			printf("calculating density contribution of species %d to species %d\n", j, i);
		      subfind_density(-(i + 1));
		    }
		}
	  t1 = second();
	  if(ThisTask == 0)
	    printf("density() of species %d took %g sec\n", j, timediff(t0, t1));

	  force_treefree();

	  /* let's sum up density contribution */
	  for(i = 0; i < NumPart; i++)
	    if((1 << P[i].Type) & (DENSITY_SPLIT_BY_TYPE))
	      if(j != P[i].Type)
		if(countall[P[i].Type] > All.DesNumNgb)
		  P[i].w.density_sum += P[i].u.DM_Density;
	}
    }


  for(i = 0; i < NumPart; i++)
    {
      P[i].u.DM_Density = P[i].w.density_sum;

      if(P[i].Type == 0)
	P[i].w.int_energy = DMAX(All.MinEgySpec,
				 SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv,
								      GAMMA_MINUS1));
      else
	P[i].w.int_energy = 0;

    }
#else
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  t0 = second();
  if(ThisTask == 0)
    printf("Tree construction.\n");

  CPU_Step[CPU_FOF] += measure_time();

  force_treebuild(NumPart, NULL);

  t1 = second();
  if(ThisTask == 0)
    printf("tree build took %g sec\n", timediff(t0, t1));


  /* let's determine the local dark matter densities */
  t0 = second();
  subfind_setup_smoothinglengths();
  subfind_density();
  t1 = second();
  if(ThisTask == 0)
    printf("dark matter density() took %g sec\n", timediff(t0, t1));

  force_treefree();
#endif /* DENSITY_SPLIT_BY_TYPE */

#ifndef SUBFIND_DENSITY_AND_POTENTIAL
#ifndef NO_SAVE_OF_HSML
  if(DumpFlag == 1)
    {
      /* let's save the densities to a file (for making images) */
      t0 = second();
      subfind_save_densities(num);
      t1 = second();
      if(ThisTask == 0)
	printf("saving densities took %g sec\n", timediff(t0, t1));
    }
#endif
#endif
#ifdef ONLY_PRODUCE_HSML_FILES
  return;
#endif

  /* count how many groups we have that should be done collectively */
#ifdef KD_CHOOSE_PSUBFIND_LIMIT
  limit = (int) (0.8 * All.PartAllocFactor * All.TotNumPart / NTask);
#else
  limit = (int) (0.6 * All.TotNumPart / NTask);
#endif

  for(i = 0, ncount = 0, ntotingrouplocal = 0; i < Ngroups; i++)
    if(Group[i].Len >= limit)
      ncount++;
    else
      ntotingrouplocal += Group[i].Len;

  MPI_Allreduce(&ncount, &Ncollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ntotingrouplocal, &nminingrouplocal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&ntotingrouplocal, &nmaxingrouplocal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nNumber of FOF halos treated with collective SubFind code = %d\n", Ncollective);
      printf("(the adopted size-limit for the collective algorithm was %d particles.)\n", limit);
      printf("the other %d FOF halos are treated in parallel with serial code\n\n", TotNgroups - Ncollective);
      printf("Unbalance in total number of particles in FOF halos is %d - %d \n\n", nminingrouplocal,
	     nmaxingrouplocal);
    }

  /*  to decide on which task a group should be:
   *  if   GrNr <= Ncollective:  collective groupfinding.
   *  the task where the group info is put is TaskNr = (GrNr - 1) % NTask
   */

  /* now we distribute the particles such that small groups are assigned in
   *  total to certain CPUs, and big groups are left where they are 
   */

  t0 = second();

  for(i = 0; i < NumPart; i++)
    {
      P[i].origintask2 = ThisTask;

      if(P[i].GrNr > Ncollective && P[i].GrNr <= TotNgroups)	/* particle is in small group */
	P[i].targettask = (P[i].GrNr - 1) % NTask;
      else
	P[i].targettask = ThisTask;
    }

  subfind_exchange();		/* distributes gas particles as well if needed */

  t1 = second();
  if(ThisTask == 0)
    printf("subfind_exchange()() took %g sec\n", timediff(t0, t1));

  subfind_distribute_groups();

#ifdef OMP_SORT
  omp_qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);
#else
  qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].GrNr > Ncollective && P[i].GrNr <= TotNgroups)
      if(((P[i].GrNr - 1) % NTask) != ThisTask)
	{
	  printf("i=%d %d task=%d type=%d\n", i, P[i].GrNr, ThisTask, P[i].Type);
	  endrun(87);
	}

  /* lets estimate the maximum number of substructures we need to store on the local CPU */
  for(i = 0, nlocid = 0; i < Ngroups; i++)
    nlocid += Group[i].Len;

  MaxNsubgroups = nlocid / All.DesLinkNgb + NTask;	/* this is a quite conservative upper limit */
  Nsubgroups = 0;
  SubGroup =
    (struct subgroup_properties *) mymalloc("SubGroup", MaxNsubgroups * sizeof(struct subgroup_properties));

  for(i = 0; i < NumPart; i++)
    P[i].SubNr = (1 << 30);	/* default */

  /* we begin by applying the collective version of subfind to distributed groups */
  t0 = second();
  for(GrNr = 1; GrNr <= Ncollective; GrNr++)
    subfind_process_group_collectively(num);
  t1 = second();
  if(ThisTask == 0)
    printf("processing of collective halos took %g sec\n", timediff(t0, t1));

#ifdef SUBFIND_COLLECTIVE_STAGE1
  if(ThisTask == 0)
    printf("stage 1 ended\n");
  endrun(0);
#endif

  for(i = 0; i < NumPart; i++)
    {
      P[i].origindex = i;
      P[i].origintask = ThisTask;
    }

  t0 = second();
#ifdef OMP_SORT
  omp_qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_GrNr_DM_Density);
#else
  qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_GrNr_DM_Density);
#endif
  t1 = second();
  if(ThisTask == 0)
    printf("sort of local particles()() took %g sec\n", timediff(t0, t1));


  /* now we have the particles of groups consecutively, but SPH particles are
     not aligned. They can however be accessed via SphP[P[i].originindex] */


  /* let's count how many local particles we have in small groups */
  for(i = 0, nlocid = 0; i < NumPart; i++)
    if(P[i].GrNr > Ncollective && P[i].GrNr <= Ngroups)	/* particle is in small group */
      nlocid++;

  if(ThisTask == 0)
    printf("contructing tree for serial subfind of local groups\n");

  subfind_loctree_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  if(ThisTask == 0)
    printf("Start to do local groups with serial subfind algorithm\n");

  t0 = second();

  /* we now apply a serial version of subfind to the local groups */
  for(gr = 0, offset = 0; gr < Ngroups; gr++)
    {
      if(Group[gr].GrNr > Ncollective)
	{
	  if(((Group[gr].GrNr - 1) % NTask) == ThisTask)
	    offset = subfind_process_group_serial(gr, offset);
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);

  t1 = second();
  if(ThisTask == 0)
    printf("\nprocessing of local groups took took %g sec\n\n", timediff(t0, t1));


  subfind_loctree_treefree();


  /* bringing back particles in original positions, such that gas particles are aligned */
  t0 = second();
#ifdef OMP_SORT
  omp_qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_origindex);
#else
  qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_origindex);
#endif
  t1 = second();
  if(ThisTask == 0)
    printf("unsorting of local particles()() took %g sec\n", timediff(t0, t1));


  GrNr = -1;			/* to ensure that domain decomposition acts normally again */

  /* now determine the remaining spherical overdensity values for the non-local groups */

  domain_free_trick();

  CPU_Step[CPU_FOF] += measure_time();


#ifdef DENSITY_SPLIT_BY_TYPE
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].origintask != ThisTask)
	printf("Task %d: Holding particle of task %d !\n", ThisTask, P[i].origintask);
      if(P[i].origindex != i)
	printf("Task %d: Particles is in wrong position (is=%d, was=%d) !\n", ThisTask, i, P[i].origindex);
    }
#endif



  t0 = second();

  for(i = 0; i < NumPart; i++)
    P[i].targettask = P[i].origintask2;

  subfind_exchange();		/* distributes gas particles as well if needed */

  t1 = second();
  if(ThisTask == 0)
    printf("subfind_exchange() (for return to original CPU)  took %g sec\n", timediff(t0, t1));



  All.DoDynamicUpdate = 0;

  domain_Decomposition(1, 0);

  force_treebuild(NumPart, NULL);


  /* compute spherical overdensities for FOF groups */
  t0 = second();

  subfind_overdensity();

  t1 = second();
  if(ThisTask == 0)
    printf("determining spherical overdensity masses took %g sec\n", timediff(t0, t1));


  /* determine which halos are contaminated by boundary particles */
  t0 = second();

  subfind_contamination();

  t1 = second();
  if(ThisTask == 0)
    printf("determining contamination of halos took %g sec\n", timediff(t0, t1));


  force_treefree();
  domain_free();

  domain_allocate_trick();

  /* now assemble final output */
  subfind_save_final(num);

  tend = second();

  if(ThisTask == 0)
    printf("\nFinished with SUBFIND.  (total time=%g sec)\n\n", timediff(tstart, tend));

  myfree(SubGroup);

#ifdef   LT_ADD_GAL_TO_SUB
#ifdef OBSERVER_FRAME
  myfree(CB07obs);
#endif
  myfree(CB07);
  myfree(tempiAS);
  myfree(Filters_Effective_Wavelenght);
#endif

  CPU_Step[CPU_FOF] += measure_time();
}




void subfind_save_final(int num)
{
  int i, j, totsubs, masterTask, groupTask, nprocgroup;
  char buf[1000];
  double t0, t1;

  /* prepare list of ids with assigned group numbers */

  t0 = second();

#ifdef KD_ALTERNATIVE_GROUP_SORT
  sort_groups();
#ifdef OMP_SORT
  omp_qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);
#else
  qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);
#endif
  t1 = second();
  if(ThisTask == 0)
    {
      printf("Alternative Global sort of Groups took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
  t0 = second();
  sort_sub_groups();
#ifdef OMP_SORT
  omp_qsort(SubGroup, Nsubgroups, sizeof(struct subgroup_properties), subfind_compare_SubGroup_GrNr_SubNr);
#else
  qsort(SubGroup, Nsubgroups, sizeof(struct subgroup_properties), subfind_compare_SubGroup_GrNr_SubNr);
#endif
  t1 = second();
  if(ThisTask == 0)
    {
      printf("Alternative Global sort of SubGroups took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

#else // KD_ALTERNATIVE_GROUP_SORT
  t0 = second();
  parallel_sort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);
  t1 = second();
  if(ThisTask == 0)
    {
      printf("Global sort of Groups took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  t0 = second();
  parallel_sort(SubGroup, Nsubgroups, sizeof(struct subgroup_properties),
		subfind_compare_SubGroup_GrNr_SubNr);
  t1 = second();
  if(ThisTask == 0)
    {
      printf("Global sort of SubGroups took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
#endif

  ID_list =
    (struct subfind_id_list *) mymalloc("ID_list",
					sizeof(struct subfind_id_list) * (int) (All.MaxPart * 1.5));

  for(i = 0, Nids = 0; i < NumPart; i++)
    {
      if(P[i].GrNr <= TotNgroups)
	{
	  ID_list[Nids].GrNr = P[i].GrNr;
	  ID_list[Nids].SubNr = P[i].SubNr;
	  ID_list[Nids].BindingEgy = P[i].v.DM_BindingEnergy;
	  ID_list[Nids].ID = P[i].ID;
#ifdef SUBFIND_SAVE_PARTICLELISTS
	  for(j = 0; j < 3; j++)
	    {
	      ID_list[Nids].Pos[j] = P[i].Pos[j];
	      ID_list[Nids].Vel[j] = P[i].Vel[j];
	    }
	  ID_list[Nids].Type = P[i].Type;
#ifdef STELLARAGE
	  ID_list[Nids].Mass = P[i].Mass;
	  if(P[i].Type == 4)
	    ID_list[Nids].StellarAge = P[i].StellarAge;
	  else
	    ID_list[Nids].StellarAge = 0;
#endif
#endif
	  Nids++;
	}
    }

  t0 = second();

#ifdef KD_ALTERNATIVE_GROUP_SORT
  sort_ids();
#ifdef OMP_SORT
  omp_qsort(ID_list, Nids, sizeof(struct subfind_id_list), subfind_compare_ID_list);
#else
  qsort(ID_list, Nids, sizeof(struct subfind_id_list), subfind_compare_ID_list);
#endif
  t1 = second();
  if(ThisTask == 0)
    {
      printf("Alternative Global sort of IDs took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
#else // KD_ALTERNATIVE_GROUP_SORT

  parallel_sort(ID_list, Nids, sizeof(struct subfind_id_list), subfind_compare_ID_list);
  t1 = second();
  if(ThisTask == 0)
    {
      printf("Global sort of IDs took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
#endif

  MPI_Allreduce(&Nsubgroups, &TotNsubgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  /* fill in the FirstSub-values */
  for(i = 0, totsubs = 0; i < Ngroups; i++)
    {
      if(i > 0)
	Group[i].FirstSub = Group[i - 1].FirstSub + Group[i - 1].Nsubs;
      else
	Group[i].FirstSub = 0;
      totsubs += Group[i].Nsubs;
    }

  MPI_Allgather(&totsubs, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  for(i = 0; i < Ngroups; i++)
    Group[i].FirstSub += Send_offset[ThisTask];


#ifdef SUBFIND_COUNT_BIG_HALOS
  for(i = 0, Nbiggroups = 0; i < Ngroups; i++)
    if(Group[i].M_TopHat200 > SUBFIND_COUNT_BIG_HALOS)
      Nbiggroups++;
  MPI_Allreduce(&Nbiggroups, &TotNbiggroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  MPI_Allgather(&Nids, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }

  t0 = second();

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	subfind_save_local_catalogue(num);
      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }

  t1 = second();

  if(ThisTask == 0)
    {
      printf("Subgroup catalogues saved. took = %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  myfree(ID_list);
}


int get_sub_entrytype_of_block(enum siofields blocknr)
{
  int type = 0;

  switch (blocknr)
    {
    case SIO_SLEN:
    case SIO_SOFF:
    case SIO_PFOF:
    case SIO_MSUB:
    case SIO_SPOS:
    case SIO_SVEL:
    case SIO_SCM:
    case SIO_SPIN:
    case SIO_DSUB:
    case SIO_VMAX:
    case SIO_RVMAX:
    case SIO_RHMS:
    case SIO_MBID:
    case SIO_GRNR:
    case SIO_SMST:
    case SIO_SLUM:
    case SIO_SLATT:
    case SIO_SLOBS:
    case SIO_DUST:
    case SIO_SAGE:
    case SIO_SZ:
    case SIO_SSFR:
      type = 1;
      break;
    case SIO_PPOS:
    case SIO_PVEL:
    case SIO_PTYP:
    case SIO_PMAS:
    case SIO_PAGE:
    case SIO_PID:
      type = 2;
      break;
    case SIO_BGPOS:
    case SIO_BGMTOP:
    case SIO_BGRTOP:
      type = 3;
      break;
    default:
      type = 0;
      break;
    }

  return type;
}

int get_values_per_sub(enum siofields blocknr)
{
  int n = 1;

  switch (blocknr)
    {
    case SIO_GPOS:
    case SIO_SPOS:
    case SIO_SVEL:
    case SIO_SCM:
    case SIO_SPIN:
    case SIO_PPOS:
    case SIO_PVEL:
    case SIO_MGAS:
    case SIO_MSTR:
    case SIO_TGAS:
    case SIO_LGAS:
    case SIO_BGPOS:
      n = 3;
      break;
    case SIO_SMST:
      n = 6;
      break;
    case SIO_SLUM:
    case SIO_SLATT:
    case SIO_SLOBS:
#ifdef LT_ADD_GAL_TO_SUB
      n = LT_ADD_GAL_TO_SUB;
#else
      n = 1;
#endif
      break;
    case SIO_DUST:
#ifdef DUSTATT
      n = DUSTATT;
#else
      n = 1;
#endif
      break;
    default:
      n = 1;
      break;
    }

  return n;
}

#ifdef HAVE_HDF5
void get_IO_Label_HDF5_sub(enum siofields blocknr, char *label)
{
  switch (blocknr)
    {
    case SIO_GLEN:
      strcpy(label, "GroupLen");
      break;
    case SIO_GOFF:
      strcpy(label, "GroupOffset");
      break;
    case SIO_MTOT:
      strcpy(label, "GroupMass");
      break;
    case SIO_GPOS:
      strcpy(label, "GroupPos");
      break;
#ifdef SO_BAR_INFO
    case SIO_MMEA:
      strcpy(label, "MVIR");
      break;
    case SIO_RMEA:
      strcpy(label, "RVIR");
      break;
    case SIO_DMEA:
      strcpy(label, "DVIR");
      break;
    case SIO_MCRI:
      strcpy(label, "M25K");
      break;
    case SIO_RCRI:
      strcpy(label, "R25K");
      break;
    case SIO_DCRI:
      strcpy(label, "D25K");
      break;
    case SIO_MTOP:
      strcpy(label, "M500");
      break;
    case SIO_RTOP:
      strcpy(label, "R500");
      break;
    case SIO_DTOP:
      strcpy(label, "D500");
      break;
#else
    case SIO_MMEA:
      strcpy(label, "Group_M_Mean200");
      break;
    case SIO_RMEA:
      strcpy(label, "Group_R_Mean200");
      break;
    case SIO_MCRI:
      strcpy(label, "Group_M_Crit200");
      break;
    case SIO_RCRI:
      strcpy(label, "Group_R_Crit200");
      break;
    case SIO_MTOP:
      strcpy(label, "Group_M_TopHat200");
      break;
    case SIO_RTOP:
      strcpy(label, "Group_R_TopHat200");
      break;
    case SIO_DMEA:
      strcpy(label, "Group_VelDisp_Mean200");
      break;
    case SIO_DCRI:
      strcpy(label, "Group_VelDisp_Crit200");
      break;
    case SIO_DTOP:
      strcpy(label, "Group_VelDisp_TopHat200");
      break;
#endif
    case SIO_MGAS:
      strcpy(label, "MGAS");
      break;
    case SIO_MSTR:
      strcpy(label, "MSTR");
      break;
    case SIO_TGAS:
      strcpy(label, "TGAS");
      break;
    case SIO_LGAS:
      strcpy(label, "LGAS");
      break;
    case SIO_NCON:
      strcpy(label, "GroupContaminationCoun");
      break;
    case SIO_MCON:
      strcpy(label, "GroupContaminationMass");
      break;
    case SIO_BGPOS:
      strcpy(label, "BGPO");
      break;
    case SIO_BGMTOP:
      strcpy(label, "BGMA");
      break;
    case SIO_BGRTOP:
      strcpy(label, "BGRA");
      break;
    case SIO_NSUB:
      strcpy(label, "GroupNsubs");
      break;
    case SIO_FSUB:
      strcpy(label, "GroupFirstSub");
      break;
    case SIO_SLEN:
      strcpy(label, "SubhaloLen");
      break;
    case SIO_SOFF:
      strcpy(label, "SubhaloOffset");
      break;
    case SIO_PFOF:
      strcpy(label, "SubhaloParent");
      break;
    case SIO_MSUB:
      strcpy(label, "SubhaloMass");
      break;
    case SIO_SPOS:
      strcpy(label, "SubhaloPos");
      break;
    case SIO_SVEL:
      strcpy(label, "SubhaloVel");
      break;
    case SIO_SCM:
      strcpy(label, "SubhaloCM");
      break;
    case SIO_SPIN:
      strcpy(label, "SubhaloSpin");
      break;
    case SIO_DSUB:
      strcpy(label, "SubhaloVelDisp");
      break;
    case SIO_VMAX:
      strcpy(label, "SubhaloVmax");
      break;
    case SIO_RVMAX:
      strcpy(label, "SubhaloVmaxRad");
      break;
    case SIO_RHMS:
      strcpy(label, "SubhaloHalfmassRad");
      break;
    case SIO_MBID:
      strcpy(label, "SubhaloIDMostbound");
      break;
    case SIO_GRNR:
      strcpy(label, "SubhaloGrNr");
      break;
    case SIO_SMST:
      strcpy(label, "SMST");
      break;
    case SIO_SLUM:
      strcpy(label, "SLUM");
      break;
    case SIO_SLATT:
      strcpy(label, "SLAT");
      break;
    case SIO_SLOBS:
      strcpy(label, "SLOB");
      break;
    case SIO_DUST:
      strcpy(label, "DUST");
      break;
    case SIO_SAGE:
      strcpy(label, "SAGE");
      break;
    case SIO_SZ:
      strcpy(label, "SZ  ");
      break;
    case SIO_SSFR:
      strcpy(label, "SSFR");
      break;
    case SIO_PPOS:
      strcpy(label, "PPOS");
      break;
    case SIO_PVEL:
      strcpy(label, "PVEL");
      break;
    case SIO_PTYP:
      strcpy(label, "PTYP");
      break;
    case SIO_PMAS:
      strcpy(label, "PMAS");
      break;
    case SIO_PAGE:
      strcpy(label, "PAGE");
      break;
    case SIO_PID:
      strcpy(label, "PID ");
      break;
    default:
      endrun(987453);
      break;
    }
}
#endif

void get_IO_Label_sub(enum siofields blocknr, char *label)
{
  switch (blocknr)
    {
    case SIO_GLEN:
      strncpy(label, "GLEN", 4);
      break;
    case SIO_GOFF:
      strncpy(label, "GOFF", 4);
      break;
    case SIO_MTOT:
      strncpy(label, "MTOT", 4);
      break;
    case SIO_GPOS:
      strncpy(label, "GPOS", 4);
      break;
#ifdef SO_BAR_INFO
    case SIO_MMEA:
      strncpy(label, "MVIR", 4);
      break;
    case SIO_RMEA:
      strncpy(label, "RVIR", 4);
      break;
    case SIO_DMEA:
      strncpy(label, "DVIR", 4);
      break;
    case SIO_MCRI:
      strncpy(label, "M25K", 4);
      break;
    case SIO_RCRI:
      strncpy(label, "R25K", 4);
      break;
    case SIO_DCRI:
      strncpy(label, "D25K", 4);
      break;
    case SIO_MTOP:
      strncpy(label, "M500", 4);
      break;
    case SIO_RTOP:
      strncpy(label, "R500", 4);
      break;
    case SIO_DTOP:
      strncpy(label, "D500", 4);
      break;
#else
    case SIO_MMEA:
      strncpy(label, "MMEA", 4);
      break;
    case SIO_RMEA:
      strncpy(label, "RMEA", 4);
      break;
    case SIO_MCRI:
      strncpy(label, "MCRI", 4);
      break;
    case SIO_RCRI:
      strncpy(label, "RCRI", 4);
      break;
    case SIO_MTOP:
      strncpy(label, "MTOP", 4);
      break;
    case SIO_RTOP:
      strncpy(label, "RTOP", 4);
      break;
    case SIO_DMEA:
      strncpy(label, "DMEA", 4);
      break;
    case SIO_DCRI:
      strncpy(label, "DCRI", 4);
      break;
    case SIO_DTOP:
      strncpy(label, "DTOP", 4);
      break;
#endif
    case SIO_MGAS:
      strncpy(label, "MGAS", 4);
      break;
    case SIO_MSTR:
      strncpy(label, "MSTR", 4);
      break;
    case SIO_TGAS:
      strncpy(label, "TGAS", 4);
      break;
    case SIO_LGAS:
      strncpy(label, "LGAS", 4);
      break;
    case SIO_NCON:
      strncpy(label, "NCON", 4);
      break;
    case SIO_MCON:
      strncpy(label, "MCON", 4);
      break;
    case SIO_BGPOS:
      strncpy(label, "BGPO", 4);
      break;
    case SIO_BGMTOP:
      strncpy(label, "BGMA", 4);
      break;
    case SIO_BGRTOP:
      strncpy(label, "BGRA", 4);
      break;
    case SIO_NSUB:
      strncpy(label, "NSUB", 4);
      break;
    case SIO_FSUB:
      strncpy(label, "FSUB", 4);
      break;
    case SIO_SLEN:
      strncpy(label, "SLEN", 4);
      break;
    case SIO_SOFF:
      strncpy(label, "SOFF", 4);
      break;
    case SIO_PFOF:
      strncpy(label, "SSUB", 4);
      break;
    case SIO_MSUB:
      strncpy(label, "MSUB", 4);
      break;
    case SIO_SPOS:
      strncpy(label, "SPOS", 4);
      break;
    case SIO_SVEL:
      strncpy(label, "SVEL", 4);
      break;
    case SIO_SCM:
      strncpy(label, "SCM ", 4);
      break;
    case SIO_SPIN:
      strncpy(label, "SPIN", 4);
      break;
    case SIO_DSUB:
      strncpy(label, "DSUB", 4);
      break;
    case SIO_VMAX:
      strncpy(label, "VMAX", 4);
      break;
    case SIO_RVMAX:
      strncpy(label, "RMAX", 4);
      break;
    case SIO_RHMS:
      strncpy(label, "RHMS", 4);
      break;
    case SIO_MBID:
      strncpy(label, "MBID", 4);
      break;
    case SIO_GRNR:
      strncpy(label, "GRNR", 4);
      break;
    case SIO_SMST:
      strncpy(label, "SMST", 4);
      break;
    case SIO_SLUM:
      strncpy(label, "SLUM", 4);
      break;
    case SIO_SLATT:
      strncpy(label, "SLAT", 4);
      break;
    case SIO_SLOBS:
      strncpy(label, "SLOB", 4);
      break;
    case SIO_DUST:
      strncpy(label, "DUST", 4);
      break;
    case SIO_SAGE:
      strncpy(label, "SAGE", 4);
      break;
    case SIO_SZ:
      strncpy(label, "SZ  ", 4);
      break;
    case SIO_SSFR:
      strncpy(label, "SSFR", 4);
      break;
    case SIO_PPOS:
      strncpy(label, "PPOS", 4);
      break;
    case SIO_PVEL:
      strncpy(label, "PVEL", 4);
      break;
    case SIO_PTYP:
      strncpy(label, "PTYP", 4);
      break;
    case SIO_PMAS:
      strncpy(label, "PMAS", 4);
      break;
    case SIO_PAGE:
      strncpy(label, "PAGE", 4);
      break;
    case SIO_PID:
      strncpy(label, "PID ", 4);
      break;
    default:
      endrun(987453);
      break;
    }
}

int get_datatype_in_sub(enum siofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case SIO_MBID:
    case SIO_PID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    case SIO_GLEN:
    case SIO_GOFF:
    case SIO_NCON:
    case SIO_NSUB:
    case SIO_FSUB:
    case SIO_SLEN:
    case SIO_SOFF:
    case SIO_PFOF:
    case SIO_GRNR:
    case SIO_PTYP:
      typekey = 0;		/* native int */
      break;

    default:
#ifdef OUTPUT_IN_DOUBLEPRECISION
      typekey = 3;
#else
      typekey = 1;		/* native MyOutputFloat */
#endif
      break;
    }

  return typekey;
}

int block_in_sub(enum siofields blocknr)
{
  int present = 0;

  switch (blocknr)
    {
    case SIO_GLEN:
    case SIO_GOFF:
    case SIO_MTOT:
    case SIO_GPOS:
    case SIO_MMEA:
    case SIO_RMEA:
    case SIO_MCRI:
    case SIO_RCRI:
    case SIO_MTOP:
    case SIO_RTOP:
#ifdef SO_VEL_DISPERSIONS
    case SIO_DMEA:
    case SIO_DCRI:
    case SIO_DTOP:
#endif
#ifdef SO_BAR_INFO
    case SIO_MGAS:
    case SIO_MSTR:
    case SIO_TGAS:
    case SIO_LGAS:
#endif
    case SIO_NCON:
    case SIO_MCON:
#ifdef SUBFIND_COUNT_BIG_HALOS
    case SIO_BGPOS:
    case SIO_BGMTOP:
    case SIO_BGRTOP:
#endif
    case SIO_NSUB:
    case SIO_FSUB:
    case SIO_SLEN:
    case SIO_SOFF:
    case SIO_PFOF:
    case SIO_MSUB:
    case SIO_SPOS:
    case SIO_SVEL:
    case SIO_SCM:
    case SIO_SPIN:
    case SIO_DSUB:
    case SIO_VMAX:
    case SIO_RVMAX:
    case SIO_RHMS:
    case SIO_MBID:
    case SIO_GRNR:
#ifdef SAVE_MASS_TAB
    case SIO_SMST:
#endif
#ifdef LT_ADD_GAL_TO_SUB
    case SIO_SLUM:
    case SIO_SAGE:
    case SIO_SZ:
    case SIO_SSFR:
#endif
#ifdef DUSTATT
    case SIO_SLATT:
    case SIO_DUST:
#endif
#ifdef OBSERVER_FRAME
    case SIO_SLOBS:
#endif
#ifdef SUBFIND_SAVE_PARTICLELISTS
    case SIO_PPOS:
    case SIO_PVEL:
    case SIO_PTYP:
    case SIO_PMAS:
#ifdef STELLARAGE
    case SIO_PAGE:
#endif
#endif
    case SIO_PID:
      present = 1;
      break;
    default:
      present = 0;
      break;
    }

  return present;
}


void subfind_save_local_catalogue(int num)
{
  FILE *fd;
  char buf[500], fname[500], label[] = "--------";

  void *IOBuffer;
  MyOutputFloat *fp;
  int *fp_int, nwrite, ndim, bytes_per_blockelement, bnr, type, datatype, bytes;
  MyIDType *fp_id;
  enum siofields blocknr;
  unsigned int blksize;

  int i, j;

#if defined(SUBFIND_COUNT_BIG_HALOS) || defined(LT_ADD_GAL_TO_SUB)
  int k;
#endif

#ifdef SUBFIND_SAVE_PARTICLELISTS
  double a3inv;
#endif

#ifdef WRITE_SUB_IN_SNAP_FORMAT
  unsigned int nextblock;
#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[3], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0;
#endif

#endif

#ifdef SUBFIND_SAVE_PARTICLELISTS
  if(All.ComovingIntegrationOn)
    a3inv = 1.0 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1.0;
#endif

#ifdef WRITE_SUB_IN_SNAP_FORMAT
  if(NTask == 1)
    sprintf(fname, "%s/groups_%03d/%s_%03d", All.OutputDir, num, "sub", num);
  else
    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "sub", num, ThisTask);
#else
  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, ThisTask);
#endif
  strcpy(buf, fname);

#if defined(HAVE_HDF5)&& defined(WRITE_SUB_IN_SNAP_FORMAT)
  if(All.SnapFormat == 3)
    {

      sprintf(buf, "%s.hdf5", fname);
      hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

      hdf5_grp[0] = H5Gcreate(hdf5_file, "/Group", 0);
      hdf5_grp[1] = H5Gcreate(hdf5_file, "/Subhalo", 0);
      hdf5_grp[2] = H5Gcreate(hdf5_file, "/IDs", 0);



    }
  else
    {
#endif
      if(!(fd = fopen(buf, "w")))

	{
	  printf("can't open file `%s`\n", buf);
	  endrun(1183);
	}
#if defined(HAVE_HDF5)&& defined(WRITE_SUB_IN_SNAP_FORMAT)
    }
#endif

#ifdef WRITE_SUB_IN_SNAP_FORMAT
  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.npartTotalHighWord[i] = 0;
      header.mass[i] = 0;
    }
  header.npart[0] = Ngroups;
  header.npartTotal[0] = TotNgroups;
  header.npart[1] = Nsubgroups;
  header.npartTotal[1] = TotNsubgroups;
  header.npart[2] = Nids;
  header.npartTotal[2] = (unsigned int) TotNids;
  header.npartTotalHighWord[2] = (unsigned int) (TotNids >> 32);
#ifdef SUBFIND_COUNT_BIG_HALOS
  header.npart[3] = Nbiggroups;
  header.npartTotal[3] = TotNbiggroups;
#endif

  header.time = All.Time;
  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

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

  header.num_files = NTask;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

#ifdef LT_ADD_GAL_TO_SUB
  for(k = 0; k < LT_ADD_GAL_TO_SUB; k++)
    {
      header.names[k][0] = bc_name[bc_use_index[k]][0];
      header.names[k][1] = bc_name[bc_use_index[k]][1];
    }
  for(k = LT_ADD_GAL_TO_SUB; k < 15; k++)
    {
      header.names[k][0] = ' ';
      header.names[k][1] = ' ';
    }
#endif

#ifdef WRITE_SUB_IN_SNAP_FORMAT
  if(All.SnapFormat == 2)
    {
      blksize = sizeof(int) + 4 * sizeof(char);
      SKIP;
      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
      nextblock = sizeof(header) + 2 * sizeof(int);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      SKIP;
    }
#ifdef HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      write_header_attributes_in_hdf5(hdf5_headergrp);
    }

  else
    {
#endif
#endif
      blksize = sizeof(header);
      SKIP;
      my_fwrite(&header, sizeof(header), 1, fd);
      SKIP;
#if defined(HAVE_HDF5) &&  defined(WRITE_SUB_IN_SNAP_FORMAT)
    }
#endif


#ifdef WRITE_SUB_IN_SNAP_FORMAT
  if(All.SnapFormat == 2)
    {
      if(!(InfoBlock = (struct info_block *) mymalloc("InfoBlock", bytes = sizeof(struct info_block) * 1000)))
	{
	  printf("failed to allocate memory for `InfoBlock' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}

      int n_info;
      for(bnr = 0, n_info = 0; bnr < 1000; bnr++)
	{
	  blocknr = (enum siofields) bnr;

	  if(blocknr == SIO_LASTENTRY)
	    break;

	  if(block_in_sub(blocknr))
	    {
	      type = get_sub_entrytype_of_block(blocknr);

	      for(i = 0; i < 6; i++)
		InfoBlock[n_info].is_present[i] = 0;
	      InfoBlock[n_info].is_present[type] = 1;

	      InfoBlock[n_info].ndim = get_values_per_sub(blocknr);

	      get_IO_Label_sub(blocknr, label);
	      for(i = 0; i < 4; i++)
		InfoBlock[n_info].label[i] = label[i];

	      datatype = get_datatype_in_sub(blocknr);
	      switch (datatype)
		{
		case 0:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "LONG    ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "LONGN   ", 8);
		  break;
		case 1:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
		  break;
		case 2:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "LLONG   ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "LLONGN  ", 8);
		  break;
		case 3:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
		  break;
		default:
		  endrun(8712519);
		  break;
		}
	      n_info++;
	    }
	}

      blksize = sizeof(int) + 4 * sizeof(char);
      SKIP;
      my_fwrite((void *) "INFO", sizeof(char), 4, fd);
      nextblock = n_info * sizeof(struct info_block) + 2 * sizeof(int);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      SKIP;
      blksize = n_info * sizeof(struct info_block);
      SKIP;
      my_fwrite(InfoBlock, sizeof(struct info_block), n_info, fd);
      SKIP;

      myfree(InfoBlock);
    }
#endif

#else
  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Nsubgroups, sizeof(int), 1, fd);
  my_fwrite(&TotNsubgroups, sizeof(int), 1, fd);
#endif


  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum siofields) bnr;

      if(blocknr == SIO_LASTENTRY)
	break;

      if(block_in_sub(blocknr))
	{
	  type = get_sub_entrytype_of_block(blocknr);
	  switch (type)
	    {
	    case 0:
	      nwrite = Ngroups;
	      break;
	    case 1:
	      nwrite = Nsubgroups;
	      break;
	    case 2:
	      nwrite = Nids;
	      break;
#ifdef SUBFIND_COUNT_BIG_HALOS
	    case 3:
	      nwrite = Nbiggroups;
	      break;
#endif
	    }

	  ndim = get_values_per_sub(blocknr);

#if defined(HAVE_HDF5) &&  defined(WRITE_SUB_IN_SNAP_FORMAT)
	  if(All.SnapFormat == 3)
	    {
	      dims[0] = nwrite;
	      dims[1] = ndim;

	      if(dims[1] == 1)
		rank = 1;
	      else
		rank = 2;

	      switch (get_datatype_in_sub(blocknr))
		{
		case 0:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
		  break;
		case 1:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
		  break;
		case 2:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
		  break;
		case 3:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
		  break;
		default:
		  endrun(8712519);
		  break;
		}
	      get_IO_Label_HDF5_sub(blocknr, buf);
	      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
	      hdf5_dataset =
		H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
	    }
#endif

	  datatype = get_datatype_in_sub(blocknr);
	  switch (datatype)
	    {
	    case 0:
	      bytes_per_blockelement = 4;
	      break;
	    case 1:
	      bytes_per_blockelement = 4;
	      break;
	    case 2:
	      bytes_per_blockelement = 8;
	      break;
	    case 3:
	      bytes_per_blockelement = 8;
	      break;
	    default:
	      endrun(8712519);
	      break;
	    }

#ifdef WRITE_SUB_IN_SNAP_FORMAT
	  get_IO_Label_sub(blocknr, label);

	  if(nwrite > 0)
	    {
	      if(All.SnapFormat == 2)
		{
		  blksize = sizeof(int) + 4 * sizeof(char);
		  SKIP;
		  my_fwrite((void *) label, sizeof(char), 4, fd);
		  nextblock = nwrite * bytes_per_blockelement * ndim + 2 * sizeof(int);
		  my_fwrite(&nextblock, sizeof(int), 1, fd);
		  SKIP;
		}
	      blksize = nwrite * bytes_per_blockelement * ndim;
	      if(All.SnapFormat != 3)
		{
		  SKIP;
		}

	    }
	  else
	    blksize = 0;
#else
	  blksize = nwrite * bytes_per_blockelement * ndim;
#endif
	  if(ThisTask == 0)
	    printf("Writing block %d (%c%c%c%c), n=%d, ptype=%d, dtype=%d, ndim=%d, bpb=%d bytes=%ud\n", bnr,
		   label[0], label[1], label[2], label[3], nwrite, type, datatype, ndim,
		   bytes_per_blockelement, blksize);

	  if(!(IOBuffer = mymalloc("IOBuffer", bytes = blksize)))
	    {
	      printf("failed to allocate memory for `IOBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	      endrun(2);
	    }

	  fp = (MyOutputFloat *) IOBuffer;
	  fp_int = (int *) IOBuffer;
	  fp_id = (MyIDType *) IOBuffer;

	  switch (blocknr)
	    {
	    case SIO_GLEN:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].Len;
	      break;
	    case SIO_GOFF:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].Offset;
	      break;
	    case SIO_MTOT:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].Mass;
	      break;
	    case SIO_GPOS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) Group[i].Pos[j];
	      break;
	    case SIO_MMEA:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].M_Mean200;
	      break;
	    case SIO_RMEA:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].R_Mean200;
	      break;
	    case SIO_MCRI:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].M_Crit200;
	      break;
	    case SIO_RCRI:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].R_Crit200;
	      break;
	    case SIO_MTOP:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].M_TopHat200;
	      break;
	    case SIO_RTOP:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].R_TopHat200;
	      break;
#ifdef SO_VEL_DISPERSIONS
	    case SIO_DMEA:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].VelDisp_Mean200;
	      break;
	    case SIO_DCRI:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].VelDisp_Crit200;
	      break;
	    case SIO_DTOP:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].VelDisp_TopHat200;
	      break;
#endif
#ifdef SO_BAR_INFO
	    case SIO_MGAS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) Group[i].gas_mass[j];
	      break;
	    case SIO_MSTR:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) Group[i].star_mass[j];
	      break;
	    case SIO_TGAS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) Group[i].temp[j];
	      break;
	    case SIO_LGAS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) Group[i].xlum[j];
	      break;
#endif
	    case SIO_NCON:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].ContaminationLen;
	      break;
	    case SIO_MCON:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].ContaminationMass;
	      break;
#ifdef SUBFIND_COUNT_BIG_HALOS
	    case SIO_BGPOS:
	      for(i = 0, k = 0; i < Ngroups; i++)
		if(Group[i].M_TopHat200 > SUBFIND_COUNT_BIG_HALOS)
		  {
		    for(j = 0; j < 3; j++)
		      fp[k * 3 + j] = (MyOutputFloat) Group[i].Pos[j];
		    k++;
		  }
	      if(k != nwrite)
		endrun(87453217);
	      break;
	    case SIO_BGMTOP:
	      for(i = 0, k = 0; i < Ngroups; i++)
		if(Group[i].M_TopHat200 > SUBFIND_COUNT_BIG_HALOS)
		  {
		    fp[k] = (MyOutputFloat) Group[i].M_TopHat200;
		    k++;
		  }
	      if(k != nwrite)
		endrun(87453218);
	      break;
	    case SIO_BGRTOP:
	      for(i = 0, k = 0; i < Ngroups; i++)
		if(Group[i].M_TopHat200 > SUBFIND_COUNT_BIG_HALOS)
		  {
		    fp[k] = (MyOutputFloat) Group[i].R_TopHat200;
		    k++;
		  }
	      if(k != nwrite)
		endrun(87453219);
	      break;
#endif
	    case SIO_NSUB:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].Nsubs;
	      break;
	    case SIO_FSUB:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].FirstSub;
	      break;
	    case SIO_SLEN:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = SubGroup[i].Len;
	      break;
	    case SIO_SOFF:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = SubGroup[i].Offset;
	      break;
	    case SIO_PFOF:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = SubGroup[i].SubParent;
	      break;
	    case SIO_MSUB:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].Mass;
	      break;
	    case SIO_SPOS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].Pos[j];
	      break;
	    case SIO_SVEL:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].Vel[j];
	      break;
	    case SIO_SCM:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].CM[j];
	      break;
	    case SIO_SPIN:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].Spin[j];
	      break;
	    case SIO_DSUB:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubVelDisp;
	      break;
	    case SIO_VMAX:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubVmax;
	      break;
	    case SIO_RVMAX:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubVmaxRad;
	      break;
	    case SIO_RHMS:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubHalfMass;
	      break;
	    case SIO_MBID:
	      for(i = 0; i < nwrite; i++)
		fp_id[i] = SubGroup[i].SubMostBoundID;
	      break;
	    case SIO_GRNR:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = SubGroup[i].GrNr;
	      break;
#ifdef SAVE_MASS_TAB
	    case SIO_SMST:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 6; j++)
		  fp[i * 6 + j] = (MyOutputFloat) SubGroup[i].MassTab[j];
	      break;
#endif
#ifdef LT_ADD_GAL_TO_SUB
	    case SIO_SLUM:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < LT_ADD_GAL_TO_SUB; j++)
		  fp[i * LT_ADD_GAL_TO_SUB + j] = (MyOutputFloat) SubGroup[i].SubLum[j];
	      break;
#ifdef DUSTATT
	    case SIO_SLATT:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < LT_ADD_GAL_TO_SUB; j++)
		  fp[i * LT_ADD_GAL_TO_SUB + j] = (MyOutputFloat) SubGroup[i].SubLumAtt[j];
	      break;
	    case SIO_DUST:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < DUSTATT; j++)
		  fp[i * DUSTATT + j] = (MyOutputFloat) SubGroup[i].DustProfile[j];
	      break;
#endif
#ifdef OBSERVER_FRAME
	    case SIO_SLOBS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < LT_ADD_GAL_TO_SUB; j++)
		  fp[i * LT_ADD_GAL_TO_SUB + j] = (MyOutputFloat) SubGroup[i].SubLumObs[j];
	      break;
#endif
	    case SIO_SAGE:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubMeanAge;
	      break;
	    case SIO_SZ:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubMeanZStar;
	      break;
	    case SIO_SSFR:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubSFR;
	      break;
#endif
#ifdef SUBFIND_SAVE_PARTICLELISTS
	    case SIO_PPOS:
#ifndef WRITE_SUB_IN_SNAP_FORMAT	/* open new file in case of old format */
	      fclose(fd);
	      sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_posvel", num, ThisTask);
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("can't open file `%s`\n", buf);
		  endrun(1184);
		}
	      my_fwrite(&Ngroups, sizeof(int), 1, fd);
	      my_fwrite(&TotNgroups, sizeof(int), 1, fd);
	      my_fwrite(&Nids, sizeof(int), 1, fd);
	      my_fwrite(&TotNids, sizeof(long long), 1, fd);
	      my_fwrite(&NTask, sizeof(int), 1, fd);
	      my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
	      my_fwrite(&All.Time, sizeof(double), 1, fd);
#endif
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) ID_list[i].Pos[j];
	      break;
	    case SIO_PVEL:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) (ID_list[i].Vel[j] * sqrt(a3inv));
	      break;
	    case SIO_PTYP:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = ID_list[i].Type;
	      break;
	    case SIO_PMAS:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) ID_list[i].Mass;
	      break;
#ifdef STELLARAGE
	    case SIO_PAGE:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) ID_list[i].StellarAge;
	      break;
#endif
#endif
	    case SIO_PID:
#ifndef WRITE_SUB_IN_SNAP_FORMAT	/* open new file in case of old format */
	      fclose(fd);
	      sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, ThisTask);
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("can't open file `%s`\n", buf);
		  endrun(1184);
		}
	      my_fwrite(&Ngroups, sizeof(int), 1, fd);
	      my_fwrite(&TotNgroups, sizeof(int), 1, fd);
	      my_fwrite(&Nids, sizeof(int), 1, fd);
	      my_fwrite(&TotNids, sizeof(long long), 1, fd);
	      my_fwrite(&NTask, sizeof(int), 1, fd);
	      my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
#endif
	      for(i = 0; i < nwrite; i++)
		fp_id[i] = ID_list[i].ID;
	      break;
	    default:
	      break;

	    }			/* closing switch(blocknr) */

	  if(blksize > 0)
	    {
#if defined(HAVE_HDF5) && defined(WRITE_SUB_IN_SNAP_FORMAT)
	      if(All.SnapFormat != 3)
		{
#endif
		  my_fwrite(IOBuffer, blksize, sizeof(char), fd);
#ifdef WRITE_SUB_IN_SNAP_FORMAT
		  SKIP;
#endif
#if defined(HAVE_HDF5) && defined(WRITE_SUB_IN_SNAP_FORMAT)
		}
	      else if(All.SnapFormat == 3)
		{
		  start[0] = 0;
		  start[1] = 0;

		  count[0] = nwrite;
		  count[1] = ndim;

		  H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

		  dims[0] = nwrite;
		  dims[1] = ndim;
		  hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

		  hdf5_status =
		    H5Dwrite(hdf5_dataset, hdf5_datatype,
			     hdf5_dataspace_memory, hdf5_dataspace_in_file, H5P_DEFAULT, IOBuffer);

		  H5Sclose(hdf5_dataspace_memory);
		}
#endif

	    }
	  myfree(IOBuffer);

	}			/* closing if block exist */

    }				/* closing for bnr loop */
#if defined(HAVE_HDF5) && defined(WRITE_SUB_IN_SNAP_FORMAT)
  if(All.SnapFormat == 3)
    {

      for(type = 0; type < 3; type++)
	H5Gclose(hdf5_grp[type]);
      H5Gclose(hdf5_headergrp);
      H5Fclose(hdf5_file);
    }
  else
    {
#endif
      fclose(fd);
#if defined(HAVE_HDF5) && defined(WRITE_SUB_IN_SNAP_FORMAT)
    }
#endif
}


#ifdef KD_ALTERNATIVE_GROUP_SORT

int GroupToTask(int i, int totngroups)
{
  int nGroupsPerTask, ret;
  double x1, x2, w1, w2, w;

  nGroupsPerTask = totngroups / NTask + 1;

  if(totngroups > 0)
    {
      x1 = (double) i;
      if(totngroups > 1)
	x2 = sqrt(x1 / (totngroups - 1)) * (totngroups - 1);
      else
	x2 = 0;

      w1 = 5 + ((double) NTask) / totngroups;
      w2 = ((double) totngroups) / NTask;
      w = w1 + w2;
      w1 /= w;
      w2 /= w;
      ret = (w1 * x1 + w2 * x2) / nGroupsPerTask;
    }
  else
    {
      ret = 0;
    }

  if(ret >= NTask)
    ret = NTask - 1;

  return ret;
}

void sort_groups()
{
  int i, nexport = 0, nimport = 0, target, ngrp, recvTask;
  group_properties *send_Group;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = GroupToTask(Group[i].GrNr, TotNgroups);
      if(target != ThisTask)
	Send_count[target]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }

  send_Group = (group_properties *) mymalloc("send_Group", nexport * sizeof(group_properties));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = GroupToTask(Group[i].GrNr, TotNgroups);
      if(target != ThisTask)
	{
	  send_Group[Send_offset[target] + Send_count[target]] = Group[i];
	  Send_count[target]++;

	  Group[i] = Group[Ngroups - 1];
	  Ngroups--;
	  i--;
	}
    }

  if(Ngroups + nimport >= MaxNgroups)
    {
      printf("Task %d: NgroupThis = %d, NgroupMax = %d, Ntotgroups=%d \n",
	     ThisTask, Ngroups + nimport, MaxNgroups, TotNgroups);
      endrun(882233);
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&send_Group[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(group_properties), MPI_BYTE,
			   recvTask, TAG_ALTSORT_GRP,
			   &Group[Ngroups + Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(group_properties), MPI_BYTE,
			   recvTask, TAG_ALTSORT_GRP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  Ngroups += nimport;

  myfree(send_Group);
}

void sort_sub_groups()
{
  int i, nexport = 0, nimport = 0, target, ngrp, recvTask;
  struct subgroup_properties *send_SubGroup;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Nsubgroups; i++)
    {
      target = GroupToTask(SubGroup[i].GrNr, TotNgroups);
      if(target != ThisTask)
	Send_count[target]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }

  send_SubGroup =
    (struct subgroup_properties *) mymalloc("send_SubGroup", nexport * sizeof(struct subgroup_properties));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Nsubgroups; i++)
    {
      target = GroupToTask(SubGroup[i].GrNr, TotNgroups);
      if(target != ThisTask)
	{
	  send_SubGroup[Send_offset[target] + Send_count[target]] = SubGroup[i];
	  Send_count[target]++;

	  SubGroup[i] = SubGroup[Nsubgroups - 1];
	  Nsubgroups--;
	  i--;
	}
    }

  if(Nsubgroups + nimport >= MaxNsubgroups)
    {
      printf("Task %d: NsubgroupThis = %d, NsubgroupMax = %d\n", ThisTask, Nsubgroups + nimport,
	     MaxNsubgroups);
      endrun(882234);
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&send_SubGroup[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct subgroup_properties), MPI_BYTE,
			   recvTask, TAG_ALTSORT_SUB,
			   &SubGroup[Nsubgroups + Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct subgroup_properties), MPI_BYTE,
			   recvTask, TAG_ALTSORT_SUB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  Nsubgroups += nimport;

  myfree(send_SubGroup);
}

void sort_ids()
{
  int i, nexport = 0, nimport = 0, target, ngrp, recvTask;
  struct subfind_id_list *send_IDs;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Nids; i++)
    {
      target = GroupToTask(ID_list[i].GrNr, TotNgroups);
      if(target != ThisTask)
	Send_count[target]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }

  send_IDs = (struct subfind_id_list *) mymalloc("send_IDs", nexport * sizeof(struct subfind_id_list));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Nids; i++)
    {
      target = GroupToTask(ID_list[i].GrNr, TotNgroups);
      if(target != ThisTask)
	{
	  send_IDs[Send_offset[target] + Send_count[target]] = ID_list[i];
	  Send_count[target]++;

	  ID_list[i] = ID_list[Nids - 1];
	  Nids--;
	  i--;
	}
    }

  if(Nids + nimport >= All.MaxPart * 1.5)
    {
      printf("Task %d: NidsThis = %d, NidsMax = %d\n", ThisTask, Nids + nimport, (int) (All.MaxPart * 1.5));
      endrun(882235);
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&send_IDs[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct subfind_id_list), MPI_BYTE,
			   recvTask, TAG_ALTSORT_IDS,
			   &ID_list[Nids + Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct subfind_id_list), MPI_BYTE,
			   recvTask, TAG_ALTSORT_IDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  Nids += nimport;

  myfree(send_IDs);
}

#endif


int subfind_compare_ID_list(const void *a, const void *b)
{
  if(((struct subfind_id_list *) a)->GrNr < ((struct subfind_id_list *) b)->GrNr)
    return -1;

  if(((struct subfind_id_list *) a)->GrNr > ((struct subfind_id_list *) b)->GrNr)
    return +1;

  if(((struct subfind_id_list *) a)->SubNr < ((struct subfind_id_list *) b)->SubNr)
    return -1;

  if(((struct subfind_id_list *) a)->SubNr > ((struct subfind_id_list *) b)->SubNr)
    return +1;

  if(((struct subfind_id_list *) a)->BindingEgy < ((struct subfind_id_list *) b)->BindingEgy)
    return -1;

  if(((struct subfind_id_list *) a)->BindingEgy > ((struct subfind_id_list *) b)->BindingEgy)
    return +1;

  return 0;
}

int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *) a)->GrNr < ((struct subgroup_properties *) b)->GrNr)
    return -1;

  if(((struct subgroup_properties *) a)->GrNr > ((struct subgroup_properties *) b)->GrNr)
    return +1;

  if(((struct subgroup_properties *) a)->SubNr < ((struct subgroup_properties *) b)->SubNr)
    return -1;

  if(((struct subgroup_properties *) a)->SubNr > ((struct subgroup_properties *) b)->SubNr)
    return +1;

  return 0;
}


int subfind_compare_P_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
    return -1;

  if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
    return +1;

  if(((struct particle_data *) a)->u.DM_Density > (((struct particle_data *) b)->u.DM_Density))
    return -1;

  if(((struct particle_data *) a)->u.DM_Density < (((struct particle_data *) b)->u.DM_Density))
    return +1;

  return 0;
}

#endif

#if defined(SUBFIND)		// || defined(WRITE_KEY_FILES)

int subfind_compare_P_origindex(const void *a, const void *b)
{
  if(((struct particle_data *) a)->origindex < (((struct particle_data *) b)->origindex))
    return -1;

  if(((struct particle_data *) a)->origindex > (((struct particle_data *) b)->origindex))
    return +1;

  return 0;
}


#endif
