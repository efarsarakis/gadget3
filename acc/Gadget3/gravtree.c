#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/sem.h>

#include "allvars.h"
#include "proto.h"


#ifdef NUM_THREADS
#include <pthread.h>
#endif

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

#ifdef NUM_THREADS
pthread_mutex_t mutex_nexport;
pthread_mutex_t mutex_workcount;
pthread_mutex_t mutex_partnodedrift;

#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#define LOCK_WORKCOUNT   pthread_mutex_lock(&mutex_workcount);
#define UNLOCK_WORKCOUNT pthread_mutex_unlock(&mutex_workcount);

#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#define LOCK_WORKCOUNT
#define UNLOCK_WORKCOUNT
#endif

//manos//
#define LOCK_PARTNODEDRIFT
#define UNLOCK_PARTNODEDRIFT

#ifdef STATICBRANDT
inline double OmegaR(int i, int mode)
{
	double r, m;
	r =
			sqrt((P[i].Pos[0] - LONG_X / 2.) * (P[i].Pos[0] - LONG_X / 2.) +
					(P[i].Pos[1] - LONG_Y / 2.) * (P[i].Pos[1] - LONG_Y / 2.));
	if(r > 0)
	{
		switch (mode)
		{
		case -1:
			m = 1;
			break;
		case 0:
			if(r < 7)
				m = 1.0;
			else
				m = abs(sin((r - 7) * M_PI / (2 * All.Alfa2H)) + 1);
			if(r > 8)
				m = 1.0;
			break;
		case 1:
			m = (BRANDT_OmegaBr) / sqrt(1 + (r / BRANDT_R0) * (r / BRANDT_R0));
			break;
		case 2:
			if(r < BRANDT_R0)
				m = BRANDT_OmegaBr;
			else
				m = BRANDT_OmegaBr * BRANDT_R0 / r;
			break;
		default:
			printf("Wrong Immplementation of Omega\n");
			exit(80085);
			break;
		};
	}
	else
		m = 0;
	return m;
};
#endif
double Ewaldcount, Costtotal;
long long N_nodesinlist;


int Ewald_iter;			/* global in file scope, for simplicity */

extern float shortrange_table[1000];

void sum_top_level_node_costfactors(void);



/*! This function computes the gravitational forces for all active particles.
 *  If needed, a new tree is constructed, otherwise the dynamically updated
 *  tree is used.  Particles are only exported to other processors when really
 *  needed, thereby allowing a good use of the communication buffer.
 */
void gravity_tree(void)
{
	long long n_exported = 0;
	int i, j, maxnumnodes, iter = 0;
	double t0, t1;
	double timeall = 0, timetree1 = 0, timetree2 = 0;
	double timetree, timewait, timecomm;
	double timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
	double sum_costtotal, ewaldtot;
	double maxt, sumt, maxt1, sumt1, maxt2, sumt2, sumcommall, sumwaitall;
	double plb, plb_max;

#ifdef FIXEDTIMEINFIRSTPHASE
	int counter;
	double min_time_first_phase, min_time_first_phase_glob;
#endif
#ifndef NOGRAVITY
	int k, ewald_max, diff, save_NextParticle;
	int ndone, ndone_flag, ngrp;
	int place;
	int recvTask;
	double tstart, tend, ax, ay, az;
	MPI_Status status;

#ifdef DISTORTIONTENSORPS
	int i1, i2;
#endif
#endif

#ifndef GRAVITY_CENTROID
	CPU_Step[CPU_MISC] += measure_time();

	/* set new softening lengths */
	if(All.ComovingIntegrationOn)
		set_softenings();

#ifndef MODGRAV
	/* construct tree if needed */
	if(TreeReconstructFlag)
	{
		if(ThisTask == 0)
			printf("Tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

		CPU_Step[CPU_MISC] += measure_time();

		force_treebuild(NumPart, NULL);

		CPU_Step[CPU_TREEBUILD] += measure_time();

		TreeReconstructFlag = 0;

		if(ThisTask == 0)
			printf("Tree construction done.\n");
	}
#endif
#endif

#ifndef NOGRAVITY

#ifdef ADAPTGRAVSOFT
	ags_density();
	ags_force_update_hmax();
#endif

	/* allocate buffers to arrange communication */
	if(ThisTask == 0)
		printf("Begin tree force.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifdef KD_BUFFER_MANAGEMENT
	size_t MyBufferSize = 0.7 * FreeBytes / (1024.0 * 1024.0);
#else
	size_t MyBufferSize = All.BufferSize;
#endif

	All.BunchSize =
			(int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					sizeof(struct gravdata_in) + sizeof(struct gravdata_out) +
					sizemax(sizeof(struct gravdata_in), sizeof(struct gravdata_out))));
	DataIndexTable =
			(struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
	DataNodeList =
			(struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

#ifdef KD_BUFFER_MANAGEMENT
	if(ThisTask == 0)
		printf("Gravity: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
				FreeBytes / (1024.0 * 1024.0));
#endif

	if(ThisTask == 0)
		printf("All.BunchSize=%d\n", All.BunchSize);

	Ewaldcount = 0;
	Costtotal = 0;
	N_nodesinlist = 0;


	CPU_Step[CPU_TREEMISC] += measure_time();
	t0 = second();

#if defined(PERIODIC) && !defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC)
	ewald_max = 1;
#else
	ewald_max = 0;
#endif

#ifdef SCF_HYBRID
	int scf_counter, max_scf_counter = 1;

	if(SCF_HYBRID == 2)
		max_scf_counter = 0;
	/*
     calculates the following forces (depending on SCF_HYBRID value)
     STAR<->STAR, STAR->DM (scf_counter=0)
     DM<->DM (scf_counter=1)
	 */

	for(scf_counter = 0; scf_counter <= max_scf_counter; scf_counter++)
	{
		/* set DM mass to zero and set gravsum to zero */
		if(scf_counter == 0)
		{
			for(i = 0; i < NumPart; i++)
			{
				if(P[i].Type == 1)	/* DM particle */
					P[i].Mass = 0.0;

				for(j = 0; j < 3; j++)
					P[i].GravAccelSum[j] = 0.0;
			}
		}
		/* set stellar mass to zero */
		if(scf_counter == 1)
		{
			for(i = 0; i < NumPart; i++)
			{
				if(P[i].Type == 2)	/* stellar particle */
					P[i].Mass = 0.0;
			}
		}

		/* particle masses changed, so reconstruct tree */
		if(ThisTask == 0)
			printf("SCF Tree construction %d\n", scf_counter);
		force_treebuild(NumPart, NULL);
		if(ThisTask == 0)
			printf("done.\n");
#endif

		if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)
		{
			/* we have a fresh tree and would like to measure gravity cost */

			/* find the closest level */
			for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - All.HighestActiveTimeBin);
					i < GRAVCOSTLEVELS; i++)
			{
				if(diff > abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin))
				{
					TakeLevel = i;
					diff = abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin);
				}
			}

			if(diff != 0)		/* we have not found a matching slot */
			{

				if(All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
				{
					/* clear levels that are out of range */
					for(i = 0; i < GRAVCOSTLEVELS; i++)
					{
						if(All.LevelToTimeBin[i] > All.HighestOccupiedTimeBin)
							All.LevelToTimeBin[i] = 0;
						if(All.LevelToTimeBin[i] < All.HighestOccupiedTimeBin - (GRAVCOSTLEVELS - 1))
							All.LevelToTimeBin[i] = 0;
					}
				}

				for(i = 0, TakeLevel = -1; i < GRAVCOSTLEVELS; i++)
				{
					if(All.LevelToTimeBin[i] == 0)
					{
						All.LevelToTimeBin[i] = All.HighestActiveTimeBin;
						TakeLevel = i;
						break;
					}
				}

				if(TakeLevel < 0 && All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
					terminate("TakeLevel < 0, even though we should have a slot");
			}
		}
		else
		{
			/* in this case we do not measure gravity cost. Check whether this time-level
	     has previously mean measured. If yes, then delete it so to make sure that it is not out of time */

			for(i = 0; i < GRAVCOSTLEVELS; i++)
				if(All.LevelToTimeBin[i] == All.HighestActiveTimeBin)
					All.LevelToTimeBin[i] = 0;

			TakeLevel = -1;
		}


		if(TakeLevel >= 0)
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(i = 0; i < NumPart; i++)
				P[i].GravCost[TakeLevel] = 0;


		for(Ewald_iter = 0; Ewald_iter <= ewald_max; Ewald_iter++)
		{

			NextParticle = FirstActiveParticle;	/* beginn with this index */

			do
			{
				iter++;
				BufferFullFlag = 0;
				Nexport = 0;
				save_NextParticle = NextParticle;

				tstart = second();

#ifdef NUM_THREADS
				pthread_t mythreads[NUM_THREADS - 1];
				int threadid[NUM_THREADS - 1];
				pthread_attr_t attr;

				pthread_attr_init(&attr);
				pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				pthread_mutex_init(&mutex_workcount, NULL);
				pthread_mutex_init(&mutex_nexport, NULL);
				pthread_mutex_init(&mutex_partnodedrift, NULL);

				TimerFlag = 0;

				for(j = 0; j < NUM_THREADS - 1; j++)
				{
					threadid[j] = j + 1;
					pthread_create(&mythreads[j], &attr, gravity_primary_loop, &threadid[j]);
				}
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
				{
#ifdef _OPENMP
					int mainthreadid = omp_get_thread_num();
#else
					int mainthreadid = 0;
#endif
					gravity_primary_loop(&mainthreadid);	/* do local particles and prepare export list */
				}

#ifdef NUM_THREADS
				for(j = 0; j < NUM_THREADS - 1; j++)
					pthread_join(mythreads[j], NULL);
#endif

				tend = second();
				timetree1 += timediff(tstart, tend);


				if(BufferFullFlag)
				{
					int last_nextparticle = NextParticle;

					NextParticle = save_NextParticle;

					while(NextParticle >= 0)
					{
						if(NextParticle == last_nextparticle)
							break;

						if(ProcessedFlag[NextParticle] != 1)
							break;

						ProcessedFlag[NextParticle] = 2;

						NextParticle = NextActiveParticle[NextParticle];
					}

					if(NextParticle == save_NextParticle)
					{
						/* in this case, the buffer is too small to process even a single particle */
						endrun(12998);
					}


					int new_export = 0;

					for(j = 0, k = 0; j < Nexport; j++)
						if(ProcessedFlag[DataIndexTable[j].Index] != 2)
						{
							if(k < j + 1)
								k = j + 1;

							for(; k < Nexport; k++)
								if(ProcessedFlag[DataIndexTable[k].Index] == 2)
								{
									int old_index = DataIndexTable[j].Index;

									DataIndexTable[j] = DataIndexTable[k];
									DataNodeList[j] = DataNodeList[k];
									DataIndexTable[j].IndexGet = j;
									new_export++;

									DataIndexTable[k].Index = old_index;
									k++;
									break;
								}
						}
						else
							new_export++;

					Nexport = new_export;

				}


				n_exported += Nexport;

				for(j = 0; j < NTask; j++)
					Send_count[j] = 0;
				for(j = 0; j < Nexport; j++)
					Send_count[DataIndexTable[j].Task]++;


				MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

				tstart = second();

				MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

				tend = second();
				timewait1 += timediff(tstart, tend);


				for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
				{
					Nimport += Recv_count[j];

					if(j > 0)
					{
						Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
						Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
					}
				}

				GravDataGet =
						(struct gravdata_in *) mymalloc("GravDataGet", Nimport * sizeof(struct gravdata_in));
				GravDataIn =
						(struct gravdata_in *) mymalloc("GravDataIn", Nexport * sizeof(struct gravdata_in));

				/* prepare particle data for export */

				for(j = 0; j < Nexport; j++)
				{
					place = DataIndexTable[j].Index;

#ifdef GRAVITY_CENTROID
					if(P[place].Type == 0)
					{
						for(k = 0; k < 3; k++)
							GravDataIn[j].Pos[k] = SphP[place].Center[k];
					}
					else
					{
						for(k = 0; k < 3; k++)
							GravDataIn[j].Pos[k] = P[place].Pos[k];
					}
#else
					for(k = 0; k < 3; k++)
						GravDataIn[j].Pos[k] = P[place].Pos[k];
#endif
#ifdef FS_TURB_ESTIM
					for(k = 0; k < 3; k++)
						GravDataIn[j].Vel[k] = P[place].Vel[k];
					if(P[place].Type == 0)
						GravDataIn[j].Vel[3] = PPP[place].Hsml;
					else
						GravDataIn[j].Vel[3] = -1;
#endif
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
					GravDataIn[j].Type = P[place].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
					if(P[place].Type == 0)
						GravDataIn[j].Soft = PPP[place].Hsml;
#endif
#endif

					GravDataIn[j].OldAcc = P[place].OldAcc;

#ifdef ADAPTGRAVSOFT
					GravDataIn[j].AGS_zeta = P[place].AGS_zeta;
					GravDataIn[j].AGS_omega = P[place].AGS_omega;
					GravDataIn[j].AGS_Hsml = P[place].AGS_Hsml;
					GravDataIn[j].Mass = P[place].Mass;
#endif

					memcpy(GravDataIn[j].NodeList,
							DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
				}


				/* exchange particle data */

				tstart = second();
				for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
				{
					recvTask = ThisTask ^ ngrp;

					if(recvTask < NTask)
					{
						if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
						{
							/* get the particles */
							MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
									Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
									recvTask, TAG_GRAV_A,
									&GravDataGet[Recv_offset[recvTask]],
									Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
									recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
						}
					}
				}
				tend = second();
				timecommsumm1 += timediff(tstart, tend);


				myfree(GravDataIn);
				GravDataResult =
						(struct gravdata_out *) mymalloc("GravDataResult", Nimport * sizeof(struct gravdata_out));
				GravDataOut =
						(struct gravdata_out *) mymalloc("GravDataOut", Nexport * sizeof(struct gravdata_out));

				report_memory_usage(&HighMark_gravtree, "GRAVTREE");

				/* now do the particles that were sent to us */
				tstart = second();

				NextJ = 0;

#ifdef NUM_THREADS
				for(j = 0; j < NUM_THREADS - 1; j++)
					pthread_create(&mythreads[j], &attr, gravity_secondary_loop, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
				{
#ifdef _OPENMP
					int mainthreadid = omp_get_thread_num();
#else
					int mainthreadid = 0;
#endif
					gravity_secondary_loop(&mainthreadid);
				}

#ifdef NUM_THREADS
				for(j = 0; j < NUM_THREADS - 1; j++)
					pthread_join(mythreads[j], NULL);

				pthread_mutex_destroy(&mutex_partnodedrift);
				pthread_mutex_destroy(&mutex_nexport);
				pthread_mutex_destroy(&mutex_workcount);
				pthread_attr_destroy(&attr);
#endif

				tend = second();
				timetree2 += timediff(tstart, tend);

				if(NextParticle < 0)
					ndone_flag = 1;
				else
					ndone_flag = 0;

				tstart = second();
				MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				tend = second();
				timewait2 += timediff(tstart, tend);

				/* get the result */
				tstart = second();
				for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
				{
					recvTask = ThisTask ^ ngrp;
					if(recvTask < NTask)
					{
						if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
						{
							/* send the results */
							MPI_Sendrecv(&GravDataResult[Recv_offset[recvTask]],
									Recv_count[recvTask] * sizeof(struct gravdata_out),
									MPI_BYTE, recvTask, TAG_GRAV_B,
									&GravDataOut[Send_offset[recvTask]],
									Send_count[recvTask] * sizeof(struct gravdata_out),
									MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);
						}
					}

				}
				tend = second();
				timecommsumm2 += timediff(tstart, tend);

				/* add the results to the local particles */
				tstart = second();
				for(j = 0; j < Nexport; j++)
				{
					place = DataIndexTable[j].Index;

					for(k = 0; k < 3; k++)
						P[place].g.dGravAccel[k] += GravDataOut[j].Acc[k];
#ifdef FS_TURB_ESTIM
					for(k = 0; k < FS_BINS; k++)
					{
						P[place].StrFnc[k] += GravDataOut[j].StrFnc[k];
						P[place].StrFnc_count[k] += GravDataOut[j].StrFnc_count[k];
					}
#endif

#ifdef DISTORTIONTENSORPS
					for(i1 = 0; i1 < 3; i1++)
						for(i2 = 0; i2 < 3; i2++)
							P[place].tidal_tensorps[i1][i2] += GravDataOut[j].tidal_tensorps[i1][i2];
#endif

#ifdef EVALPOTENTIAL
					P[place].p.dPotential += GravDataOut[j].Potential;
#endif

#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
					P[place].AGS_corr += GravDataOut[j].AGS_corr;
#endif
				}
				tend = second();
				timetree1 += timediff(tstart, tend);

				myfree(GravDataOut);
				myfree(GravDataResult);
				myfree(GravDataGet);
			}
			while(ndone < NTask);
		}			/* Ewald_iter */

#ifdef SCF_HYBRID
		/* restore particle masses */
		for(i = 0; i < NumPart; i++)
			P[i].Mass = P[i].MassBackup;


		/* add up accelerations from tree to AccelSum */
		for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
		{
			/* ignore STAR<-DM contribution */
			if(scf_counter == 1 && P[i].Type == 2)
			{
				continue;
			}
			else
			{
				for(j = 0; j < 3; j++)
					P[i].GravAccelSum[j] += P[i].g.dGravAccel[j];
			}
		}
	}				/* scf_counter */

	/* set acceleration to summed up accelerations */
	for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
		for(j = 0; j < 3; j++)
			P[i].g.dGravAccel[j] = P[i].GravAccelSum[j];
	}
#endif



#ifdef FLTROUNDOFFREDUCTION
	for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
#ifdef EVALPOTENTIAL
		P[i].p.Potential = FLT(P[i].p.dPotential);
#endif
		for(j = 0; j < 3; j++)
			P[i].g.GravAccel[j] = FLT(P[i].g.dGravAccel[j]);
	}
#endif


	myfree(DataNodeList);
	myfree(DataIndexTable);

	/* assign node cost to particles */
	if(TakeLevel >= 0)
	{
		sum_top_level_node_costfactors();

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i = 0; i < NumPart; i++)
		{
#ifdef NEUTRINOS
			if(P[i].Type != 2)
#endif
			{
				int no = Father[i];

				while(no >= 0)
				{
					if(Nodes[no].u.d.mass > 0)
						P[i].GravCost[TakeLevel] += Nodes[no].GravCost * P[i].Mass / Nodes[no].u.d.mass;

					no = Nodes[no].u.d.father;
				}
			}
		}
	}

	/* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
	if(All.ComovingIntegrationOn)
	{
		double fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

		for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
		{
			for(j = 0; j < 3; j++)
				P[i].g.GravAccel[j] += fac * P[i].Pos[j];
		}
	}
#endif
#endif

	/* to prevent that we overwrite OldAcc in the first evaluation for 2lpt ICs */
	if(!(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0))
	{
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
#ifdef PMGRID
				ax = P[i].g.GravAccel[0] + P[i].GravPM[0] / All.G;
				ay = P[i].g.GravAccel[1] + P[i].GravPM[1] / All.G;
				az = P[i].g.GravAccel[2] + P[i].GravPM[2] / All.G;
#else
				ax = P[i].g.GravAccel[0];
				ay = P[i].g.GravAccel[1];
				az = P[i].g.GravAccel[2];
#endif
				P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
			}
		}

		if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
		{
			if(!(All.Ti_Current == 0 && RestartFlag == 0))
				if(All.TypeOfOpeningCriterion == 1)
					All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */
		}
		else
		{
			if(All.TypeOfOpeningCriterion == 1)
				All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */
		}

		/*  muliply by G */
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
				for(j = 0; j < 3; j++)
					P[i].g.GravAccel[j] *= All.G;

#ifdef FS_TURB_ESTIM
				if(P[i].Type == 0)
				{
					SphP[i].Vturb = -1;
					for(k = 0; k < FS_BINS; k++)
					{
						if(P[i].StrFnc_count[k] > 0)
							P[i].StrFnc[k] = sqrt(P[i].StrFnc[k] / P[i].StrFnc_count[k]);
						else
							P[i].StrFnc[k] = 0;
						if(k > 1)
						{
							SphP[i].Lturb = (P[i].StrFnc[k] - P[i].StrFnc[k - 1]) / PPP[i].Hsml;
							if(SphP[i].Lturb < 0)
							{
								SphP[i].Vturb = (P[i].StrFnc[k - 1]);
								SphP[i].Lturb = PPP[i].Hsml * (k - 0.5);
								k = FS_BINS;
							}
						}
					}
					if(SphP[i].Vturb == -1)
					{
						SphP[i].Vturb = P[i].StrFnc[0];
						SphP[i].Lturb = PPP[i].Hsml;
					}
				}
#endif


#ifdef DISTORTIONTENSORPS
				/*
         Diaganol terms of tidal tensor need correction, because tree is running over
         all particles -> also over target particle -> extra term -> correct it
				 */
				/* 3D -> full forces */
				P[i].tidal_tensorps[0][0] +=
						P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] *
								All.ForceSoftening[P[i].Type]) * 10.666666666667;

				P[i].tidal_tensorps[1][1] +=
						P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] *
								All.ForceSoftening[P[i].Type]) * 10.666666666667;

				P[i].tidal_tensorps[2][2] +=
						P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] *
								All.ForceSoftening[P[i].Type]) * 10.666666666667;

				if(All.ComovingIntegrationOn)
				{
					P[i].tidal_tensorps[0][0] -= All.TidalCorrection / All.G;
					P[i].tidal_tensorps[1][1] -= All.TidalCorrection / All.G;
					P[i].tidal_tensorps[2][2] -= All.TidalCorrection / All.G;
				}
				/*now muliply by All.G */
				for(i1 = 0; i1 < 3; i1++)
					for(i2 = 0; i2 < 3; i2++)
						P[i].tidal_tensorps[i1][i2] *= All.G;
#endif /* DISTORTIONTENSORPS */

#ifdef EVALPOTENTIAL
				/* remove self-potential */
				P[i].p.Potential += P[i].Mass / All.SofteningTable[P[i].Type];

				if(All.ComovingIntegrationOn)
					if(All.PeriodicBoundariesOn)
						P[i].p.Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
						pow(All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G), 1.0 / 3);

				P[i].p.Potential *= All.G;

#ifdef PMGRID
				P[i].p.Potential += P[i].PM_Potential;	/* add in long-range potential */
#endif

				if(All.ComovingIntegrationOn)
				{
#ifndef PERIODIC
					double fac, r2;

					fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

					for(k = 0, r2 = 0; k < 3; k++)
						r2 += P[i].Pos[k] * P[i].Pos[k];

					P[i].p.Potential += fac * r2;
#endif
				}
				else
				{
					double fac, r2;

					fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;

					if(fac != 0)
					{
						for(k = 0, r2 = 0; k < 3; k++)
							r2 += P[i].Pos[k] * P[i].Pos[k];

						P[i].p.Potential += fac * r2;
					}
				}
#endif
			}

			/* Finally, the following factor allows a computation of a cosmological simulation
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
			if(All.ComovingIntegrationOn == 0)
			{
				double fac = All.OmegaLambda * All.Hubble * All.Hubble;

				for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
					for(j = 0; j < 3; j++)
						P[i].g.GravAccel[j] += fac * P[i].Pos[j];
			}
#endif
#endif


			if(ThisTask == 0)
				printf("tree is done.\n");

#else /* gravity is switched off */
			t0 = second();

			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
				for(j = 0; j < 3; j++)
					P[i].g.GravAccel[j] = 0;


#ifdef DISTORTIONTENSORPS
			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{
				P[i].tidal_tensorps[0][0] = 0.0;
				P[i].tidal_tensorps[0][1] = 0.0;
				P[i].tidal_tensorps[0][2] = 0.0;
				P[i].tidal_tensorps[1][0] = 0.0;
				P[i].tidal_tensorps[1][1] = 0.0;
				P[i].tidal_tensorps[1][2] = 0.0;
				P[i].tidal_tensorps[2][0] = 0.0;
				P[i].tidal_tensorps[2][1] = 0.0;
				P[i].tidal_tensorps[2][2] = 0.0;
			}
#endif
#endif /* end of NOGRAVITY */

#ifdef NOGRAVITY
			int k;
#endif

#ifdef SCFPOTENTIAL
			MyDouble xs, ys, zs;
			MyDouble pots, axs, ays, azs;

			if(ThisTask == 0)
			{
				printf("Starting SCF calculation...\n");
				fflush(stdout);
			}

			/* reset the expansion coefficients to zero */
			SCF_reset();
#ifdef SCF_HYBRID
			/*
     calculate SCF coefficients for local DM particles.
     sum them up from all processors, so every processor
     sees the same expansion coefficients 
			 */
			SCF_calc_from_particles();

			/* sum up local coefficients */
			MPI_Allreduce(sinsum, sinsum_all, (SCF_NMAX + 1) * (SCF_LMAX + 1) * (SCF_LMAX + 1), MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);
			MPI_Allreduce(cossum, cossum_all, (SCF_NMAX + 1) * (SCF_LMAX + 1) * (SCF_LMAX + 1), MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);

			/* update local coefficients to global coefficients -> every processor has now complete SCF expansion */
			SCF_collect_update();
			if(ThisTask == 0)
			{
				printf("calculated and collected coefficients.\n");
				fflush(stdout);
			}

#else
			long old_seed, global_seed_min, global_seed_max;

			/*
     resample coefficients for expansion 
     make sure that every processors sees the SAME potential, 
     i.e. has the same seed to generate coefficients  
			 */
			old_seed = scf_seed;
			SCF_calc_from_random(&scf_seed);
			/* check that all cpus have the same random seed (min max must be the same) */
			MPI_Allreduce(&scf_seed, &global_seed_max, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&scf_seed, &global_seed_min, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
			if(ThisTask == 0)
			{
				printf("sampled coefficients with old/new seed = %ld/%ld         min/max=%ld/%ld\n", old_seed, scf_seed,
						global_seed_min, global_seed_max);
				fflush(stdout);
			}
#endif


			/* get accelerations for all active particles based on current expansion */
			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{
				/* convert to unit sphere */
				to_unit(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], &xs, &ys, &zs);
				/* OR: not */
				//xs = P[i].Pos[0]; ys = P[i].Pos[1]; zs = P[i].Pos[2];

				/* evaluate potential and acceleration */
				SCF_evaluate(xs, ys, zs, &pots, &axs, &ays, &azs);

				/* scale to system size and add to acceleration */
#ifdef SCF_HYBRID
				/*
         add missing STAR<-DM force from SCF (was excluded in tree above)
				 */
				if(P[i].Type == 2 || SCF_HYBRID == 2)
				{
#endif
					/* scale */
					P[i].g.GravAccel[0] += All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * axs;
					P[i].g.GravAccel[1] += All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * ays;
					P[i].g.GravAccel[2] += All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * azs;
					/* OR: not */
					//P[i].g.GravAccel[0] += All.G * axs;
					//P[i].g.GravAccel[1] += All.G * ays;
					//P[i].g.GravAccel[2] += All.G * azs;

#ifdef DEBUG
					if(P[i].ID == 150000)
					{
						printf("SCF-ACCEL (scf)   %d  (%g|%g|%g)\n", All.NumCurrentTiStep,
								All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * axs,
								All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * ays,
								All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * azs);
						/* analyic potential of zeroth order of expansion */
						sphere_acc(xs, ys, zs, &axs, &ays, &azs);
						printf("SCF-ACCEL (exact) %d  (%g|%g|%g)\n", All.NumCurrentTiStep, All.G * axs, All.G * ays,
								All.G * azs);
					}
#endif

#ifdef SCF_HYBRID
				}
#endif
			}

			if(ThisTask == 0)
			{
				printf("done.\n");
				fflush(stdout);
			}
#endif


#ifdef SUB_TURB_DRIVING
			sub_turb_add_forces();
#endif


#ifdef STATICISO
			{
				double r, m;
				double dx, dy, dz;

				for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
				{
					dx = P[i].Pos[0];
					dy = P[i].Pos[1];
					dz = P[i].Pos[2];

					r = sqrt(dx * dx + dy * dy + dz * dz);

					if(r > ISO_R200)
						m = ISO_M200;
					else
						m = ISO_M200 * r / ISO_R200;

#ifdef ISO_FRACTION
					m *= ISO_FRACTION;
#endif
					if(r > 0)
					{
						P[i].g.GravAccel[0] += -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
						P[i].g.GravAccel[1] += -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
						P[i].g.GravAccel[2] += -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
					}
				}
			}
#endif

#ifdef DISKPOT
			gravity_tree_subfunc_diskpot();
#endif



#ifdef GROWING_DISK_POTENTIAL
			{
				double mdisk, dx, dy, dz, r, z, aR, az;

				growing_disk_init();

				mdisk = get_disk_mass(All.Time);

				for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
				{
					dx = P[i].Pos[0];
					dy = P[i].Pos[1];
					dz = P[i].Pos[2];

					r = sqrt(dx * dx + dy * dy);
					z = fabs(dz);

					get_disk_forces(r, z, &aR, &az);

					aR *= mdisk;
					az *= mdisk;

					if(r > 0)
					{
						P[i].g.GravAccel[0] += -dx / r * aR;
						P[i].g.GravAccel[1] += -dy / r * aR;
						P[i].g.GravAccel[2] += -dz / z * az;
					}
				}
			}
#endif


#ifdef STATICNFW
			double r, m;

			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{
				double dx, dy, dz;

#ifdef NFW_BOXCENTERED
				dx = P[i].Pos[0] - boxHalf_X;
				dy = P[i].Pos[1] - boxHalf_Y;
				dz = P[i].Pos[2] - boxHalf_Z;
#else
				dx = P[i].Pos[0];
				dy = P[i].Pos[1];
				dz = P[i].Pos[2];
#endif

				r = sqrt(dx * dx + dy * dy + dz * dz);

				m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
				m *= NFW_DARKFRACTION;
#endif
				if(r > 0)
				{
					P[i].g.GravAccel[0] += -All.G * m * dx / (r * r * r);
					P[i].g.GravAccel[1] += -All.G * m * dy / (r * r * r);
					P[i].g.GravAccel[2] += -All.G * m * dz / (r * r * r);

#ifdef DISTORTIONTENSORPS
					double R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
					double Rs = R200 / NFW_C;
					double K = All.G * NFW_M200 / (Rs * (log(1 + NFW_C) - NFW_C / (1 + NFW_C)));
					double r_red = r / Rs;
					double x, y, z;

					x = P[i].Pos[0];
					y = P[i].Pos[1];
					z = P[i].Pos[2];

					P[i].tidal_tensorps[0][0] +=
							-(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (1 / r - x * x / (r * r * r)) -
									K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
											2.0 * Rs * log(1 + r_red) / (r * r * r)) * x * x / (r * r));
					P[i].tidal_tensorps[0][1] +=
							-(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (0 - x * y / (r * r * r)) -
									K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
											2.0 * Rs * log(1 + r_red) / (r * r * r)) * x * y / (r * r));
					P[i].tidal_tensorps[0][2] +=
							-(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (0 - x * z / (r * r * r)) -
									K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
											2.0 * Rs * log(1 + r_red) / (r * r * r)) * x * z / (r * r));
					P[i].tidal_tensorps[1][1] +=
							-(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (1 / r - y * y / (r * r * r)) -
									K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
											2.0 * Rs * log(1 + r_red) / (r * r * r)) * y * y / (r * r));
					P[i].tidal_tensorps[1][2] +=
							-(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (0 - y * z / (r * r * r)) -
									K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
											2.0 * Rs * log(1 + r_red) / (r * r * r)) * y * z / (r * r));
					P[i].tidal_tensorps[2][2] +=
							-(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (1 / r - z * z / (r * r * r)) -
									K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
											2.0 * Rs * log(1 + r_red) / (r * r * r)) * z * z / (r * r));

					P[i].tidal_tensorps[1][0] += P[i].tidal_tensorps[0][1];
					P[i].tidal_tensorps[2][0] += P[i].tidal_tensorps[0][2];
					P[i].tidal_tensorps[2][1] += P[i].tidal_tensorps[1][2];
#endif

				}
			}
#endif

#ifdef CA_BH_ACCRETION
			double r, m;
			int l;

			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{
				r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
#ifdef CA_BH_ACCRETION_MASS_BH
				m = CA_BH_ACCRETION_MASS_BH;
#else
				m = enclosed_mass(r);
#endif
				if(r > CA_BH_ACCRETION_ACCRAD)
				{
					for(l = 0; l < 3; l++)
						P[i].g.GravAccel[l] += -All.G * m * P[i].Pos[l] / (r * r * r);
				}
				else
				{
					P[i].Mass = 0;
				}
			}
#endif



#ifdef STATICPLUMMER
			int l;
			double r;


			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{
				r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);

				for(l = 0; l < 3; l++)
					P[i].g.GravAccel[l] += -P[i].Pos[l] / pow(r * r + 1, 1.5);

#ifdef DISTORTIONTENSORPS
				double x, y, z, r2, f, f2;

				x = P[i].Pos[0];
				y = P[i].Pos[1];
				z = P[i].Pos[2];

				r2 = r * r;;
				f = pow(r2 + 1, 1.5);
				f2 = pow(r2 + 1, 2.5);


				P[i].tidal_tensorps[0][0] += -1.0 / f + 3.0 * x * x / f2;
				P[i].tidal_tensorps[0][1] += -0.0 / f + 3.0 * x * y / f2;
				P[i].tidal_tensorps[0][2] += -0.0 / f + 3.0 * x * z / f2;
				P[i].tidal_tensorps[1][1] += -1.0 / f + 3.0 * y * y / f2;
				P[i].tidal_tensorps[1][2] += -0.0 / f + 3.0 * y * z / f2;
				P[i].tidal_tensorps[2][2] += -1.0 / f + 3.0 * z * z / f2;
				P[i].tidal_tensorps[1][0] += P[i].tidal_tensorps[0][1];
				P[i].tidal_tensorps[2][0] += P[i].tidal_tensorps[0][2];
				P[i].tidal_tensorps[2][1] += P[i].tidal_tensorps[1][2];
#endif
			}
#endif



#ifdef STATICHQ
			double r, m, a;
			int l;


			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{
				r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);

				a = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C *
						sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));

				m = HQ_M200 * pow(r / (r + a), 2);
#ifdef HQ_DARKFRACTION
				m *= HQ_DARKFRACTION;
#endif
				if(r > 0)
				{
					for(l = 0; l < 3; l++)
						P[i].g.GravAccel[l] += -All.G * m * P[i].Pos[l] / (r * r * r);

#ifdef DISTORTIONTENSORPS
					double x, y, z, r2, r3, f, f2, f3;

					x = P[i].Pos[0];
					y = P[i].Pos[1];
					z = P[i].Pos[2];

					r2 = r * r;
					r3 = r * r2;
					f = r + a;
					f2 = f * f;
					f3 = f2 * f;


					P[i].tidal_tensorps[0][0] +=
							All.G * (2.0 * HQ_M200 / (r2 * f3) * x * x + HQ_M200 / (r3 * f2) * x * x - HQ_M200 / (r * f2));
					P[i].tidal_tensorps[0][1] +=
							All.G * (2.0 * HQ_M200 / (r2 * f3) * x * y + HQ_M200 / (r3 * f2) * x * y);
					P[i].tidal_tensorps[0][2] +=
							All.G * (2.0 * HQ_M200 / (r2 * f3) * x * z + HQ_M200 / (r3 * f2) * x * z);
					P[i].tidal_tensorps[1][1] +=
							All.G * (2.0 * HQ_M200 / (r2 * f3) * y * y + HQ_M200 / (r3 * f2) * y * y - HQ_M200 / (r * f2));
					P[i].tidal_tensorps[1][2] +=
							All.G * (2.0 * HQ_M200 / (r2 * f3) * y * z + HQ_M200 / (r3 * f2) * y * z);
					P[i].tidal_tensorps[2][2] +=
							All.G * (2.0 * HQ_M200 / (r2 * f3) * z * z + HQ_M200 / (r3 * f2) * z * z - HQ_M200 / (r * f2));
					P[i].tidal_tensorps[1][0] += P[i].tidal_tensorps[0][1];
					P[i].tidal_tensorps[2][0] += P[i].tidal_tensorps[0][2];
					P[i].tidal_tensorps[2][1] += P[i].tidal_tensorps[1][2];
#endif
				}
			}
#endif

#ifdef STATICBRANDT

			for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
			{

				/* note there is no acceleration in z */
				P[i].g.GravAccel[2] = 0;
				P[i].g.GravAccel[0] -= OmegaR(i, BRANDT_MODE) * OmegaR(i, BRANDT_MODE) * (P[i].Pos[0] - LONG_X / 2.);
				P[i].g.GravAccel[1] -= OmegaR(i, BRANDT_MODE) * OmegaR(i, BRANDT_MODE) * (P[i].Pos[1] - LONG_Y / 2.);
			}

#endif



			/* Now the force computation is finished */

			t1 = WallclockTime = second();
			timeall += timediff(t0, t1);

			/*  gather some diagnostic information */

			timetree = timetree1 + timetree2;
			timewait = timewait1 + timewait2;
			timecomm = timecommsumm1 + timecommsumm2;

			MPI_Reduce(&timetree, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timetree, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timetree1, &sumt1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timetree1, &maxt1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timetree2, &sumt2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timetree2, &maxt2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timewait, &sumwaitall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&timecomm, &sumcommall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&Costtotal, &sum_costtotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&Ewaldcount, &ewaldtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			sumup_longs(1, &n_exported, &n_exported);
			sumup_longs(1, &N_nodesinlist, &N_nodesinlist);

			All.TotNumOfForces += GlobNumForceUpdate;

			plb = (NumPart / ((double) All.TotNumPart)) * NTask;
			MPI_Reduce(&plb, &plb_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&Numnodestree, &maxnumnodes, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

			CPU_Step[CPU_TREEMISC] += timeall - (timetree + timewait + timecomm);
			CPU_Step[CPU_TREEWALK1] += timetree1;
			CPU_Step[CPU_TREEWALK2] += timetree2;
			CPU_Step[CPU_TREESEND] += timecommsumm1;
			CPU_Step[CPU_TREERECV] += timecommsumm2;
			CPU_Step[CPU_TREEWAIT1] += timewait1;
			CPU_Step[CPU_TREEWAIT2] += timewait2;


#ifdef FIXEDTIMEINFIRSTPHASE
			MPI_Reduce(&min_time_first_phase, &min_time_first_phase_glob, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
			if(ThisTask == 0)
			{
				printf("FIXEDTIMEINFIRSTPHASE=%g  min_time_first_phase_glob=%g\n",
						FIXEDTIMEINFIRSTPHASE, min_time_first_phase_glob);
			}
#endif

			if(ThisTask == 0)
			{
				fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
				fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g (%g) iter= %d\n",
						(int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
						(int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
						n_exported / ((double) GlobNumForceUpdate), N_nodesinlist / ((double) n_exported + 1.0e-10),
						iter);
				/* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

				fprintf(FdTimings, "work-load balance: %g (%g %g) rel1to2=%g   max=%g avg=%g\n",
						maxt / (1.0e-6 + sumt / NTask), maxt1 / (1.0e-6 + sumt1 / NTask),
						maxt2 / (1.0e-6 + sumt2 / NTask), sumt1 / (1.0e-6 + sumt1 + sumt2), maxt, sumt / NTask);
				fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
				fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
						maxnumnodes / (All.TreeAllocFactor * All.MaxPart + NTopnodes));
				fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", GlobNumForceUpdate / (sumt + 1.0e-20),
						GlobNumForceUpdate / (1.0e-6 + maxt * NTask),
						((double) (sum_costtotal)) / (1.0e-20 + GlobNumForceUpdate),
						((double) ewaldtot) / (1.0e-20 + GlobNumForceUpdate));
				fprintf(FdTimings, "\n");

				fflush(FdTimings);
			}

			CPU_Step[CPU_TREEMISC] += measure_time();

			double costtotal_new = 0, sum_costtotal_new;
			if(TakeLevel >= 0)
			{
#ifdef _OPENMP
#pragma omp parallel for reduction(+:costtotal_new)
#endif
				for(i = 0; i < NumPart; i++)
					costtotal_new += P[i].GravCost[TakeLevel];
				MPI_Reduce(&costtotal_new, &sum_costtotal_new, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				if(ThisTask == 0)
					printf("relative error in the total number of tree-gravity interactions = %g\n",
							(sum_costtotal - sum_costtotal_new) / sum_costtotal);
				/* can be non-zero if THREAD_SAFE_COSTS is not used (and due to round-off errors). */
			}
		}



		void *gravity_secondary_loop(void *p)
		{
			int j, nodesinlist, dummy, ret;

			while(1)
			{
				LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
				{
					j = NextJ;
					NextJ++;
				}
				UNLOCK_NEXPORT;

				if(j >= Nimport)
					break;

#if !defined(PMGRID)
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
				if(Ewald_iter)
				{
					int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy);

					LOCK_WORKCOUNT;
#ifdef _OPENMP
#pragma omp critical(_workcount_)
#endif
					Ewaldcount += cost;
					UNLOCK_WORKCOUNT;
				}
				else
#endif
				{
					ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy);
					LOCK_WORKCOUNT;
#ifdef _OPENMP
#pragma omp critical(_workcount_)
#endif
					{
						N_nodesinlist += nodesinlist;
						Costtotal += ret;
					}
					UNLOCK_WORKCOUNT;
				}
#else
				ret = force_treeevaluate_shortrange(j, 1, &nodesinlist, &dummy, &dummy);
				LOCK_WORKCOUNT;
#ifdef _OPENMP
#pragma omp critical(_workcount_)
#endif
				{
					N_nodesinlist += nodesinlist;
					Costtotal += ret;
				}
				UNLOCK_WORKCOUNT;
#endif
			}

			return NULL;
		}



		void sum_top_level_node_costfactors(void)
		{
			int i;

			double *costlist = (double *) mymalloc("costlist", NTopnodes * sizeof(double));
			double *costlist_all = (double *) mymalloc("costlist_all", NTopnodes * sizeof(double));

			for(i = 0; i < NTopnodes; i++)
				costlist[i] = Nodes[All.MaxPart + i].GravCost;

			MPI_Allreduce(costlist, costlist_all, NTopnodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			for(i = 0; i < NTopnodes; i++)
				Nodes[All.MaxPart + i].GravCost = costlist_all[i];

			myfree(costlist_all);
			myfree(costlist);
		}





		/*! This function sets the (comoving) softening length of all particle
		 *  types in the table All.SofteningTable[...].  We check that the physical
		 *  softening length is bounded by the Softening-MaxPhys values.
		 */
		void set_softenings(void)
		{
			int i;

			if(All.ComovingIntegrationOn)
			{
				if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
					All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
				else
					All.SofteningTable[0] = All.SofteningGas;

				if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
					All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
				else
					All.SofteningTable[1] = All.SofteningHalo;

				if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
					All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
				else
					All.SofteningTable[2] = All.SofteningDisk;

				if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
					All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
				else
					All.SofteningTable[3] = All.SofteningBulge;

				if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
					All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
				else
					All.SofteningTable[4] = All.SofteningStars;

				if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
					All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
				else
					All.SofteningTable[5] = All.SofteningBndry;
#ifdef SINKS
				All.SofteningTable[5] = All.SinkHsml / All.Time * All.HubbleParam;
#endif
			}
			else
			{
				All.SofteningTable[0] = All.SofteningGas;
				All.SofteningTable[1] = All.SofteningHalo;
				All.SofteningTable[2] = All.SofteningDisk;
				All.SofteningTable[3] = All.SofteningBulge;
				All.SofteningTable[4] = All.SofteningStars;
				All.SofteningTable[5] = All.SofteningBndry;
#ifdef SINKS
				All.SofteningTable[5] = All.SinkHsml;
#endif
			}

			for(i = 0; i < 6; i++)
				All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

			All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];

		}


		/*! This function is used as a comparison kernel in a sort routine. It is
		 *  used to group particles in the communication buffer that are going to
		 *  be sent to the same CPU.
		 */
		int data_index_compare(const void *a, const void *b)
		{
			if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
				return -1;

			if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
				return +1;

			if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
				return -1;

			if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
				return +1;

			if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet))
				return -1;

			if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet))
				return +1;

			return 0;
		}

		static void msort_dataindex_with_tmp(struct data_index *b, size_t n, struct data_index *t)
		{
			struct data_index *tmp;
			struct data_index *b1, *b2;
			size_t n1, n2;

			if(n <= 1)
				return;

			n1 = n / 2;
			n2 = n - n1;
			b1 = b;
			b2 = b + n1;

			msort_dataindex_with_tmp(b1, n1, t);
			msort_dataindex_with_tmp(b2, n2, t);

			tmp = t;

			while(n1 > 0 && n2 > 0)
			{
				if(b1->Task < b2->Task || (b1->Task == b2->Task && b1->Index <= b2->Index))
				{
					--n1;
					*tmp++ = *b1++;
				}
				else
				{
					--n2;
					*tmp++ = *b2++;
				}
			}

			if(n1 > 0)
				memcpy(tmp, b1, n1 * sizeof(struct data_index));

			memcpy(b, t, (n - n2) * sizeof(struct data_index));
		}

		void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
		{
			const size_t size = n * s;

			struct data_index *tmp = (struct data_index *) mymalloc("struct data_index *tmp", size);

			msort_dataindex_with_tmp((struct data_index *) b, n, tmp);

			myfree(tmp);
		}




		void *gravity_primary_loop(void *p)
		{
			int i, j, ret;
			int thread_id = *(int *) p;

			int *exportflag, *exportnodecount, *exportindex;

			exportflag = Exportflag + thread_id * NTask;
			exportnodecount = Exportnodecount + thread_id * NTask;
			exportindex = Exportindex + thread_id * NTask;

			/* Note: exportflag is local to each thread */
			for(j = 0; j < NTask; j++)
				exportflag[j] = -1;



			//manos variables
			int m_index, m_active_part[All.MaxPart], m_num_active_part, m_temp, m_break;
			m_num_active_part=0;




			//manos get active

			while(NextParticle>=0){
				//printf("ProcessedFlag[%d]=%d \n",NextParticle, ProcessedFlag[NextParticle]);
				LOCK_NEXPORT;
				m_index=NextParticle;
				m_active_part[m_num_active_part]=m_index;
				m_temp=ProcessedFlag[m_index];
				ProcessedFlag[m_index]=0;
				NextParticle = NextActiveParticle[m_index];
				m_num_active_part++;
				ProcessedFlag[m_index]=m_temp;
				UNLOCK_NEXPORT;
			}

			//manos//simple arrays
			MyDouble m_PPos[All.MaxPart][3], m_PMass[All.MaxPart];
			short int m_PType[All.MaxPart];
			MyFloat m_POldAcc[All.MaxPart];



			for(m_index=0; m_index<All.MaxPart; m_index++)
			{
				m_PPos[m_index][0]=P[m_index].Pos[0];
				m_PPos[m_index][1]=P[m_index].Pos[1];
				m_PPos[m_index][2]=P[m_index].Pos[2];
				m_PType[m_index]=P[m_index].Type;
				m_POldAcc[m_index]=P[m_index].OldAcc;
				m_PMass[m_index]=P[m_index].Mass;
			}

			//manos//simple arrays #2
			int m_Nextnode[All.MaxPart + NTopnodes];


			//manos//simple arrays for output
			MyLongDouble m_out_PdGravAccel[All.MaxPart + NTopnodes][3];




			for(m_index=0; m_index<(All.MaxPart + NTopnodes); m_index++)
			{
				m_Nextnode[m_index]= Nextnode[m_index];
				m_out_PdGravAccel[m_index][0] = P[m_index].g.dGravAccel[0];
				m_out_PdGravAccel[m_index][1] = P[m_index].g.dGravAccel[1];
				m_out_PdGravAccel[m_index][2] = P[m_index].g.dGravAccel[2];

			}

			//NextParticle=m_temp2;
			//manos//printf("Active parts with for loop: %d \n", m_index);
			m_break = 0;
			MyLongDouble m_acc_x;


			//manos//shortrange vars

			//input
			int *m_exportflag = exportflag;
			int *m_exportnodecount = exportnodecount;
			int *m_exportindex = exportindex;

			//private
			struct NODE *m_nop = 0;
			int m_no, m_nodesinlist, m_ptype, m_ninteractions, m_nexp, m_tabindex, m_task, m_listindex = 0;
			double m_r2, m_dx, m_dy, m_dz, m_mass, m_r, m_fac, m_u, m_h, m_h_inv, m_h3_inv;
			double m_dxx, m_dyy, m_dzz, m_pdxx, m_pdyy, m_pdzz;
			double m_pos_x, m_pos_y, m_pos_z, m_aold;
			double m_eff_dist;
			double m_rcut, m_asmth, m_asmthfac, m_rcut2, m_dist;
			MyLongDouble m_acc_y, m_acc_z;
			// cache some global vars in local vars to help compiler with alias analysis
			int m_maxPart = All.MaxPart;
			int m_bunchSize = All.BunchSize;
			int m_maxNodes = MaxNodes;
			integertime m_ti_Current = All.Ti_Current;
			double m_errTol2 = All.ErrTolTheta * All.ErrTolTheta;
			int m_exitFlag = 0;

			//used/1ST ////////////////////////////////////////////////////////////////////////////////////////////
			double m_xtmp;

			//initializations


			m_maxPart = All.MaxPart;
			m_bunchSize = All.BunchSize;
			m_maxNodes = MaxNodes;
			m_ti_Current = All.Ti_Current;
			m_errTol2 = All.ErrTolTheta * All.ErrTolTheta;

			//manos//other

			m_rcut = All.Rcut[0];
			m_asmth = All.Asmth[0];





			//manos//end shotrange vars
			//manos acc

			//#pragma acc data copy(ProcessedFlag[0:All.MaxPart],P[0:All.MaxPart], All, m_break) //create(m_acc_x, m_acc_y, m_acc_z, m_maxPart, m_bunchSize, m_maxNodes, \
			m_ti_Current, m_errTol2, m_rcut, m_asmth, m_asmthfac, m_rcut2, m_dist, m_eff_dist, m_no, m_exitFlag)
			{

				//#pragma acc parallel loop copy(ProcessedFlag[0:All.MaxPart])
#pragma acc parallel loop gang worker vector //private(m_acc_x, m_acc_y, m_acc_z, m_maxPart, m_bunchSize, m_maxNodes, \
				m_ti_Current, m_errTol2, m_rcut, m_asmth, m_asmthfac, m_rcut2, m_dist, m_eff_dist, m_index, m_break, m_no, m_exitFlag) \
				reduction(+:m_break,m_exitFlag)
				for (m_index=0; m_index<m_num_active_part; m_index++) //manos
				{
					int exitFlag = 0;
					LOCK_NEXPORT;

					//if(BufferFullFlag != 0 || m_break)
					if(m_break)
					{
						exitFlag = 1;
					}
					else
					{
						i = m_active_part[m_index];
						ProcessedFlag[i] = 0;
					}
					UNLOCK_NEXPORT;

					if(exitFlag || m_break)
					{
						m_break=1;
					}
					else
					{

						//inlining it//ret = force_treeevaluate_shortrange(i, 0, exportflag, exportnodecount, exportindex);


						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						////////////// INLINE shortrange start...   //////////////////////////////////////////

						{

							int m_target = i;



							m_nop=0;
							m_listindex=0;
							m_exitFlag = 0;


							m_acc_x = 0;
							m_acc_y = 0;
							m_acc_z = 0;
							m_ninteractions = 0;
							m_nodesinlist = 0;




							m_pos_x = m_PPos[m_target][0];//P[m_target].Pos[0];
							m_pos_y = m_PPos[m_target][1];//P[m_target].Pos[1];
							m_pos_z = m_PPos[m_target][2];//P[m_target].Pos[2];
							m_ptype = m_PType[m_target];//P[m_target].Type;
							m_aold = All.ErrTolForceAcc * m_POldAcc[m_target];//P[m_target].OldAcc;


							m_rcut2 = m_rcut * m_rcut;

							m_asmthfac = 0.5 / m_asmth * (NTAB / 3.0);

							//used/15TH ////////////////////////////////////////////////////////////////////////////////////////////
							m_h = All.ForceSoftening[m_ptype];
							m_h_inv = 1.0 / m_h;
							m_h3_inv = m_h_inv * m_h_inv * m_h_inv;

							m_no = m_maxPart;		/* root node */


							while(m_no >= 0)
							{
								while(m_no >= 0)
								{
									if(m_no < m_maxPart)
									{
										/* the index of the node is the index of the particle */
										//manos//
										//if(P[m_no].Ti_current != m_ti_Current)
										//{
										//	LOCK_PARTNODEDRIFT;
										//#pragma omp critical(_partnodedrift_)
										//	drift_particle(m_no, m_ti_Current);
										//	printf("Drift-particle()");
										//	UNLOCK_PARTNODEDRIFT;
										//}

										m_dx = m_PPos[m_no][0] - m_pos_x;
										m_dy = m_PPos[m_no][1] - m_pos_y;
										m_dz = m_PPos[m_no][2] - m_pos_z;

										m_dx = NEAREST(m_dx);
										m_dy = NEAREST(m_dy);
										m_dz = NEAREST(m_dz);

										m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

										m_mass = m_PMass[m_no];

										if(TakeLevel >= 0)
										{
											//entered!!!
											LOCK_WORKCOUNT;

											P[m_no].GravCost[TakeLevel] += 1.0;
											UNLOCK_WORKCOUNT;
										}
										m_no = m_Nextnode[m_no];//Nextnode[m_no];
									}
									else			/* we have an  internal node */
									{
										//entered!!!
										if(m_no >= m_maxPart + m_maxNodes)	/* pseudo particle */
										{

											{
												if(m_exportflag[m_task = DomainTask[m_no - (m_maxPart + m_maxNodes)]] != m_target)
												{
													m_exportflag[m_task] = m_target;
													m_exportnodecount[m_task] = NODELISTLENGTH;
												}

												if(m_exportnodecount[m_task] == NODELISTLENGTH)
												{
													//int m_exitFlag=0;

													//manos//					#pragma omp critical(_nexport_)

													if(Nexport >= m_bunchSize)
													{
														/* out of buffer space. Need to discard work for this particle and interrupt */
														//BufferFullFlag = 1;
														m_exitFlag = 1;
													}
													else
													{
														m_nexp = Nexport;
														Nexport++;
													}



													m_exportnodecount[m_task] = 0;
													m_exportindex[m_task] = m_nexp;
													DataIndexTable[m_nexp].Task = m_task;
													DataIndexTable[m_nexp].Index = m_target;
													DataIndexTable[m_nexp].IndexGet = m_nexp;

												}



												DataNodeList[m_exportindex[m_task]].NodeList[m_exportnodecount[m_task]++] =
														DomainNodeIndex[m_no - (m_maxPart + m_maxNodes)];

												if(m_exportnodecount[m_task] < NODELISTLENGTH)
													DataNodeList[m_exportindex[m_task]].NodeList[m_exportnodecount[m_task]] = -1;

											}

										} //pseudoparticle region end



										m_nop = &Nodes[m_no];


										//	if(m_nop->Ti_current != m_ti_Current)
										//	{
										//		LOCK_PARTNODEDRIFT;
										//		#pragma omp critical(_partnodedrift_)
										//		//force_drift_node(m_no, m_ti_Current);
										//		//printf("Force_drift_node()");
										//		UNLOCK_PARTNODEDRIFT;
										//	}

										m_mass = m_nop->u.d.mass;

										m_dx = m_nop->u.d.s[0] - m_pos_x;
										m_dy = m_nop->u.d.s[1] - m_pos_y;
										m_dz = m_nop->u.d.s[2] - m_pos_z;

										m_dx = NEAREST(m_dx);
										m_dy = NEAREST(m_dy);
										m_dz = NEAREST(m_dz);
										m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

										/* check whether we can stop walking along this branch */
										if(m_r2 > m_rcut2)
										{
											m_eff_dist = m_rcut + 0.5 * m_nop->len;

											m_dist = NEAREST(m_nop->center[0] - m_pos_x);
											if(m_dist < -m_eff_dist || m_dist > m_eff_dist)
											{
												m_no = m_nop->u.d.sibling;
												continue;
											}

											m_dist = NEAREST(m_nop->center[1] - m_pos_y);
											if(m_dist < -m_eff_dist || m_dist > m_eff_dist)
											{
												m_no = m_nop->u.d.sibling;
												continue;
											}

											m_dist = NEAREST(m_nop->center[2] - m_pos_z);
											if(m_dist < -m_eff_dist || m_dist > m_eff_dist)
											{
												m_no = m_nop->u.d.sibling;
												continue;
											}
										}


										if(m_errTol2)	/* check Barnes-Hut opening criterion */
										{
											if(m_nop->len * m_nop->len > m_r2 * m_errTol2)
											{
												/* open cell */
												m_no = m_nop->u.d.nextnode;
												continue;
											}
										}
										else		/* check relative opening criterion */
										{

											//used/23RD ////////////////////////////////////////////////////////////////////////////////////////////
											if(m_mass * m_nop->len * m_nop->len > m_r2 * m_r2 * m_aold)
											{
												/* open cell */
												m_no = m_nop->u.d.nextnode;
												continue;
											}

											/* check in addition whether we lie inside the cell */

											if(fabs(m_nop->center[0] - m_pos_x) < 0.60 * m_nop->len)
											{
												if(fabs(m_nop->center[1] - m_pos_y) < 0.60 * m_nop->len)
												{
													if(fabs(m_nop->center[2] - m_pos_z) < 0.60 * m_nop->len)
													{
														m_no = m_nop->u.d.nextnode;
														continue;
													}
												}
											}

										}


										if(TakeLevel >= 0)
										{
											LOCK_WORKCOUNT;
											m_nop->GravCost += 1.0;
											UNLOCK_WORKCOUNT;
										}

										m_no = m_nop->u.d.sibling;	/* ok, node can be used */


									}




									m_r = sqrt(m_r2);

									if(m_r >= m_h)
									{
										m_fac = m_mass / (m_r2 * m_r);
									}

									else
									{

										m_u = m_r * m_h_inv;
										if(m_u < 0.5)
											m_fac = m_mass * m_h3_inv * (10.666666666667 + m_u * m_u * (32.0 * m_u - 38.4));
										else
											m_fac =
													m_mass * m_h3_inv * (21.333333333333 - 48.0 * m_u +
															38.4 * m_u * m_u - 10.666666666667 * m_u * m_u * m_u - 0.066666666667 / (m_u * m_u * m_u));


									}

									m_tabindex = (int) (m_asmthfac * m_r);

									if(m_tabindex < NTAB)
									{
										m_fac *= shortrange_table[m_tabindex];

										m_acc_x += FLT(m_dx * m_fac);
										m_acc_y += FLT(m_dy * m_fac);
										m_acc_z += FLT(m_dz * m_fac);
									}

									m_ninteractions++;



								}

							}//outer while loop end


							/* store result at the proper place */
							m_out_PdGravAccel[m_target][0] = m_acc_x;
							m_out_PdGravAccel[m_target][1] = m_acc_y;
							m_out_PdGravAccel[m_target][2] = m_acc_z;



							ret = m_ninteractions;
						}


						//////////////INLINE shortrange finish... ////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


						if(ret < 0)
							m_break=1;		/* export buffer has filled up */


						if(m_break){
							continue;
						}
						else{
							LOCK_WORKCOUNT;

							Costtotal += ret;
							UNLOCK_WORKCOUNT;



							ProcessedFlag[i] = 1;	/* particle successfully finished */
						}
					}///manos//end m_break part (containing whole shortrange function)



				}//manos// end of for loop

			}//manos// end of data region


			//manos//copy computed values to cpu values
			for(m_index=0; m_index<(All.MaxPart + NTopnodes); m_index++)
			{
				m_Nextnode[m_index]= Nextnode[m_index];
				P[m_index].g.dGravAccel[0] = m_out_PdGravAccel[m_index][0];
				P[m_index].g.dGravAccel[1] = m_out_PdGravAccel[m_index][1];
				P[m_index].g.dGravAccel[2] = m_out_PdGravAccel[m_index][2];

			}

			if(m_break)BufferFullFlag = 1;



			return NULL;
		}

