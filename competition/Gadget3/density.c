#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef CS_MODEL
#include "cs_model/cs_metals.h"
#endif

#ifdef NUM_THREADS
#include <pthread.h>
#endif

#ifdef FS_ETA2_DYN
#define OUTPUT_ROTB
#endif

#ifdef NUM_THREADS
extern pthread_mutex_t mutex_nexport;

extern pthread_mutex_t mutex_partnodedrift;

#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

#ifdef LT_STELLAREVOLUTION
int LastActive, SaveFirstActiveParticle;
#endif

#if defined(LT_MV_CHEMICALDIFFUSION)
double hubble_a;
#endif


#ifdef ADAPTGRAVSOFT
 /*! Structure for communication during the computation of the adaptive gravitational softening.
  *  Holds data that is sent to other processors.
  */
static struct AGS_densdata_in
{
  MyDouble Pos[3];
  MyFloat AGS_Hsml;
  int NodeList[NODELISTLENGTH];
}
 *AGS_DensDataIn, *AGS_DensDataGet;


static struct AGS_densdata_out
{
  MyLongDouble AGS_Rho;
  MyLongDouble AGS_Ngb;
#ifdef AGS_OUTPUTNGBS
  int AGS_REALNgb;
#endif
  MyFloat AGS_omega;
  MyFloat AGS_zeta;

}
 *AGS_DensDataResult, *AGS_DensDataOut;
#endif


struct kernel_density
{
  double dx, dy, dz;
  double r;
  double dvx, dvy, dvz;
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
#ifdef WINDS
  MyFloat DelayTime;
#endif
  int NodeList[NODELISTLENGTH];
#if defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB) || defined(OUTPUT_ROTB)
  MyFloat BPred[3];
#endif
#ifdef EULERPOTENTIALS
  MyFloat EulerA, EulerB;
#endif
#ifdef CS_MODEL
  MyFloat DensityOld, Entropy;
#endif
#ifdef AB_SHOCK
  MyFloat Pressure;
#endif
#if defined(LT_STELLAREVOLUTION) || defined(KD_FRICTION_DYNAMIC)
  int Type;
#endif
}
 *DensDataIn, *DensDataGet;

static struct densdata_out
{
  MyLongDouble Rho;
  MyLongDouble DhsmlDensity;
  MyLongDouble Ngb;
#ifndef NAVIERSTOKES
  MyLongDouble Div, Rot[3];
#else
  MyFloat DV[3][3];
#endif
#ifdef JD_VTURB
  MyFloat Vbulk[3];
#endif
#if defined(JD_VTURB) || defined(HIGH_ORDER_INDUCTION)
  int TrueNGB;
#endif

#ifdef AB_SHOCK
  MyFloat Shock_N1[3];
#endif

#ifdef MAGNETIC
#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB)
  MyFloat RotB[3];
#endif
#ifdef TRACEDIVB
  MyFloat divB;
#endif
#ifdef HIGH_ORDER_INDUCTION
  MyFloat Chi[3][3];
  MyFloat Xix[3], Xiy[3], Xiz[3];
#endif
#if defined (EULERPOTENTIALS)
  MyFloat dEulerA[3];
  MyFloat dEulerB[3];
#endif
#endif

#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat Grad_ngamma[3][N_BINS];
#endif

#if defined(BLACK_HOLES)
  MyLongDouble SmoothedEntr;
  int BH_TimeBinGasNeighbor;
#endif
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif

#if defined(BLACK_HOLES) || defined(VS_TURB) || defined(AB_TURB)
  MyLongDouble GasVel[3];
#endif
#ifdef BLACK_HOLES
#ifdef KD_FRICTION
  MyFloat SurroundingVel[3];
  MyFloat SurroundingDensity;
#endif
#endif

#if defined(LT_BH_GUESSHSML) || defined(LT_STARS_GUESSHSML)
  MyFloat SmoothedHsml;
  MyFloat SmoothedRho;
#endif

#ifdef LT_MV_CHEMICALDIFFUSION
  /*int vrms_NGB; */
  MyFloat v_rms;
  float GradAbnd[LT_NMetP][3];
  float AvgAbnd[LT_NMetP];
#endif
}
 *DensDataResult, *DensDataOut;

void particle2in_density(struct densdata_in *in, int i);
void out2particle_density(struct densdata_out *out, int i, int mode);
void density_evaluate_extra_physics_gas(struct densdata_in *local, struct densdata_out *out,
					struct kernel_density *kernel, int j);

#ifdef LT_MV_CHEMICALDIFFUSION
int updatediffusion_and_normalize();
#endif

void particle2in_density(struct densdata_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];
  in->Hsml = PPP[i].Hsml;

#if defined(LT_STELLAREVOLUTION) || defined(KD_FRICTION_DYNAMIC)
  in->Type = P[i].Type;
#endif

#if defined(BLACK_HOLES)
  if(P[i].Type != 0)
    for(k = 0; k < 3; k++)
      in->Vel[k] = 0;
  else
#endif
    for(k = 0; k < 3; k++)
      in->Vel[k] = SphP[i].VelPred[k];

  if(P[i].Type == 0)
    {
#ifdef CS_MODEL
      in->DensityOld = SphP[i].DensityOld;
      in->Entropy = SphP[i].Entropy;
#endif
#ifdef AB_SHOCK
      in->Pressure = SphP[i].Pressure;
#endif
#ifdef EULERPOTENTIALS
      in->EulerA = SphP[i].EulerA;
      in->EulerB = SphP[i].EulerB;
#endif

#ifdef WINDS
      in->DelayTime = SphP[i].DelayTime;
#endif

#if defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB) || defined(OUTPUT_ROTB)
      for(k = 0; k < 3; k++)
#ifdef SFR
	in->BPred[k] = SphP[i].b2.BPred[k] * pow(1. - SphP[i].XColdCloud, 2. * POW_CC);
#else
	in->BPred[k] = SphP[i].b2.BPred[k];
#endif
#endif
    }
}

void out2particle_density(struct densdata_out *out, int i, int mode)
{
  int k;

  ASSIGN_ADD(PPP[i].n.dNumNgb, out->Ngb, mode);

  if(P[i].Type == 0)
    {
      ASSIGN_ADD(SphP[i].d.dDensity, out->Rho, mode);
      ASSIGN_ADD(SphP[i].h.dDhsmlDensityFactor, out->DhsmlDensity, mode);

#if defined(VS_TURB) || defined(AB_TURB)
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(SphP[i].SmoothedVel[k], out->GasVel[k], mode);
#endif

#ifndef NAVIERSTOKES
      ASSIGN_ADD(SphP[i].v.dDivVel, out->Div, mode);
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(SphP[i].r.dRot[k], out->Rot[k], mode);
#else
      int j;
      for(k = 0; k < 3; k++)
	for(j = 0; j < 3; j++)
	  ASSIGN_ADD(SphP[i].u.DV[k][j], out->DV[k][j], mode);
#endif

#ifdef CONDUCTION_SATURATION
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(SphP[i].GradEntr[k], out->GradEntr[k], mode);
#endif

#ifdef LT_MV_CHEMICALDIFFUSION
      /*  SphP[place].vrms_NGB += DensDataOut[j].vrm_NGB; */
      int cc;
      for(k = 0; k < LT_NMetP; k++)
	{
	  for(cc = 0; cc < 3; cc++)
	    ASSIGN_ADD(SphP[i].GradAbnd[k][cc], out->GradAbnd[k][cc], mode);
	  ASSIGN_ADD(SphP[i].AvgAbnd[k], out->AvgAbnd[k], mode);
	}
#if !defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
      ASSIGN_ADD(SphP[i].v_rms.dvel_rms, out->v_rms, mode);
#endif
#endif


#ifdef RADTRANSFER_FLUXLIMITER
      for(k = 0; k < N_BINS; k++)
	for(j = 0; j < 3; j++)
	  ASSIGN_ADD(SphP[i].Grad_ngamma[j][k], out->Grad_ngamma[j][k], mode);
#endif

#ifdef AB_SHOCK
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(SphP[i].Shock_N1[k], out->Shock_N1[k], mode);
#endif

#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB)
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(SphP[i].RotB[k], out->RotB[k], mode);
#endif

#ifdef HIGH_ORDER_INDUCTION
      for(k = 0; k < 3; k++)
	{
	  for(j = 0; j < 3; j++)
	    ASSIGN_ADD(SphP[i].Chi[j][k], out->Chi[j][k], mode);
	  ASSIGN_ADD(SphP[i].Xix[k], out->Xix[k], mode);
	  ASSIGN_ADD(SphP[i].Xiy[k], out->Xiy[k], mode);
	  ASSIGN_ADD(SphP[i].Xiz[k], out->Xiz[k], mode);
	}
#endif

#ifdef TRACEDIVB
      ASSIGN_ADD(SphP[i].divB, out->divB, mode);
#endif

#ifdef JD_VTURB
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(SphP[i].Vbulk[k], out->Vbulk[k], mode);
#endif

#if defined(EULERPOTENTIALS)
      for(k = 0; k < 3; k++)
	{
	  ASSIGN_ADD(SphP[i].b1.dEulerA[k], out->dEulerA[k], mode);
	  ASSIGN_ADD(SphP[i].b2.dEulerB[k], out->dEulerB[k], mode);
	}
#endif
    }

#if defined(JD_VTURB) || defined(HIGH_ORDER_INDUCTION)
  ASSIGN_ADD(P[i].TrueNGB, out->TrueNGB, mode);
#endif

#if (defined(RADTRANSFER) && defined(EDDINGTON_TENSOR_STARS)) || defined(SNIA_HEATING)
  if(P[i].Type == 4)
    ASSIGN_ADD(P[i].DensAroundStar, out->Rho, mode);
#endif

#ifdef LT_STELLAREVOLUTION
  if(P[i].Type == 4)
    {
#ifdef LT_DONTUSE_DENSITY_in_WEIGHT
      ASSIGN_ADD(MetP[P[i].pt.MetID].weight, out->Rho, mode);
#endif

#ifdef LT_STARS_GUESSHSML
      ASSIGN_ADD(MetP[P[i].pt.MetID].mean_hsml, out->SmoothedHsml, mode);
      ASSIGN_ADD(MetP[P[i].pt.MetID].mean_rho, out->SmoothedRho, mode);
#endif
    }
#endif


#ifdef BLACK_HOLES
  if(P[i].Type == 5)
    {
      ASSIGN_ADD(BPP(i).b1.dBH_Density, out->Rho, mode);
      ASSIGN_ADD(BPP(i).b2.dBH_Entropy, out->SmoothedEntr, mode);
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(BPP(i).b3.dBH_SurroundingGasVel[k], out->GasVel[k], mode);
#ifdef KD_FRICTION
      ASSIGN_ADD(BPP(i).BH_SurroundingDensity, out->SurroundingDensity, mode);
      for(k = 0; k < 3; k++)
	ASSIGN_ADD(BPP(i).BH_SurroundingVel[k], out->SurroundingVel[k], mode);
#endif
#ifdef LT_BH_GUESSHSML
      ASSIGN_ADD(BPP(i).mean_hsml, out->SmoothedHsml, mode);
      ASSIGN_ADD(BPP(i).mean_rho, out->SmoothedRho, mode);
#endif
      if(mode == 0)
	BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;
      else
	{
	  if(BPP(i).BH_TimeBinGasNeighbor > out->BH_TimeBinGasNeighbor)
	    BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;
	}
    }
#endif
}

/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */


/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  MyFloat *Left, *Right;
  int i, j, k, ndone, ndone_flag, npleft, iter = 0;
  int ngrp, recvTask, place;
  long long ntot;
  double fac;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported = 0;
  int redo_particle;

#ifdef HIGH_ORDER_INDUCTION
  double detChi, Chi[3][3], Xix[3], Xiy[3], Xiz[3];
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#ifdef NAVIERSTOKES
  double dvel[3][3];
  double rotx, roty, rotz;
#endif

#if defined(SOFTEREQS) || defined(MHM)
  double a3inv, afac;

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      afac = pow(All.Time, 3 * GAMMA_MINUS1);
    }
  else
    a3inv = afac = 1;
#endif

#if defined(EULERPOTENTIALS) || defined(TRACEDIVB)
  double efak;

  if(All.ComovingIntegrationOn)
    efak = 1. / All.Time / All.HubbleParam;
  else
    efak = 1;
#endif

#ifdef LT_MV_CHEMICALDIFFUSION
  int cc;
#endif

  CPU_Step[CPU_DENSMISC] += measure_time();

  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

#if defined(LT_STELLAREVOLUTION)
  unsigned int appended = 0, tot_appended, alreadyactive = 0, tot_alreadyactive, hitMaxChemSpreadL = 0;
#ifndef LT_LOCAL_IRA
  if(N_stars > 0)
#endif
    appended = append_chemicallyactive_particles(&alreadyactive);
  MPI_Reduce(&appended, &tot_appended, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&alreadyactive, &tot_alreadyactive, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0 && (tot_appended > 0 || tot_alreadyactive > 0))
    printf("%u chemically active particles queued for density calculation (%u already active)..\n",
	   tot_appended, tot_alreadyactive);
  fflush(stdout);
#endif

#ifdef LT_MV_CHEMICALDIFFUSION
  updatediffusion_and_normalize();
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
      if(density_isactive(i))
	{
	  Left[i] = Right[i] = 0;

#ifdef BLACK_HOLES
	  P[i].SwallowID = 0;
#endif
#if defined(BLACK_HOLES) && defined(FLTROUNDOFFREDUCTION)
	  if(P[i].Type == 0)
	    SphP[i].i.dInjected_BH_Energy = SphP[i].i.Injected_BH_Energy;
#endif
	}
    }

  /* allocate buffers to arrange communication */

#ifdef KD_BUFFER_MANAGEMENT
  size_t MyBufferSize = 0.7 * FreeBytes / (1024.0 * 1024.0);
#else
  size_t MyBufferSize = All.BufferSize;
#endif

  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					   sizeof(struct densdata_in) + sizeof(struct densdata_out) +
					   sizemax(sizeof(struct densdata_in), sizeof(struct densdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

#ifdef KD_BUFFER_MANAGEMENT
  if(ThisTask == 0)
    printf("Density: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	   FreeBytes / (1024.0 * 1024.0));
#endif

  t0 = second();

  desnumngb = All.DesNumNgb;
  desnumngbdev = All.MaxNumNgbDeviation;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {

      NextParticle = FirstActiveParticle;	/* beginn with this index */

      do
	{
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
	  pthread_mutex_init(&mutex_nexport, NULL);
	  pthread_mutex_init(&mutex_partnodedrift, NULL);

	  TimerFlag = 0;

	  for(j = 0; j < NUM_THREADS - 1; j++)
	    {
	      threadid[j] = j + 1;
	      pthread_create(&mythreads[j], &attr, density_evaluate_primary, &threadid[j]);
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
	    density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
	  }

#ifdef NUM_THREADS
	  for(j = 0; j < NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);
#endif

	  tend = second();
	  timecomp1 += timediff(tstart, tend);

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
		  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n", ThisTask, P[NextParticle].Type,
			 P[NextParticle].Pos[0], P[NextParticle].Pos[1], P[NextParticle].Pos[2],
			 P[NextParticle].Mass);
		  if(P[NextParticle].Type == 0)
		    printf("   rho=%g hsml=%g\n", SphP[NextParticle].d.Density, PPP[NextParticle].Hsml);

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

	  DensDataGet = (struct densdata_in *) mymalloc("DensDataGet", Nimport * sizeof(struct densdata_in));
	  DensDataIn = (struct densdata_in *) mymalloc("DensDataIn", Nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      particle2in_density(&DensDataIn[j], place);

	      memcpy(DensDataIn[j].NodeList,
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
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(DensDataIn);
	  DensDataResult =
	    (struct densdata_out *) mymalloc("DensDataResult", Nimport * sizeof(struct densdata_out));
	  DensDataOut =
	    (struct densdata_out *) mymalloc("DensDataOut", Nexport * sizeof(struct densdata_out));

	  report_memory_usage(&HighMark_sphdensity, "SPH_DENSITY");

	  /* now do the particles that were sent to us */

	  tstart = second();

	  NextJ = 0;

#ifdef NUM_THREADS
	  for(j = 0; j < NUM_THREADS - 1; j++)
	    pthread_create(&mythreads[j], &attr, density_evaluate_secondary, &threadid[j]);
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
	    density_evaluate_secondary(&mainthreadid);
	  }

#ifdef NUM_THREADS
	  for(j = 0; j < NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);

	  pthread_mutex_destroy(&mutex_partnodedrift);
	  pthread_mutex_destroy(&mutex_nexport);
	  pthread_attr_destroy(&attr);
#endif

	  tend = second();
	  timecomp2 += timediff(tstart, tend);

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
		      MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = second();
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
	      out2particle_density(&DensDataOut[j], place, 1);
	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(DensDataOut);
	  myfree(DensDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);

#ifdef FLTROUNDOFFREDUCTION
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(density_isactive(i))
	  {
	    PPP[i].n.NumNgb = FLT(PPP[i].n.dNumNgb);

	    if(P[i].Type == 0)
	      {
		SphP[i].d.Density = FLT(SphP[i].d.dDensity);
		SphP[i].h.DhsmlDensityFactor = FLT(SphP[i].h.dDhsmlDensityFactor);
		SphP[i].v.DivVel = FLT(SphP[i].v.dDivVel);
		for(j = 0; j < 3; j++)
		  SphP[i].r.Rot[j] = FLT(SphP[i].r.dRot[j]);
	      }

#ifdef BLACK_HOLES
	    if(P[i].Type == 5)
	      {
		BPP(i).b1.BH_Density = FLT(BPP(i).b1.dBH_Density);
		BPP(i).b2.BH_Entropy = FLT(BPP(i).b2.dBH_Entropy);
		for(j = 0; j < 3; j++)
		  BPP(i).b3.BH_SurroundingGasVel[j] = FLT(BPP(i).b3.dBH_SurroundingGasVel[j]);
	      }
#endif
	  }
#endif


      /* do final operations on results */
      tstart = second();

      npleft = 0;

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
      int il;
#pragma omp parallel for reduction(+:npleft) private(i,k,desnumngb,desnumngbdev,redo_particle)
      for(il = 0; il < NActivePart; il++)
	{
	  i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
#endif
	  if(density_isactive(i))
	    {
#ifdef BLACK_HOLES
	      if(P[i].Type == 0 && P[i].Mass > 0)
#else
	      if(P[i].Type == 0)
#endif
		{
#if defined(LT_STELLAREVOLUTION)
		  if(TimeBinActive[P[i].TimeBin])
		    {
#endif
		      if(SphP[i].d.Density > 0)
			{
			  SphP[i].h.DhsmlDensityFactor *= PPP[i].Hsml / (NUMDIMS * SphP[i].d.Density);
			  if(SphP[i].h.DhsmlDensityFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
			    SphP[i].h.DhsmlDensityFactor = 1 / (1 + SphP[i].h.DhsmlDensityFactor);
			  else
			    SphP[i].h.DhsmlDensityFactor = 1;

#if defined(WC4_BIAS_CORRECTION) && defined(WENDLAND_C4_KERNEL)	/* Dehnen & Aly 2012, eq (15), (16) */
			  struct kernel_density kernel;
			  kernel_hinv(SphP[i].Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
			  kernel_main(0, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);

			  /* subtract self-contribution of the kernel */
			  SphP[i].d.Density -=
			    0.0132 * pow(P[i].TrueNGB * 0.01, -1.579) * P[i].Mass * kernel.wk;
#endif

#if defined(WC6_BIAS_CORRECTION) && defined(WENDLAND_C6_KERNEL)	/* Dehnen & Aly 2012, eq (15), (16) */
			  struct kernel_density kernel;
			  kernel_hinv(SphP[i].Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
			  kernel_main(0, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);

			  /* subtract self-contribution of the kernel */
			  SphP[i].d.Density -=
			    0.0116 * pow(P[i].TrueNGB * 0.01, -2.236) * P[i].Mass * kernel.wk;
#endif

#ifndef NAVIERSTOKES
#if defined (VS_TURB) || defined (AB_TURB) || defined(OUTPUT_VORTICITY)
			  SphP[i].Vorticity[0] = SphP[i].r.Rot[0] / SphP[i].d.Density;
			  SphP[i].Vorticity[1] = SphP[i].r.Rot[1] / SphP[i].d.Density;
			  SphP[i].Vorticity[2] = SphP[i].r.Rot[2] / SphP[i].d.Density;
#endif
			  SphP[i].r.CurlVel = sqrt(SphP[i].r.Rot[0] * SphP[i].r.Rot[0] +
						   SphP[i].r.Rot[1] * SphP[i].r.Rot[1] +
						   SphP[i].r.Rot[2] * SphP[i].r.Rot[2]) / SphP[i].d.Density;

			  SphP[i].v.DivVel /= SphP[i].d.Density;
#else
			  for(k = 0; k < 3; k++)
			    {
			      dvel[k][0] = SphP[i].u.DV[k][0] / SphP[i].d.Density;
			      dvel[k][1] = SphP[i].u.DV[k][1] / SphP[i].d.Density;
			      dvel[k][2] = SphP[i].u.DV[k][2] / SphP[i].d.Density;
			    }
			  SphP[i].u.s.DivVel = dvel[0][0] + dvel[1][1] + dvel[2][2];

			  SphP[i].u.s.StressDiag[0] = 2 * dvel[0][0] - 2.0 / 3 * SphP[i].u.s.DivVel;
			  SphP[i].u.s.StressDiag[1] = 2 * dvel[1][1] - 2.0 / 3 * SphP[i].u.s.DivVel;
			  SphP[i].u.s.StressDiag[2] = 2 * dvel[2][2] - 2.0 / 3 * SphP[i].u.s.DivVel;

			  SphP[i].u.s.StressOffDiag[0] = dvel[0][1] + dvel[1][0];	/* xy */
			  SphP[i].u.s.StressOffDiag[1] = dvel[0][2] + dvel[2][0];	/* xz */
			  SphP[i].u.s.StressOffDiag[2] = dvel[1][2] + dvel[2][1];	/* yz */

#ifdef NAVIERSTOKES_BULK
			  SphP[i].u.s.StressBulk = All.NavierStokes_BulkViscosity * SphP[i].u.s.DivVel;
#endif
			  rotx = dvel[1][2] - dvel[2][1];
			  roty = dvel[2][0] - dvel[0][2];
			  rotz = dvel[0][1] - dvel[1][0];
			  SphP[i].u.s.CurlVel = sqrt(rotx * rotx + roty * roty + rotz * rotz);
#endif

#ifdef CONDUCTION_SATURATION
			  SphP[i].GradEntr[0] /= SphP[i].d.Density;
			  SphP[i].GradEntr[1] /= SphP[i].d.Density;
			  SphP[i].GradEntr[2] /= SphP[i].d.Density;
#endif

#if defined(LT_MV_CHEMICALDIFFUSION)
#if !defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
			  SphP[i].v_rms.vel_rms /= P[i].n.NumNgb;
#endif
			  for(k = 0; k < LT_NMetP; k++)
			    {
			      for(cc = 0; cc < 3; cc++)
				SphP[i].GradAbnd[k][cc] /= SphP[i].d.Density;
			      if(SphP[i].AvgAbnd[k] /= SphP[i].d.Density < 1e-30)
				SphP[i].AvgAbnd[k] = 0;
			    }
#endif


#if defined(VS_TURB) || defined(AB_TURB)
			  SphP[i].SmoothedVel[0] /= SphP[i].d.Density;
			  SphP[i].SmoothedVel[1] /= SphP[i].d.Density;
			  SphP[i].SmoothedVel[2] /= SphP[i].d.Density;
#endif

#ifdef RADTRANSFER_FLUXLIMITER
			  for(k = 0; k < N_BINS; k++)
			    {
			      SphP[i].Grad_ngamma[0][k] /= SphP[i].d.Density;
			      SphP[i].Grad_ngamma[1][k] /= SphP[i].d.Density;
			      SphP[i].Grad_ngamma[2][k] /= SphP[i].d.Density;
			    }
#endif

#ifdef AB_SHOCK
			  /* Division with density drops out when normalizing */
			  double Shock_Abs1 = -1.0 / sqrt(SphP[i].Shock_N1[0] * SphP[i].Shock_N1[0] +
							  SphP[i].Shock_N1[1] * SphP[i].Shock_N1[1] +
							  SphP[i].Shock_N1[2] * SphP[i].Shock_N1[2]);
			  SphP[i].Shock_N1[0] *= Shock_Abs1;
			  SphP[i].Shock_N1[1] *= Shock_Abs1;
			  SphP[i].Shock_N1[2] *= Shock_Abs1;
#ifdef AB_SHOCK_VELDIV
			  double Shock_W2[3];
			  Shock_W2[0] = 0;
			  Shock_W2[1] = -SphP[i].Shock_N1[2];
			  Shock_W2[2] = SphP[i].Shock_N1[1];

			  double Shock_Scalar1 =
			    SphP[i].Shock_N1[0] * Shock_W2[0] + SphP[i].Shock_N1[1] * Shock_W2[1] +
			    SphP[i].Shock_N1[2] * Shock_W2[2];
			  SphP[i].Shock_N2[0] = Shock_W2[0] - Shock_Scalar1 * SphP[i].Shock_N1[0];
			  SphP[i].Shock_N2[1] = Shock_W2[1] - Shock_Scalar1 * SphP[i].Shock_N1[1];
			  SphP[i].Shock_N2[2] = Shock_W2[2] - Shock_Scalar1 * SphP[i].Shock_N1[2];
			  double Shock_Abs2 = 1.0 / sqrt(SphP[i].Shock_N2[0] * SphP[i].Shock_N2[0] +
							 SphP[i].Shock_N2[1] * SphP[i].Shock_N2[1] +
							 SphP[i].Shock_N2[2] * SphP[i].Shock_N2[2]);
			  SphP[i].Shock_N2[0] *= Shock_Abs2;
			  SphP[i].Shock_N2[1] *= Shock_Abs2;
			  SphP[i].Shock_N2[2] *= Shock_Abs2;

			  SphP[i].Shock_N3[0] =
			    SphP[i].Shock_N1[1] * SphP[i].Shock_N2[2] -
			    SphP[i].Shock_N1[2] * SphP[i].Shock_N2[1];
			  SphP[i].Shock_N3[1] =
			    SphP[i].Shock_N1[2] * SphP[i].Shock_N2[0] -
			    SphP[i].Shock_N1[0] * SphP[i].Shock_N2[2];
			  SphP[i].Shock_N3[2] =
			    SphP[i].Shock_N1[0] * SphP[i].Shock_N2[1] -
			    SphP[i].Shock_N1[1] * SphP[i].Shock_N2[0];
			  double Shock_Abs3 =
			    1.0 / sqrt(SphP[i].Shock_N3[0] * SphP[i].Shock_N3[0] +
				       SphP[i].Shock_N3[1] * SphP[i].Shock_N3[1] +
				       SphP[i].Shock_N3[2] * SphP[i].Shock_N3[2]);
			  SphP[i].Shock_N3[0] *= Shock_Abs3;
			  SphP[i].Shock_N3[1] *= Shock_Abs3;
			  SphP[i].Shock_N3[2] *= Shock_Abs3;
#endif
#endif /* AB_SHOCK */

#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB)
			  SphP[i].RotB[0] /= SphP[i].d.Density;
			  SphP[i].RotB[1] /= SphP[i].d.Density;
			  SphP[i].RotB[2] /= SphP[i].d.Density;
#endif
#ifdef TRACEDIVB
			  SphP[i].divB /= SphP[i].d.Density;
#endif
#ifdef FS_ALFA2_TURB
			  double rr[3], final_mom = 0.0, norm = 0.0;
			  final_mom = 0.0;
			  norm = 0.0;
			  for(k = 0; k < 3; k++)
			    rr[k] = (P[i].Pos[k] - SysState.DensestPart[k]);
			  for(k = 0; k < 3; k++)
			    final_mom += rr[k] * SysState.AngMomentumComp[0][k];
			  for(k = 0; k < 3; k++)
			    norm += pow(rr[k] * SysState.AngMomentumComp[0][k], 2);
			  norm = sqrt(norm);
#ifdef FS_TURB_ESTIM
			  if(P[i].TimeBin == 0)
			    norm = 1;
			  SphP[i].alfa2 = final_mom / norm * SphP[i].Lturb * SphP[i].Vturb / PPP[i].Hsml * 1E-1;	// ( 10. * SphP[i].Ne);
#endif
			  if(fabs(SphP[i].alfa2) < 1E-4)
			    SphP[i].alfa2 = SphP[i].alfa2 > 0 ? 1E-4 : -1E-4;
			  if(fabs(SphP[i].alfa2) > 2.0)
			    SphP[i].alfa2 = SphP[i].alfa2 > 0 ? 2. : -2.;
#endif
#ifdef FS_ALFA2_LOCAL_QUENCH
			  if(SphP[i].Alfa2 != 0)
			    SphP[i].alfa2 *= All.Alfa0Quench /
			      ((SphP[i].b2.BPred[0] * SphP[i].b2.BPred[0] +
				SphP[i].b2.BPred[1] * SphP[i].b2.BPred[1] +
				SphP[i].b2.BPred[2] * SphP[i].b2.BPred[2]) / All.AlfaB2eq + 1.);
#endif
#ifdef FS_ETA2_TURB
			  SphP[i].eta2 = 1. / 3. * PPP[i].Hsml * PPP[i].Hsml * SphP[i].r.CurlVel;
#ifdef FS_TURB_ESTIM
			  SphP[i].eta2 = 1. / 3. * SphP[i].Lturb * SphP[i].Vturb * 1E-2;
#endif
			  if(SphP[i].eta2 > 5.0)
			    SphP[i].eta2 = 5.0;
#endif
#ifdef JD_VTURB
			  SphP[i].Vbulk[0] /= P[i].TrueNGB;
			  SphP[i].Vbulk[1] /= P[i].TrueNGB;
			  SphP[i].Vbulk[2] /= P[i].TrueNGB;
#endif


#ifdef EULERPOTENTIALS
			  double bx = SphP[i].b1.dEulerA[1] * SphP[i].b2.dEulerB[2] -
			    SphP[i].b1.dEulerA[2] * SphP[i].b2.dEulerB[1];
			  double by = SphP[i].b1.dEulerA[2] * SphP[i].b2.dEulerB[0] -
			    SphP[i].b1.dEulerA[0] * SphP[i].b2.dEulerB[2];
			  double bz = SphP[i].b1.dEulerA[0] * SphP[i].b2.dEulerB[1] -
			    SphP[i].b1.dEulerA[1] * SphP[i].b2.dEulerB[0];
			  SphP[i].b2.BPred[0] = bx / (SphP[i].d.Density * SphP[i].d.Density);
			  SphP[i].b2.BPred[1] = by / (SphP[i].d.Density * SphP[i].d.Density);
			  SphP[i].b2.BPred[2] = bz / (SphP[i].d.Density * SphP[i].d.Density);
#endif

#ifdef HIGH_ORDER_INDUCTION
			  int ii;
			  detChi = SphP[i].Chi[0][0] * SphP[i].Chi[1][1] * SphP[i].Chi[2][2] +
			    SphP[i].Chi[0][1] * SphP[i].Chi[1][2] * SphP[i].Chi[2][0] +
			    SphP[i].Chi[0][2] * SphP[i].Chi[1][0] * SphP[i].Chi[2][1] -
			    SphP[i].Chi[0][2] * SphP[i].Chi[1][1] * SphP[i].Chi[2][0] -
			    SphP[i].Chi[0][1] * SphP[i].Chi[1][0] * SphP[i].Chi[2][2] -
			    SphP[i].Chi[0][0] * SphP[i].Chi[1][2] * SphP[i].Chi[2][1];
			  if(fabs(detChi) > 0)
			    {
			      Chi[0][0] =
				(SphP[i].Chi[1][1] * SphP[i].Chi[2][2] -
				 SphP[i].Chi[1][2] * SphP[i].Chi[2][1]) / detChi;
			      Chi[0][1] =
				(SphP[i].Chi[0][2] * SphP[i].Chi[2][1] -
				 SphP[i].Chi[0][1] * SphP[i].Chi[2][2]) / detChi;
			      Chi[0][2] =
				(SphP[i].Chi[0][1] * SphP[i].Chi[1][2] -
				 SphP[i].Chi[0][2] * SphP[i].Chi[1][1]) / detChi;
			      Chi[1][0] =
				(SphP[i].Chi[1][2] * SphP[i].Chi[2][0] -
				 SphP[i].Chi[1][0] * SphP[i].Chi[2][2]) / detChi;
			      Chi[1][1] =
				(SphP[i].Chi[0][0] * SphP[i].Chi[2][2] -
				 SphP[i].Chi[0][2] * SphP[i].Chi[2][0]) / detChi;
			      Chi[1][2] =
				(SphP[i].Chi[0][2] * SphP[i].Chi[1][0] -
				 SphP[i].Chi[0][0] * SphP[i].Chi[1][2]) / detChi;
			      Chi[2][0] =
				(SphP[i].Chi[1][0] * SphP[i].Chi[2][1] -
				 SphP[i].Chi[1][1] * SphP[i].Chi[2][0]) / detChi;
			      Chi[2][1] =
				(SphP[i].Chi[0][1] * SphP[i].Chi[2][0] -
				 SphP[i].Chi[0][0] * SphP[i].Chi[2][1]) / detChi;
			      Chi[2][2] =
				(SphP[i].Chi[0][0] * SphP[i].Chi[1][1] -
				 SphP[i].Chi[0][1] * SphP[i].Chi[1][0]) / detChi;

			      Xix[0] =
				Chi[0][0] * SphP[i].Xix[0] + Chi[1][0] * SphP[i].Xix[1] +
				Chi[2][0] * SphP[i].Xix[2];
			      Xix[1] =
				Chi[0][1] * SphP[i].Xix[0] + Chi[1][1] * SphP[i].Xix[1] +
				Chi[2][1] * SphP[i].Xix[2];
			      Xix[2] =
				Chi[0][2] * SphP[i].Xix[0] + Chi[1][2] * SphP[i].Xix[1] +
				Chi[2][2] * SphP[i].Xix[2];

			      Xiy[0] =
				Chi[0][0] * SphP[i].Xiy[0] + Chi[1][0] * SphP[i].Xiy[1] +
				Chi[2][0] * SphP[i].Xiy[2];
			      Xiy[1] =
				Chi[0][1] * SphP[i].Xiy[0] + Chi[1][1] * SphP[i].Xiy[1] +
				Chi[2][1] * SphP[i].Xiy[2];
			      Xiy[2] =
				Chi[0][2] * SphP[i].Xiy[0] + Chi[1][2] * SphP[i].Xiy[1] +
				Chi[2][2] * SphP[i].Xiy[2];

			      Xiz[0] =
				Chi[0][0] * SphP[i].Xiz[0] + Chi[1][0] * SphP[i].Xiz[1] +
				Chi[2][0] * SphP[i].Xiz[2];
			      Xiz[1] =
				Chi[0][1] * SphP[i].Xiz[0] + Chi[1][1] * SphP[i].Xiz[1] +
				Chi[2][1] * SphP[i].Xiz[2];
			      Xiz[2] =
				Chi[0][2] * SphP[i].Xiz[0] + Chi[1][2] * SphP[i].Xiz[1] +
				Chi[2][2] * SphP[i].Xiz[2];

			      for(ii = 0; ii < 3; ii++)
				{
				  SphP[i].Xix[ii] = Xix[ii] / P[i].TrueNGB;
				  SphP[i].Xiy[ii] = Xiy[ii] / P[i].TrueNGB;
				  SphP[i].Xiz[ii] = Xiz[ii] / P[i].TrueNGB;
				}
			    }
			  else
			    {
			      for(ii = 0; ii < 3; ii++)
				{
				  SphP[i].Xix[ii] = 0.0;
				  SphP[i].Xiy[ii] = 0.0;
				  SphP[i].Xiz[ii] = 0.0;
				}
			    }
#endif

			}

		      SphP[i].Pressure = get_pressure(i);


#if defined(LT_STELLAREVOLUTION)
		    }
#endif
		}

#ifdef BLACK_HOLES
	      if(P[i].Type == 5)
		{
		  if(BPP(i).b1.BH_Density > 0)
		    {
		      BPP(i).b2.BH_Entropy /= BPP(i).b1.BH_Density;
		      BPP(i).b3.BH_SurroundingGasVel[0] /= BPP(i).b1.BH_Density;
		      BPP(i).b3.BH_SurroundingGasVel[1] /= BPP(i).b1.BH_Density;
		      BPP(i).b3.BH_SurroundingGasVel[2] /= BPP(i).b1.BH_Density;
		    }
#ifdef KD_FRICTION
		  if(BPP(i).BH_SurroundingDensity > 0)
		    {
		      BPP(i).BH_SurroundingVel[0] /= BPP(i).BH_SurroundingDensity;
		      BPP(i).BH_SurroundingVel[1] /= BPP(i).BH_SurroundingDensity;
		      BPP(i).BH_SurroundingVel[2] /= BPP(i).BH_SurroundingDensity;
		    }
#endif
#ifdef LT_BH_GUESSHSML
		  if(BPP(i).mean_rho > 0)
		    BPP(i).mean_hsml /= BPP(i).mean_rho;
#endif
		}
#endif

#ifdef LT_STARS_GUESSHSML
#ifdef LT_STELLAREVOLUTION
	      if(P[i].Type == 4)
		if(MetP[P[i].pt.MetID].mean_rho > 0)
		  MetP[P[i].pt.MetID].mean_hsml /= MetP[P[i].pt.MetID].mean_rho;
#endif
#endif

	      /* now check whether we had enough neighbours */

	      desnumngb = All.DesNumNgb;
	      desnumngbdev = All.MaxNumNgbDeviation;

#ifdef BLACK_HOLES
	      if(P[i].Type == 5)
		{
		  desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
#ifdef LT_STELLAREVOLUTION
		  desnumngbdev = All.SpreadNumNgbDev * All.BlackHoleNgbFactor;
#else
		  desnumngbdev = 4 * All.BlackHoleNgbFactor;
#endif
		}
#endif

#if defined(RADTRANSFER) && defined(EDDINGTON_TENSOR_STARS)
	      if(P[i].Type == 4)
		desnumngb = 64;	//NORM_COEFF * KERNEL_COEFF_1;   /* will assign the stellar luminosity to very few (one actually) gas particles */
#endif

#if defined(LT_STELLAREVOLUTION)
	      if(P[i].Type == 4)
		{
		  desnumngb = All.DesNumNgbSN;
		  desnumngbdev = All.SpreadNumNgbDev;
		}
#endif

	      redo_particle = 0;

	      if(PPP[i].n.NumNgb < (desnumngb - desnumngbdev) ||
		 (PPP[i].n.NumNgb > (desnumngb + desnumngbdev) && PPP[i].Hsml > (1.01 * All.MinGasHsml)))
		redo_particle = 1;

#ifdef KD_RESTRICT_NEIGHBOURS
	      if(P[i].TrueNGB > desnumngb * 3 || P[i].TrueNGB < desnumngb * 0.4)
		redo_particle = 1;
#endif

#ifdef MAXHSML
	      if(PPP[i].n.NumNgb < (desnumngb - desnumngbdev) && PPP[i].Hsml >= All.MaxHsml)
		{
		  PPP[i].Hsml = All.MaxHsml;
		  redo_particle = 0;
		}
#endif

	      if(redo_particle)
		{
		  /* need to redo this particle */
		  npleft++;

#ifdef LT_STARS_GUESSHSML
		  if(P[i].Type == 4)
		    if(Right[i] == 0 && Left[i] == 0)
		      if(MetP[P[i].pt.MetID].mean_hsml > 0)
			{
			  if(PPP[i].n.NumNgb < (desnumngb - desnumngbdev)
			     && MetP[P[i].pt.MetID].mean_hsml > PPP[i].Hsml)
			    {
			      Left[i] = PPP[i].Hsml;
			      PPP[i].Hsml = MetP[P[i].pt.MetID].mean_hsml;
			      if(PPP[i].Hsml > All.MaxChemSpreadL)
				PPP[i].Hsml = Left[i] = Right[i] = All.MaxChemSpreadL;
			      continue;
			    }
			  if(PPP[i].n.NumNgb > (desnumngb + desnumngbdev)
			     && MetP[P[i].pt.MetID].mean_hsml < PPP[i].Hsml)
			    {
			      Right[i] = PPP[i].Hsml;
			      PPP[i].Hsml = MetP[P[i].pt.MetID].mean_hsml;
			      continue;
			    }
			}
#endif

#ifdef LT_BH_GUESSHSML
		  if(P[i].Type == 5)
		    if(Right[i] == 0 && Left[i] == 0)
		      if(BPP(i).mean_hsml > 0)
			{
			  if(PPP[i].n.NumNgb < (desnumngb - desnumngbdev) && BPP(i).mean_hsml > PPP[i].Hsml)
			    {
			      Left[i] = PPP[i].Hsml;
			      PPP[i].Hsml = BPP(i).mean_hsml;
			      if(PPP[i].Hsml > All.BlackHoleMaxAccretionRadius)
				PPP[i].Hsml = Left[i] = Right[i] = All.BlackHoleMaxAccretionRadius;
			      continue;
			    }
			  if(PPP[i].n.NumNgb > (desnumngb + desnumngbdev) && BPP(i).mean_hsml < PPP[i].Hsml)
			    {
			      Right[i] = PPP[i].Hsml;
			      PPP[i].Hsml = BPP(i).mean_hsml;
			      continue;
			    }
			}
#endif

		  if(Left[i] > 0 && Right[i] > 0)
		    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
		      {
			/* this one should be ok */
			npleft--;
			P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			continue;
		      }
#ifdef KD_RESTRICT_NEIGHBOURS
		  if(PPP[i].n.NumNgb < (desnumngb - desnumngbdev) && P[i].TrueNGB > desnumngb * 2)
		    {
		      if(P[i].TrueNGB > desnumngb * 3)
			{
			  /* Try to do better */
			  if(Right[i] > 0)
			    Right[i] = DMIN(PPP[i].Hsml, Right[i]);
			  else
			    Right[i] = PPP[i].Hsml;
			  if(Left[i] > 0)
			    PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
			  else
			    PPP[i].Hsml *= 0.9;
			}
		      else
			{
			  /* Stop here ! */
			  npleft--;
			  P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			  continue;
			}
		    }

		  if((PPP[i].n.NumNgb < (desnumngb + desnumngbdev)) &&
		     (PPP[i].n.NumNgb > (desnumngb - desnumngbdev)) && P[i].TrueNGB > desnumngb * 3)
		    {
		      /* Try to do better */
		      if(Right[i] > 0)
			Right[i] = DMIN(PPP[i].Hsml, Right[i]);
		      else
			Right[i] = PPP[i].Hsml;
		      if(Left[i] > 0)
			PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		      else
			PPP[i].Hsml *= 0.9;
		    }

		  if(PPP[i].n.NumNgb > (desnumngb + desnumngbdev) && P[i].TrueNGB < desnumngb * 0.5)
		    {
		      if(P[i].TrueNGB < desnumngb * 0.4)
			{
			  /* Try to do better */
			  if(Left[i] > 0)
			    Left[i] = DMAX(PPP[i].Hsml, Left[i]);
			  else
			    Left[i] = PPP[i].Hsml;
			  if(Right[i] > 0)
			    PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
			  else
			    PPP[i].Hsml *= 1.1;
			}
		      else
			{
			  /* Stop here ! */
			  npleft--;
			  P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			  continue;
			}
		    }
		  if((PPP[i].n.NumNgb < (desnumngb + desnumngbdev)) &&
		     (PPP[i].n.NumNgb > (desnumngb - desnumngbdev)) && P[i].TrueNGB < desnumngb * 0.4)
		    {
		      /* Try to do better */
		      if(Left[i] > 0)
			Left[i] = DMAX(PPP[i].Hsml, Left[i]);
		      else
			Left[i] = PPP[i].Hsml;
		      if(Right[i] > 0)
			PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		      else
			PPP[i].Hsml *= 1.1;
		    }
#endif

		  if(PPP[i].n.NumNgb < (desnumngb - desnumngbdev))
		    Left[i] = DMAX(PPP[i].Hsml, Left[i]);
		  else
		    {
		      if(Right[i] != 0)
			{
			  if(PPP[i].Hsml < Right[i])
			    Right[i] = PPP[i].Hsml;
			}
		      else
			Right[i] = PPP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
#ifdef KD_RESTRICT_NEIGHBOURS
		      printf
			("i=%d task=%d ID=%llu Hsml=%g Left=%g Right=%g Ngbs=%g Tngbs=%d Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (unsigned long long) P[i].ID, PPP[i].Hsml, Left[i], Right[i],
			 (float) PPP[i].n.NumNgb, P[i].TrueNGB, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
#else
		      printf
			("i=%d task=%d ID=%llu Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (unsigned long long) P[i].ID, PPP[i].Hsml, Left[i], Right[i],
			 (float) PPP[i].n.NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
#endif
		      fflush(stdout);
		    }

		  if(Right[i] > 0 && Left[i] > 0)
		    PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		  else
		    {
		      if(Right[i] == 0 && Left[i] == 0)
			{
			  char buf[1000];
			  sprintf(buf, "Right[i] == 0 && Left[i] == 0 PPP[i].Hsml=%g\n", PPP[i].Hsml);
			  terminate(buf);
			}

		      if(Right[i] == 0 && Left[i] > 0)
			{
			  if(P[i].Type == 0 && fabs(PPP[i].n.NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (PPP[i].n.NumNgb -
					 desnumngb) / (NUMDIMS * PPP[i].n.NumNgb) *
				SphP[i].h.DhsmlDensityFactor;

			      if(fac < 1.26)
				PPP[i].Hsml *= fac;
			      else
				PPP[i].Hsml *= 1.26;
			    }
			  else
			    PPP[i].Hsml *= 1.26;
			}

		      if(Right[i] > 0 && Left[i] == 0)
			{
			  if(P[i].Type == 0 && fabs(PPP[i].n.NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (PPP[i].n.NumNgb -
					 desnumngb) / (NUMDIMS * PPP[i].n.NumNgb) *
				SphP[i].h.DhsmlDensityFactor;
#ifdef KD_RESTRICT_NEIGHBOURS
			      if((PPP[i].n.NumNgb < (desnumngb + desnumngbdev)) &&
				 (PPP[i].n.NumNgb > (desnumngb - desnumngbdev)) &&
				 P[i].TrueNGB > desnumngb * 3)
				fac = 1 / 1.26;
#endif

			      if(fac > 1 / 1.26)
				PPP[i].Hsml *= fac;
			      else
				PPP[i].Hsml /= 1.26;
			    }
			  else
			    PPP[i].Hsml /= 1.26;
			}
		    }

#if defined(LT_STELLAREVOLUTION)
		  if(P[i].Type != 4)
#endif
		    if(PPP[i].Hsml < All.MinGasHsml)
		      PPP[i].Hsml = All.MinGasHsml;

#ifdef BLACK_HOLES
		  if(P[i].Type == 5)
		    if(Left[i] > All.BlackHoleMaxAccretionRadius)
		      {
			/* this will stop the search for a new BH smoothing length in the next iteration */
			PPP[i].Hsml = Left[i] = Right[i] = All.BlackHoleMaxAccretionRadius;
		      }
#endif

#if defined(LT_STELLAREVOLUTION)
		  if(P[i].Type == 4)
		    {
		      /*                       if(PPP[i].Hsml < All.MinChemSpreadL) */
		      /*                         PPP[i].Hsml = Left[i] = Right[i] = All.MinChemSpreadL; */
		      if(PPP[i].Hsml > All.MaxChemSpreadL)
			{
			  PPP[i].Hsml = Left[i] = Right[i] = All.MaxChemSpreadL;
			  /*
			     #ifndef LT_STRAS_GUESSHSML
			     if(PPP[i].n.NumNgb == 0)
			     {
			     printf
			     ("Task %d : particle %d (ID %llu) hits Maximum spreading length = %g with 0 Neighbours!\n",
			     ThisTask, i, (unsigned long long) P[i].ID, All.MaxChemSpreadL);
			     endrun(778899);
			     }
			     #endif
			   */
			  hitMaxChemSpreadL++;
			}
		      if(PPP[i].Hsml < All.MinChemSpreadL)
			PPP[i].Hsml = Left[i] = Right[i] = All.MinChemSpreadL;
		    }
#endif
		}
	      else
		P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
    }
  while(ntot > 0);

#if defined(LT_STELLAREVOLUTION)
  unsigned int tot_hitMaxChemSpreadL = 0;

  MPI_Reduce(&hitMaxChemSpreadL, &tot_hitMaxChemSpreadL, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0 && tot_hitMaxChemSpreadL > 0)
    printf("%u particles have hit on the maximum allowed spreading hsml (%g)\n", tot_hitMaxChemSpreadL,
	   All.MaxChemSpreadL);
#endif

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Right);
  myfree(Left);
  myfree(Ngblist);

#if defined(CS_MODEL) && defined(CS_FEEDBACK)
  double xhyd, yhel, ne, mu, energy, temp, a3inv;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0 && P[i].EnergySN < 0)	/* only gas particles should enter here */
	{

	  /* printf("prom 2i=%d type=%d energySN%g\n", i, P[i].Type, P[i].EnergySN);
	     fflush(0); */

	  P[i].EnergySN = 0;

	}

      if(All.ComovingIntegrationOn)
	a3inv = 1 / (All.Time * All.Time * All.Time);
      else
	a3inv = 1;

      if(P[i].Type == 0 && (SphP[i].TempPromotion > 0 || SphP[i].DensPromotion > 0))
	{
	  xhyd = P[i].Zm[6] / P[i].Mass;
	  yhel = (1 - xhyd) / (4. * xhyd);
	  ne = SphP[i].Ne;
	  mu = (1 + 4 * yhel) / (1 + yhel + ne);
	  energy = SphP[i].Entropy * P[i].Mass / GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);	/* Total Energy */
	  temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[i].Mass * PROTONMASS * mu;
	  temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */
#ifndef LONGIDS
	  fprintf(FdPromotion, "%g %u %g %g %g %g %g %g\n", All.Time, P[i].ID, SphP[i].DensPromotion,
		  SphP[i].TempPromotion, SphP[i].d.Density, temp, SphP[i].da.DensityAvg,
		  SphP[i].ea.EntropyAvg);
#else
	  fprintf(FdPromotion, "%g %llu %g %g %g %g %g %g\n", All.Time, P[i].ID, SphP[i].DensPromotion,
		  SphP[i].TempPromotion, SphP[i].d.Density, temp, SphP[i].da.DensityAvg,
		  SphP[i].ea.EntropyAvg);
#endif
	  fflush(FdPromotion);
	}
    }
#endif

  /* mark as active again */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].TimeBin < 0)
	P[i].TimeBin = -P[i].TimeBin - 1;
    }

#if defined(LT_STELLAREVOLUTION)
#ifndef LT_LOCAL_IRA
  if(N_stars > 0)
#endif
    drop_chemicallyactive_particles();
#endif

#if defined(LT_MV_CHEMICALDIFFUSION) && defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
  if(All.NumCurrentTiStep == 0 && All.ConstantD == 0)
    {
      double Avg[4], AvgSum[4];
      Avg[0] = Avg[1] = Avg[2] = 0;
      Avg[3] = (double) N_gas;
      AvgSum[0] = AvgSum[1] = AvgSum[2] = 0;
      for(i = 0; i < N_gas; i++)
	{
	  Avg[0] += SphP[i].d.Density;	/* coefficient D  */
	  Avg[1] += pow(P[i].Mass / SphP[i].d.Density, 0.33333333);	/* hsml */
	  Avg[2] += sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density);	/* sound velocity */
	}

      MPI_Allreduce(Avg, AvgSum, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      All.ConstantD = (AvgSum[1] / AvgSum[3] * AvgSum[2] / AvgSum[3]);
#if !defined(LT_MV_CHEMICALDIFFUSION_NO_DENSITY_IN_D)
      All.ConstantD *= AvgSum[0] / AvgSum[3];
#endif
      if(ThisTask == 0)
	{
	  printf("The constant Diffusion Coefficient is %g\n", All.ConstantD);
	  fflush(stdout);
	  printf("Density avg %g Smooth Lenght avg %g Sound velocity avg %g \n", AvgSum[0] / AvgSum[3],
		 AvgSum[1] / AvgSum[3], AvgSum[2] / AvgSum[3]);
	  fflush(stdout);
	}
    }
#endif


  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_DENSCOMPUTE] += timecomp;
  CPU_Step[CPU_DENSWAIT] += timewait;
  CPU_Step[CPU_DENSCOMM] += timecomm;
  CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		     int *ngblist)
{
  int j, n;
  int startnode, numngb_inbox, listindex = 0;
  double r2, h2, u, mass_j;

  struct kernel_density kernel;
  struct densdata_in local;
  struct densdata_out out;
  memset(&out, 0, sizeof(struct densdata_out));

#ifdef BLACK_HOLES
  out.BH_TimeBinGasNeighbor = TIMEBINS;
#endif

  if(mode == 0)
    particle2in_density(&local, target);
  else
    local = DensDataGet[target];

  h2 = local.Hsml * local.Hsml;

  kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifdef CS_MODEL
	  /* CECILIA CHECK: shouldn't we use with or without threads according to input options? */
	  numngb_inbox =
	    cs_ngb_treefind_variable_decoupling_threads(&local.Pos[0], local.Hsml, target, &startnode,
							local.densityOld, local.entropy, &local.Vel[0], mode,
							exportflag, exportnodecount, exportindex, ngblist);
#else
#ifndef KD_FRICTION
	  numngb_inbox =
	    ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
					  exportnodecount, exportindex, ngblist);
#else
	  numngb_inbox =
	    ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
					  exportnodecount, exportindex, ngblist, local.Type);
#endif
#endif

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = ngblist[n];
#ifdef WINDS
#ifdef KD_FRICTION
	      if(P[j].Type == 0)
#endif
		if(SphP[j].DelayTime > 0)	/* partner is a wind particle */
		  if(!(local.DelayTime > 0))	/* if I'm not wind, then ignore the wind particle */
		    continue;
#endif
#if defined(BLACK_HOLES) || defined(CA_BH_ACCRETION)
	      if(P[j].Mass == 0)
		continue;
#endif

	      kernel.dx = local.Pos[0] - P[j].Pos[0];
	      kernel.dy = local.Pos[1] - P[j].Pos[1];
	      kernel.dz = local.Pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      kernel.dx = NEAREST_X(kernel.dx);
	      kernel.dy = NEAREST_Y(kernel.dy);
	      kernel.dz = NEAREST_Z(kernel.dz);
#endif
	      r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;

	      if(r2 < h2)
		{
#if defined(JD_VTURB) || defined(HIGH_ORDER_INDUCTION)
		  out.TrueNGB++;
#endif
		  kernel.r = sqrt(r2);
		  u = kernel.r * kernel.hinv;

		  kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);

		  mass_j = P[j].Mass;
#ifdef KD_FRICTION
		  if(P[j].Type == 0)
		    {
#endif
		      kernel.mj_wk = FLT(mass_j * kernel.wk);

		      out.Rho += kernel.mj_wk;

		      out.Ngb += FLT(NORM_COEFF * kernel.wk / kernel.hinv3);	/* 4.0/3 * PI = 4.188790204786 */

		      out.DhsmlDensity += FLT(-mass_j * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk));

		      if(kernel.r > 0)
			{
			  kernel.mj_dwk_r = mass_j * kernel.dwk / kernel.r;

			  kernel.dvx = local.Vel[0] - SphP[j].VelPred[0];
			  kernel.dvy = local.Vel[1] - SphP[j].VelPred[1];
			  kernel.dvz = local.Vel[2] - SphP[j].VelPred[2];
#ifdef FS_SHEAR_BOX
			  kernel.dvz =
			    (local.Pos[0] - P[j].Pos[0]) > boxHalf_X ?
			    kernel.dvz + FS_SHEAR_Q * FS_SHEAR_OMEGA * boxSize_X :
			    (local.Pos[0] - P[j].Pos[0]) < -boxHalf_X ?
			    kernel.dvz - FS_SHEAR_Q * FS_SHEAR_OMEGA * boxSize_X : kernel.dvz;
#endif

#ifndef NAVIERSTOKES
			  out.Div +=
			    FLT(-kernel.mj_dwk_r *
				(kernel.dx * kernel.dvx + kernel.dy * kernel.dvy + kernel.dz * kernel.dvz));

			  out.Rot[0] +=
			    FLT(kernel.mj_dwk_r * (kernel.dz * kernel.dvy - kernel.dy * kernel.dvz));
			  out.Rot[1] +=
			    FLT(kernel.mj_dwk_r * (kernel.dx * kernel.dvz - kernel.dz * kernel.dvx));
			  out.Rot[2] +=
			    FLT(kernel.mj_dwk_r * (kernel.dy * kernel.dvx - kernel.dx * kernel.dvy));
#endif
			}

		      density_evaluate_extra_physics_gas(&local, &out, &kernel, j);

#ifdef KD_FRICTION
		    }
		  else
		    {
		      if(P[j].Type == 1 || P[j].Type == 4)
			{
			  out.SurroundingDensity += kernel.mj_wk;
			  out.SurroundingVel[0] += kernel.mj_wk * P[j].Vel[0];
			  out.SurroundingVel[1] += kernel.mj_wk * P[j].Vel[1];
			  out.SurroundingVel[2] += kernel.mj_wk * P[j].Vel[2];
			}
		    }
#endif
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    out2particle_density(&out, target, 0);
  else
    DensDataResult[target] = out;

  return 0;
}


void *density_evaluate_primary(void *p)
{
  int thread_id = *(int *) p;
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	if(BufferFullFlag != 0 || NextParticle < 0)
	  {
	    exitFlag = 1;
	  }
	else
	  {
	    i = NextParticle;
	    ProcessedFlag[i] = 0;
	    NextParticle = NextActiveParticle[NextParticle];
	  }
      }
      UNLOCK_NEXPORT;
      if(exitFlag)
	break;

      if(density_isactive(i))
	{
	  if(density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *density_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;


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

      density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}


int density_isactive(int n)
{
  if(P[n].TimeBin < 0)
    return 0;

#if (defined(RADTRANSFER) && defined(EDDINGTON_TENSOR_STARS))|| defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION)
  if(P[n].Type == 4)
    return 1;
#endif

#ifdef BLACK_HOLES
  if(P[n].Type == 5)
    return 1;
#endif

#ifdef BLACK_HOLES
  if(P[n].Type == 0 && P[n].Mass > 0)
#else
  if(P[n].Type == 0)
#endif
    return 1;

  return 0;
}


void density_evaluate_extra_physics_gas(struct densdata_in *local, struct densdata_out *out,
					struct kernel_density *kernel, int j)
{
#ifdef JD_VTURB
  out->Vbulk[0] += SphP[j].VelPred[0];
  out->Vbulk[1] += SphP[j].VelPred[1];
  out->Vbulk[2] += SphP[j].VelPred[2];
#endif

#ifdef BLACK_HOLES
  out->SmoothedEntr += FLT(kernel->mj_wk * SphP[j].Entropy);
  if(out->BH_TimeBinGasNeighbor > P[j].TimeBin)
    out->BH_TimeBinGasNeighbor = P[j].TimeBin;
#endif
#if defined(BLACK_HOLES) || defined(VS_TURB) || defined(AB_TURB)
  out->GasVel[0] += FLT(kernel->mj_wk * SphP[j].VelPred[0]);
  out->GasVel[1] += FLT(kernel->mj_wk * SphP[j].VelPred[1]);
  out->GasVel[2] += FLT(kernel->mj_wk * SphP[j].VelPred[2]);
#endif

#if defined(LT_BH_GUESSHSML) || defined(LT_STARS_GUESSHSML)
  out->SmoothedHsml += kernel->mj_wk * PPP[j].Hsml;
  out->SmoothedRho += kernel->mj_wk;
#endif

  if(kernel->r > 0)
    {
#ifdef CONDUCTION_SATURATION
      out->GradEntr[0] += kernel->mj_dwk_r * kernel->dx * SphP[j].Entropy;
      out->GradEntr[1] += kernel->mj_dwk_r * kernel->dy * SphP[j].Entropy;
      out->GradEntr[2] += kernel->mj_dwk_r * kernel->dz * SphP[j].Entropy;
#endif

#ifdef RADTRANSFER_FLUXLIMITER
      int k;
      for(k = 0; k < N_BINS; k++)
	{
	  out->Grad_ngamma[0][k] += kernel->mj_dwk_r * kernel->dx * SphP[j].n_gamma[k];
	  out->Grad_ngamma[1][k] += kernel->mj_dwk_r * kernel->dy * SphP[j].n_gamma[k];
	  out->Grad_ngamma[2][k] += kernel->mj_dwk_r * kernel->dz * SphP[j].n_gamma[k];
	}
#endif

#ifdef NAVIERSTOKES
      out->dvel[0][0] -= kernel->mj_dwk_r * kernel->dx * kernel->dvx;
      out->dvel[0][1] -= kernel->mj_dwk_r * kernel->dx * kernel->dvy;
      out->dvel[0][2] -= kernel->mj_dwk_r * kernel->dx * kernel->dvz;
      out->dvel[1][0] -= kernel->mj_dwk_r * kernel->dy * kernel->dvx;
      out->dvel[1][1] -= kernel->mj_dwk_r * kernel->dy * kernel->dvy;
      out->dvel[1][2] -= kernel->mj_dwk_r * kernel->dy * kernel->dvz;
      out->dvel[2][0] -= kernel->mj_dwk_r * kernel->dz * kernel->dvx;
      out->dvel[2][1] -= kernel->mj_dwk_r * kernel->dz * kernel->dvy;
      out->dvel[2][2] -= kernel->mj_dwk_r * kernel->dz * kernel->dvz;
#endif

#ifdef AB_SHOCK
      out->Shock_N1[0] += kernel->mj_dwk_r * kernel->dx * (SphP[j].Pressure - local->Pressure);
      out->Shock_N1[1] += kernel->mj_dwk_r * kernel->dy * (SphP[j].Pressure - local->Pressure);
      out->Shock_N1[2] += kernel->mj_dwk_r * kernel->dz * (SphP[j].Pressure - local->Pressure);
#endif

#if defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB) || defined(OUTPUT_ROTB)
      double bpred_j[3];
#ifdef SFR
      double dbx = local->BPred[0] - SphP[j].b2.BPred[0] * pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
      double dby = local->BPred[1] - SphP[j].b2.BPred[1] * pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
      double dbz = local->BPred[2] - SphP[j].b2.BPred[2] * pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
      bpred_j[0] = SphP[j].b2.BPred[0] * pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
      bpred_j[1] = SphP[j].b2.BPred[1] * pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
      bpred_j[2] = SphP[j].b2.BPred[2] * pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
#else
      double dbx = local->BPred[0] - SphP[j].b2.BPred[0];
      double dby = local->BPred[1] - SphP[j].b2.BPred[1];
      double dbz = local->BPred[2] - SphP[j].b2.BPred[2];
      bpred_j[0] = SphP[j].b2.BPred[0];
      bpred_j[1] = SphP[j].b2.BPred[1];
      bpred_j[2] = SphP[j].b2.BPred[2];
#endif
#endif

#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB)
      out->RotB[0] += FLT(kernel->mj_dwk_r * (kernel->dz * dby - kernel->dy * dbz));
      out->RotB[1] += FLT(kernel->mj_dwk_r * (kernel->dx * dbz - kernel->dz * dbx));
      out->RotB[2] += FLT(kernel->mj_dwk_r * (kernel->dy * dbx - kernel->dx * dby));
#endif
#ifdef HIGH_ORDER_INDUCTION
      out->Chi[0][0] += kernel->mj_dwk_r * kernel->dx * kernel->dx;
      out->Chi[0][1] += kernel->mj_dwk_r * kernel->dy * kernel->dx;
      out->Chi[0][2] += kernel->mj_dwk_r * kernel->dz * kernel->dx;
      out->Chi[1][0] += kernel->mj_dwk_r * kernel->dx * kernel->dy;
      out->Chi[1][1] += kernel->mj_dwk_r * kernel->dy * kernel->dy;
      out->Chi[1][2] += kernel->mj_dwk_r * kernel->dz * kernel->dy;
      out->Chi[2][0] += kernel->mj_dwk_r * kernel->dx * kernel->dz;
      out->Chi[2][1] += kernel->mj_dwk_r * kernel->dy * kernel->dz;
      out->Chi[2][2] += kernel->mj_dwk_r * kernel->dz * kernel->dz;

      out->Xix[0] += kernel->mj_dwk_r * kernel->dvx * kernel->dx;
      out->Xiy[0] += kernel->mj_dwk_r * kernel->dvy * kernel->dx;
      out->Xiz[0] += kernel->mj_dwk_r * kernel->dvz * kernel->dx;
      out->Xix[1] += kernel->mj_dwk_r * kernel->dvx * kernel->dy;
      out->Xiy[1] += kernel->mj_dwk_r * kernel->dvy * kernel->dy;
      out->Xiz[1] += kernel->mj_dwk_r * kernel->dvz * kernel->dy;
      out->Xix[2] += kernel->mj_dwk_r * kernel->dvx * kernel->dz;
      out->Xiy[2] += kernel->mj_dwk_r * kernel->dvy * kernel->dz;
      out->Xiz[2] += kernel->mj_dwk_r * kernel->dvz * kernel->dz;
#endif

#ifdef TRACEDIVB
#ifdef RELDIVB
      double babs_i =
	sqrt(local->BPred[0] * local->BPred[0] + local->BPred[1] * local->BPred[1] +
	     local->BPred[2] * local->BPred[2]);
      double babs_j = sqrt(bpred_j[0] * bpred_j[0] + bpred_j[1] * bpred_j[1] + bpred_j[2] * bpred_j[2]);
      double divbfac = (local->Hsml + PPP[j].Hsml) / (babs_i + babs_j) * babs_i / local->Hsml;
      out->divB +=
	FLT(-divbfac * kernel->mj_dwk_r * (dbx * kernel->dx + dby * kernel->dy + dbz * kernel->dz));
#else
      out->divB += FLT(-kernel->mj_dwk_r * (dbx * kernel->dx + dby * kernel->dy + dbz * kernel->dz));
#endif
#endif

#ifdef EULERPOTENTIALS
      double dea = kernel->mj_dwk_r * (local->EulerA - SphP[j].EulerA);
      double deb = kernel->mj_dwk_r * (local->EulerB - SphP[j].EulerB);
#ifdef EULER_VORTEX
      deb = NEAREST_Z(deb);
#endif
      out->dEulerA[0] -= kernel->dx * dea;
      out->dEulerA[1] -= kernel->dy * dea;
      out->dEulerA[2] -= kernel->dz * dea;
      out->dEulerB[0] -= kernel->dx * deb;
      out->dEulerB[1] -= kernel->dy * deb;
      out->dEulerB[2] -= kernel->dz * deb;
#endif
    }

#if defined(LT_MV_CHEMICALDIFFUSION)
  int el;

  for(el = 0; el < LT_NMetP; el++)
    {
      out->GradAbnd[el][0] += kernel->mj_dwk_r * kernel->dx * SphP[j].Abundance[el];
      out->GradAbnd[el][1] += kernel->mj_dwk_r * kernel->dy * SphP[j].Abundance[el];
      out->GradAbnd[el][2] += kernel->mj_dwk_r * kernel->dz * SphP[j].Abundance[el];
      out->AvgAbnd[el] += kernel->mj_wk * SphP[j].Abundance[el];
    }
#if !defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
  out->v_rms += (SphP[j].VelPred[0] - local->Vel[0]) * (SphP[j].VelPred[0] - local->Vel[0]) +
    (SphP[j].VelPred[1] - local->Vel[1]) * (SphP[j].VelPred[1] - local->Vel[1]) +
    (SphP[j].VelPred[2] - local->Vel[2]) * (SphP[j].VelPred[2] - local->Vel[2]);
#endif
#endif
}

#ifdef NAVIERSTOKES
double get_shear_viscosity(int i)
{
  return All.NavierStokes_ShearViscosity;
}
#endif

#if defined(LT_STELLAREVOLUTION)
int append_chemicallyactive_particles(unsigned int *already_active)
{
  int i, n, prev, appended_particles, alreadyactive_particles;

  appended_particles = alreadyactive_particles = 0;

  LastActive = prev = FirstActiveParticle;
  while(LastActive >= 0)
    {
      prev = LastActive;
      LastActive = NextActiveParticle[LastActive];
    }

  LastActive = prev;

  for(i = 0; i < NumPart; i++)
    {
#ifdef BLACK_HOLES
      if(P[i].Mass == 0)
	continue;
#endif
      n = P[i].TimeBin;
      if(TimeBinActive[n])
	{
	  if((P[i].Type & 15) == 4)
	    alreadyactive_particles++;
	  continue;
	}

#ifdef LT_LOCAL_IRA
      if(P[i].Type == 0)
	n = SphP[i].ChemTimeBin;
      else if(P[i].Type == 4)
#else
      if((P[i].Type & 15) == 4)
#endif
	n = MetP[P[i].pt.MetID].ChemTimeBin;

      if(TimeBinActive[n])
	{
	  if(is_chemically_active(i) == 0)
#ifndef LONGIDS
	    printf("\n\t task %d @%d has found a mismatch A for particle %d, ID = %u \n",
		   ThisTask, All.NumCurrentTiStep, i, P[i].ID);
#else
	    printf("\n\t task %d @%d has found a mismatch A for particle %d, ID = %llu \n",
		   ThisTask, All.NumCurrentTiStep, i, P[i].ID);
#endif
	  /*             printf("\n\t task %d @%d has found a mismatch A for particle %d, ID = %llu :: %g %g %d %g\n", */
	  /*                    ThisTask, All.NumCurrentTiStep, i, (unsigned long long)P[i].ID, All.Time_Age, MetP[P[i].pt.MetID].NextChemTime, n, get_age(P[i].StellarAge) - All.Time_Age, get_age(P[i].StellarAge)-MetP[P[i].pt.MetID].NextChemTime); fflush(stdout); */

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
	  ActiveParticleList[NActivePart + appended_particles] = i;
#endif
	  appended_particles++;

	  if(prev == -1)
	    FirstActiveParticle = i;

	  if(prev >= 0)
	    NextActiveParticle[prev] = i;

	  prev = i;
	}
    }

  if(prev >= 0)
    NextActiveParticle[prev] = -1;

  /* KD: There would be a more optimized way to do it !! */
#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  NActivePartMerk = NActivePart;
  NActivePart += appended_particles;
#endif

  *already_active = alreadyactive_particles;

  return appended_particles;
}


void drop_chemicallyactive_particles(void)
{
  if(LastActive >= 0)
    NextActiveParticle[LastActive] = -1;
  else
    FirstActiveParticle = -1;

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  NActivePart = NActivePartMerk;
#endif

}
#endif

#ifdef ADAPTGRAVSOFT
#include"AdaptGravSoft/ads_density.c"
#endif


#ifdef LT_MV_CHEMICALDIFFUSION

int updatediffusion_and_normalize()
/* perform the update of diffusion
     accumulate the resulting total mass of metals
*/
{
  double hubble_a;
  double expfact, dt;
  double prevMetalsAmount[LT_NMet], prevMetalsAmount_tot[LT_NMet];
  double newMetalsAmount[LT_NMet], newMetalsAmount_tot[LT_NMet];
  double factor[LT_NMet];

  int tistep, i, j;

  if(All.ComovingIntegrationOn)
    hubble_a = hubble_function(All.Time);
  else
    hubble_a = 1;


  /* initialize quantities for statitics
   */

  if((All.NumCurrentTiStep > 0))
    {
#if defined(LT_MV_CHEMICALDIFFUSION_FULLSTEP)
      if(Flag_FullStep)
#endif
	{
	  memset(prevMetalsAmount, 0, LT_NMetP * sizeof(double));
	  memset(newMetalsAmount, 0, LT_NMetP * sizeof(double));
	}

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(SphP[i].A != 0)
	  {
	    tistep = P[i].Ti_current - P[i].Ti_begstep;
	    //tistep = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
	    dt = tistep * All.Timebase_interval / hubble_a;
	    expfact = exp(SphP[i].A * dt);

	    for(j = 0; j < LT_NMetP; j++)
	      {
		if(Flag_FullStep)
		  prevMetalsAmount[j] += SphP[i].Metals[j] + SphP[i].deltaMetals[j];
		SphP[i].Abundance[j] =
		  SphP[i].Abundance[j] * expfact + SphP[i].B[j] / SphP[i].A * (1 - expfact);
		if(Flag_FullStep)
		  newMetalsAmount[j] += SphP[i].Abundance[j] * SphP[i].prevMass;
	      }
	  }

#if defined(LT_MV_CHEMICALDIFFUSION_FULLSTEP)
      if(Flag_FullStep)
#endif
	/* now re-normalize the amount of metals of each specie accross
	   the particles that participated to diffusion
	 */
	{

	  MPI_Allreduce(prevMetalsAmount, prevMetalsAmount_tot, LT_NMetP, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	  MPI_Allreduce(newMetalsAmount, newMetalsAmount_tot, LT_NMetP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  for(j = 0; j < LT_NMetP; j++)
	    if(newMetalsAmount_tot[j] > 0)
	      factor[j] = prevMetalsAmount_tot[j] / newMetalsAmount_tot[j];

	  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	    if(SphP[i].A != 0)
	      for(j = 0; j < LT_NMetP; j++)
		{
		  SphP[i].Metals[j] = SphP[i].Abundance[j] * factor[j] + SphP[i].deltaMetals[j];
		  SphP[i].Abundance[j] = SphP[i].Metals[j] / P[i].Mass;
		  SphP[i].deltaMetals[j] = 0;
		}
	}
      /* NOTE: the following code just perform a consistency check in a case in which the amount of metals */
      /* is strictly costant. */

      /* if(Flag_FullStep) */
      /*   { */
      /*     memset(newMetalsAmount, 0, LT_NMetP * sizeof(double)); */
      /*     for(i = 0; i < N_gas; i++) */
      /*       for(j = 0; j < LT_NMetP; j++) */
      /*        newMetalsAmount[j] += SphP[i].Abundance[j] * P[i].Mass; */

      /*     /\* check whether some significant variation occurred */
      /*      *\/ */
      /*     MPI_Reduce(newMetalsAmount, newMetalsAmount_tot, LT_NMetP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
      /*     if(ThisTask == 0) */
      /*       for(j = 0; j < LT_NMetP; j++) */
      /*      if(fabs(All.MetalsAmount[j] - newMetalsAmount_tot[j]) / All.MetalsAmount[j] > 0.001) */
      /*        printf("\t\t>>> warning: step %d mass of element %d has not been conserved by %g percent", */
      /*               All.NumCurrentTiStep, j, fabs(All.MetalsAmount[j] - newMetalsAmount_tot[j]) / All.MetalsAmount[j]); */
      /*   } */

      MPI_Barrier(MPI_COMM_WORLD);
    }

  return 0;
}
#endif
