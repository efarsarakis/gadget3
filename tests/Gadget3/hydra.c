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
#ifdef MACHNUM
#include "machfinder.h"
#endif
#ifdef CS_MODEL
#include "cs_model/cs_metals.h"
#endif

#ifdef JD_DPP
#include "cr_electrons.h"
#endif

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

#ifdef NUM_THREADS
#include <pthread.h>
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

/*! \file hydra.c
*  \brief Computation of SPH forces and rate of entropy generation
*
*  This file contains the "second SPH loop", where the SPH forces are
*  computed, and where the rate of change of entropy due to the shock heating
*  (via artificial viscosity) is computed.
*/

#ifdef MACHNUM
double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;
#else
static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv;
#endif

struct kernel_hydra
{
  double dx, dy, dz;
  double r, vsig, sound_i, sound_j;
  double dvx, dvy, dvz, vdotr2;
  double wk_i, wk_j, dwk_i, dwk_j;
  double h_i, h_j, dwk_ij, rho_ij_inv;
#ifdef MAGNETIC
  double mj_r;
#ifdef MAGFORCE
  double mf_i, mf_j;
#endif
#if defined(MAGNETIC_SIGNALVEL)
  double alfven2_i, alfven2_j;
#endif
#if defined(MAGFORCE) || defined(MAGNETIC_SIGNALVEL)
  double b2_i, b2_j;
#endif
#if defined(EULER_DISSIPATION) || defined(MAGNETIC_DISSIPATION) || defined(FS_ETA2_DYN) || defined(MAGNETIC_DIFFUSION)
  double mf_dissInd;
#endif
#if defined(EULER_DISSIPATION_HEAT) || defined(MAGNETIC_DISSIPATION) || defined(MAGNETIC_DIFFUSION_HEAT)
  double mf_dissEnt;
#endif
#ifndef EULERPOTENTIALS
  double mf_Ind;
#endif
#endif				/* MAGNETIC */
#ifdef LT_MV_CHEMICALDIFFUSION_WAKEUP
  double Dvsig;
#endif
};

struct hydrodata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat Density;
  MyFloat Pressure;
  MyFloat F1;
  MyFloat DhsmlDensityFactor;
  int Timestep;

#ifdef CS_MODEL
  MyFloat DensityNow;
  MyFloat Entropy;
#endif

#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

#ifdef JD_VTURB
  MyFloat Vbulk[3];
#endif

#ifdef AB_SHOCK
  MyFloat Shock_N1[3];
#ifdef AB_SHOCK_VELDIV
  MyFloat Shock_N2[3];
  MyFloat Shock_N3[3];
#endif
#endif

#ifdef BP_REAL_CRs
  MyFloat CRpPressure;
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  MyFloat CRpE[BP_REAL_CRs];
  MyFloat CRpN[BP_REAL_CRs];
#endif
#endif

#ifdef MAGNETIC
#if !(!defined(MAGFORCE) && defined(EULERPOTENTIALS))
  MyFloat BPred[3];
#endif
#ifdef EULER_DISSIPATION
  MyFloat EulerA, EulerB;
#endif
#ifdef HIGH_ORDER_INDUCTION
  MyFloat Xix[3], Xiy[3], Xiz[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat PhiPred;
#endif
#ifdef FS_ALFA2_DYN
  MyFloat alfa2;
#endif
#ifdef FS_ETA2_DYN
  MyFloat eta2;
  MyFloat RotB[3];
#endif
#endif				// Magnetic

#ifdef TIME_DEP_ART_VISC
  MyFloat alpha;
#endif
#ifdef NAVIERSTOKES
  MyFloat Entropy;
  MyFloat stressoffdiag[3];
  MyFloat stressdiag[3];
  MyFloat shear_viscosity;
#endif

#ifdef NAVIERSTOKES_BULK
  MyFloat divvel;
#endif

#ifdef EOS_DEGENERATE
  MyFloat dpdr;
#endif

#ifdef LT_MV_CHEMICALDIFFUSION_WAKEUP
  MyFloat Dvsig;
#endif

#ifndef DONOTUSENODELIST
  int NodeList[NODELISTLENGTH];
#endif
}
 *HydroDataIn, *HydroDataGet;


struct hydrodata_out
{
  MyLongDouble Acc[3];
  MyLongDouble DtEntropy;
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;
#endif

#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  MyFloat DtCRpE[BP_REAL_CRs];
  MyFloat DtCRpN[BP_REAL_CRs];
#endif

#ifdef AB_SHOCK
  MyFloat Shock_Up_V1[3];
  MyFloat Shock_Down_V1[3];
  MyFloat Shock_Up_Signal;
  MyFloat Shock_Down_Signal;
  MyFloat Shock_Up_Weight;
  MyFloat Shock_Down_Weight;
  MyFloat Shock_Up_Rho;
  MyFloat Shock_Down_Rho;
  MyFloat Shock_Up_Pressure;
  MyFloat Shock_Down_Pressure;
#ifdef AB_SHOCK_VELDIV
  MyFloat Shock_Up_V2[3];
  MyFloat Shock_Down_V2[3];
  MyFloat Shock_Up_V3[3];
  MyFloat Shock_Down_V3[3];
  MyFloat Shock_Up_Weight23[2];
  MyFloat Shock_Down_Weight23[2];
#endif
#endif				/* AB_SHOCK */

#ifdef JD_VTURB
  MyFloat Vrms;
#endif

#if defined(JD_VTURB) && defined(JD_DECOMPOSE_VTURB)
  MyFloat Vrad;
  MyFloat Vtan;
#endif

#if defined(MAGNETIC)
#if !defined(EULERPOTENTIALS)
  MyFloat DtB[3];
#endif
#ifdef DIVBFORCE3
  MyFloat magacc[3];
  MyFloat magcorr[3];
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(FS_DEDNER_DIRECT)
  MyFloat GradPhi[3];
#endif
#endif
#if defined(EULERPOTENTIALS) && defined(EULER_DISSIPATION)
  MyFloat DtEulerA, DtEulerB;
#endif
#if  defined(CR_SHOCK)
  MyFloat CR_EnergyChange[NUMCRPOP];
  MyFloat CR_BaryonFractionChange[NUMCRPOP];
#endif
}
 *HydroDataResult, *HydroDataOut;


static inline void particle2in_hydra(struct hydrodata_in *in, int i);
static inline void out2particle_hydra(struct hydrodata_out *out, int i, int mode);

#ifdef AB_SHOCK
static inline void hydra_evaluate_shockfinder(struct hydrodata_in *local, struct hydrodata_out *out,
					      struct kernel_hydra *kernel, int j);
#endif

static inline void particle2in_hydra(struct hydrodata_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    {
      in->Pos[k] = P[i].Pos[k];
      in->Vel[k] = SphP[i].VelPred[k];
    }
  in->Hsml = PPP[i].Hsml;
  in->Mass = P[i].Mass;
  in->DhsmlDensityFactor = SphP[i].h.DhsmlDensityFactor;
  in->Density = SphP[i].d.Density;
  in->Pressure = SphP[i].Pressure;
  in->Timestep = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);

#ifdef EOS_DEGENERATE
  in->dpdr = SphP[i].dpdr;
#endif

#ifdef BP_REAL_CRs
  in->CRpPressure = SphP[i].CRpPressure;
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  int Nbin;
  for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
    {
      in->CRpE[Nbin] = SphP[i].CRpE[Nbin];
      in->CRpN[Nbin] = SphP[i].CRpN[Nbin];
    }
#endif
#endif

#ifndef NO_SHEAR_VISCOSITY_LIMITER
#ifndef ALTVISCOSITY
#ifndef EOS_DEGENERATE
  double sound = sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density);
#else
  double sound = sqrt(SphP[i].dpdr);
#endif
#ifndef NAVIERSTOKES
  in->F1 = fabs(SphP[i].v.DivVel) / (fabs(SphP[i].v.DivVel) + SphP[i].r.CurlVel +
				     0.0001 * sound / PPP[i].Hsml / fac_mu);
#else
  in->F1 = fabs(SphP[i].v.DivVel) / (fabs(SphP[i].v.DivVel) + SphP[i].u.s.CurlVel +
				     0.0001 * sound / PPP[i].Hsml / fac_mu);
#endif
#else
  in->F1 = SphP[i].v.DivVel;
#endif
#endif

#ifdef CS_MODEL
  in->DensityNow = SphP[i].d.Density;
  in->Entropy = SphP[i].Entropy;
#endif

#ifdef JD_VTURB
  in->Vbulk[0] = SphP[i].Vbulk[0];
  in->Vbulk[1] = SphP[i].Vbulk[1];
  in->Vbulk[2] = SphP[i].Vbulk[2];
#endif

#ifdef AB_SHOCK
  for(k = 0; k < 3; k++)
    {
      in->Shock_N1[k] = SphP[i].Shock_N1[k];
#ifdef AB_SHOCK_VELDIV
      in->Shock_N2[k] = SphP[i].Shock_N2[k];
      in->Shock_N3[k] = SphP[i].Shock_N3[k];
#endif
    }
#endif

#ifdef MAGNETIC
  for(k = 0; k < 3; k++)
    {
#ifdef HIGH_ORDER_INDUCTION
      in->Xix[k] = SphP[i].Xix[k];
      in->Xiy[k] = SphP[i].Xiy[k];
      in->Xiz[k] = SphP[i].Xiz[k];
#endif
#if !(!defined(MAGFORCE) && defined(EULERPOTENTIALS))
      in->BPred[k] = SphP[i].b2.BPred[k];
#if defined(SFR)
      in->BPred[k] *= pow(1. - SphP[i].XColdCloud, 2. * POW_CC);
#endif
#endif
    }
#ifdef FS_ALFA2_DYN
#ifdef FS_DYN_DEF
  if((P[i].Pos[0] - LONG_Y / 2.) * (P[i].Pos[0] - LONG_Y / 2.) +
     (P[i].Pos[1] - LONG_X / 2) * (P[i].Pos[1] - LONG_X / 2.) < 16)
    {
      SphP[i].alfa2 =
	All.Alfa0Quench / (SysState.EnergyMag / All.AlfaB2eq +
			   1.) * sin((fabs(P[i].Pos[2] - LONG_Z / 2.) -
				      LONG_Z / 2.) * M_PI / (All.Alfa2H)) * OmegaR(i, BRANDT_MODE);
      SphP[i].alfa2 = (fabs(P[i].Pos[2] - LONG_Z / 2.) > All.Alfa2H) ? SphP[i].alfa2 =
	0. : (P[i].Pos[2] - LONG_Z / 2.) < 0.0 ? -fabs(SphP[i].alfa2) : fabs(SphP[i].alfa2);
    }
  else
    SphP[i].alfa2 = 0;
#ifdef FS_FIX_ALPHA
  SphP[i].alfa2 = FS_FIX_ALPHA;
#endif
#endif
  in->alfa2 = SphP[i].alfa2;
#endif
#ifdef FS_ETA2_DYN
  in->eta2 = SphP[i].eta2;

  in->RotB[0] = SphP[i].RotB[0];
  in->RotB[1] = SphP[i].RotB[1];
  in->RotB[2] = SphP[i].RotB[2];
#endif
#ifdef TIME_DEP_MAGN_DISP
  in->Balpha = SphP[i].Balpha;
#endif
#if defined(EULERPOTENTIALS) && defined(EULER_DISSIPATION)
  in->EulerA = SphP[i].EulerA;
  in->EulerB = SphP[i].EulerB;
#endif
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
  in->PhiPred = SphP[i].SmoothPhi;
#else
  in->PhiPred = SphP[i].PhiPred;
#endif
#endif
#endif

#ifdef TIME_DEP_ART_VISC
  in->alpha = SphP[i].alpha;
#endif

#ifdef PARTICLE_DEBUG
  in->ID = P[i].ID;
#endif

#ifdef NAVIERSTOKES
  in->Entropy = SphP[i].Entropy;
  for(k = 0; k < 3; k++)
    {
      in->stressdiag[k] = SphP[i].u.s.StressDiag[k];
      in->stressoffdiag[k] = SphP[i].u.s.StressOffDiag[k];
    }
  in->shear_viscosity = get_shear_viscosity(i);
#ifdef NAVIERSTOKES_BULK
  in->divvel = SphP[i].u.s.DivVel;
#endif
#endif

#ifdef LT_MV_CHEMICALDIFFUSION_WAKEUP
  double myD, maxgrad;
#ifdef LT_MV_CHEMICALDIFFUSION_CONSTANT_D
  myD = All.ConstantD;
#else
  myD = SphP[place].D;
#endif
  maxgrad = 0;
  for(k = 0; k < LT_NMetP; k++)
    if(maxgrad < fabs(SphP[i].GradAbnd[k]))
      maxgrad = fabs(SphP[i].GradAbnd[k]);
  in->Dvsig = myD * maxgrad;
#endif
}

static inline void out2particle_hydra(struct hydrodata_out *out, int i, int mode)
{
  int k;

  for(k = 0; k < 3; k++)
    ASSIGN_ADD(SphP[i].a.dHydroAccel[k], out->Acc[k], mode);
  ASSIGN_ADD(SphP[i].e.dDtEntropy, out->DtEntropy, mode);

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  if(mode == 0)
    SphP[i].MinViscousDt = out->MinViscousDt;
  else if(SphP[i].MinViscousDt > out->MinViscousDt)
    SphP[i].MinViscousDt = out->MinViscousDt;
#else
  if(mode == 0)
    SphP[i].MaxSignalVel = out->MaxSignalVel;
  else if(SphP[i].MaxSignalVel < out->MaxSignalVel)
    SphP[i].MaxSignalVel = out->MaxSignalVel;
#endif

#ifdef JD_VTURB
  ASSIGN_ADD(SphP[i].Vrms, out->Vrms, mode);
#endif

#if defined(JD_VTURB) && defined(JD_DECOMPOSE_VTURB)
  ASSIGN_ADD(SphP[i].Vtan, out->Vtan, mode);
  ASSIGN_ADD(SphP[i].Vrad, out->Vrad, mode);
#endif

#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  int Nbin;
  for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
    {
      ASSIGN_ADD(SphP[i].DtCRpE[Nbin], out->DtCRpE[Nbin], mode);
      ASSIGN_ADD(SphP[i].DtCRpN[Nbin], out->DtCRpN[Nbin], mode);
    }
#endif

#ifdef AB_SHOCK
  for(k = 0; k < 3; k++)
    {
      ASSIGN_ADD(SphP[i].Shock_Up_V1[k], out->Shock_Up_V1[k], mode);
      ASSIGN_ADD(SphP[i].Shock_Down_V1[k], out->Shock_Down_V1[k], mode);
    }
  ASSIGN_ADD(SphP[i].Shock_Up_Signal, out->Shock_Up_Signal, mode);
  ASSIGN_ADD(SphP[i].Shock_Down_Signal, out->Shock_Down_Signal, mode);
  ASSIGN_ADD(SphP[i].Shock_Up_Weight, out->Shock_Up_Weight, mode);
  ASSIGN_ADD(SphP[i].Shock_Down_Weight, out->Shock_Down_Weight, mode);
  ASSIGN_ADD(SphP[i].Shock_Up_Rho, out->Shock_Up_Rho, mode);
  ASSIGN_ADD(SphP[i].Shock_Down_Rho, out->Shock_Down_Rho, mode);
  ASSIGN_ADD(SphP[i].Shock_Up_Pressure, out->Shock_Up_Pressure, mode);
  ASSIGN_ADD(SphP[i].Shock_Down_Pressure, out->Shock_Down_Pressure, mode);
#ifdef AB_SHOCK_VELDIV
  for(k = 0; k < 3; k++)
    {
      ASSIGN_ADD(SphP[i].Shock_Up_V2[k], out->Shock_Up_V2[k], mode);
      ASSIGN_ADD(SphP[i].Shock_Down_V2[k], out->Shock_Down_V2[k], mode);
      ASSIGN_ADD(SphP[i].Shock_Up_V3[k], out->Shock_Up_V3[k], mode);
      ASSIGN_ADD(SphP[i].Shock_Down_V3[k], out->Shock_Down_V3[k], mode);
    }
  ASSIGN_ADD(SphP[i].Shock_Up_Weight23[0], out->Shock_Up_Weight23[0], mode);
  ASSIGN_ADD(SphP[i].Shock_Down_Weight23[0], out->Shock_Down_Weight23[0], mode);
  ASSIGN_ADD(SphP[i].Shock_Up_Weight23[1], out->Shock_Up_Weight23[1], mode);
  ASSIGN_ADD(SphP[i].Shock_Down_Weight23[1], out->Shock_Down_Weight23[1], mode);
#endif
#endif /* AB_SHOCK */

#if defined(MAGNETIC)
  for(k = 0; k < 3; k++)
    {
#ifndef EULERPOTENTIALS
      ASSIGN_ADD(SphP[i].DtB[k], out->DtB[k], mode);
#endif
#ifdef DIVBFORCE3
      ASSIGN_ADD(SphP[i].magacc[k], out->magacc[k], mode);
      ASSIGN_ADD(SphP[i].magcorr[k], out->magcorr[k], mode);
#endif
#if defined(DIVBCLEANING_DEDNER) && !defined(FS_DEDNER_DIRECT)
      ASSIGN_ADD(SphP[i].GradPhi[k], out->GradPhi[k], mode);
#endif
    }
#endif
#if defined(EULERPOTENTIALS) && defined(EULER_DISSIPATION)
  ASSIGN_ADD(SphP[i].DtEulerA, out->DtEulerA, mode);
  ASSIGN_ADD(SphP[i].DtEulerB, out->DtEulerB, mode);
#endif
}

/*! This function is the driver routine for the calculation of hydrodynamical
*  force and rate of change of entropy due to shock heating for all active
*  particles .
*/
void hydro_force(void)
{
  int i, j, k, ngrp, ndone, ndone_flag;
  int recvTask, place;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0, timenetwork = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;

  int save_NextParticle;

  long long n_exported = 0;

#ifdef NAVIERSTOKES
  double fac;
#endif

#if (!defined(COOLING) && !defined(CR_SHOCK) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION)))
  double utherm;
  double dt;
  int CRpop;
#endif

#if defined(CR_SHOCK)
  double rShockEnergy;
  double rNonRethermalizedEnergy;

#ifndef COOLING
  double utherm, CRpop;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
  double f, cs_h;
#endif

#if defined(HEALPIX)
  double r_new, t[3];
  long ipix;
  int count = 0;
  int count2 = 0;
  int total_count = 0;
#endif
#if defined(HEALPIX) && defined(DIVBCLEANING_DEDNER)
  double ded_heal_fac = 0;
#endif

#ifdef NUCLEAR_NETWORK
  double dedt_nuc;
  int nuc_particles = 0;
  int nuc_particles_sum;
#endif

#if defined(WAKEUP) || defined(LT_MV_CHEMICALDIFFUSION_WAKEUP)
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	SphP[i].wakeup = 0;
    }
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = hubble_function(All.Time);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
#ifdef MACHNUM
      fac_egy = pow(All.Time, 3 * (GAMMA - 1));
#endif
      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    {
      hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = 1.0;
#ifdef MACHNUM
      fac_egy = 1.0;
#endif
    }

  /* allocate buffers to arrange communication */

  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct hydrodata_in) +
					     sizeof(struct hydrodata_out) +
					     sizemax(sizeof(struct hydrodata_in),
						     sizeof(struct hydrodata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = second();

  NextParticle = FirstActiveParticle;	/* beginn with this index */

  do
    {

      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
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
	  pthread_create(&mythreads[j], &attr, hydro_evaluate_primary, &threadid[j]);
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
	hydro_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
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

      HydroDataGet = (struct hydrodata_in *) mymalloc("HydroDataGet", Nimport * sizeof(struct hydrodata_in));
      HydroDataIn = (struct hydrodata_in *) mymalloc("HydroDataIn", Nexport * sizeof(struct hydrodata_in));

      /* prepare particle data for export */

      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  particle2in_hydra(&HydroDataIn[j], place);
#ifndef DONOTUSENODELIST
	  memcpy(HydroDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif

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
		  MPI_Sendrecv(&HydroDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &HydroDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(HydroDataIn);
      HydroDataResult =
	(struct hydrodata_out *) mymalloc("HydroDataResult", Nimport * sizeof(struct hydrodata_out));
      HydroDataOut =
	(struct hydrodata_out *) mymalloc("HydroDataOut", Nexport * sizeof(struct hydrodata_out));


      report_memory_usage(&HighMark_sphhydro, "SPH_HYDRO");

      /* now do the particles that were sent to us */

      tstart = second();

      NextJ = 0;

#ifdef NUM_THREADS
      for(j = 0; j < NUM_THREADS - 1; j++)
	pthread_create(&mythreads[j], &attr, hydro_evaluate_secondary, &threadid[j]);
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
	hydro_evaluate_secondary(&mainthreadid);
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
		  MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B,
			       &HydroDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	  out2particle_hydra(&HydroDataOut[j], place, 1);
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(HydroDataOut);
      myfree(HydroDataResult);
      myfree(HydroDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


  /* do final operations on results */


#ifdef FLTROUNDOFFREDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	SphP[i].e.DtEntropy = FLT(SphP[i].e.dDtEntropy);

	for(j = 0; j < 3; j++)
	  SphP[i].a.HydroAccel[j] = FLT(SphP[i].a.dHydroAccel[j]);
      }
#endif



#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  int il;
#pragma omp parallel for private(i,k)
  for(il = 0; il < NActivePart; il++)
    {
      i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
#ifdef BLACK_HOLES
      if(P[i].Type == 0 && P[i].Mass > 0)
#else
      if(P[i].Type == 0)
#endif
	{
#ifdef CR_SHOCK
	  /* state right here:
	   *
	   * _c denotes comoving quantities
	   * _p denotes physical quantities
	   *
	   *
	   * Delta u_p = rho_p^(gamma-1)/(gamma-1) Delta A
	   *
	   * Delta A = dA/dloga * Delta loga
	   *
	   * dA/dloga = DtE * (gamma-1) / ( H(a) a^2 rho_c^(gamma-1)
	   *
	   * => Delta u_p = DtE * dloga / ( H(a) a^2 a^(3(gamma-1)) )
	   */

	  if(SphP[i].e.DtEntropy > 0.0)
	    {
	      rShockEnergy = SphP[i].e.DtEntropy *
		(P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a2 / fac_egy;
	    }
	  else
	    {
	      rShockEnergy = 0.0;
	    }

#endif /* CR_SHOCK */

#if !defined(EOS_DEGENERATE)

#ifndef TRADITIONAL_SPH_FORMULATION
	  /* Translate energy change rate into entropy change rate */
	  SphP[i].e.DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1));
#endif

#else
	  /* DtEntropy stores the energy change rate in internal units */
	  SphP[i].e.DtEntropy *= All.UnitEnergy_in_cgs / All.UnitTime_in_s;
#endif

#ifdef MACHNUM

	  /* Estimates the Mach number of particle i for non-radiative runs,
	   * or the Mach number, density jump and specific energy jump
	   * in case of cosmic rays!
	   */
#if (CR_SHOCK == 2)
	  GetMachNumberCR(SphP + i);
#else

#ifndef CS_MODEL
	  GetMachNumber(SphP + i);
#else
	  GetMachNumber(SphP + i, P + i);
#endif /* CS_MODEL */
#endif /* COSMIC_RAYS */
#endif /* MACHNUM */
#ifdef MACHSTATISTIC
	  GetShock_DtEnergy(SphP + i);
#endif

#ifdef CR_SHOCK
	  if(rShockEnergy > 0.0)
	    {
	      /* Feed fraction "All.CR_ShockEfficiency" into CR and see what
	       * amount of energy instantly gets rethermalized
	       *
	       * for this, we need the physical time step, which is
	       * Delta t_p = Delta t_c / hubble_a
	       */

	      /* The  CR_find_alpha_InjectTo induces an error in the density jump since it can set
	       *  Particle->Shock_DensityJump = 1.0 + 1.0e-6 which is used in ShockInject as the input DensityJump
	       *  if (NUMCRPOP > 1)
	       *  {
	       *  #if ( CR_SHOCK == 1 )
	       *  InjPopulation = CR_Find_Alpha_to_InjectTo(All.CR_ShockAlpha);
	       *  #else
	       *  InjPopulation = CR_find_alpha_InjectTo(SphP + i);
	       *  #endif
	       *  }
	       *  else
	       *  InjPopulation = 0;
	       */

	      rNonRethermalizedEnergy =
		CR_Particle_ShockInject(SphP + i,
					rShockEnergy,
					(P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval /
					hubble_a);

	      /* Fraction of total energy that went and remained in CR is
	       * rNonRethermalizedEnergy / rShockEnergy,
	       * hence, we conserve energy if we do:
	       */
#ifndef CR_NO_CHANGE
	      SphP[i].e.DtEntropy *= (1.0 - rNonRethermalizedEnergy / rShockEnergy);
#endif /* CR_NO_CHANGE */

	      assert(rNonRethermalizedEnergy >= 0.0);

	      assert(rNonRethermalizedEnergy <= (rShockEnergy * All.CR_ShockEfficiency));


#if (!defined(COOLING) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION)))
	      utherm = 0.0;
	      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
		utherm +=
		  CR_Particle_ThermalizeAndDissipate(SphP + i,
						     (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) *
						     All.Timebase_interval / hubble_a, CRpop);

	      /* we need to add this thermalized energy to the internal energy */

	      SphP[i].e.DtEntropy += GAMMA_MINUS1 * utherm * fac_egy / pow(SphP[i].d.Density, GAMMA_MINUS1) /
		((P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval);
#endif

	    }
#endif /* CR_SHOCK */


#if (!defined(COOLING) && !defined(CR_SHOCK) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION)))
	  double utherm;
	  double dt;
	  int CRpop;

	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a;

	  if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
	    {
	      if(dt > 0)
		{
		  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
		    {
		      utherm = CR_Particle_ThermalizeAndDissipate(SphP + i, dt, CRpop);

		      SphP[i].e.DtEntropy +=
			GAMMA_MINUS1 * utherm * fac_egy / pow(SphP[i].d.Density,
							      GAMMA_MINUS1) / (dt * hubble_a);
		    }
		}
	    }
#endif

#ifdef NAVIERSTOKES
	  /* sigma_ab * sigma_ab */
	  for(k = 0, fac = 0; k < 3; k++)
	    {
	      fac += SphP[i].u.s.StressDiag[k] * SphP[i].u.s.StressDiag[k] +
		2 * SphP[i].u.s.StressOffDiag[k] * SphP[i].u.s.StressOffDiag[k];
	    }

#ifndef NAVIERSTOKES_CONSTANT	/*entropy increase due to the shear viscosity */
#ifdef NS_TIMESTEP
	  SphP[i].ViscEntropyChange = 0.5 * GAMMA_MINUS1 /
	    (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	    get_shear_viscosity(i) / SphP[i].d.Density * fac *
	    pow((SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);

	  SphP[i].e.DtEntropy += SphP[i].ViscEntropyChange;
#else
	  SphP[i].e.DtEntropy += 0.5 * GAMMA_MINUS1 /
	    (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	    get_shear_viscosity(i) / SphP[i].d.Density * fac *
	    pow((SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);
#endif

#else
	  SphP[i].e.DtEntropy += 0.5 * GAMMA_MINUS1 /
	    (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	    get_shear_viscosity(i) / SphP[i].d.Density * fac;

#ifdef NS_TIMESTEP
	  SphP[i].ViscEntropyChange = 0.5 * GAMMA_MINUS1 /
	    (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	    get_shear_viscosity(i) / SphP[i].d.Density * fac;
#endif

#endif

#ifdef NAVIERSTOKES_BULK	/*entropy increase due to the bulk viscosity */
	  SphP[i].e.DtEntropy += GAMMA_MINUS1 /
	    (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	    All.NavierStokes_BulkViscosity / SphP[i].d.Density * pow(SphP[i].u.s.a4.DivVel, 2);

#ifdef NS_TIMESTEP
	  SphP[i].ViscEntropyChange = GAMMA_MINUS1 /
	    (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	    All.NavierStokes_BulkViscosity / SphP[i].d.Density * pow(SphP[i].u.s.a4.DivVel, 2);
#endif

#endif

#endif /* these entropy increases directly follow from the general heat transfer equation */


#ifdef JD_VTURB
	  SphP[i].Vrms += (SphP[i].VelPred[0] - SphP[i].Vbulk[0]) * (SphP[i].VelPred[0] - SphP[i].Vbulk[0])
	    + (SphP[i].VelPred[1] - SphP[i].Vbulk[1]) * (SphP[i].VelPred[1] - SphP[i].Vbulk[1])
	    + (SphP[i].VelPred[2] - SphP[i].Vbulk[2]) * (SphP[i].VelPred[2] - SphP[i].Vbulk[2]);
	  SphP[i].Vrms = sqrt(SphP[i].Vrms / P[i].TrueNGB);
#endif

#if defined(JD_VTURB) && defined(JD_DECOMPOSE_VTURB)
	  /* particle i has only a tangential component by definition */
	  SphP[i].Vtan += (SphP[i].VelPred[0] - SphP[i].Vbulk[0]) * (SphP[i].VelPred[0] - SphP[i].Vbulk[0])
	    + (SphP[i].VelPred[1] - SphP[i].Vbulk[1]) * (SphP[i].VelPred[1] - SphP[i].Vbulk[1])
	    + (SphP[i].VelPred[2] - SphP[i].Vbulk[2]) * (SphP[i].VelPred[2] - SphP[i].Vbulk[2]);

	  /* do the root mean of the added squares */
	  SphP[i].Vtan = sqrt(SphP[i].Vtan / P[i].TrueNGB);
	  SphP[i].Vrad = sqrt(SphP[i].Vrad / P[i].TrueNGB);

#endif

#if defined(JD_DPP) && !defined(JD_DPPONSNAPSHOTONLY)
	  compute_Dpp(i);
#endif

#ifdef WINDS
	  /* if we have winds, we decouple particles briefly if delaytime>0 */

	  if(SphP[i].DelayTime > 0)
	    {

#ifdef EULER_DISSIPATION
	      SphP[i].DtEulerA = 0;
	      SphP[i].DtEulerB = 0;
#endif

	      for(k = 0; k < 3; k++)
		SphP[i].a.HydroAccel[k] = 0;

	      SphP[i].e.DtEntropy = 0;

#ifdef NOWINDTIMESTEPPING
	      SphP[i].MaxSignalVel = 2 * sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density);
#else
	      double windspeed, hsml_c;
#if !defined(LT_WIND_VELOCITY) && !defined(LT_STELLAREVOLUTION)
	      windspeed = sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			       All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency) * All.Time;
#else
#ifdef LT_WIND_VELOCITY
	      windspeed = LT_WIND_VELOCITY * All.Time;
#else
	      int IMFi, SFi;
	      get_SF_index(i, &SFi, &IMFi);
	      windspeed = sqrt(2 * SFs[SFi].WindEnergyFraction * SFs[SFi].totFactorSN *
			       SFs[SFi].EgySpecSN / (1 - SFs[SFi].totFactorSN) / SFs[SFi].WindEfficiency) *
		All.Time;
#endif
#endif
	      windspeed *= fac_mu;
#ifdef QUICK_LYALPHA
	      hsml_c = pow(All.WindFreeTravelDensFac * All.OverDensThresh / SphP[i].d.Density, (1. / 3.));
#else
	      hsml_c = pow(All.WindFreeTravelDensFac * All.PhysDensThresh /
			   (SphP[i].d.Density * a3inv), (1. / 3.));
#endif
	      SphP[i].MaxSignalVel = hsml_c * DMAX((2 * windspeed), SphP[i].MaxSignalVel);
#endif
	    }
#endif

#if defined(HEALPIX)
	  r_new = 0;
	  ded_heal_fac = 1.;
	  for(k = 0; k < 3; k++)
	    {
	      t[k] = P[i].Pos[k] - SysState.CenterOfMassComp[0][k];
	      r_new = r_new + t[k] * t[k];
	    }
	  r_new = sqrt(r_new);
	  vec2pix_nest((long) All.Nside, t, &ipix);
	  if(r_new > All.healpixmap[ipix] * HEALPIX)
	    {
	      SphP[i].e.DtEntropy = 0;
	      for(k = 0; k < 3; k++)
		{
		  SphP[i].a.HydroAccel[k] = 0;
		  SphP[i].VelPred[k] = P[i].Vel[k] = 0.0;
#ifdef DIVBFORCE3
		  SphP[i].magcorr[k] = 0;
		  SphP[i].magacc[k] = 0;
#endif
		}
#ifdef FS_ETA2_TURB
	      SphP[i].eta2 = 0;
#endif
#if defined(DIVBCLEANING_DEDNER)
	      ded_heal_fac = 2.;
#endif
	      SphP[i].v.DivVel = 0.0;
	      SphP[i].r.CurlVel = 0.0;
	      count++;
	      if(r_new > All.healpixmap[ipix] * HEALPIX * 1.5)
		{
		  count2++;
		}
	    }
#endif

#ifdef TIME_DEP_ART_VISC
#if !defined(EOS_DEGENERATE)
	  cs_h = sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density) / PPP[i].Hsml;
#else
	  cs_h = sqrt(SphP[i].dpdr) / PPP[i].Hsml;
#endif
	  f = fabs(SphP[i].v.DivVel) / (fabs(SphP[i].v.DivVel) + SphP[i].r.CurlVel + 0.0001 * cs_h / fac_mu);
	  SphP[i].Dtalpha = -(SphP[i].alpha - All.AlphaMin) * All.DecayTime *
#ifdef JD_TIME_DEP_ART_VISC
	    0.5 * cs_h / (PPP[i].Hsml * fac_mu)
#else
	    0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu)
#endif
	    + f * All.ViscSource * DMAX(0.0, -SphP[i].v.DivVel);
	  if(All.ComovingIntegrationOn)
	    SphP[i].Dtalpha /= (hubble_a * All.Time * All.Time);
#endif
#ifdef MAGNETIC
#ifdef TIME_DEP_MAGN_DISP
	  SphP[i].DtBalpha = -(SphP[i].Balpha - All.ArtMagDispMin) * All.ArtMagDispTime *
	    0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu)
#ifndef ROT_IN_MAG_DIS
	    + All.ArtMagDispSource * fabs(SphP[i].divB) / sqrt(MU0 * SphP[i].d.Density);
#else
#ifdef SMOOTH_ROTB
	    + All.ArtMagDispSource / sqrt(MU0 * SphP[i].d.Density) *
	    DMAX(fabs(SphP[i].divB), fabs(sqrt(SphP[i].SmoothedRotB[0] * SphP[i].SmoothedRotB[0] +
					       SphP[i].SmoothedRotB[1] * SphP[i].SmoothedRotB[1] +
					       SphP[i].SmoothedRotB[2] * SphP[i].SmoothedRotB[2])));
#else
	    + All.ArtMagDispSource / sqrt(MU0 * SphP[i].d.Density) *
	    DMAX(fabs(SphP[i].divB), fabs(sqrt(SphP[i].RotB[0] * SphP[i].RotB[0] +
					       SphP[i].RotB[1] * SphP[i].RotB[1] +
					       SphP[i].RotB[2] * SphP[i].RotB[2])));
#endif /* End SMOOTH_ROTB        */
#endif /* End ROT_IN_MAG_DIS     */
#endif /* End TIME_DEP_MAGN_DISP */

#if defined(DIVBFORCE3) || defined(DIVBCLEANING_DEDNER)
	  double tmpb, phiphi;
#endif
#ifdef DIVBFORCE3
	  phiphi =
	    sqrt(pow(SphP[i].magcorr[0], 2.) + pow(SphP[i].magcorr[1], 2.) + pow(SphP[i].magcorr[2], 2.));
	  tmpb = sqrt(pow(SphP[i].magacc[0], 2.) + pow(SphP[i].magacc[1], 2.) + pow(SphP[i].magacc[2], 2.));

	  if(phiphi > DIVBFORCE3 * tmpb)
	    for(k = 0; k < 3; k++)
	      SphP[i].magcorr[k] *= DIVBFORCE3 * tmpb / phiphi;

	  for(k = 0; k < 3; k++)
	    SphP[i].a.HydroAccel[k] += (SphP[i].magacc[k] - SphP[i].magcorr[k]);

#endif

#ifdef DIVBCLEANING_DEDNER
#if defined(FS_DEDNER_DIRECT)
	  tmpb = sqrt(SphP[i].divB * SphP[i].divB * MU0_1 / SphP[i].d.Density) * PPP[i].Hsml;
#else
	  tmpb = 0.5 * SphP[i].MaxSignalVel;
#endif
	  phiphi = tmpb * tmpb * All.DivBcleanHyperbolicSigma
#if defined(HEALPIX)
	    / ded_heal_fac
#endif
#if defined( SMOOTH_DIVB) && !defined(BSMOOTH_TIME)
	    * SphP[i].SmoothDivB;
#else
	    * SphP[i].divB;
#endif
#ifdef SMOOTH_PHI
	  phiphi += SphP[i].SmoothPhi *
#else
	  phiphi += SphP[i].PhiPred *
#endif
#if defined(FS_DEDNER_DIRECT)
	    sqrt((pow(SphP[i].b2.BPred[0], 2.) + pow(SphP[i].b2.BPred[1], 2.) +
		  pow(SphP[i].b2.BPred[2], 2.)) * MU0_1 / SphP[i].d.Density) / PPP[i].Hsml *
#else
	    0.5 * SphP[i].MaxSignalVel / PPP[i].Hsml *
#endif
#if defined(HEALPIX)
	    ded_heal_fac *
#endif
	    All.DivBcleanParabolicSigma;

	  SphP[i].DtPhi = -phiphi * atime * atime;	/* Compensate for the 1/Ha^2 in dt_mag */
#ifndef FS_DEDNER_DIRECT
	  phiphi =
	    sqrt(pow(SphP[i].GradPhi[0], 2.) + pow(SphP[i].GradPhi[1], 2.) + pow(SphP[i].GradPhi[2], 2.));
	  tmpb = sqrt(pow(SphP[i].DtB[0], 2.) + pow(SphP[i].DtB[1], 2.) + pow(SphP[i].DtB[2], 2.));
	  if(phiphi > All.DivBcleanQ * tmpb && tmpb != 0)
	    for(k = 0; k < 3; k++)
	      SphP[i].GradPhi[k] *= tmpb * All.DivBcleanQ / phiphi;
	  SphP[i].e.DtEntropy -=
	    (SphP[i].b2.BPred[0] * SphP[i].GradPhi[0] + SphP[i].b2.BPred[1] * SphP[i].GradPhi[1] +
	     SphP[i].b2.BPred[2] * SphP[i].GradPhi[2])
	    * GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) * MU0_1 / SphP[i].d.Density;

	  SphP[i].DtB[0] += SphP[i].GradPhi[0];
	  SphP[i].DtB[1] += SphP[i].GradPhi[1];
	  SphP[i].DtB[2] += SphP[i].GradPhi[2];
#endif

#endif /* End DEDNER */
#endif /* End Magnetic */

#ifdef AB_SHOCK
	  SphP[i].Shock_Up_Rho /= SphP[i].Shock_Up_Weight;
	  SphP[i].Shock_Down_Rho /= SphP[i].Shock_Down_Weight;
	  SphP[i].Shock_Up_Pressure /= SphP[i].Shock_Up_Weight;
	  SphP[i].Shock_Down_Pressure /= SphP[i].Shock_Down_Weight;
	  SphP[i].Shock_Up_Signal /= SphP[i].Shock_Up_Weight;
	  SphP[i].Shock_Down_Signal /= SphP[i].Shock_Down_Weight;

	  double Shock_Up_V1 = (SphP[i].Shock_Up_V1[0] * SphP[i].Shock_N1[0]
				+ SphP[i].Shock_Up_V1[1] * SphP[i].Shock_N1[1]
				+ SphP[i].Shock_Up_V1[2] * SphP[i].Shock_N1[2]) / SphP[i].Shock_Up_Weight;

	  double Shock_Down_V1 = (SphP[i].Shock_Down_V1[0] * SphP[i].Shock_N1[0]
				  + SphP[i].Shock_Down_V1[1] * SphP[i].Shock_N1[1]
				  + SphP[i].Shock_Down_V1[2] * SphP[i].Shock_N1[2])
	    / SphP[i].Shock_Down_Weight;

	  double Shock_Veldiv1 = Shock_Down_V1 - Shock_Up_V1;

	  SphP[i].Shock_Compress = SphP[i].Shock_Down_Rho / SphP[i].Shock_Up_Rho;
	  SphP[i].Shock_Mach = 0;
	  SphP[i].Shock_Mach1 = 0;
	  SphP[i].Shock_Mach2 = 0;
	  SphP[i].Shock_Mach3 = 0;

	  if(Shock_Veldiv1 > 0)
	    {
#ifdef AB_SHOCK_VELDIV
	      double Shock_Up_V2 = (SphP[i].Shock_Up_V2[0] * SphP[i].Shock_N2[0]
				    + SphP[i].Shock_Up_V2[1] * SphP[i].Shock_N2[1]
				    + SphP[i].Shock_Up_V2[2] * SphP[i].Shock_N2[2])
		/ SphP[i].Shock_Up_Weight23[0];

	      double Shock_Down_V2 = (SphP[i].Shock_Down_V2[0] * SphP[i].Shock_N2[0]
				      + SphP[i].Shock_Down_V2[1] * SphP[i].Shock_N2[1]
				      + SphP[i].Shock_Down_V2[2] * SphP[i].Shock_N2[2])
		/ SphP[i].Shock_Down_Weight23[0];

	      double Shock_Up_V3 = (SphP[i].Shock_Up_V3[0] * SphP[i].Shock_N3[0]
				    + SphP[i].Shock_Up_V3[1] * SphP[i].Shock_N3[1]
				    + SphP[i].Shock_Up_V3[2] * SphP[i].Shock_N3[2])
		/ SphP[i].Shock_Up_Weight23[1];

	      double Shock_Down_V3 = (SphP[i].Shock_Down_V3[0] * SphP[i].Shock_N3[0]
				      + SphP[i].Shock_Down_V3[1] * SphP[i].Shock_N3[1]
				      + SphP[i].Shock_Down_V3[2] * SphP[i].Shock_N3[2])
		/ SphP[i].Shock_Down_Weight23[1];

	      double Shock_Veldiv2 = fabs(Shock_Down_V2 - Shock_Up_V2);
	      double Shock_Veldiv3 = fabs(Shock_Down_V3 - Shock_Up_V3);
	      if(Shock_Veldiv1 > 2.0 * Shock_Veldiv2)
		if(Shock_Veldiv1 > 2.0 * Shock_Veldiv3)
#endif
		  {
		    double Shock_Compress = fabs(1.0 - 1.0 / SphP[i].Shock_Compress);

		    double Shock_dpr = SphP[i].Shock_Down_Pressure / SphP[i].Shock_Up_Pressure;
		    SphP[i].Shock_Mach2 =
		      sqrt(fabs(Shock_dpr * (GAMMA + 1) / 2 / GAMMA + (GAMMA - 1) / 2 / GAMMA));
		    SphP[i].Shock_Mach3 =
		      sqrt(fabs(1 / ((GAMMA + 1) / SphP[i].Shock_Compress / 2 - (GAMMA - 1) / 2)));

		    if(Shock_Compress > SHOCK_NOISELEVEL)
		      {
			double Shock_vdif_fraction =
			  Shock_Veldiv1 / (atime * atime * SphP[i].Shock_Up_Signal);

			SphP[i].Shock_Mach1 = Shock_vdif_fraction / Shock_Compress;

			if(Shock_vdif_fraction > SHOCK_MINIMAL_VDIF_FRACTION &&
			   SphP[i].Shock_Down_Pressure > SHOCK_GRADIENT_TOLERANCE * SphP[i].Shock_Up_Pressure
			   && SphP[i].Shock_Down_Pressure > SphP[i].Pressure
			   && SphP[i].Pressure > SphP[i].Shock_Up_Pressure
			   && SphP[i].Shock_Down_Rho > SHOCK_GRADIENT_TOLERANCE * SphP[i].Shock_Up_Rho
			   && SphP[i].Shock_Down_Rho > SphP[i].d.Density
			   && SphP[i].d.Density > SphP[i].Shock_Up_Rho)
			  SphP[i].Shock_Mach = SphP[i].Shock_Mach1;
		      }
		  }
	    }

	  SphP[i].Shock_Speed = SphP[i].Shock_Mach * SphP[i].Shock_Up_Signal;
#endif /* AB_SHOCK */

#ifdef SPH_BND_PARTICLES
	  if(P[i].ID == 0)
	    {
	      SphP[i].e.DtEntropy = 0;
#ifdef NS_TIMESTEP
	      SphP[i].ViscEntropyChange = 0;
#endif

#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
#ifdef SPH_BND_DTB
	      for(k = 0; k < 3; k++)
		SphP[i].DtB[k] = 0;
#endif
#endif
#ifdef SPH_BND_BFLD
	      for(k = 0; k < 3; k++)
		SphP[i].b1.B[k] = 0.0;
#endif
#ifdef SPH_BND_MRI
	      SphP[i].a.HydroAccel[0] -= SphP[i].a.HydroAccel[0] * (P[i].Pos[0] - 10.)
		/ sqrt(pow(P[i].Pos[1] - 10., 2.) + pow(P[i].Pos[0] - 10., 2.));
	      SphP[i].a.HydroAccel[1] -= SphP[i].a.HydroAccel[1] * (P[i].Pos[1] - 10.)
		/ sqrt(pow(P[i].Pos[1] - 10., 2.) + pow(P[i].Pos[0] - 10., 2.));
	      SphP[i].a.HydroAccel[2] = 0;
#else
	      for(k = 0; k < 3; k++)
		SphP[i].a.HydroAccel[k] = 0;
#endif
	    }
#endif
	}
    }


#if defined(CS_MODEL) && defined(CS_FEEDBACK)
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0 && (SphP[i].TempPromotion > 0 || SphP[i].DensPromotion > 0))
      {
	SphP[i].TempPromotion = 0;
	SphP[i].DensPromotion = 0;
      }
#endif
#if defined(HEALPIX)
  MPI_Allreduce(&count, &total_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  count = 0;
  MPI_Allreduce(&count2, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(total_count > 0)
    {
      if(ThisTask == 0)
	printf(" hey %i (%i) particles where freeezed and limit is %f \n", total_count, count,
	       (float) All.TotN_gas / 1000.0);
      if(total_count * 1000.0 > All.TotN_gas)	/*//for normal resolution ~100 */
	{
	  if(ThisTask == 0)
	    printf(" Next calculation of Healpix\n");
	  healpix_halo_cond(All.healpixmap);

	}
      total_count = 0;
      count2 = 0;
      fflush(stdout);
    }
#endif

#ifdef NUCLEAR_NETWORK
  if(ThisTask == 0)
    {
      printf("Doing nuclear network.\n");
    }
  MPI_Barrier(MPI_COMM_WORLD);

  tstart = second();

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	/* evaluate network here, but do it only for high enough temperatures */
	if(SphP[i].temp > All.NetworkTempThreshold)
	  {
	    nuc_particles++;
	    network_integrate(SphP[i].temp, SphP[i].d.Density * All.UnitDensity_in_cgs, SphP[i].xnuc,
			      SphP[i].dxnuc,
			      (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval *
			      All.UnitTime_in_s, &dedt_nuc, NULL, &All.nd, &All.nw);
	    SphP[i].e.DtEntropy += dedt_nuc * All.UnitEnergy_in_cgs / All.UnitTime_in_s;
	  }
	else
	  {
	    for(k = 0; k < EOS_NSPECIES; k++)
	      {
		SphP[i].dxnuc[k] = 0;
	      }
	  }
      }

  tend = second();
  timenetwork += timediff(tstart, tend);

  MPI_Allreduce(&nuc_particles, &nuc_particles_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      printf("Nuclear network done for %d particles.\n", nuc_particles_sum);
    }

  timewait1 += timediff(tend, second());
#endif

#ifdef RT_RAD_PRESSURE
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;

	if(All.Time != All.TimeBegin)
	  for(k = 0; k < 3; k++)
	    {
	      SphP[i].RadAccel[k] = 0.0;

	      for(j = 0; j < N_BINS; j++)
		SphP[i].RadAccel[k] += SphP[i].n_gamma[j] * nu[j];

	      SphP[i].RadAccel[k] *= SphP[i].n[k] / P[i].Mass * ELECTRONVOLT_IN_ERGS /
		All.UnitEnergy_in_cgs * All.HubbleParam /
		(C / All.UnitVelocity_in_cm_per_s) / dt / SphP[i].d.Density;
	    }
      }
#endif

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDNETWORK] += timenetwork;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}




/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*/
int hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		   int *ngblist)
{
  int startnode, numngb, listindex = 0;
  int j, k, n, l;


  double hinv, hinv3, hinv4, r2, u, visc;

  struct kernel_hydra kernel;
  struct hydrodata_in local;
  struct hydrodata_out out;
  memset(&out, 0, sizeof(struct hydrodata_out));

  if(mode == 0)
    particle2in_hydra(&local, target);
  else
    local = HydroDataGet[target];

  double p_over_rho2_i = (double) local.Pressure / ((double) local.Density * (double) local.Density);

#ifndef EOS_DEGENERATE
  kernel.sound_i = sqrt(GAMMA * p_over_rho2_i * local.Density);
#else
  kernel.sound_i = sqrt(local.dpdr);
#endif

  kernel.h_i = local.Hsml;

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  out.MinViscousDt = 1.0e32;
#else
  out.MaxSignalVel = kernel.sound_i;
#endif

#ifdef MAGNETIC
#ifndef EULERPOTENTIALS
  double mf_Ind = 1.0 / local.Density;
#endif
#ifdef MAGFORCE
  double mf_i = MU0_1 * pow(atime, 3 * GAMMA - 2) / (local.Density * local.Density);
  double mf_j = MU0_1 * pow(atime, 3 * GAMMA - 2);
#endif
#ifdef DIVBFORCE4
  double magacc[3], magcorr[3], norm_acc, norm_corr;
#endif
#if defined(EULER_DISSIPATION) || defined(MAGNETIC_DISSIPATION) || defined(MAGNETIC_DIFFUSION) || defined(FS_ETA2_DYN)
  double mf_dissInd = local.Density * atime * atime;
#endif
#if defined(EULER_DISSIPATION_HEAT) || defined(MAGNETIC_DISSIPATION) || defined(MAGNETIC_DIFFUSION_HEAT)
  double mf_dissEnt = 0.5 * MU0_1 * atime * atime;
#endif
#if defined(MAGFORCE) || defined(MAGNETIC_SIGNALVEL)
  kernel.b2_i =
    local.BPred[0] * local.BPred[0] + local.BPred[1] * local.BPred[1] + local.BPred[2] * local.BPred[2];
#endif
#ifdef MAGFORCE
  double mm_i[3][3], mm_j[3][3];
  for(k = 0; k < 3; k++)
    {
      for(l = 0; l < 3; l++)
	mm_i[k][l] = local.BPred[k] * local.BPred[l];
    }
  for(k = 0; k < 3; k++)
    mm_i[k][k] -= 0.5 * kernel.b2_i;
#endif
#if defined(MAGNETIC_SIGNALVEL)
  kernel.alfven2_i = kernel.b2_i * MU0_1 / local.Density;
#ifdef ALFVEN_VEL_LIMITER
  kernel.alfven2_i = DMIN(kernel.alfven2_i, ALFVEN_VEL_LIMITER * kernel.sound_i * kernel.sound_i);
#endif
#endif
#ifdef MAGNETIC_SIGNALVEL
  double vcsa2_i = kernel.sound_i * kernel.sound_i + kernel.alfven2_i;
#endif
#endif /* MAGNETIC */

  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
#ifndef DONOTUSENODELIST
      startnode = HydroDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
#else
      startnode = All.MaxPart;	/* root node */
#endif
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifdef CS_MODEL
	  numngb =
	    cs_ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, local.DensityNow,
					  local.Entropy, &local.Vel[0], mode, exportflag, exportnodecount,
					  exportindex, ngblist);
#else
	  numngb =
	    ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
				       exportnodecount, exportindex, ngblist);
#endif

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];

#if defined(BLACK_HOLES) || defined(CA_BH_ACCRETION)
	      if(P[j].Mass == 0)
		continue;
#endif

#ifdef NOWINDTIMESTEPPING
#ifdef WINDS
	      if(P[j].Type == 0)
		if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
		  continue;
#endif
#endif

	      kernel.dx = local.Pos[0] - P[j].Pos[0];
	      kernel.dy = local.Pos[1] - P[j].Pos[1];
	      kernel.dz = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      kernel.dx = NEAREST_X(kernel.dx);
	      kernel.dy = NEAREST_Y(kernel.dy);
	      kernel.dz = NEAREST_Z(kernel.dz);
#endif

#ifdef AB_SHOCK
	      hydra_evaluate_shockfinder(&local, &out, &kernel, j);
#endif

	      r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;
	      kernel.h_j = PPP[j].Hsml;

	      if(r2 < kernel.h_i * kernel.h_i || r2 < kernel.h_j * kernel.h_j)
		{
		  kernel.r = sqrt(r2);
		  if(kernel.r > 0)
		    {
		      double p_over_rho2_j =
			(double) SphP[j].Pressure / ((double) SphP[j].d.Density * (double) SphP[j].d.Density);
#ifndef EOS_DEGENERATE
		      kernel.sound_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].d.Density);
#else
		      kernel.sound_j = sqrt(SphP[j].dpdr);
#endif
		      kernel.dvx = local.Vel[0] - SphP[j].VelPred[0];
		      kernel.dvy = local.Vel[1] - SphP[j].VelPred[1];
		      kernel.dvz = local.Vel[2] - SphP[j].VelPred[2];
		      kernel.vdotr2 =
			kernel.dx * kernel.dvx + kernel.dy * kernel.dvy + kernel.dz * kernel.dvz;
		      kernel.rho_ij_inv = 2.0 / (local.Density + SphP[j].d.Density);

		      if(All.ComovingIntegrationOn)
			kernel.vdotr2 += hubble_a2 * r2;

		      if(kernel.r < kernel.h_i)
			{
			  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
			  u = kernel.r * hinv;
			  kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, 1);
			}
		      else
			{
			  kernel.dwk_i = 0;
			  kernel.wk_i = 0;
			}

		      if(kernel.r < kernel.h_j)
			{
			  kernel_hinv(kernel.h_j, &hinv, &hinv3, &hinv4);
			  u = kernel.r * hinv;
			  kernel_main(u, hinv3, hinv4, &kernel.wk_j, &kernel.dwk_j, 1);
			}
		      else
			{
			  kernel.dwk_j = 0;
			  kernel.wk_j = 0;
			}

		      kernel.dwk_ij = 0.5 * (kernel.dwk_i + kernel.dwk_j);

#ifndef TRADITIONAL_SPH_FORMULATION
		      kernel.dwk_i *= local.DhsmlDensityFactor;
		      kernel.dwk_j *= SphP[j].h.DhsmlDensityFactor;
#endif

#ifdef JD_VTURB
		      if(kernel.h_i >= kernel.r)	/* Make sure j is inside targets hsml */
			out.Vrms +=
			  (SphP[j].VelPred[0] - local.Vbulk[0]) * (SphP[j].VelPred[0] - local.Vbulk[0]) +
			  (SphP[j].VelPred[1] - local.Vbulk[1]) * (SphP[j].VelPred[1] - local.Vbulk[1]) +
			  (SphP[j].VelPred[2] - local.Vbulk[2]) * (SphP[j].VelPred[2] - local.Vbulk[2]);
#endif

#ifdef MAGNETIC_SN_SEEDING
		      if(kernel.h_j >= kernel.r)
			{
			  double ms_r =
			    kernel.dx * SphP[j].MagSeed[0] + kernel.dy * SphP[j].MagSeed[1] +
			    kernel.dz * SphP[j].MagSeed[2];
			  double soft_r = All.SnSeedSoftening * kernel.h_j;
			  out.DtB[0] +=
			    (3.0 * ms_r * kernel.dx / pow(kernel.r, 2) -
			     SphP[j].MagSeed[0]) / (pow(kernel.r, 3) + pow(soft_r, 3));
			  out.DtB[1] +=
			    (3.0 * ms_r * kernel.dy / pow(kernel.r, 2) -
			     SphP[j].MagSeed[1]) / (pow(kernel.r, 3) + pow(soft_r, 3));
			  out.DtB[2] +=
			    (3.0 * ms_r * kernel.dz / pow(kernel.r, 2) -
			     SphP[j].MagSeed[2]) / (pow(kernel.r, 3) + pow(soft_r, 3));
			}
#endif

#if defined(JD_VTURB) && defined(JD_DECOMPOSE_VTURB)
		      MyFloat vturb[3], e_rad[3], l_rad, vrad[3], vtan[3];

		      if(kernel.h_i >= kernel.r)
			{	/* Make sure j is inside targets hsml */

			  /* remove laminar velocity of j's velocity */
			  vturb[0] = SphP[j].VelPred[0] - local.Vbulk[0];
			  vturb[1] = SphP[j].VelPred[1] - local.Vbulk[1];
			  vturb[2] = SphP[j].VelPred[2] - local.Vbulk[2];

			  /* unit vector from i to j */
			  e_rad[0] = kernel.dx / kernel.r;
			  e_rad[1] = kernel.dy / kernel.r;
			  e_rad[2] = kernel.dz / kernel.r;

			  /* project vturb onto radial unit vector */
			  l_rad = vturb[0] * e_rad[0] + vturb[1] * e_rad[1] + vturb[2] * e_rad[2];

			  /* radial component of vturb */
			  vrad[0] = e_rad[0] * l_rad;
			  vrad[1] = e_rad[1] * l_rad;
			  vrad[2] = e_rad[2] * l_rad;

			  /* tangential component of vturb */
			  vtan[0] = vturb[0] - vrad[0];
			  vtan[1] = vturb[1] - vrad[1];
			  vtan[2] = vturb[2] - vrad[2];

			  /* We need dispersion in the end */
			  out.Vrad = l_rad;
			  out.Vtan = vtan[0] * vtan[0] + vtan[1] * vtan[1] + vtan[2] * vtan[2];
			}
#endif /* JD_VTURB && JD_DECOMPOSE_VTURB */

#ifdef MAGNETIC
		      kernel.mj_r = P[j].Mass / kernel.r;
#ifndef EULERPOTENTIALS
#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
			  kernel.mf_Ind = 0;
			else
#endif
			  kernel.mf_Ind = mf_Ind * kernel.mj_r * kernel.dwk_i;
#endif
#ifdef MAGFORCE
		      kernel.mf_i = mf_i * kernel.mj_r * kernel.dwk_i;
		      kernel.mf_j =
			mf_j * kernel.mj_r * kernel.dwk_j / (SphP[j].d.Density * SphP[j].d.Density);
#endif
#if defined(EULER_DISSIPATION) || defined(MAGNETIC_DISSIPATION) || defined(MAGNETIC_DIFFUSION) || defined(FS_ETA2_DYN)
		      double dissfac = kernel.mj_r * kernel.dwk_ij * kernel.rho_ij_inv * kernel.rho_ij_inv;
#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
			  kernel.mf_dissInd = 0;
			else
#endif
			  kernel.mf_dissInd = mf_dissInd * dissfac;

#endif
#if defined(EULER_DISSIPATION_HEAT) || defined(MAGNETIC_DISSIPATION) || defined(MAGNETIC_DIFFUSION_HEAT)
		      kernel.mf_dissEnt = mf_dissEnt * dissfac;
#endif

#if (!defined(EULERPOTENTIALS))
#ifdef HIGH_ORDER_INDUCTION
		      double ddvx =
			kernel.dx * local.Xix[0] + kernel.dy * local.Xix[1] + kernel.dz * local.Xix[2];
		      double ddvy =
			kernel.dx * local.Xiy[0] + kernel.dy * local.Xiy[1] + kernel.dz * local.Xiy[2];
		      double ddvz =
			kernel.dx * local.Xiz[0] + kernel.dy * local.Xiz[1] + kernel.dz * local.Xiz[2];

		      out.DtB[0] +=
			kernel.mf_Ind *
			((local.BPred[0] * (kernel.dvy + ddvy) -
			  local.BPred[1] * (kernel.dvx + ddvx)) * kernel.dy + (local.BPred[0] * (kernel.dvz +
												 ddvz) -
									       local.BPred[2] * (kernel.dvx +
												 ddvx)) *
			 kernel.dz);
		      out.DtB[1] +=
			kernel.mf_Ind *
			((local.BPred[1] * (kernel.dvz + ddvz) -
			  local.BPred[2] * (kernel.dvy + ddvy)) * kernel.dz + (local.BPred[1] * (kernel.dvx +
												 ddvx) -
									       local.BPred[0] * (kernel.dvy +
												 ddvy)) *
			 kernel.dx);
		      out.DtB[2] +=
			kernel.mf_Ind *
			((local.BPred[2] * (kernel.dvx + ddvx) -
			  local.BPred[0] * (kernel.dvz + ddvz)) * kernel.dx + (local.BPred[2] * (kernel.dvy +
												 ddvy) -
									       local.BPred[1] * (kernel.dvz +
												 ddvz)) *
			 kernel.dy);
#else
		      out.DtB[0] +=
			kernel.mf_Ind * ((local.BPred[0] * kernel.dvy - local.BPred[1] * kernel.dvx) *
					 kernel.dy + (local.BPred[0] * kernel.dvz -
						      local.BPred[2] * kernel.dvx) * kernel.dz);
		      out.DtB[1] +=
			kernel.mf_Ind * ((local.BPred[1] * kernel.dvz - local.BPred[2] * kernel.dvy) *
					 kernel.dz + (local.BPred[1] * kernel.dvx -
						      local.BPred[0] * kernel.dvy) * kernel.dx);
		      out.DtB[2] +=
			kernel.mf_Ind * ((local.BPred[2] * kernel.dvx - local.BPred[0] * kernel.dvz) *
					 kernel.dx + (local.BPred[2] * kernel.dvy -
						      local.BPred[1] * kernel.dvz) * kernel.dy);
#endif /* HIGH_ORDER */
#endif /* EULER */

#ifdef DIVBCLEANING_DEDNER
		      double phifac = kernel.mf_Ind * atime * atime;
#ifdef SMOOTH_PHI
		      phifac *= (local.PhiPred - SphP[j].SmoothPhi);
#else
		      phifac *= (local.PhiPred - SphP[j].PhiPred);
#endif
#ifndef FS_DEDNER_DIRECT
		      out.GradPhi[0] += phifac * kernel.dx;
		      out.GradPhi[1] += phifac * kernel.dy;
		      out.GradPhi[2] += phifac * kernel.dz;
#else
		      out.DtB[0] += phifac * kernel.dx;
		      out.DtB[1] += phifac * kernel.dy;
		      out.DtB[2] += phifac * kernel.dz;
		      out.DtEntropy -=
			phifac * (local.BPred[0] * kernel.dx + local.BPred[1] * kernel.dy +
				  local.BPred[2] * kernel.dz) * MU0_1 / local.Density;
#endif
#endif

#if defined(MAGFORCE) || defined(MAGNETIC_DISSIPATION) || defined(MAGNETIC_DIFFUSION) || defined(FS_ALFA2_DYN) || defined(MAGNETIC_SIGNALVEL)
		      double bpred_j[3];
		      bpred_j[0] = SphP[j].b2.BPred[0];
		      bpred_j[1] = SphP[j].b2.BPred[1];
		      bpred_j[2] = SphP[j].b2.BPred[2];
#ifdef SFR
		      bpred_j[0] *= pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
		      bpred_j[1] *= pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
		      bpred_j[2] *= pow(1. - SphP[j].XColdCloud, 2. * POW_CC);
#endif
		      double dBx = local.BPred[0] - bpred_j[0];
		      double dBy = local.BPred[1] - bpred_j[1];
		      double dBz = local.BPred[2] - bpred_j[2];
#endif
#if defined(FS_ALFA2_DYN) || defined(FS_ETA2_DYN)
		      double dr[3], Alfa_i[3], Alfa_j[3], dB[3];
		      dr[0] = kernel.dx;
		      dr[1] = kernel.dy;
		      dr[2] = kernel.dz;
		      dB[0] = dBx;
		      dB[1] = dBy;
		      dB[2] = dBz;
#if defined(FS_ALFA2_DYN)
#ifdef FS_DYN_DEF
		      if((P[j].Pos[0] - LONG_Y / 2.) * (P[j].Pos[0] - LONG_Y / 2.) +
			 (P[j].Pos[1] - LONG_X / 2) * (P[j].Pos[1] - LONG_X / 2.) < 16)
			{
			  SphP[j].alfa2 =
			    All.Alfa0Quench / (SysState.EnergyMag / All.AlfaB2eq +
					       1.) * sin((fabs(P[j].Pos[2] - LONG_Z / 2.) -
							  LONG_Z / 2.) * M_PI / All.Alfa2H) * OmegaR(j,
												     BRANDT_MODE);
			  SphP[j].alfa2 = (fabs(P[j].Pos[2] - LONG_Z / 2.) > All.Alfa2H) ? SphP[j].alfa2 =
			    0. : (P[j].Pos[2] - LONG_Z / 2.) <
			    0.0 ? -fabs(SphP[j].alfa2) : fabs(SphP[j].alfa2);
			}
		      else
			SphP[j].alfa2 = 0;
#ifdef FS_FIX_ALPHA
		      SphP[j].alfa2 = FS_FIX_ALPHA;
#endif
#endif
		      for(k = 0; k < 3; k++)
			{
			  Alfa_i[k] = local.alfa2;
			  Alfa_j[k] = SphP[j].alfa2;
			}
#ifdef FS_ALFA2_DIAG
		      Alfa_i[2] = 0.0;
		      Alfa_j[2] = 0.0;
#endif
		      for(k = 0; k < 3; k++)
			out.DtB[k] +=
			  kernel.mf_Ind * atime * atime * ((local.BPred[(k + 1) % 3] * Alfa_i[(k + 1) % 3] -
							    bpred_j[(k + 1) % 3] * Alfa_j[(k +
											   1) % 3]) * dr[(k +
													  2) %
													 3] -
							   (local.BPred[(k + 2) % 3] * Alfa_i[(k + 2) % 3] -
							    bpred_j[(k + 2) % 3] * Alfa_j[(k +
											   2) % 3]) * dr[(k +
													  1) %
													 3]);
#endif
#endif
#if defined(MAGFORCE) || defined(MAGNETIC_SIGNALVEL)
		      kernel.b2_j =
			bpred_j[0] * bpred_j[0] + bpred_j[1] * bpred_j[1] + bpred_j[2] * bpred_j[2];
#endif
#ifdef MAGFORCE
		      for(k = 0; k < 3; k++)
			{
			  for(l = 0; l < 3; l++)
			    mm_j[k][l] = bpred_j[k] * bpred_j[l];
			}
		      for(k = 0; k < 3; k++)
			mm_j[k][k] -= 0.5 * kernel.b2_j;

		      for(k = 0; k < 3; k++)
#ifdef DIVBFORCE3
			out.magacc[k] +=
#else
#ifdef DIVBFORCE4
			magacc[k] =
#else /* DIVBFORCE or normal behavior */
			out.Acc[k] +=
#endif
#endif
			  (mm_i[k][0] * kernel.mf_i + mm_j[k][0] * kernel.mf_j) * kernel.dx +
			  (mm_i[k][1] * kernel.mf_i + mm_j[k][1] * kernel.mf_j) * kernel.dy +
			  (mm_i[k][2] * kernel.mf_i + mm_j[k][2] * kernel.mf_j) * kernel.dz;
#if defined(DIVBFORCE) || defined(DIVBFORCE3) || defined(DIVBFORCE4)
		      for(k = 0; k < 3; k++)
#if defined(DIVBFORCE)
			out.Acc[k] -=
#else
#if defined(DIVBFORCE3)
			out.magcorr[k] +=
#else /* DIVBFORCE4 */
			magcorr[k] =
#endif
#endif
			  local.BPred[k] *
			  ((local.BPred[0] * kernel.mf_i + bpred_j[0] * kernel.mf_j) * kernel.dx +
			   (local.BPred[1] * kernel.mf_i + bpred_j[1] * kernel.mf_j) * kernel.dy +
			   (local.BPred[2] * kernel.mf_i + bpred_j[2] * kernel.mf_j) * kernel.dz);
#endif
#ifdef DIVBFORCE4
		      norm_acc =
			10. * sqrt(magacc[0] * magacc[0] + magacc[1] * magacc[1] + magacc[2] * magacc[2]);
		      norm_corr =
			sqrt(magcorr[0] * magcorr[0] + magcorr[1] * magcorr[1] + magcorr[2] * magcorr[2]);
		      norm_corr = norm_corr < norm_acc ? 1.0 : norm_acc / norm_corr;
		      for(k = 0; k < 3; k++)
			out.Acc[k] += (magacc[k] - norm_corr * magcorr[k]);
#endif
#endif /* end MAG FORCE   */
#endif /* end of MAGNETIC */


#if defined(MAGNETIC_SIGNALVEL)
		      kernel.alfven2_j = kernel.b2_j * MU0_1 / SphP[j].d.Density;
#ifdef ALFVEN_VEL_LIMITER
		      kernel.alfven2_j =
			DMIN(kernel.alfven2_j, ALFVEN_VEL_LIMITER * kernel.sound_j * kernel.sound_j);
#endif
#else
		      kernel.vsig = kernel.sound_i + kernel.sound_j;
#endif
#ifdef MAGNETIC_SIGNALVEL
		      double vcsa2_j = kernel.sound_j * kernel.sound_j + kernel.alfven2_j;

		      double Bpro2_j =
			(bpred_j[0] * kernel.dx + bpred_j[1] * kernel.dy + bpred_j[2] * kernel.dz) / kernel.r;

		      Bpro2_j *= Bpro2_j;
		      double magneticspeed_j = sqrt(0.5 * (vcsa2_j + sqrt(DMAX((vcsa2_j * vcsa2_j -
										4 * kernel.sound_j *
										kernel.sound_j * Bpro2_j *
										MU0_1 / SphP[j].d.Density),
									       0))));
		      double Bpro2_i =
			(local.BPred[0] * kernel.dx + local.BPred[1] * kernel.dy +
			 local.BPred[2] * kernel.dz) / kernel.r;
		      Bpro2_i *= Bpro2_i;
		      double magneticspeed_i = sqrt(0.5 * (vcsa2_i + sqrt(DMAX((vcsa2_i * vcsa2_i -
										4 * kernel.sound_i *
										kernel.sound_i * Bpro2_i *
										MU0_1 / local.Density), 0))));
		      kernel.vsig = magneticspeed_i + magneticspeed_j;
#endif

#ifndef ALTERNATIVE_VISCOUS_TIMESTEP
		      if(kernel.vsig > out.MaxSignalVel)
			out.MaxSignalVel = kernel.vsig;
#endif

		      if(kernel.vdotr2 < 0)	/* ... artificial viscosity */
			{
#ifndef ALTVISCOSITY
#ifndef CONVENTIONAL_VISCOSITY
			  double mu_ij = fac_mu * kernel.vdotr2 / kernel.r;	/* note: this is negative! */
#else
			  double c_ij = 0.5 * (kernel.sound_i + kernel.sound_j);
			  double h_ij = 0.5 * (kernel.h_i + kernel.h_j);
			  double mu_ij = fac_mu * h_ij * kernel.vdotr2 / (r2 + 0.0001 * h_ij * h_ij);
#endif
#ifdef MAGNETIC
			  kernel.vsig -= 1.5 * mu_ij;
#else
			  kernel.vsig -= 3 * mu_ij;
#endif


#ifndef ALTERNATIVE_VISCOUS_TIMESTEP
			  if(kernel.vsig > out.MaxSignalVel)
			    out.MaxSignalVel = kernel.vsig;
#endif

#ifdef NO_SHEAR_VISCOSITY_LIMITER
			  double f1, f2;
			  f1 = f2 = 1;
#else
			  double f1 = local.F1;
#ifndef NAVIERSTOKES
			  double f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].r.CurlVel +
						      0.0001 * kernel.sound_j / fac_mu / kernel.h_j);
#else
			  double f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].u.s.CurlVel +
						      0.0001 * kernel.sound_j / fac_mu / kernel.h_j);
#endif
#endif

#ifdef TIME_DEP_ART_VISC
			  double BulkVisc_ij = 0.5 * (local.alpha + SphP[j].alpha);
#else
			  double BulkVisc_ij = All.ArtBulkViscConst;
#endif
#ifndef CONVENTIONAL_VISCOSITY
			  visc = 0.25 * BulkVisc_ij * kernel.vsig * (-mu_ij) * kernel.rho_ij_inv * (f1 + f2);
#else
			  visc =
			    (-BulkVisc_ij * mu_ij * c_ij + 2 * BulkVisc_ij * mu_ij * mu_ij) *
			    kernel.rho_ij_inv * (f1 + f2) * 0.5;
#endif

#ifdef JD_TIME_DEP_ART_VISC
			  visc *= (kernel.sound_i + kernel.sound_j - 3 * mu_ij) / kernel.vsig;
#endif

#else /* start of ALTVISCOSITY block */
			  double mu_i;
			  if(f1 < 0)
			    mu_i = h_i * fabs(f1);	/* local.F1 hold here the velocity divergence of particle i */
			  else
			    mu_i = 0;

			  double mu_j;
			  if(SphP[j].v.DivVel < 0)
			    mu_j = h_j * fabs(SphP[j].v.DivVel);
			  else
			    mu_j = 0;
			  double visc =
			    All.ArtBulkViscConst * ((kernel.sound_i + mu_i) * mu_i / local.Density +
						    (kernel.sound_j + mu_j) * mu_j / SphP[j].d.Density);
#endif /* end of ALTVISCOSITY block */


			  /* .... end artificial viscosity evaluation */
			  /* now make sure that viscous acceleration is not too large */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
			  if(visc > 0)
			    {
			      double dt = fac_vsic_fix * kernel.vdotr2 /
				((local.Mass + P[j].Mass) * kernel.dwk_ij * kernel.r * visc);

			      dt /= hubble_a;

			      if(dt < out.MinViscousDt)
				out.MinViscousDt = dt;
			    }
#endif

#ifndef NOVISCOSITYLIMITER
			  double dt = 2 * IMAX(local.Timestep,
					       (P[j].TimeBin ? (1 << P[j].TimeBin) : 0)) *
			    All.Timebase_interval;
			  if(dt > 0 && kernel.dwk_ij < 0)
			    {
#ifdef BLACK_HOLES
			      if((local.Mass + P[j].Mass) > 0)
#endif
				visc = DMIN(visc, 0.5 * fac_vsic_fix * kernel.vdotr2 /
					    ((local.Mass + P[j].Mass) * kernel.dwk_ij * kernel.r * dt));
			    }
#endif
			}
		      else
			{
			  visc = 0;
			}
#ifdef TRADITIONAL_SPH_FORMULATION
		      double hfc_visc = P[j].Mass * visc * kernel.dwk_ij / kernel.r;

		      double hfc = hfc_visc +
			P[j].Mass * kernel.dwk_ij / kernel.r * (p_over_rho2_i + p_over_rho2_j);

		      /* hfc_egy = 0.5 * P[j].Mass * (kernel.dwk_i + kernel.dwk_j) / kernel.r * (p_over_rho2_i + p_over_rho2_j); */
		      double hfc_egy = P[j].Mass * kernel.dwk_ij / kernel.r * (p_over_rho2_i);
#else
		      double hfc_visc = P[j].Mass * visc * kernel.dwk_ij / kernel.r;
		      /* Formulation derived from the Lagrangian */
		      double hfc =
			hfc_visc + P[j].Mass * (p_over_rho2_i * kernel.dwk_i +
						p_over_rho2_j * kernel.dwk_j) / kernel.r;
#endif

#ifdef EOS_DEGENERATE
		      double hfc_egy = P[j].Mass * kernel.dwk_ij / kernel.r * (p_over_rho2_i);
#endif

#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No force by wind particles */
			  {
			    hfc = hfc_visc = 0;
			  }
#endif

#ifndef NOACCEL
		      out.Acc[0] += FLT(-hfc * kernel.dx);
		      out.Acc[1] += FLT(-hfc * kernel.dy);
		      out.Acc[2] += FLT(-hfc * kernel.dz);
#endif

#if !defined(EOS_DEGENERATE) && !defined(TRADITIONAL_SPH_FORMULATION)
		      out.DtEntropy += FLT(0.5 * hfc_visc * kernel.vdotr2);
#else

#ifdef TRADITIONAL_SPH_FORMULATION
		      out.DtEntropy += FLT(0.5 * (hfc_visc + hfc_egy) * kernel.vdotr2);
#else
		      /* EOS_DEGENERATE */
		      out.DtEntropy += FLT((0.5 * hfc_visc + hfc_egy) * kernel.vdotr2);
#endif
#endif

#ifdef NAVIERSTOKES
		      double faci =
			local.Mass * local.shear_viscosity / (local.Density * local.Density) * kernel.dwk_i /
			kernel.r;

#ifndef NAVIERSTOKES_CONSTANT
		      faci *= pow((local.Entropy * pow(local.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);	/*multiplied by E^5/2 */
#endif
		      double facj = P[j].Mass * get_shear_viscosity(j) /
			(SphP[j].d.Density * SphP[j].d.Density) * kernel.dwk_j / kernel.r;

#ifndef NAVIERSTOKES_CONSTANT
		      facj *= pow((SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);	/*multiplied by E^5/2 */
#endif

#ifdef NAVIERSTOKES_BULK
		      double facbi =
			local.Mass * All.NavierStokes_BulkViscosity / (local.Density * local.Density) *
			kernel.dwk_i / kernel.r;
		      double facbj =
			P[j].Mass * All.NavierStokes_BulkViscosity / (SphP[j].d.Density * SphP[j].d.Density) *
			kernel.dwk_j / kernel.r;
#endif

#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No visc for wind particles */
			  {
			    faci = facj = 0;
#ifdef NAVIERSTOKES_BULK
			    facbi = facbj = 0;
#endif
			  }
#endif

#ifdef VISCOSITY_SATURATION
		      double IonMeanFreePath_i = All.IonMeanFreePath * pow((local.Entropy * pow(local.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / local.Density;	/* u^2/local.Density */
		      double IonMeanFreePath_j = All.IonMeanFreePath * pow((SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / SphP[j].d.Density;	/* u^2/local.Density */
		      double VelLengthScale_i, VelLengthScale_j;

		      for(k = 0, VelLengthScale_i = 0, VelLengthScale_j = 0; k < 3; k++)
			{
			  if(fabs(local.stressdiag[k]) > 0)
			    {
			      VelLengthScale_i = 2 * kernel.sound_i / fabs(local.stressdiag[k]);

			      if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
				{
				  local.stressdiag[k] =
				    local.stressdiag[k] * (VelLengthScale_i / IonMeanFreePath_i);

				}
			    }
			  if(fabs(SphP[j].u.s.StressDiag[k]) > 0)
			    {
			      VelLengthScale_j = 2 * kernel.sound_j / fabs(SphP[j].u.s.StressDiag[k]);

			      if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
				{
				  SphP[j].u.s.StressDiag[k] = SphP[j].u.s.StressDiag[k] *
				    (VelLengthScale_j / IonMeanFreePath_j);

				}
			    }
			  if(fabs(local.stressoffdiag[k]) > 0)
			    {
			      VelLengthScale_i = 2 * kernel.sound_i / fabs(local.stressoffdiag[k]);

			      if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
				{
				  local.stressoffdiag[k] =
				    local.stressoffdiag[k] * (VelLengthScale_i / IonMeanFreePath_i);
				}
			    }
			  if(fabs(SphP[j].u.s.StressOffDiag[k]) > 0)
			    {
			      VelLengthScale_j = 2 * kernel.sound_j / fabs(SphP[j].u.s.StressOffDiag[k]);

			      if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
				{
				  SphP[j].u.s.StressOffDiag[k] = SphP[j].u.s.StressOffDiag[k] *
				    (VelLengthScale_j / IonMeanFreePath_j);
				}
			    }
			}
#endif

		      /* Acceleration due to the shear viscosity */
		      out.Acc[0] +=
			faci * (local.stressdiag[0] * kernel.dx + local.stressoffdiag[0] * kernel.dy +
				local.stressoffdiag[1] * kernel.dz) +
			facj * (SphP[j].u.s.StressDiag[0] * kernel.dx +
				SphP[j].u.s.StressOffDiag[0] * kernel.dy +
				SphP[j].u.s.StressOffDiag[1] * kernel.dz);

		      out.Acc[1] +=
			faci * (local.stressoffdiag[0] * kernel.dx + local.stressdiag[1] * kernel.dy +
				local.stressoffdiag[2] * kernel.dz) +
			facj * (SphP[j].u.s.StressOffDiag[0] * kernel.dx +
				SphP[j].u.s.StressDiag[1] * kernel.dy +
				SphP[j].u.s.StressOffDiag[2] * kernel.dz);

		      out.Acc[2] +=
			faci * (local.stressoffdiag[1] * kernel.dx + local.stressoffdiag[2] * kernel.dy +
				local.stressdiag[2] * kernel.dz) +
			facj * (SphP[j].u.s.StressOffDiag[1] * kernel.dx +
				SphP[j].u.s.StressOffDiag[2] * kernel.dy +
				SphP[j].u.s.StressDiag[2] * kernel.dz);

		      /*Acceleration due to the bulk viscosity */
#ifdef NAVIERSTOKES_BULK
#ifdef VISCOSITY_SATURATION
		      VelLengthScale_i = 0;
		      VelLengthScale_j = 0;

		      if(fabs(local.divvel) > 0)
			{
			  VelLengthScale_i = 3 * kernel.sound_i / fabs(local.divvel);

			  if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
			    {
			      local.divvel = local.divvel * (VelLengthScale_i / IonMeanFreePath_i);
			    }
			}

		      if(fabs(SphP[j].u.s.a4.DivVel) > 0)
			{
			  VelLengthScale_j = 3 * kernel.sound_j / fabs(SphP[j].u.s.a4.DivVel);

			  if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
			    {
			      SphP[j].u.s.a4.DivVel = SphP[j].u.s.a4.DivVel *
				(VelLengthScale_j / IonMeanFreePath_j);

			    }
			}
#endif

		      out.Acc[0] +=
			facbi * local.divvel * kernel.dx + facbj * SphP[j].u.s.a4.DivVel * kernel.dx;
		      out.Acc[1] +=
			facbi * local.divvel * kernel.dy + facbj * SphP[j].u.s.a4.DivVel * kernel.dy;
		      out.Acc[2] +=
			facbi * local.divvel * kernel.dz + facbj * SphP[j].u.s.a4.DivVel * kernel.dz;
#endif
#endif /* end NAVIERSTOKES */

#ifdef ARTIFICIAL_CONDUCTIVITY
#ifdef BP_REAL_CRs
		      double vsigu =
			sqrt(fabs(local.Pressure - SphP[j].Pressure - local.CRpPressure + SphP[j].CRpPressure)
			     * kernel.rho_ij_inv);
		      double u_i = (local.Pressure - local.CRpPressure) / (GAMMA_MINUS1 * local.Density);
		      double u_j =
			(SphP[j].Pressure - SphP[j].CRpPressure) / (GAMMA_MINUS1 * SphP[j].d.Density);
		      out.DtEntropy +=
			P[j].Mass * All.ArtCondConstant * vsigu * (u_i -
								   u_j) * kernel.rho_ij_inv * kernel.dwk_ij;
#else
		      double vsigu = sqrt(fabs(local.Pressure - SphP[j].Pressure) * kernel.rho_ij_inv);
		      double u_i = local.Pressure / (GAMMA_MINUS1 * local.Density);
		      double u_j = SphP[j].Pressure / (GAMMA_MINUS1 * SphP[j].d.Density);
		      out.DtEntropy +=
			P[j].Mass * All.ArtCondConstant * vsigu * (u_i -
								   u_j) * kernel.rho_ij_inv * kernel.dwk_ij;
#endif
#endif

#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
		      int Nbin;
		      double vsigu_cr =
			sqrt(fabs(local.CRpPressure - SphP[j].CRpPressure) * kernel.rho_ij_inv);
		      for(Nbin = 0; Nbin < BP_REAL_CRs; Nbin++)
			{
			  out.DtCRpE[Nbin] +=
			    P[j].Mass * All.CRsArtCondConstant * vsigu_cr * (local.CRpE[Nbin] -
									     SphP[j].CRpE[Nbin]) *
			    kernel.rho_ij_inv * kernel.dwk_ij;
			  out.DtCRpN[Nbin] +=
			    P[j].Mass * All.CRsArtCondConstant * vsigu_cr * (local.CRpN[Nbin] -
									     SphP[j].CRpN[Nbin]) *
			    kernel.rho_ij_inv * kernel.dwk_ij;
			}
#endif

#ifdef MAGNETIC
#if defined(MAGNETIC_DIFFUSION)
		      out.DtB[0] += 2.0 * All.MagneticEta * kernel.mf_dissInd * dBx;
		      out.DtB[1] += 2.0 * All.MagneticEta * kernel.mf_dissInd * dBy;
		      out.DtB[2] += 2.0 * All.MagneticEta * kernel.mf_dissInd * dBz;
#ifdef MAGNETIC_DIFFUSION_HEAT
		      out.DtEntropy -=
			2.0 * All.MagneticEta * kernel.mf_dissEnt * (dBx * dBx + dBy * dBy + dBz * dBz);
#endif
#endif
#if defined(FS_ETA2_DYN)
		      for(k = 0; k < 3; k++)
			out.DtB[k] += kernel.mf_dissInd * ((local.eta2 + SphP[j].eta2) * dB[k]
							   - (local.eta2 -
							      SphP[j].eta2) * (dr[(k + 1) % 3] * (dr[k] *
												  dB[(k +
												      1) %
												     3] -
												  dr[(k +
												      1) %
												     3] *
												  dB[k]) +
									       dr[(k +
										   2) % 3] * (dr[k] * dB[(k +
													  2) %
													 3] -
											      dr[(k +
												  2) % 3] *
											      dB[k])));
#endif
#ifdef EULER_DISSIPATION
		      out.DtEulerA +=
			All.ArtMagDispConst * 0.5 * kernel.vsig * (local.EulerA -
								   SphP[j].EulerA) * kernel.mf_dissInd *
			kernel.r / hubble_a2;
		      out.DtEulerB +=
			All.ArtMagDispConst * 0.5 * kernel.vsig * (local.EulerB -
								   SphP[j].EulerB) * kernel.mf_dissInd *
			kernel.r / hubble_a2;
#endif
#ifdef MAGNETIC_DISSIPATION
		      double vsigb = 0.5 * sqrt(kernel.alfven2_i + kernel.alfven2_j);
#ifdef TIME_DEP_MAGN_DISP
		      double eta = 0.5 * (local.Balpha + SphP[j].Balpha) * vsigb * kernel.r;
#else
		      double eta = All.ArtMagDispConst * vsigb * kernel.r;
#endif
		      out.DtEntropy -= eta * kernel.mf_dissEnt * (dBx * dBx + dBy * dBy + dBz * dBz);
		      out.DtB[0] += eta * kernel.mf_dissInd * dBx;
		      out.DtB[1] += eta * kernel.mf_dissInd * dBy;
		      out.DtB[2] += eta * kernel.mf_dissInd * dBz;
#endif
#endif

#ifdef LT_MV_CHEMICALDIFFUSION_WAKEUP
		      double myD, maxgrad;
#ifdef LT_MV_CHEMICALDIFFUSION_CONSTANT_D
		      myD = All.ConstantD;
#else
		      myD = SphP[j].D;
#endif
		      maxgrad = 0;
		      for(k = 0; k < LT_NMetP; k++)
			if(maxgrad < fabs(SphP[j].GradAbnd[k]))
			  maxgrad = fabs(SphP[j].GradAbnd[k]);
		      Dvsig = myD * maxgrad + local.Dvsig;
#endif

#ifdef WAKEUP
		      if(kernel.vsig > WAKEUP * SphP[j].MaxSignalVel)
			SphP[j].wakeup = 1;
#endif
#ifdef LT_MV_CHEMICALDIFFUSION_WAKEUP
		      if(Dvsig > WAKEUP * SphP[j].MaxSignalVel)
			SphP[j].wakeup = 1;
#endif
		    }
		}
	    }
	}

#ifndef DONOTUSENODELIST
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = HydroDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
#endif
    }


  /* Now collect the result at the right place */
  if(mode == 0)
    out2particle_hydra(&out, target, 0);
  else
    HydroDataResult[target] = out;

  return 0;
}

void *hydro_evaluate_primary(void *p)
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

#ifdef BLACK_HOLES
      if(P[i].Type == 0 && P[i].Mass > 0)
#else
      if(P[i].Type == 0)
#endif
	{
	  if(hydro_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *hydro_evaluate_secondary(void *p)
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

      hydro_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}


#ifdef AB_SHOCK
static inline void hydra_evaluate_shockfinder(struct hydrodata_in *local, struct hydrodata_out *out,
					      struct kernel_hydra *kernel, int j)
{
  double mj_wki, u, hinv, hinv3, hinv4;

  /* Parallel to shock normal */
  double Shock_X1_r =
    kernel->dx * local->Shock_N1[0] + kernel->dy * local->Shock_N1[1] + kernel->dz * local->Shock_N1[2];
  double r2 = kernel->dx * kernel->dx + kernel->dy * kernel->dy + kernel->dz * kernel->dz;
  kernel_hinv(local->Hsml, &hinv, &hinv3, &hinv4);

  u = fabs(1.0 - fabs(Shock_X1_r * hinv));
  if(u < 1.0)
    {
      kernel_main(u, hinv3, hinv4, &kernel->wk_i, &kernel->dwk_i, -1);
      mj_wki = P[j].Mass * kernel->wk_i;

      u = sqrt(r2 - Shock_X1_r * Shock_X1_r) * hinv;
      kernel_main(u, hinv3, hinv4, &kernel->wk_i, &kernel->dwk_i, -1);
      mj_wki *= P[j].Mass * kernel->wk_i;

      if(Shock_X1_r < 0)
	{
	  out->Shock_Up_V1[0] += mj_wki * SphP[j].VelPred[0];
	  out->Shock_Up_V1[1] += mj_wki * SphP[j].VelPred[1];
	  out->Shock_Up_V1[2] += mj_wki * SphP[j].VelPred[2];
	  out->Shock_Up_Rho += mj_wki * SphP[j].d.Density;
	  out->Shock_Up_Pressure += mj_wki * SphP[j].Pressure;
	  out->Shock_Up_Weight += mj_wki;

#ifndef BP_REAL_CRs
	  double shock_up_signal2 = GAMMA * SphP[j].Pressure / SphP[j].d.Density;
#else
	  double shock_up_signal2 =
	    (GAMMA * SphP[j].Pressure + (All.CR_Gamma - GAMMA) * SphP[j].CRpPressure) / SphP[j].d.Density;
#endif

#if defined(AB_SHOCK_MAGNETIC) && defined(MAGNETIC)
	  double shockb2j;
	  shockb2j =
	    SphP[j].b2.BPred[0] * SphP[j].b2.BPred[0] + SphP[j].b2.BPred[1] * SphP[j].b2.BPred[1] +
	    SphP[j].b2.BPred[2] * SphP[j].b2.BPred[2];
	  out->Shock_Up_Signal +=
	    mj_wki * sqrt(0.5 * (shockb2j * MU0_1 / SphP[j].d.Density + shock_up_signal2));
#else
	  out->Shock_Up_Signal += mj_wki * sqrt(shock_up_signal2);
#endif

	}
      else
	{
	  out->Shock_Down_V1[0] += mj_wki * SphP[j].VelPred[0];
	  out->Shock_Down_V1[1] += mj_wki * SphP[j].VelPred[1];
	  out->Shock_Down_V1[2] += mj_wki * SphP[j].VelPred[2];
	  out->Shock_Down_Rho += mj_wki * SphP[j].d.Density;
	  out->Shock_Down_Pressure += mj_wki * SphP[j].Pressure;
	  out->Shock_Down_Weight += mj_wki;

	  double shock_down_signal2 = GAMMA * SphP[j].Pressure / SphP[j].d.Density;

#if defined(AB_SHOCK_MAGNETIC) && defined(MAGNETIC)
	  double shockb2j;
	  shockb2j =
	    SphP[j].b2.BPred[0] * SphP[j].b2.BPred[0] + SphP[j].b2.BPred[1] * SphP[j].b2.BPred[1] +
	    SphP[j].b2.BPred[2] * SphP[j].b2.BPred[2];
	  out->Shock_Down_Signal +=
	    mj_wki * sqrt(0.5 * (shockb2j * MU0_1 / SphP[j].d.Density + shock_down_signal2));
#else
	  out->Shock_Down_Signal += mj_wki * sqrt(shock_down_signal2);
#endif

	}
    }

#ifdef AB_SHOCK_VELDIV
  /* Perpendicular to shock normal */
  double Shock_X2_r =
    kernel->dx * local->Shock_N2[0] + kernel->dy * local->Shock_N2[1] + kernel->dz * local->Shock_N2[2];

  u = fabs(1.0 - fabs(Shock_X2_r * hinv));
  if(u < 1.0)
    {
      kernel_main(u, hinv3, hinv4, &kernel->wk_i, &kernel->dwk_i, -1);
      mj_wki = P[j].Mass * kernel->wk_i;

      u = sqrt(r2 - Shock_X2_r * Shock_X2_r) * hinv;
      kernel_main(u, hinv3, hinv4, &kernel->wk_i, &kernel->dwk_i, -1);
      mj_wki *= P[j].Mass * kernel->wk_i;

      if(Shock_X2_r < 0)
	{
	  out->Shock_Up_V2[0] += mj_wki * SphP[j].VelPred[0];
	  out->Shock_Up_V2[1] += mj_wki * SphP[j].VelPred[1];
	  out->Shock_Up_V2[2] += mj_wki * SphP[j].VelPred[2];
	  out->Shock_Up_Weight23[0] += mj_wki;
	}
      else
	{
	  out->Shock_Down_V2[0] += mj_wki * SphP[j].VelPred[0];
	  out->Shock_Down_V2[1] += mj_wki * SphP[j].VelPred[1];
	  out->Shock_Down_V2[2] += mj_wki * SphP[j].VelPred[2];
	  out->Shock_Down_Weight23[0] += mj_wki;
	}
    }

  /* Perpendicular to shock normal */
  double Shock_X3_r =
    kernel->dx * local->Shock_N3[0] + kernel->dy * local->Shock_N3[1] + kernel->dz * local->Shock_N3[2];

  u = fabs(1.0 - fabs(Shock_X3_r * hinv));
  if(u < 1.0)
    {
      kernel_main(u, hinv3, hinv4, &kernel->wk_i, &kernel->dwk_i, -1);
      mj_wki = P[j].Mass * kernel->wk_i;

      u = sqrt(r2 - Shock_X3_r * Shock_X3_r) * hinv;
      kernel_main(u, hinv3, hinv4, &kernel->wk_i, &kernel->dwk_i, -1);
      mj_wki *= P[j].Mass * kernel->wk_i;

      if(Shock_X3_r < 0)
	{
	  out->Shock_Up_V3[0] += mj_wki * SphP[j].VelPred[0];
	  out->Shock_Up_V3[1] += mj_wki * SphP[j].VelPred[1];
	  out->Shock_Up_V3[2] += mj_wki * SphP[j].VelPred[2];
	  out->Shock_Up_Weight23[1] += mj_wki;
	}
      else
	{
	  out->Shock_Down_V3[0] += mj_wki * SphP[j].VelPred[0];
	  out->Shock_Down_V3[1] += mj_wki * SphP[j].VelPred[1];
	  out->Shock_Down_V3[2] += mj_wki * SphP[j].VelPred[2];
	  out->Shock_Down_Weight23[1] += mj_wki;
	}
    }
#endif
}
#endif /* AB_SHOCK */
