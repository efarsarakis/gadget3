
/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#ifdef _OPENMP
#include <omp.h>
#endif


#include "gadgetconfig.h"


#ifdef CHEMCOOL
#include "chemcool_consts.h"
#include "f2c.h"
#endif

#ifdef MPISENDRECV_CHECKSUM
#define MPI_Sendrecv MPI_Check_Sendrecv
#endif

#ifdef MPISENDRECV_SIZELIMIT
#define MPI_Sendrecv MPI_Sizelimited_Sendrecv
#endif

#include "tags.h"
#if defined(CHEMISTRY) || defined (UM_CHEMISTRY)
#include "chemistry.h"
#endif
#include "assert.h"

#ifdef EOS_DEGENERATE
#include "helm_eos.h"
#endif

#ifdef NUCLEAR_NETWORK
#include "network.h"
#include "integrate.h"
#endif

#ifdef MYSORT
#ifdef OMP_MYSORT
#define MYSORT_DATAINDEX serial_sort_omp
#else
#define MYSORT_DATAINDEX mysort_dataindex
#endif
#else // MYSORT
#ifdef OMP_SORT
#define MYSORT_DATAINDEX omp_qsort
#else
#define MYSORT_DATAINDEX qsort
#endif
#endif

// compiler specific data alignment hints
// XLC compiler
#if defined(__xlC__)
#define ALIGN(n) __attribute__((__aligned__(n)))
// GNU compiler 
#elif defined(__GNUC__)
#define ALIGN(n) __attribute__((__aligned__(n)))
// Intel Compiler
#elif defined(__INTEL_COMPILER)
// GNU Intel Compiler
#define ALIGN(n) __declspec(align(n))
// Unknown Compiler
#else
#define ALIGN(n) 
#endif

#ifdef PERIODIC
#ifdef POWER6
#define NEAREST_X(x) ((x)  - boxSize_X * __frin ( (x) * inverse_boxSize_X))
#define NEAREST_Y(y) ((y)  - boxSize_Y * __frin ( (y) * inverse_boxSize_Y))
#define NEAREST_Z(z) ((z)  - boxSize_Z * __frin ( (z) * inverse_boxSize_Z))
#define NEAREST(x) ((x)  - boxSize * __frin ( (x) * inverse_boxSize))
#else
#define NEAREST_X(x) (((x)>boxHalf_X)?((x)-boxSize_X):(((x)<-boxHalf_X)?((x)+boxSize_X):(x)))
#define NEAREST_Y(y) (((y)>boxHalf_Y)?((y)-boxSize_Y):(((y)<-boxHalf_Y)?((y)+boxSize_Y):(y)))
#define NEAREST_Z(z) (((z)>boxHalf_Z)?((z)-boxSize_Z):(((z)<-boxHalf_Z)?((z)+boxSize_Z):(z)))
#define NEAREST(x) (((x)>boxHalf)?((x)-boxSize):(((x)<-boxHalf)?((x)+boxSize):(x)))
#define __fsel(crit,age,alt) (((crit) >= 0.0) ? (age) : (alt))
#endif
#else
#define NEAREST_X(x) (x)
#define NEAREST_Y(x) (x)
#define NEAREST_Z(x) (x)
#define NEAREST(x) (x)
#endif

#define ASSIGN_ADD(x,y,mode) (mode == 0 ? (x=y) : (x+=y))

#define  GADGETVERSION   "3.0"	/*!< code version string */

#ifndef  GENERATIONS
#define  GENERATIONS     2	/*!< Number of star particles that may be created per gas particle */
#endif


#ifdef   ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef  long long integertime;
#define  TIMEBINS        60
#define  TIMEBASE        (((long long)1)<<TIMEBINS)     /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                                         *   where TIMESPAN needs to be a power of 2.
                                                         */
#else
typedef  int integertime;
#define  TIMEBINS        29
#define  TIMEBASE        (1<<TIMEBINS)  /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                         *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
                                         *   to 2^29
                                         */
#endif

#ifdef RADTRANSFER
#ifdef UM_CHEMISTRY
#define N_BINS N_nu
#else
#ifndef RT_MULTI_FREQUENCY
#define N_BINS 1
#else
#define N_BINS 4
#endif
#endif
#endif

#ifndef  MULTIPLEDOMAINS
#define  MULTIPLEDOMAINS     1
#endif

#ifdef ONEDIM
#define DIMS 1
#else
#ifdef TWODIMS    /* will only be compiled in 2D case */
#define DIMS 2
#else
#define DIMS 3
#endif
#endif

#ifndef  TOPNODEFACTOR
#define  TOPNODEFACTOR       2.5
#endif

#ifndef  GRAVCOSTLEVELS
#define  GRAVCOSTLEVELS      6
#endif

#define  NUMBER_OF_MEASUREMENTS_TO_RECORD  6  /* this is the number of past executions of a timebin that the reported average CPU-times average over */

#define  NODELISTLENGTH      8

typedef unsigned long long peanokey;


#define  BITS_PER_DIMENSION 21	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))

#define  BITS_PER_DIMENSION_SAVE_KEYS 10
#define  PEANOCELLS_SAVE_KEYS (((peanokey)1)<<(3*BITS_PER_DIMENSION_SAVE_KEYS))

#ifdef LT_STELLAREVOLUTION
#include "lt_sfr/lt.h"
#endif

#ifdef LT_ADD_GAL_TO_SUB
extern float *tempiAS, *CB07,*Filters_Effective_Wavelenght;
#ifdef OBSERVER_FRAME
extern float *CB07obs;
#endif
#endif

#define  check_particles()          check_particles_info( __FUNCTION__, __FILE__, __LINE__)

#define  terminate(x) {char termbuf[2000]; sprintf(termbuf, "code termination on task=%d, function '%s()', file '%s', line %d: '%s'\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); printf("%s", termbuf); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, 1); exit(0);}

#ifndef DISABLE_MEMORY_MANAGER
#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)
#else
#define  mymalloc(x, y)            malloc(y)
#define  mymalloc_movable(x, y, z) malloc(z)

#define  myrealloc(x, y)           realloc(x, y)
#define  myrealloc_movable(x, y)   realloc(x, y)

#define  myfree(x)                 free(x)
#define  myfree_movable(x)         free(x)

#define  report_memory_usage(x, y) printf("Memory manager disabled.\n")
#endif

#ifdef LT_STELLAREVOLUTION
#define endrun(x) EndRun(x, __FUNCTION__, __FILE__, __LINE__)
#endif

#ifdef AB_SHOCK
#define SHOCK_GRADIENT_TOLERANCE 1.2
#define SHOCK_MINIMAL_VDIF_FRACTION 0.0001
#define SHOCK_NOISELEVEL 0.05
#endif

#ifndef  GAMMA
#define  GAMMA         (5.0/3.0)	/*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_MINUS1  (GAMMA-1)
#define  GAMMA_MINUS1_INV  (1./(GAMMA-1))

#ifndef  HYDROGEN_ONLY
#define  HYDROGEN_MASSFRAC 0.76 /*!< mass fraction of hydrogen, relevant only for radiative cooling */
#else
#define  HYDROGEN_MASSFRAC 1 /*!< mass fraction of hydrogen, relevant only for radiative cooling */
#endif

#define  METAL_YIELD       0.02	/*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 8192

/* ... often used physical constants (cgs units) */

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.38066e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  LYMAN_ALPHA      1215.6e-8	/* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8	/* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615
#define  ELECTRONVOLT_IN_ERGS      1.60217733e-12


#ifdef MAGNETIC
#ifdef SFR 
#define POW_CC 1./3.
#endif

#ifdef MU0_UNITY
#define MU0 1.0
#define MU0_1 1.0
#else
#define MU0  12.566370614  /* 4pi */
#define MU0_1 0.079577472  /* 1/4pi */
#endif
#endif

#ifdef NAVIERSTOKES
#define  LOG_LAMBDA      37.8	/* logarithmic Coulomb factor */
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#define  T_CMB0      2.728	/* present-day CMB temperature */
#endif

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

/*Determines the maximum size of arrays related to the number of CR populations */
#ifndef NUMCRPOP   /*!< Number of CR populations pressent in parameter file */
#define NUMCRPOP 1
#endif



#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif


/* some flags for the field "flag_ic_info" in the file header */
#define FLAG_ZELDOVICH_ICS     1
#define FLAG_SECOND_ORDER_ICS  2
#define FLAG_EVOLVED_ZELDOVICH 3
#define FLAG_EVOLVED_2LPT      4
#define FLAG_NORMALICS_2LPT    5


#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5
#endif

#define COND_TIMESTEP_PARAMETER 0.25
#define VISC_TIMESTEP_PARAMETER 0.25

#define MAXLEN_OUTPUTLIST 1400	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000	/*!< length of the lookup table used to hold the drift and kick factors */


#define MAXITER 150

#ifndef LINKLENGTH
#define LINKLENGTH 0.2
#endif

#ifndef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_LEN 32
#endif

#define MINRESTFAC 0.05


#ifdef SUBFIND_DENSITY_AND_POTENTIAL  /*!< activate needed options */
 #define ONLY_PRODUCE_HSML_FILES
 #define COMPUTE_POTENTIAL_ENERGY
 #define SUBFIND_RESHUFFLE_AND_POTENTIAL
 #define SUBFIND_RESHUFFLE_CATALOGUE
#endif
#ifdef ADAPTGRAVSOFT
#ifndef AGS_MAXSOFT
/*!< AGS_MAXSOFT gives the maximum allowed gravitational softening when using the TreePM method.
 *  The quantity is given in units of the scale used for the force split (ASMTH) */
#define AGS_MAXSOFT 0.5
#endif
#ifndef AGS_MINSOFT
/*!< AGS_MINSOFT gives the minimum allowed gravitational softening when using the TreePM method.
 *  The quantity is given in units of the scale used for the force split (ASMTH) */
#define AGS_MINSOFT 0.001
#endif
#endif

#ifndef GDE_TYPES 
#define GDE_TYPES 2
#endif

#ifdef SORT_FROM_L3
typedef int  sort_type;
#define sizesort sizeof(sort_type)
#endif

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
#else
#if (DOUBLEPRECISION+0) == 2 
typedef float   MyFloat;
typedef double  MyDouble;
#else                        /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
#endif
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif

struct unbind_data
{
  int index;
};


#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif

#ifdef FLTROUNDOFFREDUCTION
#define FLT(x) ((MyFloat)(x))
#ifdef SOFTDOUBLEDOUBLE      /* this requires a C++ compilation */
#include "dd.h"
typedef dd MyLongDouble;
#else
typedef long double MyLongDouble;
#endif
#else  /* not enabled */
#define FLT(x) (x)
typedef MyFloat MyLongDouble;
#endif  /* end FLTROUNDOFFREDUCTION */


#ifdef GDE_BIGFLOAT
#include "bigfloat.h"
typedef BigFloat MyBigFloat;
extern int xErrNo;
#else
typedef MyDouble MyBigFloat;
#endif

#ifdef GDE_BIGFLOAT
#define GDE_ABS(x) ((x).abs())
#define GDE_SQRT(x) ((x).sqrt())
#define GDE_LOG(x) ((x).log())
#else
#define GDE_ABS(x) (fabs(x))
#define GDE_SQRT(x) (sqrt(x))
#define GDE_LOG(x) (log(x))
#endif



#define CPU_ALL            0
#define CPU_TREEWALK1      1
#define CPU_TREEWALK2      2
#define CPU_TREEWAIT1      3
#define CPU_TREEWAIT2      4
#define CPU_TREESEND       5
#define CPU_TREERECV       6
#define CPU_TREEMISC       7
#define CPU_TREEBUILD      8
#define CPU_TREEUPDATE     9
#define CPU_TREEHMAXUPDATE 10
#define CPU_DOMAIN         11
#define CPU_DENSCOMPUTE    12
#define CPU_DENSWAIT       13
#define CPU_DENSCOMM       14
#define CPU_DENSMISC       15
#define CPU_HYDCOMPUTE     16
#define CPU_HYDWAIT        17
#define CPU_HYDCOMM        18
#define CPU_HYDMISC        19
#define CPU_DRIFT          20
#define CPU_TIMELINE       21
#define CPU_POTENTIAL      22
#define CPU_MESH           23
#define CPU_PEANO          24
#define CPU_COOLINGSFR     25
#define CPU_SNAPSHOT       26
#define CPU_FOF            27
#define CPU_BLACKHOLES     28
#define CPU_MISC           29
#define CPU_SMTHCOMPUTE    30
#define CPU_SMTHWAIT       31
#define CPU_SMTHCOMM       32
#define CPU_SMTHMISC       33
#define CPU_HOTNGBS        34
#define CPU_WEIGHTS_HOT    35
#define CPU_ENRICH_HOT     36
#define CPU_WEIGHTS_COLD   37
#define CPU_ENRICH_COLD    38
#define CPU_CSMISC         39
#define CPU_HYDNETWORK     40
#define CPU_AGSDENSCOMPUTE 41
#define CPU_AGSDENSWAIT    42
#define CPU_AGSDENSCOMM    43
#define CPU_AGSDENSMISC    44
#define CPU_AGSTREEHMAXUPD 45
#define CPU_MG_CIC         46
#define CPU_MG_FIELDSOLVE  47
#define CPU_MG_EFF_MASS    48
#define CPU_PARTS          49  /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 120

#if !defined(QUINTIC_KERNEL)
#if !defined(TWODIMS) && !defined(ONEDIM)
#define  NUMDIMS 3		/*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786	/*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#ifndef  ONEDIM
#define  NUMDIMS 2		/*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI	/*!< Coefficient for kernel normalization. */
#else
#define  NUMDIMS 1             /*!< For 1D-normalized kernel */
#define  KERNEL_COEFF_1  (4.0/3)
#define  KERNEL_COEFF_2  (8.0)
#define  KERNEL_COEFF_3  (24.0)
#define  KERNEL_COEFF_4  (16.0)
#define  KERNEL_COEFF_5  (8.0/3)
#define  KERNEL_COEFF_6  (-8.0)
#define  NORM_COEFF      2.0
#endif
#endif
#else /* here comes the QUINTIC kernel */
#if !defined(TWODIMS) && !defined(ONEDIM)
#define  NUMDIMS 3
#define  NORM_COEFF      4.188790204786
#else
#ifndef  ONEDIM
#define  NUMDIMS 2
#define  NORM_COEFF      M_PI
#else
#define  NUMDIMS 1 
#define  NORM_COEFF      2.0
#endif
#endif
#endif /* end of !QUINTIC */


#if defined (BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION)
#define PPP P
#else
#define PPP SphP
#endif


#ifdef PERIODIC
extern MyDouble boxSize, boxHalf, inverse_boxSize;
#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X, inverse_boxSize_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#define inverse_boxSize_X inverse_boxSize
#endif
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y, inverse_boxSize_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#define inverse_boxSize_Y inverse_boxSize
#endif
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z, inverse_boxSize_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#define inverse_boxSize_Z inverse_boxSize
#endif
#endif

#ifdef PERIODIC
#ifndef POWER6
#define NGB_PERIODIC_LONG(x,box,hbox) ((fabs(x)>hbox)?(box-fabs(x)):fabs(x))
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp)
#else
#ifdef DOUBLEPRECISION
#define NGB_PERIODIC_LONG(x,box,hbox) __fsel(hbox-fabs(x),fabs(x),box-fabs(x))
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),__fsel(boxHalf_X-xtmp,xtmp,boxSize_X-xtmp))
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),__fsel(boxHalf_Y-xtmp,xtmp,boxSize_Y-xtmp))
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),__fsel(boxHalf_Z-xtmp,xtmp,boxSize_Z-xtmp))
#else
#define NGB_PERIODIC_LONG(x,box,hbox) __fsel(hbox-fabsf(x),fabsf(x),box-fabsf(x))
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabsf(x),__fsels(boxHalf_X-xtmp,xtmp,boxSize_X-xtmp))
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabsf(x),__fsels(boxHalf_Y-xtmp,xtmp,boxSize_Y-xtmp))
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabsf(x),__fsels(boxHalf_Z-xtmp,xtmp,boxSize_Z-xtmp))
#endif
#endif
#else
#define NGB_PERIODIC_LONG(x,box,hbox) fabs(x)
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#endif

#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACT2 0.86602540        /* FACT2 = 0.5 * sqrt(3) */

#ifdef FS_TURB_ESTIM
#define FS_BINS 30
#endif


/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

#ifdef UM_SEN_H2_SOLOMON_PROCESS
#define JLWn 23
extern float JLW[JLWn];
extern float z_JLW[JLWn];
#endif

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
extern int NActivePart;
extern int *ActiveParticleList;
#if defined(LT_STELLAREVOLUTION)
extern int NActivePartMerk;
#endif
#endif

extern int FirstActiveParticle;
extern int *NextActiveParticle;
extern unsigned char *ProcessedFlag;

extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinActive[TIMEBINS];

extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern int *NextInTimeBin;
extern int *PrevInTimeBin;

#ifdef SFR
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
extern double TimeBin_BH_Medd[TIMEBINS];
#endif

extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */

#ifdef INVARIANCETEST
extern int World_ThisTask;
extern int World_NTask;
extern int Color;
extern MPI_Comm MPI_CommLocal;
#ifndef DO_NOT_REDEFINE_MPI_COMM_WORLD
#undef  MPI_COMM_WORLD
#define MPI_COMM_WORLD MPI_CommLocal
#endif
#endif


extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;

extern int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */

extern int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
extern int RestartSnapNum;
extern int SelRnd;

extern int TakeLevel;

extern int *Exportflag;	        /*!< Buffer used for flagging whether a particle needs to be exported to another process */
extern int *Exportnodecount;
extern int *Exportindex;

extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset;

#ifdef VORONOI
extern int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;
#endif

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN + 1];

extern double WallclockTime;    /*!< This holds the last wallclock time measurement for timings measurements */

extern int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */

extern size_t HighMark_run,  HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic,  HighMark_sphdensity, HighMark_sphhydro, HighMark_addSPH;

#ifdef VORONOI
extern size_t HighMark_voronoi;
#endif

#ifdef CS_MODEL
#if defined(CS_SNI) || defined(CS_SNII)
extern size_t HighMark_cs_enrich;
#endif
extern size_t HighMark_cs_hotngbs;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
extern size_t HighMark_turbpower;
#endif

#ifdef ADAPTGRAVSOFT
extern size_t HighMark_agsdensity;
#endif

extern int TreeReconstructFlag;
extern int GlobFlag;

extern char DumpFlag;

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
#ifdef LT_STELLAREVOLUTION
extern int N_stars;             /*!< number of star particles in the LOCAL processor */
#endif
#if defined(BLACK_HOLES) && defined(DETACH_BLACK_HOLES)
extern int N_BHs;
#endif

#ifdef SINKS
extern int NumSinks;
#endif

extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */

extern gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif

#ifdef MODGRAV
extern int NumSavedAMRCells; /*!< number of AMR cells from old tree saved on LOCAL task */
extern int NumAMRCells; /*!< number of AMR cells on LOCAL task */

extern int NumSavedMultigridCells;  /*!< number of multigird cells from old tree saved on LOCAL task */
extern int SavedMultigridMinLevel;  /*!< minimum tree level for which multigrid cells have been saved */
extern int SavedMultigridMaxLevel;  /*!< maximum tree level for which multigrid cells have been saved */

extern int FirstAMRCell; /*!< index of first AMR node */
extern int *NextAMRCell, *NextAMRCell_base; /*!< defined for all nodes, for AMR nodes it gives the index of the next AMR node */

extern int phi_guess_available;

extern struct amr_cell_data /*!< used to store AMR data of old tree */
{
  MyFloat center[3];
  int level;
  int bitflags;
  MyFloat phi;
  
  int target_task;
} *saved_amr_cells, *saved_multigrid_cells;

extern double *mg_time_multi_grid_nbg_search;
extern double *mg_time_multi_grid_update_cells;
extern double *mg_time_multi_grid_update_cells_comm_imbal;
extern int *mg_num_it_multi_grid;
extern long long *mg_num_cells_multi_grid;
#endif

extern double TimeOfLastTreeConstruction;	/*!< holds what it says */

extern int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */

extern double *R2ngblist;

extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;

extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;

extern peanokey *Key, *KeySorted;

#ifdef RADTRANSFER
double lum[N_BINS];
#ifdef RT_POPIII
double lum_popIII[N_BINS];
#endif
double rt_sigma_HI[N_BINS];
double rt_sigma_HeI[N_BINS];
double rt_sigma_HeII[N_BINS];
double G_HI[N_BINS];
double G_HeI[N_BINS];
double G_HeII[N_BINS];
#ifndef UM_CHEMISTRY
double nu[N_BINS];
#endif
#endif

extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  MyFloat GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;

extern double RndTable[RNDTABLE];


#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif


#ifdef WRITE_KEY_FILES
extern peanokey *KeyIndex;
extern int *NPartPerKey, *PartKeyOffset;
extern int NKeys[6];
extern long long NKeysTotal[6];
#endif


/* variables for input/output , usually only used on process 0 */


extern char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,		/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU,			/*!< file handle for cpu.txt log-file. */
 *FdTimebin;

#ifdef SCFPOTENTIAL
extern FILE *FdSCF;
#endif

#ifdef SFR
extern FILE *FdSfr;		/*!< file handle for sfr.txt log-file. */
#endif

#if defined (VS_TURB) || defined (AB_TURB)
extern FILE *FdTurb;
#endif

#ifdef RADTRANSFER
extern FILE *FdRad;		/*!< file handle for radtransfer.txt log-file. */
extern FILE *FdStar;		/*!< file handle for lum_star.txt log-file. */
#endif

#ifdef DISTORTIONTENSORPS
#ifdef PMGRID
extern FILE *FdTidaltensor;     /*!< file handle for tidaltensor.txt log-file. */
#endif
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;	/*!< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
#endif


#ifdef FORCETEST
extern FILE *FdForceTest;	/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef MODGRAV
extern FILE *FdMOG;  /*!< file handle for modgrav.txt log-file. */
extern FILE *FdMultiGridTimes;  /*!< file handle for multi_grid_times.txt log-file */
#endif

#ifdef DARKENERGY
extern FILE *FdDE;  /*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
extern FILE *FdXXL;		/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
extern double MeanB;

#ifdef TRACEDIVB
extern double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
extern double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];

#ifdef MAGNETIC
/*! table for the cosmological kick factor for induction equation */
extern double MagKickTable[DRIFT_TABLE_LENGTH];
#endif


#ifdef VORONOI
extern struct individual_data
{
  double AllocFacNdp;
  double AllocFacNdt;
  double AllocFacNvf;
  double AllocFacNinlist;
  double AllocFacN_DP_Buffer;
}
Indi;
#endif


extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef LT_STELLAREVOLUTION
  long long TotN_stars;         /*!<  total star particle number (global value) */
#endif

#ifdef NEUTRINOS
  long long TotNumNeutrinos;
#endif

#ifdef BLACK_HOLES
  int TotBHs;
#if defined(DETACH_BLACK_HOLES)
  int MaxPartBH;
  double BHfactor;
#endif
#endif

#ifdef SIDM
  double CrossSectionPerMass_in_cgs;
  double PeakSigma;
#ifdef SIDM_INELASTIC
  double Splitting;
#endif
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				   processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				   processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int DoDynamicUpdate;

  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when
					   writing/reading restart-files, or when writing snapshot files */

  double BufferSize;		/*!< size of communication buffer in MB */
  int BunchSize;     	        /*!< number of particles fitting into the buffer in the parallel tree algorithm  */


  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

#ifdef SCALARFIELD
  double ScalarBeta;
  double ScalarScreeningLength;
#endif

#ifdef WINDTUNNEL
  double WindGrid;
  double WindVel;
  double WindDmean;
  double WindEntr;
  double WindDens;
  double WindPmass;
  double WindCurrentX;
  double WindSizeOfInjectionRegion;
#endif

  /* some SPH parameters */

  int DesNumNgb;		/*!< Desired number of SPH neighbours */
#ifdef SUBFIND
  int DesLinkNgb;
  double ErrTolThetaSubfind;
#endif

  double MaxNumNgbDeviation;	/*!< Maximum allowed deviation neighbour number */
#ifdef START_WITH_EXTRA_NGBDEV
  double MaxNumNgbDeviationStart;    /*!< Maximum allowed deviation neighbour number to start with*/
#endif

  double ArtBulkViscConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;		/*!< the minimum allowed temperature expressed as energy per unit mass */
#ifdef ARTIFICIAL_CONDUCTIVITY
  double ArtCondConstant;
#endif

#ifdef KSPACE_NEUTRINOS
  int KspaceNeutrinoSeed;
  int Nsample;
  int SphereMode;
  char	KspaceDirWithTransferfunctions[500];
  char	KspaceBaseNameTransferfunctions[500];
  double PrimordialIndex;
  double Sigma8;
  double InputSpectrum_UnitLength_in_cm;
  double OmegaNu;
#endif


  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;	/*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */

  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a;   /* various cosmological factors that are only a function of the current scale factor, and in Newtonian runs are set to 1 */

  /* system of units  */

  double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,		/*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/*!< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,	/*!< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/*!< Gravity-constant in internal units */
  double UnitDensity_in_Gev_per_cm3; /*!< factor to convert internal density unit to GeV/c^2 / cm^3 */
  /* Cosmology */

  double Hubble;		/*!< Hubble-constant in internal units */
  double Omega0,		/*!< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/*!< baryon density in units of the critical density (at z=0) */
    HubbleParam;		/*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;		/*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;	/*!< flags that periodic boundaries are enabled */
  int ResubmitOn;		/*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;	/*!< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;		/*!< flags that output times are listed in a specified file */
  int CoolingOn;		/*!< flags that cooling is enabled */
  int StarformationOn;		/*!< flags that star formation is enabled */

  int HighestActiveTimeBin;
  int HighestOccupiedTimeBin;

  /* parameters determining output frequency */

  int SnapshotFileCount;	/*!< number of snapshot that is written next */
  double TimeBetSnapshot,	/*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/*!< cpu-time when last restart-file was written */
    TimeBetStatistics,		/*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;		/*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,			/*!< current time of the simulation */
    TimeBegin,			/*!< time of initial conditions of the simulation */
    TimeStep,			/*!< difference between current times of previous and current timestep */
    TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;	/*!< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;		/*!< current time on integer timeline */
  integertime Previous_Ti_Current;
  integertime Ti_nextoutput;		/*!< next output time on integer timeline */
  integertime Ti_lastoutput;

#ifdef PMGRID
  integertime PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  double Epsilon;
#endif

  integertime Ti_nextlineofsight;
#ifdef OUTPUTLINEOFSIGHT
  double TimeFirstLineOfSight;
#endif

  int    CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];

  int LevelToTimeBin[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];    /*!< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,	/*!< minimum allowed timestep. Normally, the simulation terminates if the
				   timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;		/*!< maximum allowed timestep */

  double MaxRMSDisplacementFac;	/*!< this determines a global timestep criterion for cosmological simulations
				   in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement
				   is a fraction of the mean particle separation (determined from the
				   particle mass and the cosmological parameters). This parameter specifies
				   this fraction. */

  int MaxMemSize;

  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional,	/*!< minimum allowed SPH smoothing length in units of SPH gravitational
				   softening length */
    MinGasHsml;			/*!< minimum allowed SPH smoothing length */

#ifdef MODGRAV
  int MaxNumAMRCells; /*!< Maximum number of AMR cells per task */
  int MaxNumMultigridCells; /*!< Maximum number of AMR cells per task */
#endif

#ifdef MAXHSML
  double MaxHsml;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  double RefDensity;
  double RefEntropy;
  double IsoSoundSpeed;
  double TurbInjectedEnergy;
  double TurbDissipatedEnergy;
  double TimeBetTurbSpectrum;
  double TimeNextTurbSpectrum;
  int FileNumberTurbSpectrum;
  double SetLastTime;
#endif

  double SofteningGas,		/*!< for type 0 */
    SofteningHalo,		/*!< for type 1 */
    SofteningDisk,		/*!< for type 2 */
    SofteningBulge,		/*!< for type 3 */
    SofteningStars,		/*!< for type 4 */
    SofteningBndry;		/*!< for type 5 */

  double SofteningGasMaxPhys,	/*!< for type 0 */
    SofteningHaloMaxPhys,	/*!< for type 1 */
    SofteningDiskMaxPhys,	/*!< for type 2 */
    SofteningBulgeMaxPhys,	/*!< for type 3 */
    SofteningStarsMaxPhys,	/*!< for type 4 */
    SofteningBndryMaxPhys;	/*!< for type 5 */

  double SofteningTable[6];	/*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];	/*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], TimingsFile[100], TimebinFile[100], RestartFile[100], ResubmitCommand[100], OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength;		/*!< number of times stored in table of desired output times */



#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
  double ReferenceGasMass;
#endif

#ifdef VORONOI_MESHRELAX
  double MeanMass;
  double MeanPressure;
#endif

#ifdef RADTRANSFER
  double IonizingLumPerSolarMass;
  double IonizingLumPerSolarMassPopIII;
  double IonizingLumPerSFR;
  double star_Teff;
  double star_Teff_popIII;
  integertime Radiation_Ti_begstep;
  integertime Radiation_Ti_endstep;
#endif

#ifdef ADAPTIVE_FORCE_ACC
  double ErrTolForceAccParam;
#endif


#ifdef DISTORTIONTENSORPS
  /* present day velocity dispersion of DM particle in cm/s (e.g. Neutralino = 0.03 cm/s) */
  double DM_velocity_dispersion;
  double TidalCorrection;
#ifdef GDE_LEAN
  double GDEInitStreamDensity;
#endif
#endif

#ifdef SFR		/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
#ifdef LT_STELLAREVOLUTION
  double OrigGasMass;
#endif
  double EgySpecSN;
  double FactorSN;
  double EgySpecCold;
  double FactorEVP;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double FactorForSofterEQS;
#ifdef WINDS
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelMaxTimeFactor;  /* maximum free travel time in units of the Hubble time at the current simulation redshift */
  double WindFreeTravelDensFac;
#ifdef VARIABLE_WINDS
  double VariableWindVelFactor;  /* wind velocity in units of the halo escape velcoity */
  double VariableWindSpecMomentum;  /* momentum available for wind per unit mass of stars formed, in internal velocity units */
  double HaloConcentrationNorm;  /* concentration c0 of a halo of unit mass */
  double HaloConcentrationSlope;  /* slope n of mass concentration relation, namely c = c0 * M_200,crit^n */
#endif
#endif
#endif

#ifdef CS_MODEL
  double FactorSFR;
  double DecouplingParam;
  double MinTlifeSNI;
  double MaxTlifeSNI;
  double TlifeSNII;
  int    Raiteri_TlifeSNII;
  double RateSNI;
  double SN_Energy_cgs;
  double Tcrit_Phase;
  double DensFrac_Phase;
  double SN_Energy_frac_cold;
  double MaxHotHsmlParam;
  double InitialHotHsmlFactor;
  double DensityTailThreshold;
  double MaxNumHotNgbDeviation;	/*!< Maximum allowed deviation HOT neighbour number */
#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#ifdef TIMEDEPGRAV
  double Gini;
#endif
#endif
#endif

#ifdef MODGRAV
  char ModifiedGravityFile[100];
  int MaxAMRLevel;
  int MaxLevelToptree;
  int MaxLevelFullTree;
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif

#if defined(SNIA_HEATING)
  double SnIaHeatingRate;
#endif

#ifdef VORONOI_SHAPESCHEME
  double VoronoiStiffNess;
  double VoronoiRoundNess;
#endif

#ifdef TIME_DEP_ART_VISC
  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution */
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */
#endif

#if defined(LT_MV_CHEMICALDIFFUSION)
  int StandardCourantDiffusion;
  double MaxDiffusionTimeStep;
  double DiffusionCourantFac;
  double DiffusionCourantMinFac;
  double MetalsAmount[LT_NMet];
#if defined(LT_MV_CHEMICALDIFFUSION_CONSTANT_D)
  double ConstantD;
#endif
#endif
  
#ifdef CONDUCTION
  double ConductionCoeff;	/*!< Thermal Conductivity */
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;	/*!< Factor to get electron mean free path */
#endif

  integertime Conduction_Ti_begstep, Conduction_Ti_endstep;
  double MaxSizeConductionStep;
#endif

#if defined(HEALPIX)
  //change this to read in the Parameterfile
  int Nside;
#define NSIDE2NPIX(nside)  (12*nside*nside)
  float *healpixmap;
  double Minmass,Maxmass;
#endif

#ifdef VSMOOTH
  double VSmoothScale;
#endif

#ifdef MAGNETIC
#ifdef BINISET
  double BiniX, BiniY, BiniZ;	/*!< Initial values for B */
#endif

#if defined(BSMOOTH) || defined(BSMOOTH_TIME)
  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;
#ifdef SETMAINTIMESTEPCOUNT
  int MainTimestepCountIni;
#endif
#endif /* BSMOOTH || BSMOOTH_TIME */

#if defined(MAGNETIC_DISSIPATION) || defined(EULER_DISSIPATION)
#if !defined(EULER_DISSIPATION)
#define MAGNETIC_SIGNALVEL      /*!< Make explicit the dependence with MAGNETIC_SIGNALVEL as described in Dolag2009 */
#endif
  double ArtMagDispConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial magnetic disipation */
#ifdef TIME_DEP_MAGN_DISP
  double ArtMagDispMin;
  double ArtMagDispSource;
  double ArtMagDispTime;
#endif
#endif

#ifdef MAGNETIC_SN_SEEDING
  double SnSeedRadius;
  double SnSeedBubble;
  double SnSeedField;
  double SnSeedSoftening;
#endif

#ifdef DIVBCLEANING_DEDNER
  double DivBcleanParabolicSigma;
  double DivBcleanHyperbolicSigma;
  double DivBcleanQ;
#endif

#ifdef FS_ALFA2_QUENCH
  double Alfa0Quench;
  double AlfaB2eq;
#endif
#ifdef FS_DYN_DEF
  double Eta2Min;
  double Eta2Max;
  double Alfa2H;
#endif

#ifdef MAGNETIC_DIFFUSION
  double MagneticEta;
#endif
#endif /* MAGNETIC */

#if defined(BLACK_HOLES) || defined(VARIABLE_WINDS)
  double TimeNextOnTheFlyFoF;
  double TimeBetOnTheFlyFoF;
#endif

#ifdef BLACK_HOLES
  double BlackHoleAccretionFactor;	/*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;	/*!< Fraction of the black luminosity feed into thermal feedback */
  double SeedBlackHoleMass;	/*!< Seed black hole mass */
  double MinFoFMassForNewSeed;	/*!< Halo mass required before new seed is put in */
  double BlackHoleNgbFactor;	/*!< Factor by which the normal SPH neighbour should be increased/decreased */
  double BlackHoleMaxAccretionRadius;
  double BlackHoleEddingtonFactor;	/*!< Factor above Eddington */
  double BlackHoleRadiativeEfficiency;  /**< Radiative efficiency determined by the spin value, default value is 0.1 */
#ifdef FOF
  double massDMpart;
#endif
#ifdef MODIFIEDBONDI
  double BlackHoleRefDensity;
  double BlackHoleRefSoundspeed;
#endif
#ifdef LT_BH_ACCRETE_SLICES
  int    NBHslices;
#endif
#ifdef LT_DF_BH_BHAR_SWITCH
  double BH_radio_treshold;
#endif
#endif

#ifdef COSMIC_RAYS
  double CR_Alpha[NUMCRPOP];	/*!< Cosmic ray spectral index [2..3] */
  double CR_SNEff;		/*!< SN injection efficiency [0..1] */
  double CR_SNAlpha;		/*!< SN injection spectral index [2..3] */
  int bDebugFlag;		/*!< enables debug outputs after triggered */

#if defined(CR_DIFFUSION)
  double CR_DiffusionCoeff;	/*!< (temporary) fixed value for CR diffusivity */

  double CR_DiffusionDensScaling;	/*!< grade of density dependence of diffusivity */
  double CR_DiffusionDensZero;	/*!< Reference point density for diffusivity */

  double CR_DiffusionEntropyScaling;	/*!< grade of specific energy dependence of diffusivity */

  double CR_DiffusionEntropyZero;	/*!< Reference Entropic function for diffusivity */

  double CR_DiffusionMaxSizeTimestep;
  integertime CR_Diffusion_Ti_begstep, CR_Diffusion_Ti_endstep;
#endif				/* CR_DIFFUSION */

#if defined(CR_SHOCK)
#if (CR_SHOCK == 1)
  double CR_ShockAlpha;		/*!< spectral index to be used in shock injection */
#else
  double CR_ShockCutoff;	/*!< Cutoff factor x_inj for CR accel */
#endif
  double CR_ShockEfficiency;	/*!< energy fraction of shock energy fed into CR */
#endif				/* CR_SHOCK */

#ifdef FIX_QINJ
  double Shock_Fix_Qinj;	/*!< inject only CRps with threshold cutoff Shock_Fix_Qinj */
#endif

#ifdef CR_BUBBLES
  double CR_AGNEff;               /*!< AGN injection efficiency [0..1] */
#endif
#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  double Shock_Length;		/*!< length scale on which the shock is smoothed out */
  double Shock_DeltaDecayTimeMax;	/*!< maximum time interval (Dloga) for which the
					   Mach number is kept at its maximum */
#endif

#ifdef REIONIZATION
  int not_yet_reionized;	/*!< flag that makes sure that there is only one reionization */
#endif



#ifdef BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double FirstBubbleRedshift;
#ifdef FOF
  int BiggestGroupLen;
  float BiggestGroupCM[3];
  double BiggestGroupMass;
#endif
#endif

#ifdef BH_BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleEnergy;
  double BlackHoleRadioTriggeringFactor;
  double DefaultICMDensity;
  double RadioFeedbackFactor;
#ifdef UNIFIED_FEEDBACK
  double RadioThreshold;
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
#ifndef BLACK_HOLES
  double MinFoFMassForNewSeed;	/*!< Halo mass required before new seed is put in */
  double massDMpart;
#endif
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double ClusterMass200;
  double FirstBubbleRedshift;
#endif

#ifdef NAVIERSTOKES
  double NavierStokes_ShearViscosity;
  double FractionSpitzerViscosity;
  double ShearViscosityTemperature;
#endif
#ifdef NAVIERSTOKES_BULK
  double NavierStokes_BulkViscosity;
#endif
#ifdef VISCOSITY_SATURATION
  double IonMeanFreePath;
#endif

#ifdef EOS_DEGENERATE
  char EosTable[100];
  char EosSpecies[100];
#endif

#ifdef SINKS
  int TotNumSinks;
  double SinkHsml;
  double SinkDensThresh;
#endif

#ifdef NUCLEAR_NETWORK
  char NetworkRates[100];
  char NetworkPartFunc[100];
  char NetworkMasses[100];
  char NetworkWeakrates[100];
  struct network_data nd;
  struct network_workspace nw;
  double NetworkTempThreshold;
#endif

#ifdef RELAXOBJECT
  double RelaxBaseFac;
  double RelaxFac;
#endif

#ifdef BP_REAL_CRs
  double ecr_min, ecr_max;		/*!< min and max momenta of cosmic rays */
  double ecr_bound[BP_REAL_CRs+1];	/*!< boundaries of cosmic rays momentum bins */
  double CR_Gamma;			/*!< adiabatic index of CRs  */
#ifdef BP_SEED_CRs
  double pSlope_init, eSlope_init;      /*!< init slopes of spectra of protons and electrons */
#endif
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  double CRsArtCondConstant;
#endif
#endif

#ifdef LT_STELLAREVOLUTION
  int           MaxPartMet;           /*!< This gives the maxmimum number of STAR particles that can be stored on one
                                        processor. */

  double        Time_Age;             /*!< current cosmic time in Gyrs */

  int           TestSuite, StarBits;
  double        SFfactor;             /*!< the expected maximum factor of proportionality between TotN_gas and TotN_star */
  int           Generations;

  int           NeighInfNum;          /*!< minimum number of neighbours for metal and egy spreading */
  int           DesNumNgbSN,          /*!< desired number of neighbours for sn spreading */
                SpreadNumNgbDev,      /*!< maximum deviation in number of neighbours for sn spreading */
                LeftNumNgbSN,
                RightNumNgbSN;        /*!< range of neighbours for sn spreading*/
  double        SpreadNeighCoeff;     /*!< store DesNumNgbSN / DesNumNgb */

  double        MinChemSpreadL;
  double        MinChemTimeStep;
  double        Enrich_SFGas_Th;

  double        LocalSpreadFactor;
  char          SFfilename[300], IMFfilename[300], SnIaDataFile[300], SnIIDataFile[300], AGBDataFile[300];
  int           Ia_Nset_ofYields, II_Nset_ofYields, AGB_Nset_ofYields;

  double        Mup,        /*!< SnII threshold (e.g. 8 Msun) */
                MBm,        /*!< min bin. system mass */
                MBM,        /*!< max bin. system mass */
                BinFrac,    /*!< bin. system fraction */
                MBms;       /*!< min star mass in SnIa binary systems */
  double        SnIaRemn;   /*!< SnIa Remn (1.4Msun ?) */
                            /* lifetimes */
  double        inf_lifetime, mean_lifetime, sup_lifetime;
                            /* energy provided by sn explosions */
  double        SnIaEgy, SnIIEgy;
                            /* define IRA range */
  double        metIRA_ThMass, egyIRA_ThMass;
  double        SnII_Step_Prec, LLv_Step_Prec;
  double        referenceZ_toset_SF_DensTh;
  int           SFTh_Zdep, referenceZbin_SFTh;

#ifdef LT_POPIII
  double PopIII_Zlimit;
  int PopIII_IMF_idx;
#endif
  
#ifdef LT_HOT_EJECTA
  double EgySpecEjecta;
#endif

#ifdef LT_STARBURSTS
  int StarBurstCondition;
  double SB_Density_Thresh;
  double SB_DEntropy_Thresh;
#endif

#ifdef LT_DF_BH
  double BH_Radiative_Efficiency;
#endif

#ifdef LT_STOP_COOL_BELOW_Z
  double Below_this_redshift_stop_cooling;
#endif

#ifdef LT_SMOOTH_Z
#if defined(LT_SMOOTH_SIZE) && !defined(LT_SMOOTH_NGB)
  double SmoothRegionSize, SmoothRegionSizeMax;
#endif
#if defined(LT_SMOOTH_NGB) && !defined(LT_SMOOTH_SIZE)
  int DesNumNgbSmooth;
#endif
#endif

  double MaxChemSpreadL;

#endif

#ifdef LT_METAL_COOLING_WAL
  char WalCool_CoolTables_path[200];
#endif

#ifdef LT_ADD_GAL_TO_SUB
  char BC_SED_Path[200];
#endif

#ifdef CHEMCOOL
  int NeedAbundancesForOutput;
  double H2RefDustEff;
  double OxyAbund;
  double CarbAbund;
  double SiAbund;
  double DeutAbund;
  double MgAbund;
  double UVField;
  double PhiPAH;
  double InitDustTemp;
  double DustToGasRatio;
  double AVConversionFactor;
  double CosmicRayIonRate;
  double InitRedshift;
  double ExternalDustExtinction;
  double H2FormEx;
  double H2FormKin;
  int PhotochemApprox;
  int ChemistryNetwork;
  int ADRateFlag;
  int MNRateFlag;
  int AtomicFlag;
  int ThreeBodyFlagA;
  int ThreeBodyFlagB;
  int H3PlusRateFlag;
  int DMAFlag;
  int RadHeatFlag;
  double InitMolHydroAbund;
  double InitHPlusAbund;
  double InitDIIAbund;
  double InitHDAbund;
  double InitHeIIAbund;
  double InitHeIIIAbund;
  double InitCIIAbund;
  double InitSiIIAbund;
  double InitOIIAbund;
  double InitCOAbund;
  double InitC2Abund;
  double InitOHAbund;
  double InitH2OAbund;
  double InitO2Abund;
  double InitHCOPlusAbund;
  double InitCHAbund;
  double InitCH2Abund;
  double InitSiIIIAbund;
  double InitCH3PlusAbund;
  double InitMgPlusAbund;
#endif

#ifdef GENERATE_GAS_IN_ICS
#ifdef GENERATE_GAS_TG
  int GenGasRefFac;
#endif
#endif

#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
          /* used if read initial composition from the parameter file */
  double Startelec;
  double StartHI, StartHII, StartHM;
  double StartHeI, StartHeII, StartHeIII;
  double StartH2I, StartH2II;
  double StartHD, StartDI, StartDII;
  double StartHeHII;
#endif

#ifdef MOL_CLOUDS
  unsigned int MOL_CLOUDS_NumMCs;
  unsigned int MOL_CLOUDS_NumStars;
  double MOL_CLOUDS_MassMCs;
  double MOL_CLOUDS_MassStars;

  MyFloat MOL_CLOUDS_m1_tilde, MOL_CLOUDS_m2_tilde;
#endif

#ifdef ADAPTGRAVSOFT
  int AGS_DesNumNgb;
  double AGS_MaxNumNgbDeviation;
#ifdef PMGRID
  double AGS_MaxSoft, AGS_MinSoft;
#endif
#endif

#if defined(AB_TURB) || defined(TURB_DRIVING)
  double StDecay;
  double StEnergy;
  double StDtFreq;
  double StKmin;
  double StKmax;
  double StSolWeight;
  double StAmplFac;

  int StSpectForm;
  int StSeed;
#endif

#ifdef ADJ_BOX_POWERSPEC 
  double BoxWidth;     /*lenght of the transformation box side*/
  double BoxCenter_x;  /*x coordinate of the box center*/
  double BoxCenter_y;  /*y coordinate of the box center*/
  double BoxCenter_z;  /*z coordinate of the box center*/
  int FourierGrid;     /*dimension of the Fourier transform (actual size is FourierGrid^3)*/
#endif

}
All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern ALIGN(32) struct particle_data
{
  short int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  short int TimeBin;
  MyIDType ID;

  integertime Ti_begstep;		/*!< marks start of current timestep of particle on integer timeline */
  integertime Ti_current;		/*!< current time of the particle */

  ALIGN(32) MyDouble Pos[3];   /*!< particle position at its current time */
  MyDouble Mass;     /*!< particle mass */

  MyDouble Vel[3];   /*!< particle velocity at its current time */
  MyDouble dp[3];

  union
  {
    MyFloat       GravAccel[3];		/*!< particle acceleration due to gravity */
    MyLongDouble dGravAccel[3];
  } g;
#ifdef PMGRID
  MyFloat GravPM[3];		/*!< particle acceleration due to long-range PM gravity force */
#endif
#ifdef FORCETEST
  MyFloat GravAccelDirect[3];	/*!< particle acceleration calculated by direct summation */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union
  {
    MyFloat       Potential;		/*!< gravitational potential */
    MyLongDouble dPotential;
  } p;
#endif

#ifdef FS_TURB_ESTIM
  MyFloat StrFnc[FS_BINS];
  long int StrFnc_count[FS_BINS];
#endif
#ifdef MODGRAV
  MyDouble MassOrig;  /* the mass that particles would have without mass correction */ 
  MyDouble mg_phi;    /* value of the scalar field fluctuation */
  MyDouble mg_grad_phi[3];    /* value of the scalar field fluctuation */
  MyDouble mg_beta;   /* local scalar coupling */
  MyDouble ModGravAccel[3];  /* acceleration due to scalar field fluctuations */
  int mg_treelevel;   /* level of parent AMR node */
#endif

#ifdef DISTORTIONTENSORPS
  MyBigFloat distortion_tensorps[6][6];               /*!< phase space distortion tensor */
  MyBigFloat last_determinant;                        /*!< last real space distortion tensor determinant */
  MyBigFloat stream_density;                          /*!< physical stream density that is going to be integrated */
  double tidal_tensorps[3][3];                        /*!< tidal tensor (=second derivatives of grav. potential) */
  float caustic_counter;                              /*!< caustic counter */
#ifndef GDE_LEAN
  MyBigFloat annihilation;                            /*!< integrated annihilation rate */
  MyBigFloat analytic_annihilation;                   /*!< analytically integrated annihilation rate */
  MyBigFloat rho_normed_cutoff_current;               /*!< current and last normed_cutoff density in rho_max/rho_init * sqrt(sigma) */
  MyBigFloat rho_normed_cutoff_last;
  MyBigFloat s_1_current, s_2_current, s_3_current;   /*! < current and last stretching factor */
  MyBigFloat s_1_last, s_2_last, s_3_last;
  MyBigFloat second_deriv_current;                    /*! < current and last second derivative */
  MyBigFloat second_deriv_last;
  double V_matrix[3][3];                              /*!< initial orientation of CDM sheet the particle is embedded in */
  float init_density;                                 /*!< initial stream density */
  float analytic_caustics;                            /*!< number of caustics that were integrated analytically */
  float a0;
#endif

#ifdef OUTPUT_LAST_CAUSTIC
  MyFloat lc_Time;                                  /*!< time of caustic passage */
  MyFloat lc_Pos[3];                                /*!< position of caustic */
  MyFloat lc_Vel[3];                                /*!< particle velocity when passing through caustic */
  MyFloat lc_rho_normed_cutoff;                     /*!< normed_cutoff density at caustic */
  MyFloat lc_Dir_x[3];                              /*!< principal axis frame of smear out */
  MyFloat lc_Dir_y[3];
  MyFloat lc_Dir_z[3];
  MyFloat lc_smear_x;                               /*!< smear out length */
  MyFloat lc_smear_y;
  MyFloat lc_smear_z;
#endif
#ifdef PMGRID
  double tidal_tensorpsPM[3][3];	            /*!< for TreePM simulations, long range tidal field */
#endif
#endif

  MyFloat OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                          criterion */
#if defined(EVALPOTENTIAL) && defined(PMGRID)
  MyFloat PM_Potential;
#endif

#ifdef STELLARAGE
  MyFloat StellarAge;		/*!< formation time of star particle */
#endif
#ifdef METALS
  MyFloat Metallicity;		/*!< metallicity of gas or star particle */
#endif				/* closes METALS */

#if defined (BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION) || defined(LT_MV_CHEMICALDIFFUSION) || defined(UM_CHEMISTRY)
  MyFloat Hsml;

  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#if defined(RADTRANSFER) || defined(SNIA_HEATING)
  MyFloat DensAroundStar;
#endif
#endif


#if defined(BLACK_HOLES)
  MyIDType SwallowID;
#if !defined(DETACH_BLACK_HOLES)
#ifdef BH_COUNTPROGS
  int BH_CountProgs;
#endif
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
  int     BH_TimeBinGasNeighbor;
#ifdef BH_BUBBLES
  MyFloat BH_Mass_bubbles;
  MyFloat BH_Mass_ini;
#ifdef UNIFIED_FEEDBACK
  MyFloat BH_Mass_radio;
#endif
#endif
  union
  {
    MyFloat BH_Density;
    MyLongDouble dBH_Density;
  } b1;
  union
  {
    MyFloat BH_Entropy;
    MyLongDouble dBH_Entropy;
  } b2;
  union
  {
    MyFloat BH_SurroundingGasVel[3];
    MyLongDouble dBH_SurroundingGasVel[3];
  } b3;
  union
  {
    MyFloat BH_accreted_Mass;
    MyLongDouble dBH_accreted_Mass;
  } b4;
  union
  {
    MyFloat BH_accreted_BHMass;
    MyLongDouble dBH_accreted_BHMass;
  } b5;
  union
  {
    MyFloat BH_accreted_momentum[3];
    MyLongDouble dBH_accreted_momentum[3];
  } b6;
#ifdef BH_BUBBLES
  union
  {
    MyFloat BH_accreted_BHMass_bubbles;
    MyLongDouble dBH_accreted_BHMass_bubbles;
  } b7;
#ifdef UNIFIED_FEEDBACK
  union
  {
    MyFloat BH_accreted_BHMass_radio;
    MyLongDouble dBH_accreted_BHMass_radio;
  } b8;
#endif
#endif
#ifdef KD_FRICTION
  MyFloat BH_SurroundingVel[3];
  MyFloat BH_SurroundingDensity;
#endif
#ifdef KD_FRICTION_DYNAMIC
  MyFloat BH_sigma;
  MyFloat BH_bmax;
#endif
#ifdef KD_TAKE_CENTER_OF_MASS_FOR_BH_MERGER
  MyFloat BH_SwallowPos[3];
#endif
#ifdef LT_DF_BH_BHAR_SWITCH
  MyFloat BlackHoleFeedbackFactor;
#endif
#ifdef LT_BH_GUESSHSML
  MyFloat mean_hsml;
  MyFloat mean_rho;
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  MyFloat ActiveTime;
  MyFloat ActiveEnergy;
#endif
#endif  /* if !defined(DETACH_BLACK_HOLES) */
#endif  /* if deffined(BLACK_HOLES) */ 

#ifdef SIDM
  int ShouldScatterInStep;
  int DidScatterInStep;
  MyIDType ScatterID;
  MyDouble RandX; 
  MyDouble sidm_PSum;
  MyDouble sidm_NumTotalScatter;
  MyDouble sidm_Hsml;
  MyDouble sidm_Density;
  MyDouble sidm_VelDisp;
  MyDouble sidm_NumNgb;
#endif

  
#if defined(DISKPOT) && !defined(SUBFIND) && !defined(ORDER_SNAPSHOTS_BY_ID)
	int GrNr;
#endif
  
#if defined(SUBFIND) 
  int GrNr;
  int SubNr;
  int DM_NumNgb;
  unsigned short targettask, origintask2;
  int origintask, submark, origindex;
  MyFloat DM_Hsml;
  union
  {
    MyFloat DM_Density;
    MyFloat DM_Potential;
  } u;
  union
  {
    MyFloat DM_VelDisp;
    MyFloat DM_BindingEnergy;
  } v;
#ifdef DENSITY_SPLIT_BY_TYPE
  union
  {
    MyFloat int_energy;
    MyFloat density_sum;
  } w;
#endif

#ifdef SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI
  MyFloat DM_Hsml_V;
  MyFloat DM_Density_V;
#endif

#ifdef SAVE_HSML_IN_IC_ORDER
  MyIDType ID_ic_order;
#endif
#ifdef SUBFIND_ALTERNATIVE_COLLECTIVE
  peanokey Key;
#endif
#endif

#if defined(ORDER_SNAPSHOTS_BY_ID) && !defined(SUBFIND)
  int     GrNr;
  int     SubNr;
#endif


#ifdef CS_MODEL
  MyFloat Zm[12];
  MyFloat ZmReservoir[12];
#ifdef CS_FEEDBACK
  MyFloat EnergySN;
  MyFloat EnergySNCold;
#endif
#endif

  float GravCost[GRAVCOSTLEVELS];   /*!< weight factor used for balancing the work-load */

#ifdef WAKEUP
  int dt_step;
#endif

#if defined(LT_STELLAREVOLUTION) || defined(DETACH_BLACK_HOLES)
  union
  {
    unsigned int BHID;
    unsigned int MetID;
  } pt;
#endif

#ifdef SCF_HYBRID
  MyDouble GravAccelSum[3];
  MyFloat MassBackup;
#endif

#ifdef MOL_CLOUDS
  MyFloat MOL_CLOUDS_TimeBorn;
  MyFloat MOL_CLOUDS_LifeTime;

  unsigned int MOL_CLOUDS_index;
#endif

#if defined(JD_VTURB) || defined(HIGH_ORDER_INDUCTION) || defined(KD_RESTRICT_NEIGHBOURS)
  int TrueNGB;			/*!< Number of neighbours inside hsml */
#endif

#ifdef ADAPTGRAVSOFT
  MyFloat AGS_Density;		/* !< (mass/number) density of particle */
#ifdef AGS_OUTPUTNGBS
  int AGS_REALNumNgb;
#endif
  MyFloat AGS_NumNgb;
  MyFloat AGS_Hsml;
  MyFloat AGS_zeta;           /*!< factor in the correction term */
  MyFloat AGS_omega;          /*!< factor in the correction term */
#ifdef AGS_OUTPUTCORR
  MyFloat AGS_corr;          /*!< factor in the correction term */
#endif
#endif
}
 *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */


#ifndef GDE_LEAN
#define GDE_TIMEBEGIN(i) (P[i].a0)
#define GDE_VMATRIX(i, a, b) (P[i].V_matrix[a][b])
#define GDE_INITDENSITY(i) (P[i].init_density)
#else
#define GDE_TIMEBEGIN(i) (All.TimeBegin)
#define GDE_VMATRIX(i, a, b) (0.0)
#define GDE_INITDENSITY(i) (All.GDEInitStreamDensity)
#endif

#if defined(BLACK_HOLES)
#if defined(DETACH_BLACK_HOLES)

#define BPP(i) BHP[P[(i)].pt.BHID]

extern struct bh_particle_data
{
  unsigned int PID;
#ifdef BH_COUNTPROGS
  int BH_CountProgs;
#endif
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
  int     BH_TimeBinGasNeighbor;
#ifdef BH_BUBBLES
  MyFloat BH_Mass_bubbles;
  MyFloat BH_Mass_ini;
#ifdef UNIFIED_FEEDBACK
  MyFloat BH_Mass_radio;
#endif
#endif
  union
  {
    MyFloat BH_Density;
    MyLongDouble dBH_Density;
  } b1;
  union
  {
    MyFloat BH_Entropy;
    MyLongDouble dBH_Entropy;
  } b2;
  union
  {
    MyFloat BH_SurroundingGasVel[3];
    MyLongDouble dBH_SurroundingGasVel[3];
  } b3;
  union
  {
    MyFloat BH_accreted_Mass;
    MyLongDouble dBH_accreted_Mass;
  } b4;
  union
  {
    MyFloat BH_accreted_BHMass;
    MyLongDouble dBH_accreted_BHMass;
  } b5;
  union
  {
    MyFloat BH_accreted_momentum[3];
    MyLongDouble dBH_accreted_momentum[3];
  } b6;
#ifdef BH_BUBBLES
  union
  {
    MyFloat BH_accreted_BHMass_bubbles;
    MyLongDouble dBH_accreted_BHMass_bubbles;
  } b7;
#ifdef UNIFIED_FEEDBACK
  union
  {
    MyFloat BH_accreted_BHMass_radio;
    MyLongDouble dBH_accreted_BHMass_radio;
  } b8;
#endif
#endif
#ifdef KD_FRICTION
  MyFloat BH_SurroundingVel[3];
  MyFloat BH_SurroundingDensity;
#endif
#ifdef KD_FRICTION_DYNAMIC
  MyFloat BH_sigma;
  MyFloat BH_bmax;
#endif
#ifdef KD_TAKE_CENTER_OF_MASS_FOR_BH_MERGER
  MyFloat BH_SwallowPos[3];
#endif
#ifdef LT_DF_BH_BHAR_SWITCH
  MyFloat BlackHoleFeedbackFactor;
#endif
#ifdef LT_BH_GUESSHSML
  MyFloat mean_hsml;
  MyFloat mean_rho;
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  MyFloat ActiveTime;
  MyFloat ActiveEnergy;
#endif
}
  *BHP,
  *DomainBHBuf;
#else

#define BPP(i) P[(i)]

#endif
#endif
                                               /* [----------- start LT block ------------- ]*/
#ifdef LT_STELLAREVOLUTION                     /* [LT] define the structure which hosts the stellar data] */

extern struct met_particle_data
{
  float        LastChemTime;                   /*!< the last and next next time of evolution for Ia and II */
  float        iMass;                          /*!< initial mass of this SSP                               */
  float        Metals[LT_NMetP];               /*!< the metal array (H is not stored here)                 */
  double       weight;                         /*!< used in spreading                                      */
  unsigned int PID;
#ifdef LT_TRACK_CONTRIBUTES
  Contrib      contrib;
#endif
#ifdef LT_ZAGE
  float        ZAge;
#endif
#ifdef LT_ZAGE_LLV
  float        ZAge_llv;
#endif
  int          ChemTimeBin;
#ifdef LT_STARS_GUESSHSML
  MyFloat      mean_hsml;
  MyFloat      mean_rho;
#endif
#ifdef LT_SMOOTH_Z
  float        Zsmooth;
#endif

}
 *MetP,                                 /*!< holds metal particle data on local processor */
 *DomainMetBuf; 			/*!< buffer for metal data in domain decomposition */
/* [----------- end LT block ------------- ]*/
#endif


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  MyDouble Entropy;		/*!< entropy (actually entropic function) of particle */
  MyFloat  EntropyPred;         /*!< predicted value of the entropy at the current time */
  
  MyFloat  Pressure;		/*!< current pressure */
  MyFloat  VelPred[3];		/*!< predicted SPH particle velocity at the current time */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;           /*!< maximum signal velocity */
#endif

#ifdef VORONOI
  MyFloat MaxDelaunayRadius;
  MyFloat Volume;
  MyFloat Center[3];
#ifdef VORONOI_SHAPESCHEME
  MyFloat W;
#endif
#endif

  union
  {
    MyFloat       Density;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensity;
  } d;
  union
  {
    MyFloat       DtEntropy;		/*!< rate of change of entropy */
    MyLongDouble dDtEntropy;
  } e;
  union
  {
    MyFloat       HydroAccel[3];	/*!< acceleration due to hydrodynamical force */
    MyLongDouble dHydroAccel[3];
  } a;
  union
  {
    MyFloat       DhsmlDensityFactor;	/*!< correction factor needed in entropy formulation of SPH */
    MyLongDouble dDhsmlDensityFactor;
  } h;
  union
  {
    MyFloat       DivVel;		/*!< local velocity divergence */
    MyLongDouble dDivVel;
  } v;
#ifndef NAVIERSTOKES
  union
  {
    MyFloat CurlVel;     	        /*!< local velocity curl */
    MyFloat       Rot[3];		/*!< local velocity curl */
    MyLongDouble dRot[3];
  } r;
#else
  union
  {
    MyFloat DV[3][3];
    struct
    {
      MyFloat DivVel;
      MyFloat CurlVel;
      MyFloat StressDiag[3];
      MyFloat StressOffDiag[3];
#ifdef NAVIERSTOKES_BULK
      MyFloat StressBulk;
#endif
    } s;
  } u;
#endif

#if defined(AB_TURB) || defined(VS_TURB) || defined(OUTPUT_VORTICITY) 
   MyFloat Vorticity[3];
   MyFloat SmoothedVel[3];
#endif

#if !(defined(BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING) || defined(LT_STELLAREVOLUTION)) || defined(UM_CHEMISTRY)
  MyFloat Hsml;			/*!< current smoothing length */
  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
  union
  {
    MyFloat       Injected_BH_Energy;
    MyLongDouble dInjected_BH_Energy;
  } i;
#endif

#ifdef COOLING
#if !defined(UM_CHEMISTRY)
  MyFloat Ne;  /*!< electron fraction, expressed as local electron number
		    density normalized to the hydrogen number density. Gives
		    indirectly ionization state and mean molecular weight. */
#endif
#endif
#ifdef SFR
  MyFloat Sfr;
#endif
#ifdef WINDS
  MyFloat DelayTime;		/*!< remaining maximum decoupling time of wind particle */
#ifdef VARIABLE_WINDS
  MyFloat HostHaloMass;
#endif
#endif

#ifdef JD_VTURB
  MyFloat Vrms;		    /*!< RMS velocity inside kernel around Vbulk */
  MyFloat Vbulk[3];	    /*!< Mean velocity inside kernel */
#ifdef JD_DECOMPOSE_VTURB 
  MyFloat Vtan;         /*!< RMS of tangential component of velocity inside kernel around Vbulk */
  MyFloat Vrad;         /*!< RMS of radial component of velocity inside kernel around Vbulk */
#endif
#ifdef JD_DPP
  MyFloat Dpp;			/*!< Reacceleration Coefficient as (Cassano+ '04) */
#endif
#endif
#ifdef JD_RELAX_CLUSTERS
    MyFloat U;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  MyDouble DuDt_diss;
  MyDouble DuDt_drive;
  MyDouble EgyDiss;
  MyDouble EgyDrive;
  MyDouble TurbAccel[3];
#endif
#if defined(TURB_DRIVING)
  MyDouble TurbAccel[3];
#endif

#ifdef AB_SHOCK
  MyFloat Shock_N1[3];           /*!< Parallel to shock >*/
  MyFloat Shock_Up_V1[3];        /*!< Upstream parallel velocity >*/
  MyFloat Shock_Down_V1[3];      /*!< Downstream parallel velocity >*/
  MyFloat Shock_Up_Signal;       /*!< Upstream signal velocity >*/
  MyFloat Shock_Down_Signal;     /*!< Downstream signal velocity >*/
  MyFloat Shock_Mach;            /*!< Machnumber >*/
  MyFloat Shock_Mach1;            /*!< Machnumber >*/
  MyFloat Shock_Mach2;            /*!< Machnumber >*/
  MyFloat Shock_Mach3;            /*!< Machnumber >*/
  MyFloat Shock_Speed;           /*!< Shock speed >*/
  MyFloat Shock_Compress;        /*!< Compression ratio >*/
  MyFloat Shock_Up_Weight;       /*!< Upstream weight >*/
  MyFloat Shock_Down_Weight;     /*!< Downstream weight >*/
  MyFloat Shock_Up_Rho;          /*!< Upstream density >*/
  MyFloat Shock_Down_Rho;        /*!< Downstream density >*/
  MyFloat Shock_Up_Pressure;     /*!< Upstream pressure >*/
  MyFloat Shock_Down_Pressure;   /*!< Downstream pressure >*/
#ifdef AB_SHOCK_VELDIV
  MyFloat Shock_N2[3];           /*!< Perpendicular to shock >*/
  MyFloat Shock_N3[3];           /*!< Perpendicular to shock >*/
  MyFloat Shock_Up_V2[3];        /*!< Upstream perpendicular velocity >*/
  MyFloat Shock_Down_V2[3];      /*!< Downstream perpendicular velocity >*/
  MyFloat Shock_Up_V3[3];        /*!< Upstream perpendicular velocity >*/
  MyFloat Shock_Down_V3[3];      /*!< Downstream perpendicular velocity >*/
  MyFloat Shock_Up_Weight23[2];  /*!< Upstream perpendicular weights >*/
  MyFloat Shock_Down_Weight23[2];/*!< Downstream perpendicular weights >*/
#endif
#endif /* AB_SHOCK */

#ifdef FS_TURB_ESTIM
  MyFloat Vturb;
  MyFloat Lturb;
#endif
#ifdef MAGNETIC
  union
  {
    MyFloat B[3];
#if defined(EULERPOTENTIALS)
    MyFloat dEulerA[3];
#endif
  } b1;

  union
  {  
    MyFloat BPred[3];
#if defined(EULERPOTENTIALS)
    MyFloat dEulerB[3];
#endif
  } b2;

#ifdef MAGNETIC_SN_SEEDING
  MyFloat MagSeed[3];
#endif

#ifdef HIGH_ORDER_INDUCTION
  MyFloat Chi[3][3];
  MyFloat Xix[3], Xiy[3], Xiz[3];
#endif
#ifdef FS_ALFA2_DYN
  MyFloat alfa2;
#endif
#ifdef FS_ETA2_DYN
  MyFloat eta2;
#endif
#ifdef DIVBFORCE3
  MyFloat magacc[3];
  MyFloat magcorr[3];
#endif
#ifdef EULERPOTENTIALS
  MyFloat EulerA,EulerB;
#ifdef EULER_DISSIPATION
  MyFloat DtEulerA,DtEulerB;
#endif
#endif
#if !defined(EULERPOTENTIALS)
  MyFloat DtB[3];
#endif
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP) || defined(DIVBCLEANING_DEDNER)
  MyFloat divB;
#endif
#if defined(BSMOOTH) || defined(BFROMROTA) || defined(BSMOOTH_TIME)
  MyFloat BSmooth[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha, DtBalpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat Phi, PhiPred, DtPhi;
  MyFloat GradPhi[3];
#ifdef SMOOTH_PHI
  MyFloat SmoothPhi;
#endif
#endif
#if defined(SMOOTH_DIVB) || defined(DIVBCLEANING_DEDNER)
  MyFloat SmoothDivB;
#endif

#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB) || defined(FS_ETA2_DYN)
  MyFloat RotB[3];
#ifdef SMOOTH_ROTB
  MyFloat SmoothedRotB[3];
#endif
#endif

#endif /* MAGNETIC */

#ifdef VSMOOTH
  MyFloat VSmooth[3];
  int SmoothNgb;
#endif

#if (defined(SMOOTH_DIVB) || defined(DIVBCLEANING_DEDNER) || defined(BSMOOTH) || defined(SMOOTH_ROTB) || ((defined(LT_STELLAREVOLUTION) && !defined(LT_DONTUSE_DENSITY_in_WEIGHT)) || defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_XCLD) || defined(LT_TRACK_WINDS))) || defined(VSMOOTH)
  MyFloat DensityNorm;
#endif

#ifdef TIME_DEP_ART_VISC
  MyFloat alpha, Dtalpha;
#endif
#ifdef NS_TIMESTEP
  MyFloat ViscEntropyChange;
#endif
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif

#ifdef MHM
  MyFloat FeedbackEnergy;
#endif

#ifdef COSMIC_RAYS
  MyFloat CR_C0[NUMCRPOP];			/*!< Cosmic ray amplitude adiabatic invariable */
  MyFloat CR_q0[NUMCRPOP];			/*!< Cosmic ray cutoff adiabatic invariable */
  MyFloat CR_E0[NUMCRPOP];			/*!< Specific Energy at Rho0 */
  MyFloat CR_n0[NUMCRPOP];			/*!< baryon fraction in cosmic rays */

  MyFloat CR_DeltaE[NUMCRPOP];		/*!< Specific Energy growth during timestep */
  MyFloat CR_DeltaN[NUMCRPOP];		/*!< baryon fraction growth during timestep */
#ifdef MACHNUM
  MyFloat CR_Gamma0[NUMCRPOP];
#endif

#ifdef CR_OUTPUT_INJECTION
  MyFloat CR_Specific_SupernovaHeatingRate;
#endif
#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  MyFloat Shock_MachNumber;	/*!< Mach number */
  MyFloat Shock_DecayTime;	/*!< Shock decay time */
#ifdef COSMIC_RAYS
  MyFloat Shock_DensityJump;	/*!< Density jump at the shock */
  MyFloat Shock_EnergyJump;	/*!< Energy jump at the shock */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
  MyFloat PreShock_PhysicalEnergy;	/*!< Density in the preshock regime */
  MyFloat PreShock_XCR;		/*!< XCR = PCR / Pth in the preshock regime */
#endif
#ifdef MACHSTATISTIC
  MyFloat Shock_DtEnergy;		/*!< Change of thermal specific energy at Shocks */
#endif
#ifdef OUTPUT_PRESHOCK_CSND
  MyFloat PreShock_PhysicalSoundSpeed;	/*!< Sound speed in the preshock regime */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
#endif
#endif				/* Mach number estimate */


#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  MyFloat elec;
  MyFloat HI;
  MyFloat HII;

  MyFloat HeI;
  MyFloat HeII;
  MyFloat HeIII;

  MyFloat H2I;
  MyFloat H2II;

  MyFloat HM;

  MyFloat Gamma;
  MyFloat t_elec, t_cool;

#ifdef UM_CHEMISTRY
  MyFloat Um_MeanMolecularWeight;
  MyFloat HD;
  MyFloat DI;
  MyFloat DII;
  MyFloat HeHII;
#endif

#endif

#if defined(RADTRANSFER)
  MyFloat ET[6];                /* eddington tensor - symmetric -> only 6 elements needed */
  MyFloat Je[N_BINS];           /* emmisivity */

#ifndef UM_CHEMISTRY
  MyFloat HI;                  /* HI fraction */
  MyFloat HII;                 /* HII fraction */
  MyFloat elec;               /* electron fraction */
#ifdef RT_INCLUDE_HE
  MyFloat HeI;                 /* HeI fraction */
  MyFloat HeII;                 /* HeII fraction */
  MyFloat HeIII;                 /* HeIII fraction */
#endif
#endif

  MyFloat n_gamma[N_BINS];

#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat Grad_ngamma[3][N_BINS];
#endif

#ifdef RT_RAD_PRESSURE
  MyFloat n[3];
  MyFloat RadAccel[3];
#endif

#ifdef SFR
  MyDouble DensitySfr;
  MyDouble HsmlSfr;
  MyDouble DhsmlDensityFactorSfr;
  MyDouble NgbSfr;
#endif

#endif //radtransfer 

#if defined CS_MODEL
  MyFloat DensityOld;
#ifdef CS_FEEDBACK
  union
  {
    MyFloat       DensityAvg;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensityAvg;
  } da;
  union
  {
    MyFloat       EntropyAvg;		/*!< current baryonic mass density of particle */
    MyLongDouble dEntropyAvg;
  } ea;
  MyFloat HotHsml;
  int     HotNgbNum;
  MyFloat DensPromotion;
  MyFloat TempPromotion;
#endif
#endif

#ifdef EOS_DEGENERATE
  MyFloat u;                            /* internal energy density */
  MyFloat temp;                         /* temperature */
  MyFloat dpdr;							/* derivative of pressure with respect to density at constant entropy */
  MyFloat xnuc[EOS_NSPECIES];           /* nuclear mass fractions */
  MyFloat dxnuc[EOS_NSPECIES];          /* change of nuclear mass fractions */
  MyFloat xnucPred[EOS_NSPECIES];
#endif

#ifdef WAKEUP
  short int wakeup;             /*!< flag to wake up particle */
#endif

#ifdef BP_REAL_CRs
  MyFloat CRpNorm[BP_REAL_CRs];         /*!< normalization of CR protons spectrum */
  MyFloat CRpSlope[BP_REAL_CRs];        /*!< slope of CR protons spectrum */
  MyFloat CRpCut;                       /*!< cutoff of CR protons spectrum  */
  MyFloat CRpN[BP_REAL_CRs];            /*!< number of CR p */
  MyFloat CRpE[BP_REAL_CRs];            /*!< energy of CR p */
  MyFloat CRpPressure;                  /*!< pressure of CR p */
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
  MyFloat DtCRpE[BP_REAL_CRs];		/*!< time derivative of CR p energy */
  MyFloat DtCRpN[BP_REAL_CRs];		/*!< time derivative of CR p number */
#endif
  MyFloat CReNorm[BP_REAL_CRs];         /*!< normalization of CR electrons spectrum */
  MyFloat CReSlope[BP_REAL_CRs];        /*!< slope of CR electrons spectrum */
  MyFloat CReCut;                       /*!< cutoff of CR electrons spectrum  */
  MyFloat CReN[BP_REAL_CRs];            /*!< number of CR e */
  MyFloat CReE[BP_REAL_CRs];            /*!< energy of CR e */
  MyFloat CRePressure;                  /*!< pressure of CR e */
  MyFloat DensityOld;
#endif

#if (defined(SFR) && defined(MAGNETIC)) || defined(LT_STELLAREVOLUTION)
  float XColdCloud;
#endif
#if defined(LT_STELLAREVOLUTION)
  float Temperature;
#endif

#ifdef LT_STELLAREVOLUTION
  MyFloat  MassRes;
  MyDouble EgyRes;                      /*!< external (Sn) energy resorvoir */
  float    Metals[LT_NMetP];            /*!< H is not stored here */
#ifndef LT_LOCAL_IRA
  double   mstar;
#endif
#ifdef LT_TRACK_CONTRIBUTES
  Contrib  contrib;
#endif
#ifdef LT_ZAGE
  MyFloat  ZAge, ZAgeW;
#endif
#ifdef LT_ZAGE_LLV
  MyFloat  ZAge_llv, ZAgeW_llv;
#endif
#ifdef LT_TRACK_WINDS
  MyFloat  AvgHsml;
#endif

#endif

#ifdef LT_MV_CHEMICALDIFFUSION
  
  /* int DiffusionIntSteps; */
  /* double DiffusionIntTime; */
  
  MyFloat A, D, B[LT_NMetP];               /*!< The two coefficients in the Greif et al. approximation */
  MyFloat GradAbnd[LT_NMetP][3];
  MyFloat Abundance[LT_NMetP];             /*!< The quantity to be diffused */
  MyFloat deltaMetals[LT_NMetP];           /*!< needed to renormalize  */
  MyFloat AvgAbnd[LT_NMetP];               /*!< only for testing purposes  */

  MyDouble prevMass;
  union
  {
    MyFloat       vel_rms;
    MyDouble      dvel_rms;
  } v_rms;           /*!< rms velocity dispersion within its smoothing length */
  int NeighboringTo, TaskNeighTo, NeighboredBy, TaskNeighBy;
  int ChemFof;
  /*  int vrms_NGB;*/
#endif

  
#if defined(LT_SMOOTH_Z)
#if defined(LT_SMOOTH_ALLMETALS)
  MyFloat Zsmooth[LT_NMetP];
#else
  MyFloat Zsmooth;
  MyFloat Zsmooth_a;
  MyFloat Zsmooth_b;
#if defined(LT_SMOOTH_SIZE) || defined(LT_SMOOTH_NGB)
  float SmoothDens;
  float SmoothDens_b;
  int SmoothNgb;
#endif
#if defined(LT_SMOOTH_NGB)
  float SmoothHsml;
#endif
#endif                                 /* closes LT_METAL_COOLING_WAL  */
#endif                                 /* closes LT_SMOOTH_Z */

#ifdef LT_SMOOTH_XCLD                  /* smooth the cloud fraction */
  float XCLDsmooth;
#endif

#ifdef CHEMCOOL
  double TracAbund[TRAC_NUM];
#endif

#if defined(BLACK_HOLES) && defined(LT_BH_ACCRETE_SLICES)
  int NSlicesSwallowed;
#endif
}
  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system
*/
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
#if defined(MAGNETIC_STATISTICS) || defined(FS_ALFA2_QUENCH)
    EnergyMag,
#endif
#if defined(FS_ALFA2_TURB) 
    DensestPart[4],
#endif
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6],
    EnergyTotComp[6],
    MomentumComp[6][4],
    AngMomentumComp[6][4],
    CenterOfMassComp[6][4];
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
    double DivBerr;
#endif
#if defined(FS_ALFA2_QUENCH)
    double MaxAlfa2;
#endif
}
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

extern struct data_index
{
  int Task;
  int Index;
  int IndexGet;
}
 *DataIndexTable;		/*!< the particles to be exported are grouped
                                  by task-number. This table allows the
                                  results to be disentangled again and to be
                                  assigned to the correct particle */

extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;

extern struct gravdata_in
{
  MyFloat Pos[3];
#ifdef FS_TURB_ESTIM
  MyFloat Vel[4];
#endif
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD) || defined(MODGRAV)
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat Soft;
#endif
#endif
  MyFloat OldAcc;
  int NodeList[NODELISTLENGTH];
#ifdef MODGRAV
  MyFloat beta;
#endif

#ifdef ADAPTGRAVSOFT
  MyFloat AGS_zeta, AGS_omega, AGS_Hsml, Mass;
#endif
}
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


extern struct gravdata_out
{
  MyLongDouble Acc[3];
#ifdef FS_TURB_ESTIM
  MyFloat StrFnc[FS_BINS];
  int StrFnc_count[FS_BINS];
#endif
#ifdef EVALPOTENTIAL
  MyLongDouble Potential;
#endif
#ifdef DISTORTIONTENSORPS
  MyLongDouble tidal_tensorps[3][3];
#endif
#ifdef MODGRAV
  MyLongDouble mg_grad_phi_tree[3];
  MyLongDouble mg_phi_tree;
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
  MyFloat AGS_corr;
#endif
}
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct potdata_out
{
  MyLongDouble Potential;
}
 *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct info_block
{
  char label[4];
  char type[8];
  int ndim;
  int is_present[6];
}
*InfoBlock;


/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
#ifdef COSMIC_RAYS
  double SpectralIndex_CR_Pop[NUMCRPOP]; /*!< spectral indices of cosmic ray populations */
#endif
  double time;			/*!< time of snapshot file */
  double redshift;		/*!< redshift of snapshot file */
  int flag_sfr;			/*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/*!< flags whether cooling was included  */
  int num_files;		/*!< number of files in multi-file snapshot */
  double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

#ifdef COSMIC_RAYS
  char fill[18-8*NUMCRPOP];	/*!< fills to 256 Bytes */
#else
  char fill[18];		/*!< fills to 256 Bytes */
#endif

  char names[15][2];
}
header;				/*!< holds header for snapshot files */


enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_SECONDORDERMASS,
  IO_U,
  IO_RHO,
  IO_NE,
  IO_NH,
  IO_HSML,
  IO_VALPHA,
  IO_SFR,
  IO_AGE,
  IO_HSMS,
  IO_ACRS,
  IO_Z,
  IO_BHMASS,
  IO_BHMDOT,
  IO_BHPROGS,
  IO_BHMBUB,
  IO_BHMINI,
  IO_BHMRAD,
  IO_ACRB,
  IO_POT,
  IO_ACCEL,
  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,
  IO_HD,
  IO_DI,
  IO_DII,
  IO_HeHII,
  IO_DTENTR,
  IO_STRESSDIAG,
  IO_STRESSOFFDIAG,
  IO_STRESSBULK,
  IO_SHEARCOEFF,
  IO_TSTP,
  IO_BFLD,
  IO_BSMTH,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_VTURB,
  IO_LTURB,
  IO_ALFA2_DYN,
  IO_ETA2_DYN,
  IO_PHI,
  IO_XPHI,
  IO_GRADPHI,
  IO_ROTB,
  IO_SROTB,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_DENN,
  IO_EGYPROM,
  IO_EGYCOLD,
  IO_MACH,
  IO_MACH1,
  IO_MACH2,
  IO_MACH3,
  IO_RUP,
  IO_RDOWN,
  IO_PUP,
  IO_PDOWN,
  IO_SHSPEED,
  IO_SHCOMPRESS,
  IO_SHNORMAL,
  IO_VUP,
  IO_VDOWN,
  IO_DTENERGY,
  IO_PRESHOCK_CSND,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,
  IO_CRINJECT,
  IO_TIDALTENSORPS,
  IO_DISTORTIONTENSORPS,
  IO_EULERA,
  IO_EULERB,
  IO_FLOW_DETERMINANT,
  IO_PHASE_SPACE_DETERMINANT,
  IO_ANNIHILATION_RADIATION,
  IO_STREAM_DENSITY,
  IO_EOSTEMP,
  IO_EOSXNUC,
  IO_PRESSURE,
  IO_RADGAMMA,
  IO_RAD_ACCEL,
  IO_EDDINGTON_TENSOR,
  IO_LAST_CAUSTIC,
  IO_SHEET_ORIENTATION,
  IO_INIT_DENSITY,
  IO_CAUSTIC_COUNTER,
  IO_DMHSML,                    /* for 'SUBFIND_RESHUFFLE_CATALOGUE' option */
  IO_DMDENSITY,
  IO_DMVELDISP,
  IO_DMHSML_V,                 /* for 'SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI' option */
  IO_DMDENSITY_V,
  IO_VRMS,
  IO_VBULK,
  IO_VRAD,
  IO_VTAN,
  IO_TRUENGB,
  IO_VDIV,
  IO_VROT,
  IO_VORT,
  IO_DPP,
  IO_BPCR_pNORM,
  IO_BPCR_eNORM,
  IO_BPCR_pSLOPE,
  IO_BPCR_eSLOPE,
  IO_BPCR_pE,
  IO_BPCR_pN,
  IO_BPCR_ePRESSURE,
  IO_BPCR_pPRESSURE,

  IO_iMass,
  IO_Zs,
  IO_ZAGE,
  IO_ZAGE_LLV,
  IO_CLDX,
  IO_HTEMP,
  IO_TEMP,
  IO_CONTRIB,
  IO_ZSMOOTH,
  IO_allZSMOOTH,
  IO_CHEM,
  IO_DELAYTIME,

  IO_ABUNDANCE,
  IO_ABND_GRAD,
  IO_ABND_SPH,
  IO_DIFFUSING_CB,
  IO_DIFFUSING_CA,
  IO_DIFFUSING_CD,
  IO_DIFF_TSTEP,
  
  IO_PSUM,
  IO_SIDMNUMNGB, 
  IO_NUMTOTALSCATTER,
  IO_SIDMHSML,
  IO_SIDMDENSITY,
  IO_SIDMVELDISP,

  IO_AGS_SOFT,
  IO_AGS_DENS,
  IO_AGS_ZETA,
  IO_AGS_OMEGA,
  IO_AGS_CORR,
  IO_AGS_NGBS,

  IO_VSTURB_DISS,
  IO_VSTURB_DRIVE,
  
  IO_MG_PHI,
  IO_MG_ACCEL,
  
  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};

enum siofields
{ SIO_GLEN,
  SIO_GOFF,
  SIO_MTOT,
  SIO_GPOS,
  SIO_MMEA,
  SIO_RMEA,
  SIO_MCRI,
  SIO_RCRI,
  SIO_MTOP,
  SIO_RTOP,
  SIO_DMEA,
  SIO_DCRI,
  SIO_DTOP,
  SIO_MGAS,
  SIO_MSTR,
  SIO_TGAS,
  SIO_LGAS,
  SIO_NCON,
  SIO_MCON,
  SIO_BGPOS,
  SIO_BGMTOP,
  SIO_BGRTOP,
  SIO_NSUB,
  SIO_FSUB,
  SIO_SLEN,
  SIO_SOFF,
  SIO_PFOF,
  SIO_MSUB,
  SIO_SPOS,
  SIO_SVEL,
  SIO_SCM,
  SIO_SPIN,
  SIO_DSUB,
  SIO_VMAX,
  SIO_RVMAX,
  SIO_RHMS,
  SIO_MBID,
  SIO_GRNR,
  SIO_SMST,
  SIO_SLUM,
  SIO_SLATT,
  SIO_SLOBS,
  SIO_DUST,
  SIO_SAGE,
  SIO_SZ,
  SIO_SSFR,
  SIO_PPOS,
  SIO_PVEL,
  SIO_PTYP,
  SIO_PMAS,
  SIO_PAGE,
  SIO_PID,

  SIO_LASTENTRY
};

/*
 * Variables for Tree
 * ------------------
 */

extern int Nexport, Nimport;
extern int BufferFullFlag;
extern int NextParticle;
extern int NextJ;
extern int TimerFlag;

// note, the ALIGN(32) directive will effectively pad the structure size
// to a multiple of 32 bytes
extern ALIGN(32) struct NODE
{
  MyFloat center[3];		/*!< geometrical center of node */
  MyFloat len;			/*!< sidelength of treenode */

#ifdef MODGRAV
  struct
#else
  union
#endif
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      MyFloat s[3];		/*!< center of mass of node */
      MyFloat mass;		/*!< mass of node */
      unsigned int bitflags;	/*!< flags certain node properties */
      int sibling;		/*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;		/*!< this gives the next node in case the current node needs to be opened */
      int father;		/*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;

  double GravCost;
  integertime Ti_current;

#ifdef RADTRANSFER
  MyFloat stellar_s[3];         /*!< enter of mass for the stars in the node*/
  MyFloat stellar_mass;         /*!< mass in stars in the node*/
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTGRAVSOFT)
  MyFloat maxsoft;		/*!< hold the maximum gravitational softening of particle in the
				             node if either the ADAPTIVE_GRAVSOFT_FORGAS or the ADAPTGRAVSOFT option is selected */
#endif

#ifdef MODGRAV
  MyFloat mass_cic;

  MyFloat phi;
  MyFloat grad_phi[3];
  MyFloat s_eff[3];
  MyFloat eff_mass;
  
  MyFloat phi_mapped_up;   /*!< phi value that was mapped up from the finer level */
  MyFloat fsource;   /*!< source term for linear equation that is iteratively solved */
  MyFloat Lphi;   /*!< value of L phi */
  int iterate;   /*!< whether this cell should be iterated */
  
  int level;      /*!< starting with one for largest node */
  int amr_node;   /*!< 1 for lowest level nodes that have 7 siblings */
  int num_daughters;   /*!< number of regular daughter nodes  */
  
  MyFloat len_orig;
  MyFloat volume;
  
  int nxm1, nxp1, nym1, nyp1, nzm1, nzp1;   /*!< neighbor nodes on same or coarser grid level */
  int red;   /*!< flag for red-black sweep */
  
  int always_open;
#endif
  
#ifdef SCALARFIELD
  MyFloat s_dm[3];
  MyFloat mass_dm;
#endif
  
#ifdef ADAPTGRAVSOFT
 int Flag;
#endif
}
 *Nodes_base,			/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
				   gives the first allocated node */


extern struct extNODE
{
  MyLongDouble dp[3];
#ifdef GRAVITY_CENTROID
  int suns[8];
#endif
#ifdef SCALARFIELD
  MyLongDouble dp_dm[3];
  MyFloat vs_dm[3];
#endif
#ifdef FLTROUNDOFFREDUCTION
  MyFloat s_base[3];
  MyFloat len_base;
#ifdef SCALARFIELD
  MyFloat s_dm_base[3];
#endif
#endif
  MyFloat vs[3];
  MyFloat vmax;
  MyFloat hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
#ifndef KD_HMAX_ESTIMATE
  MyFloat divVmax;
#endif
  integertime Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;		/*!< gives next node in tree walk  (nodes array) */
extern int *Father;		/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif

#if defined  (UM_METAL_COOLING)
/* --==[ link with LT_ stuffs]==-- */
extern float *um_ZsPoint, um_FillEl_mu, um_mass;
#endif


#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
extern double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
extern double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
  k10a[N_T], k11a[N_T];
extern double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
  k20a[N_T], k21a[N_T];
extern double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
extern double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
extern double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu], sigma30[N_nu], sigma31[N_nu];
#endif
#endif


#ifdef UM_CHEMISTRY
extern double kHeHII1a[N_T],kHeHII2a[N_T],kHeHII3a[N_T];
extern double kHD1a[N_T],kHD2a[N_T],kHD3a[N_T],kHD4a[N_T],kHD5a[N_T],kHD6a[N_T],kHD7a[N_T];
#endif


#if defined (UM_CHEMISTRY) && defined (CHEMISTRY)
#error you cannot define both UM_CHEMISTRY and CHEMISTRY !
#endif


#ifdef UM_CHEMISTRY
#define T_SUP_INTERPOL_LIMIT        1.e4
#endif


#ifdef CHEMCOOL
extern struct{
  double temptab[NMD];
  double cltab[NMD][NCLTAB];
  double chtab[NMD][NCHTAB];
  double dtcltab[NMD][NCLTAB];
  double dtchtab[NMD][NCHTAB];
  double crtab[NCRTAB];
  double crphot[NCRPHOT];
  double phtab[NPHTAB];
  double cst[NCONST];
  double dtlog;
  double tdust;
  double tmax;
  double tmin;
  double deff;
  double abundc;
  double abundo;
  double abundsi;
  double abundD;
  double abundmg;
  double G0;
  double f_rsc;
  double phi_pah;
  double dust_to_gas_ratio;
  double AV_conversion_factor;
  double cosmic_ray_ion_rate;
  double redshift;
  double AV_ext;
  double pdv_term;
  double h2_form_ex;
  double h2_form_kin;
  double dm_density;
}COOLR;

extern struct{
  int iphoto;
  int iflag_mn;
  int iflag_ad;
  int iflag_atom;
  int iflag_3bh2a;
  int iflag_3bh2b;
  int iflag_h3pra;
  int iflag_h2opc;
  int id_current;
  int index_current;
  int idma_mass_option;
  int no_chem;
  int irad_heat;
}COOLI;

#endif

#ifdef SCFPOTENTIAL
extern long scf_seed;
extern MyDouble *Anltilde, *coeflm, *twoalpha, *c1, *c2, *c3;
extern MyDouble *cosmphi, *sinmphi;
extern MyDouble *ultrasp, *ultraspt, *ultrasp1;
extern MyDouble *dblfact, *plm, *dplm;
extern MyDouble *sinsum, *cossum;
extern MyDouble *sinsum_all, *cossum_all;
#ifdef SCF_SCALEFAC
extern float *scalefac_nlm;
#endif
#endif


extern int maxThreads;

#if defined(AB_TURB) || defined(TURB_DRIVING)
//parameters
extern double StDecay;
extern double StEnergy;
extern double StDtFreq;
extern double StKmin;
extern double StKmax;
extern double StSolWeight;
extern double StAmplFac;
 
//Ornstein-Uhlenbeck variables
extern int StSeed;
extern double StOUVar;
extern double* StOUPhases;


//forcing field in fourie space
extern double* StAmpl;
extern double* StAka; //phases (real part)
extern double* StAkb; //phases (imag part)
extern double* StMode;
extern int StNModes;

extern integertime StTPrev;
extern double StSolWeightNorm;

extern int StSpectForm;
#endif
#ifdef AB_TURB
extern gsl_rng* StRng;
#endif
#ifdef TURB_DRIVING
extern int FB_Seed;
#endif


#ifdef KD_ALTERNATIVE_GROUP_SORT
extern int MaxNgroups;
#endif

#endif  /* ALLVARS_H  - please do not put anything below this line */








