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

#include "allvars.h"



#ifdef PERIODIC
MyDouble boxSize, boxHalf, inverse_boxSize;

#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X, inverse_boxSize_X;
#else
#endif
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y, inverse_boxSize_Y;
#else
#endif
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z, inverse_boxSize_Z;
#else
#endif
#endif

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/



int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
int PTask;			/*!< note: NTask = 2^PTask */

#ifdef INVARIANCETEST
int World_ThisTask;
int World_NTask;
int Color;
MPI_Comm MPI_CommLocal;
#endif

double CPUThisRun;		/*!< Sums CPU time of current process */

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */
long long GlobNumForceUpdate;
int NumSphUpdate;		/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
int RestartSnapNum;
int SelRnd;

int *Exportflag;		/*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount;

#ifdef VORONOI
int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;
#endif

int TakeLevel;

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
int NActivePart;
int *ActiveParticleList;
#if defined(LT_STELLAREVOLUTION)
int NActivePartMerk;
#endif
#endif

int FirstActiveParticle;
int *NextActiveParticle;
unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountSph[TIMEBINS];
int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];
int LastInTimeBin[TIMEBINS];
int *NextInTimeBin;
int *PrevInTimeBin;

size_t HighMark_run, HighMark_domain, HighMark_gravtree,
  HighMark_pmperiodic, HighMark_pmnonperiodic, HighMark_sphdensity, HighMark_sphhydro, HighMark_addSPH;

#ifdef VORONOI
size_t HighMark_voronoi;
#endif

#ifdef CS_MODEL
#if defined(CS_SNI) || defined(CS_SNII)
size_t HighMark_cs_enrich;
#endif
size_t HighMark_cs_hotngbs;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
size_t HighMark_turbpower;
#endif

#ifdef ADAPTGRAVSOFT
size_t HighMark_agsdensity;
#endif

#ifdef SFR
double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
double TimeBin_BH_mass[TIMEBINS];
double TimeBin_BH_dynamicalmass[TIMEBINS];
double TimeBin_BH_Mdot[TIMEBINS];
double TimeBin_BH_Medd[TIMEBINS];
#endif

#ifdef RADTRANSFER
double lum[N_BINS];
#ifdef RT_POPIII
double lum_popIII[N_BINS];
#endif
double rt_sigma_HI[N_BINS];
double rt_sigma_HeI[N_BINS];
double rt_sigma_HeII[N_BINS];
#ifndef UM_CHEMISTRY
double nu[N_BINS];
#endif
#endif


char DumpFlag = 1;

size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;

double CPU_Step[CPU_PARTS];
char CPU_Symbol[CPU_PARTS] =
  { '-', '*', '=', ';', '<', '[', '^', ':', '.', '~', '|', '+', '"', '/', '`', ',', '>', '@', '#', '&', '$',
  ']', '(', '?', ')', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '\\', '?', '{', '}'
};
char CPU_SymbolImbalance[CPU_PARTS] =
  { 'a', 't', 'u', 'v', 'b', 'w', 'd', 'r', 'h', 'm', 'n', 'l', 'o', 'p', 's', 'f', 'i', 'g', 'c', 'e', 'x',
  'y', 'z', 'A', 'I', 'W', 'T', 'V', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'
};
char CPU_String[CPU_STRING_LEN + 1];

double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */

int Flag_FullStep;		/*!< Flag used to signal that the current step involves all particles */


int TreeReconstructFlag;
int GlobFlag;


int NumPart;			/*!< number of particles on the LOCAL processor */
int N_gas;			/*!< number of gas particles on the LOCAL processor  */
#ifdef LT_STELLAREVOLUTION
int N_stars;			/*!< number of star particles in the LOCAL processor */
#endif
#if defined(BLACK_HOLES) && defined(DETACH_BLACK_HOLES)
int N_BHs;
#endif

#ifdef SINKS
int NumSinks;
#endif

long long Ntype[6];		/*!< total number of particles of each type */
int NtypeLocal[6];		/*!< local number of particles of each type */

gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef SFR
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif

#ifdef MODGRAV
int NumSavedAMRCells; /*!< number of AMR cells form old tree saved on LOCAL task */
int NumAMRCells; /*!< number of AMR cells on LOCAL processor */

int NumSavedMultigridCells;  /*!< number of multigird cells from old tree saved on LOCAL task */
int SavedMultigridMinLevel;  /*!< minimum tree level for which multigrid cells have been saved */
int SavedMultigridMaxLevel;  /*!< maximum tree level for which multigrid cells have been saved */

int FirstAMRCell; /*!< index of first AMR node */
int *NextAMRCell, *NextAMRCell_base; /*!< defined for all nodes, for AMR nodes it gives the index of the next AMR node */

int phi_guess_available; 

struct amr_cell_data *saved_amr_cells, *saved_multigrid_cells;

double *mg_time_multi_grid_nbg_search;
double *mg_time_multi_grid_update_cells;
double *mg_time_multi_grid_update_cells_comm_imbal;
int *mg_num_it_multi_grid;
long long *mg_num_cells_multi_grid;
#endif

double TimeOfLastTreeConstruction;	/*!< holds what it says */

int *Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */
double *R2ngblist;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int *DomainStartList, *DomainEndList;



double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;
int *DomainNodeIndex;
int *DomainList, DomainNumChanged;

peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;


#ifdef VORONOI
struct individual_data Indi;
#endif

double RndTable[RNDTABLE];

#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif


#ifdef WRITE_KEY_FILES
peanokey *KeyIndex;
int *NPartPerKey, *PartKeyOffset;
int NKeys[6];
long long NKeysTotal[6];
#endif

#ifdef LT_ADD_GAL_TO_SUB
float *tempiAS, *CB07, *Filters_Effective_Wavelenght;
#ifdef OBSERVER_FRAME
float *CB07obs;
#endif
#endif


/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,			/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU,			/*!< file handle for cpu.txt log-file. */
 *FdTimebin;

#ifdef SCFPOTENTIAL
FILE *FdSCF;
#endif

#ifdef SFR
FILE *FdSfr;			/*!< file handle for sfr.txt log-file. */
#endif

#if defined (VS_TURB) || defined (AB_TURB)
FILE *FdTurb;
#endif

#ifdef RADTRANSFER
FILE *FdRad;			/*!< file handle for radtransfer.txt log-file. */
FILE *FdStar;			/*!< file handle for lum_star.txt log-file. */
#endif

#ifdef DISTORTIONTENSORPS
#ifdef PMGRID
FILE *FdTidaltensor;		/*!< file handle for Tidaltensor.txt log-file. */
#endif
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;		/*!< file handle for blackholes.txt log-file. */
FILE *FdBlackHolesDetails;
#endif


#ifdef FORCETEST
FILE *FdForceTest;		/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef MODGRAV
FILE *FdMOG;			/*!< file handle for modgrav.txt log-file. */
FILE *FdMultiGridTimes;  /*!< file handle for multi_grid_times.txt log-file */
#endif

#ifdef DARKENERGY
FILE *FdDE;			/*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
FILE *FdXXL;			/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
double MeanB;

#ifdef TRACEDIVB
double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
double HydroKickTable[DRIFT_TABLE_LENGTH];

#ifdef MAGNETIC
/*! table for the cosmological kick factor for induction equation */
double MagKickTable[DRIFT_TABLE_LENGTH];
#endif


void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP,	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */

#ifdef LT_STELLAREVOLUTION
struct met_particle_data *MetP, *DomainMetBuf;
#endif

#if defined(BLACK_HOLES) && defined(DETACH_BLACK_HOLES)
struct bh_particle_data *BHP;
#endif

peanokey *DomainKeyBuf;

/* global state of system
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

struct data_index *DataIndexTable;	/*!< the particles to be exported are grouped
					   by task-number. This table allows the
					   results to be disentangled again and to be
					   assigned to the correct particle */

struct data_nodelist *DataNodeList;

struct gravdata_in *GravDataIn,	/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct potdata_out *PotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct info_block *InfoBlock;

/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */





/*
 * Variables for Tree
 * ------------------
 */

int Nexport, Nimport;
int BufferFullFlag;
int NextParticle;
int NextJ;
int TimerFlag;

struct NODE *Nodes_base,	/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
				   gives the first allocated node */


struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/*!< maximum allowed number of internal nodes */
int Numnodestree;		/*!< number of (internal) nodes in each tree */


int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */
int *Father;			/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif

#if defined (UM_CHEMISTRY) || defined (UM_METAL_COOLING)
/* --==[ link with LT_ stuffs]==-- */
float *um_ZsPoint, um_FillEl_mu, um_mass;

/* char *PT_Symbols; */
/* double *PT_Masses; */
/* int NPT; */

/* double **II_AvgFillNDens, **IIShLv_AvgFillNDens, **Ia_AvgFillNDens, **AGB_AvgFillNDens; */
#endif

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
  k10a[N_T], k11a[N_T];
double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
  k20a[N_T], k21a[N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif
#endif

#if defined (UM_CHEMISTRY)
double kHeHII1a[N_T], kHeHII2a[N_T], kHeHII3a[N_T];
double kHD1a[N_T], kHD2a[N_T], kHD3a[N_T], kHD4a[N_T], kHD5a[N_T], kHD6a[N_T], kHD7a[N_T];
#endif

#ifdef LT_STELLAREVOLUTION
#include "lt_sfr/lt.c"
#endif



#ifdef SCFPOTENTIAL
long scf_seed;
MyDouble *Anltilde, *coeflm, *twoalpha, *c1, *c2, *c3;
MyDouble *cosmphi, *sinmphi;
MyDouble *ultrasp, *ultraspt, *ultrasp1;
MyDouble *dblfact, *plm, *dplm;
MyDouble *sinsum, *cossum;
MyDouble *sinsum_all, *cossum_all;
#ifdef SCF_SCALEFAC
float *scalefac_nlm;
#endif
#endif

#if defined(NUM_THREADS)
int maxThreads = NUM_THREADS;
#else
int maxThreads = 1;
#endif

#ifdef KD_ALTERNATIVE_GROUP_SORT
int MaxNgroups;
#endif

