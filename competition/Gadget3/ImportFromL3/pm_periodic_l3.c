#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../allvars.h"
#include "../proto.h"

/*! \file pm_periodic.c
 *  \brief routines for periodic PM-force computation
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef NUM_THREADS
#include <pthread.h>
#endif

#ifdef PM_FROM_L3
#ifdef PMGRID
#ifdef PERIODIC

#define CONCAT(prefix, name) prefix ## name
#ifdef DOUBLEPRECISION_FFTW
#define FFTW(x) CONCAT(fftw_, x)
#else
#define FFTW(x) CONCAT(fftwf_, x)
#endif

#ifdef FFTW3
#include <fftw3.h>
#include <fftw3-mpi.h>
#ifdef DOUBLEPRECISION_FFTW
typedef double fftw_real;
#else
typedef float fftw_real;
#endif
#else
#ifdef NOTYPEPREFIX_FFTW
#include      <rfftw_mpi.h>
#ifdef NUM_THREADS
#include      <rfftw_threads.h>
#endif
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#ifdef NUM_THREADS
#include     <drfftw_threads.h>
#endif
#else
#include     <srfftw_mpi.h>
#ifdef NUM_THREADS
#include     <srfftw_threads.h>
#endif
#endif
#endif
#endif

#define  PMGRID2 (2*(PMGRID/2 + 1))

#if (PMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#ifdef FFTW3
static FFTW(plan) fft_forward_plan, fft_inverse_plan;
     ptrdiff_t slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;
     ptrdiff_t fftsize, maxfftsize;
     static FFTW(complex) * fft_of_rhogrid;
#else
static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;
static int fftsize, maxfftsize;
static fftw_complex *fft_of_rhogrid;
#endif

static fftw_real *rhogrid, *forcegrid, *workspace;
void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch);
void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch);
static int slab_to_task[PMGRID];
static int *slabs_per_task;
static int *first_slab_of_task;

static MyFloat to_slab_fac;

struct partbuf
{
  MyFloat Pos[3];
  MyFloat Mass;
};

struct data_for_thread
{
  int startindex, lastindex, id;
  union
  {
    int *send_count_thread;
    struct partbuf *particle_buffer;
  } u;
};

void *pmforce_binning_loop(void *p);
void pmforce_binning(int nimport, struct partbuf *partin);


#ifndef FFTW3
void rfftwnd_mpi_threads(rfftwnd_mpi_plan p,
			 int n_fields, fftw_real * local_data, fftw_real * work,
			 fftwnd_mpi_output_order output_order);
#endif


void create_plans(void)
{
#ifndef FFTW3
  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
                                             FFTW_REAL_TO_COMPLEX,
                                             FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
                                             FFTW_COMPLEX_TO_REAL,
                                             FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);
#endif
}

/*! This routines generates the FFTW-plans to carry out the parallel FFTs
 *  later on. Some auxiliary variables are also initialized.
 */
void pm_init_periodic(void)
{
  int i;
  int slab_to_task_local[PMGRID];

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0] = RCUT * All.Asmth[0];


#ifdef NUM_THREADS
  if(PMGRID % NUM_THREADS)
    terminate("PMGRID must be a multiple of NUM_THREADS");
#endif

  /* Set up the FFTW plan files. */

#ifdef FFTW3

  FFTW(mpi_init) ();
#ifdef FFTW3_THREADS
  FFTW(init_threads) ();
#endif
  fftsize = FFTW(mpi_local_size_3d_transposed) (PMGRID, PMGRID, PMGRID, MPI_COMM_WORLD,
						&nslab_x, &slabstart_x, &nslab_y, &slabstart_y);


  fftsize *= 2;
//  fft_forward_plan = FFTW(mpi_plan_dft_r2c_3d)( PMGRID, PMGRID, PMGRID, rhogrid, fft_of_rhogrid, MPI_COMM_WORLD, FFTW_ESTIMATE);
//  fft_inverse_plan = FFTW(mpi_plan_dft_c2r_3d)( PMGRID, PMGRID, PMGRID, fft_of_rhogrid, rhogrid, MPI_COMM_WORLD, FFTW_ESTIMATE);

#else

  //  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
  //					     FFTW_REAL_TO_COMPLEX,
  //					     FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  //  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
  //					     FFTW_COMPLEX_TO_REAL,
  //					     FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  /* Workspace out the ranges on each processor. */
  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);
#endif

  if(ThisTask == 0)
    printf("fftsize=%d nslab_y=%d slabstart_y=%d, nslab_y=%d slabstart_y=%d \n", fftsize, nslab_x,
	   slabstart_x, nslab_y, slabstart_y);

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  slabs_per_task = (int *) mymalloc("slabs_per_task", NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc("first_slab_of_task", NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = (MyFloat) (PMGRID / All.BoxSize);
#ifdef REPLICATE
  int GlassTileFac = REPLICATE;
  to_slab_fac /= GlassTileFac;
  All.Asmth[0] *= GlassTileFac;
  All.Rcut[0] *= GlassTileFac;
#endif

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}


void pmforce_binning(int nimport, struct partbuf *partin)
{
  struct data_for_thread main_thread;

#ifdef NUM_THREADS
  pthread_t mythreads[NUM_THREADS - 1];
  pthread_attr_t attr;
  struct data_for_thread thread_data[NUM_THREADS - 1];
  int j;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for(j = 0; j < NUM_THREADS - 1; j++)
    {
      thread_data[j].u.particle_buffer = partin;
      thread_data[j].id = j + 1;
      thread_data[j].startindex = 0;
      thread_data[j].lastindex = nimport - 1;
      pthread_create(&mythreads[j], &attr, pmforce_binning_loop, &thread_data[j]);
    }

#endif

  main_thread.u.particle_buffer = partin;
  main_thread.startindex = 0;
  main_thread.lastindex = nimport - 1;
  main_thread.id = 0;

  pmforce_binning_loop(&main_thread);

#ifdef NUM_THREADS
  for(j = 0; j < NUM_THREADS - 1; j++)
    pthread_join(mythreads[j], NULL);

  pthread_attr_destroy(&attr);
#endif
}



void *pmforce_binning_loop(void *p)
{
  int i, thread_id, id_y, id_yy;
  int flag_slab_x, flag_slab_xx;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  double dx, dy, dz;
  struct data_for_thread *thread_data = (struct data_for_thread *) p;
  struct partbuf *partin;

  partin = thread_data->u.particle_buffer;
  thread_id = thread_data->id;

  for(i = thread_data->startindex; i <= thread_data->lastindex; i++)
    {
      slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
      dy = to_slab_fac * partin[i].Pos[1] - slab_y;

      if(slab_y >= PMGRID)
	slab_y -= PMGRID;

      slab_yy = slab_y + 1;
      if(slab_yy >= PMGRID)
	slab_yy -= PMGRID;

#ifdef NUM_THREADS
      id_y = slab_y / (PMGRID / NUM_THREADS);
      id_yy = slab_yy / (PMGRID / NUM_THREADS);
#else
      id_y = id_yy = 0;
#endif

      if(id_y == thread_id || id_yy == thread_id)
	{
	  slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
	  slab_z = (int) (to_slab_fac * partin[i].Pos[2]);

	  dx = to_slab_fac * partin[i].Pos[0] - slab_x;
	  dz = to_slab_fac * partin[i].Pos[2] - slab_z;

	  if(slab_x >= PMGRID)
	    slab_x -= PMGRID;
	  if(slab_z >= PMGRID)
	    slab_z -= PMGRID;

	  slab_xx = slab_x + 1;
	  slab_zz = slab_z + 1;

	  if(slab_xx >= PMGRID)
	    slab_xx -= PMGRID;
	  if(slab_zz >= PMGRID)
	    slab_zz -= PMGRID;

	  if(slab_to_task[slab_x] == ThisTask)
	    {
	      slab_x -= first_slab_of_task[ThisTask];
	      flag_slab_x = 1;
	    }
	  else
	    flag_slab_x = 0;

	  if(slab_to_task[slab_xx] == ThisTask)
	    {
	      slab_xx -= first_slab_of_task[ThisTask];
	      flag_slab_xx = 1;
    }
	  else
	    flag_slab_xx = 0;

	  if(flag_slab_x)
	    {
	      if(id_y == thread_id)
		{
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_y) + slab_z] +=
		    (fftw_real) (partin[i].Mass * ((1.0 - dx) * (1.0 - dy) * (1.0 - dz)));
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_y) + slab_zz] +=
		    (fftw_real) (partin[i].Mass * ((1.0 - dx) * (1.0 - dy) * (dz)));
		}

	      if(id_yy == thread_id)
		{
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_yy) + slab_z] +=
		    (fftw_real) (partin[i].Mass * ((1.0 - dx) * (dy) * (1.0 - dz)));
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_yy) + slab_zz] +=
		    (fftw_real) (partin[i].Mass * ((1.0 - dx) * (dy) * (dz)));
		}
	    }


	  if(flag_slab_xx)
	    {
	      if(id_y == thread_id)
		{
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_y) + slab_z] +=
		    (fftw_real) (partin[i].Mass * ((dx) * (1.0 - dy) * (1.0 - dz)));
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_y) + slab_zz] +=
		    (fftw_real) (partin[i].Mass * ((dx) * (1.0 - dy) * (dz)));
		}

	      if(id_yy == thread_id)
		{
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_z] +=
		    (fftw_real) (partin[i].Mass * ((dx) * (dy) * (1.0 - dz)));
		  rhogrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz] +=
		    (fftw_real) (partin[i].Mass * ((dx) * (dy) * (dz)));
		}
	    }
	}
    }

  return NULL;
}


/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The potential is finite differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved. Note that the particle distribution is not in the slab
 *  decomposition that is used for the FFT. Instead, overlapping patches
 *  between local domains and FFT slabs are communicated as needed.
 *
 *  For mode=0, normal force calculation, mode=1, only PS calculation.
 */

void pmforce_periodic(int mode, int *typelist)
{
  double smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, acc_dim;
  int i, j, sendTask, recvTask;
  int k2, kx, ky, kz, x, y, z, yl, zl, yr, zr, yll, zll, yrr, zrr, ip, dim;
#ifdef AVOID_BGQ_COMPILER_BUG_IN_PM
  int iif, iir1, iir2, iir3, iir4;
#endif
  int task0, task1, nimport, nexport;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int ngrp, xx;
  int *send_count_bak;
  MyFloat *flistin, *flistout;
  struct partbuf *partin, *partout;
  /*
     double time0, time1;

     time0 = second();
   */


  if(ThisTask == 0)
    {
      printf("Starting periodic PM calculation.  (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  /*
  if(mode == 1)
    {
      foldonitself(typelist);
    }
  */

  if(mode == 0 || mode == 2)
    {
       asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
      asmth2 *= asmth2;
      /*
      fac = All.G * All.PartMass / (M_PI * All.BoxSize);
      */
      fac = All.G / (M_PI * All.BoxSize);        /* to get potential */
      fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */

    }
  else
    fac = asmth2 = 0;


  send_count_bak = (int *) mymalloc("send_count_bak", sizeof(int) * NTask);

  /* determine the slabs each particles accesses */
  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < NumPart; i++)
    {
      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;

      slab_xx = slab_x + 1;
      if(slab_xx >= PMGRID)
	slab_xx -= PMGRID;

      task0 = slab_to_task[slab_x];
      task1 = slab_to_task[slab_xx];

      Send_count[task0]++;
      if(task0 != task1)
	Send_count[task1]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  /* allocate import and export buffer */

  partin = (struct partbuf *) mymalloc("partin", nimport * sizeof(struct partbuf));
  partout = (struct partbuf *) mymalloc("partout", nexport * sizeof(struct partbuf));

  /* fill export buffer */
  for(j = 0; j < NTask; j++)
    {
      send_count_bak[j] = Send_count[j];
      Send_count[j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;

      slab_xx = slab_x + 1;
      if(slab_xx >= PMGRID)
	slab_xx -= PMGRID;

      task0 = slab_to_task[slab_x];
      task1 = slab_to_task[slab_xx];

      for(j = 0; j < 3; j++)
	partout[Send_offset[task0] + Send_count[task0]].Pos[j] = P[i].Pos[j];
      partout[Send_offset[task0] + Send_count[task0]].Mass = P[i].Mass;
      Send_count[task0]++;

      if(task0 != task1)
	{
	  for(j = 0; j < 3; j++)
	    partout[Send_offset[task1] + Send_count[task1]].Pos[j] = P[i].Pos[j];
	  partout[Send_offset[task1] + Send_count[task1]].Mass = P[i].Mass;
	  Send_count[task1]++;
	}
    }

  /* make a check */
  for(j = 0; j < NTask; j++)
    if(send_count_bak[j] != Send_count[j])
      {
	char buffer[500];

	sprintf(buffer, "Inconsitency: Task=%d: j=%d send_count_bak[j]=%d  Send_count[j]=%d", ThisTask, j,
		send_count_bak[j], Send_count[j]);
	terminate(buffer);
      }

  /* exchange particle data */

  memcpy(&partin[Recv_offset[ThisTask]],
	 &partout[Send_offset[ThisTask]], Recv_count[ThisTask] * sizeof(struct partbuf));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the particles */
	      MPI_Sendrecv(&partout[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct partbuf), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &partin[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct partbuf), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(partout);


  /* allocate density field */
  rhogrid = (fftw_real *) mymalloc_movable(&rhogrid, "rhogrid", maxfftsize * sizeof(fftw_real));

  memset(rhogrid, 0, fftsize * sizeof(fftw_real));	/* clear local FFT-mesh density field */

  /* the memset() above is a slightly faster alternative for the above loop:
     for(i = 0; i < fftsize; i++)
     rhogrid[i] = 0;
   */

  pmforce_binning(nimport, partin);	/* bin particle data onto mesh */

  if(mode == 3)
    {
#ifdef PICTURE_ON_OUTPUT
      picture();
#else
      terminate("how can this be?");
#endif
    }
  else
    {

      forcegrid = (fftw_real *) mymalloc("forcegrid", maxfftsize * sizeof(fftw_real));
      workspace = forcegrid;

      /* Do the FFT of the density field */

      report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC");

#ifdef FFTW3
#ifdef FFTW3_THREADS
      FFTW(plan_with_nthreads) (FFTW3_THREADS);
#endif
      fft_of_rhogrid = (FFTW(complex) *) & rhogrid[0];
      fft_forward_plan =
	FFTW(mpi_plan_dft_r2c_3d) (PMGRID, PMGRID, PMGRID, rhogrid, fft_of_rhogrid, MPI_COMM_WORLD,
				   FFTW_ESTIMATE | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_OUT);
      FFTW(execute) (fft_forward_plan);

      FFTW(destroy_plan) (fft_forward_plan);
#else
      fft_of_rhogrid = (fftw_complex *) & rhogrid[0];
#ifdef FFTW2_THREADS
      rfftwnd_mpi_threads(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else
      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#endif
#endif

      /*
      if(mode == 1)
	{
	  powerspec(1, typelist);
	}
      */

      if(mode == 0 || mode == 2)
	{
	  /* multiply with Green's function for the potential */

	  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
#pragma omp parallel for private(x,z,kx,ky,kz,k2,fx,fy,fz,ff,smth,ip)
	    for(x = 0; x < PMGRID; x++)
	      for(z = 0; z < PMGRID / 2 + 1; z++)
		{
		  if(x > PMGRID / 2)
		    kx = x - PMGRID;
		  else
		    kx = x;
		  if(y > PMGRID / 2)
		    ky = y - PMGRID;
		  else
		    ky = y;
		  if(z > PMGRID / 2)
		    kz = z - PMGRID;
		  else
		    kz = z;

		  k2 = kx * kx + ky * ky + kz * kz;

		  if(k2 > 0)
		    {
		      smth = -exp(-k2 * asmth2) / k2;

		      /* do deconvolution */

		      fx = fy = fz = 1;
		      if(kx != 0)
			{
			  fx = (M_PI * kx) / PMGRID;
			  fx = sin(fx) / fx;
			}
		      if(ky != 0)
			{
			  fy = (M_PI * ky) / PMGRID;
			  fy = sin(fy) / fy;
			}
		      if(kz != 0)
			{
			  fz = (M_PI * kz) / PMGRID;
			  fz = sin(fz) / fz;
			}
		      ff = 1 / (fx * fy * fz);
		      smth *= ff * ff * ff * ff;

		      /* end deconvolution */

		      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
#ifdef FFTW3
		      fft_of_rhogrid[ip][0] *= (fftw_real) smth;
		      fft_of_rhogrid[ip][1] *= (fftw_real) smth;
#else
		      fft_of_rhogrid[ip].re *= (fftw_real) smth;
		      fft_of_rhogrid[ip].im *= (fftw_real) smth;
#endif
		    }
		}

	  if(slabstart_y == 0)
#ifdef FFTW3
	    fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#else
	    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;
#endif

	  /* Do the inverse FFT to get the potential */
	  if(ThisTask == 0)
	    printf("Do the inverse FFT to get the potential\n");
#ifdef FFTW3
#ifdef FFTW3_THREADS
	  FFTW(plan_with_nthreads) (FFTW3_THREADS);
#endif
	  fft_inverse_plan =
	    FFTW(mpi_plan_dft_c2r_3d) (PMGRID, PMGRID, PMGRID, fft_of_rhogrid, rhogrid, MPI_COMM_WORLD,
				       FFTW_ESTIMATE | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_IN);
	  FFTW(execute) (fft_inverse_plan);
	  FFTW(destroy_plan) (fft_inverse_plan);
#else
#ifdef FFTW2_THREADS
	  rfftwnd_mpi_threads(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else
	  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#endif
#endif
	  /* Now rhogrid holds the potential */

	  /* get the force components by finite differencing the potential for each dimension,
	     and send back the results to the right CPUs */

	  for(dim = 2; dim >= 0; dim--)	/* Calculate each component of the force. */
	    {			/* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */
	      if(dim == 0)
		pm_periodic_transposeA(rhogrid, forcegrid);	/* compute the transpose of the potential field */

	      for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
#pragma omp parallel for private(z,  x, yrr, yll, yr, yl, zrr, zll, zr, zl)
		for(y = 0; y < PMGRID; y++)
		  for(z = 0; z < PMGRID; z++)
		    {
		      x = xx - slabstart_x;

		      yrr = yll = yr = yl = y;
		      zrr = zll = zr = zl = z;

		      switch (dim)
			{
			case 0:	/* note: for the x-direction, we difference the transposed direction (y) */
			case 1:
			  yr = y + 1;
			  yl = y - 1;
			  yrr = y + 2;
			  yll = y - 2;
			  if(yr >= PMGRID)
			    yr -= PMGRID;
			  if(yrr >= PMGRID)
			    yrr -= PMGRID;
			  if(yl < 0)
			    yl += PMGRID;
			  if(yll < 0)
			    yll += PMGRID;
			  break;
			case 2:
			  zr = z + 1;
			  zl = z - 1;
			  zrr = z + 2;
			  zll = z - 2;
			  if(zr >= PMGRID)
			    zr -= PMGRID;
			  if(zrr >= PMGRID)
			    zrr -= PMGRID;
			  if(zl < 0)
			    zl += PMGRID;
			  if(zll < 0)
			    zll += PMGRID;
			  break;
			}

		      if(dim == 0)
			{
			  forcegrid[PMGRID * (x + y * nslab_x) + z]
			    = (fftw_real) (fac * ((4.0 / 3) *
						  (rhogrid[PMGRID * (x + yl * nslab_x) + zl] -
						   rhogrid[PMGRID * (x + yr * nslab_x) + zr]) -
						  (1.0 / 6) * (rhogrid[PMGRID * (x + yll * nslab_x) + zll] -
							       rhogrid[PMGRID * (x + yrr * nslab_x) + zrr])));
			}
		      else
#ifdef AVOID_BGQ_COMPILER_BUG_IN_PM
			{
			  iif = PMGRID2 * (PMGRID * x + y) + z;
			  iir1 = PMGRID2 * (PMGRID * x + yl) + zl;
			  iir2 = PMGRID2 * (PMGRID * x + yr) + zr;
			  iir3 = PMGRID2 * (PMGRID * x + yll) + zll;
			  iir4 = PMGRID2 * (PMGRID * x + yrr) + zrr;
			  if(iif < 0 || iif >= maxfftsize)
			    printf("Task %d: Something is wrong: %d %d (x=%d,y=%d,z=%d\n",ThisTask,iif,maxfftsize,x,y,z);
			  if(iir1 < 0 || iir1 >= maxfftsize)
			    printf("Task %d: Something is wrong: %d %d (x=%d,y=%d,z=%d\n",ThisTask,iir1,maxfftsize,x,yl,zl);
			  if(iir2 < 0 || iir2 >= maxfftsize)
			    printf("Task %d: Something is wrong: %d %d (x=%d,y=%d,z=%d\n",ThisTask,iir2,maxfftsize,x,yr,zr);
			  if(iir3 < 0 || iir3 >= maxfftsize)
			    printf("Task %d: Something is wrong: %d %d (x=%d,y=%d,z=%d\n",ThisTask,iir3,maxfftsize,x,yll,zll);
			  if(iir4 < 0 || iir4 >= maxfftsize)
			    printf("Task %d: Something is wrong: %d %d (x=%d,y=%d,z=%d\n",ThisTask,iir4,maxfftsize,x,yrr,zrr);
                          forcegrid[iif]
                            = (fftw_real) (fac * ((4.0 / 3) *
                                                  (rhogrid[iir1] - rhogrid[iir2]) -
                                                  (1.0 / 6) * (rhogrid[iir3] - rhogrid[iir4])));
			}
#else
			forcegrid[PMGRID2 * (PMGRID * x + y) + z]
			  = (fftw_real) (fac * ((4.0 / 3) *
						(rhogrid[PMGRID2 * (PMGRID * x + yl) + zl] -
						 rhogrid[PMGRID2 * (PMGRID * x + yr) + zr]) -
						(1.0 / 6) * (rhogrid[PMGRID2 * (PMGRID * x + yll) + zll] -
							     rhogrid[PMGRID2 * (PMGRID * x + yrr) + zrr])));
#endif
		    }

	      if(dim == 0)
		pm_periodic_transposeB(forcegrid, rhogrid);	/* compute the transpose of the potential field */

	      /* read out the force components */

	      flistin = (MyFloat *) mymalloc("flistin", nimport * sizeof(MyFloat));
	      flistout = (MyFloat *) mymalloc("flistout", nexport * sizeof(MyFloat));

	      report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC");

#pragma omp parallel for private(slab_x, slab_y, slab_z, slab_xx, slab_yy, slab_zz, dx, dy, dz)
	      for(i = 0; i < nimport; i++)
		{
		  flistin[i] = 0;

		  slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
		  slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
		  slab_z = (int) (to_slab_fac * partin[i].Pos[2]);

		  dx = to_slab_fac * partin[i].Pos[0] - slab_x;
		  dy = to_slab_fac * partin[i].Pos[1] - slab_y;
		  dz = to_slab_fac * partin[i].Pos[2] - slab_z;

		  if(slab_x >= PMGRID)
		    slab_x -= PMGRID;
		  if(slab_y >= PMGRID)
		    slab_y -= PMGRID;
		  if(slab_z >= PMGRID)
		    slab_z -= PMGRID;

		  slab_xx = slab_x + 1;
		  slab_yy = slab_y + 1;
		  slab_zz = slab_z + 1;

		  if(slab_xx >= PMGRID)
		    slab_xx -= PMGRID;
		  if(slab_yy >= PMGRID)
		    slab_yy -= PMGRID;
		  if(slab_zz >= PMGRID)
		    slab_zz -= PMGRID;

		  if(slab_to_task[slab_x] == ThisTask)
		    {
		      slab_x -= first_slab_of_task[ThisTask];

		      flistin[i] +=
			(MyFloat) (forcegrid
				   [((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_y) +
				    slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
				   forcegrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_y) +
					     slab_zz] * (1.0 - dx) * (1.0 - dy) * (dz) +
				   forcegrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_yy) +
					     slab_z] * (1.0 - dx) * (dy) * (1.0 - dz) +
				   forcegrid[((large_array_offset) PMGRID2) * (PMGRID * slab_x + slab_yy) +
					     slab_zz] * (1.0 - dx) * (dy) * (dz));
		    }

		  if(slab_to_task[slab_xx] == ThisTask)
		    {
		      slab_xx -= first_slab_of_task[ThisTask];

		      flistin[i] +=
			(MyFloat) (forcegrid
				   [((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_y) +
				    slab_z] * (dx) * (1.0 - dy) * (1.0 - dz) +
				   forcegrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_y) +
					     slab_zz] * (dx) * (1.0 - dy) * (dz) +
				   forcegrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) +
					     slab_z] * (dx) * (dy) * (1.0 - dz) +
				   forcegrid[((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) +
					     slab_zz] * (dx) * (dy) * (dz));
		    }
		}

	      /* exchange the force component data */

	      memcpy(&flistout[Send_offset[ThisTask]],
		     &flistin[Recv_offset[ThisTask]], Recv_count[ThisTask] * sizeof(MyFloat));

	      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
		{
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
			{
			  /* get the force components */
			  MPI_Sendrecv(&flistin[Recv_offset[recvTask]],
				       Recv_count[recvTask] * sizeof(MyFloat), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &flistout[Send_offset[recvTask]],
				       Send_count[recvTask] * sizeof(MyFloat), MPI_BYTE,
				       recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    }
		}


	      /* now assign them to the correct particles */

	      /* fill export buffer */
	      for(j = 0; j < NTask; j++)
		Send_count[j] = 0;

	      for(i = 0; i < NumPart; i++)
		{
		  slab_x = (int) (to_slab_fac * P[i].Pos[0]);
		  if(slab_x >= PMGRID)
		    slab_x = PMGRID - 1;

		  slab_xx = slab_x + 1;
		  if(slab_xx >= PMGRID)
		    slab_xx -= PMGRID;

		  task0 = slab_to_task[slab_x];
		  task1 = slab_to_task[slab_xx];

		  acc_dim = 0;

		  acc_dim += flistout[Send_offset[task0] + Send_count[task0]++];

		  if(task0 != task1)
		    acc_dim += flistout[Send_offset[task1] + Send_count[task1]++];
		  /*
		  P[i].Vel[dim] += (MyFloat) (fac2 * acc_dim);
		  */
                  P[i].GravPM[dim] += (MyFloat) (acc_dim);

#ifdef RELATIVE_OPENING_CRITERION
		  acc_dim /= All.G;
		  set_OldAcc(&P[i], get_OldAcc(&P[i]) + (MyFloat) (acc_dim * acc_dim));
#endif
		}

	      /* make a check */
	      for(j = 0; j < NTask; j++)
		if(send_count_bak[j] != Send_count[j])
		  {
		    char buffer[500];

		    sprintf(buffer, "Inconsitency: Task=%d: j=%d send_count_bak[j]=%d  Send_count[j]=%d",
			    ThisTask, j, send_count_bak[j], Send_count[j]);
		    terminate(buffer);
		  }

	      myfree(flistout);
	      myfree(flistin);
	    }
	}			/* end of if(mode==0 || mode==2) block */

      myfree(forcegrid);
    }
  myfree_movable(rhogrid);
  myfree(partin);

  myfree(send_count_bak);

#ifdef RELATIVE_OPENING_CRITERION
  if(mode == 0 || mode == 2)
    {
#pragma omp parallel for
      for(i = 0; i < NumPart; i++)
	set_OldAcc(&P[i], sqrt(get_OldAcc(&P[i])));
    }
#endif

  if(ThisTask == 0)
    {
      printf("done PM.\n");
      fflush(stdout);
    }
}




void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#pragma omp parallel for private(y,z,task)
  for(x = 0; x < nslab_x; x++)
    for(task = 0; task < NTask; task++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	{
	  for(z = 0; z < PMGRID; z++)
	    {
	      scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
				x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z] =
		field[PMGRID2 * (PMGRID * x + y) + z];
	    }
	}


#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif
}



void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }


  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif

#pragma omp parallel for private(y,z,task)
  for(x = 0; x < nslab_x; x++)
    for(task = 0; task < NTask; task++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	{
	  for(z = 0; z < PMGRID; z++)
	    {
	      field[PMGRID2 * (PMGRID * x + y) + z] =
		scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
				  x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z];
	    }
	}
}


#endif
#endif
#endif
