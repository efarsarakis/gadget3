#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif


#ifndef FFTW3			/* the following is only for FFTW2 */
#ifdef  FFTW2_THREADS

#ifdef NOTYPEPREFIX_FFTW
#include      <rfftw_mpi.h>
#include      <rfftw_threads.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include <drfftw_mpi.h>
#ifdef NUM_THREADS
#include <drfftw_threads.h>
#endif
#else
#include <srfftw_mpi.h>
#ifdef NUM_THREADS
#include <srfftw_threads.h>
#endif
#endif
#endif


#include "../allvars.h"
#include "../proto.h"

/* the routines in this file are rewrites of stuff in the FFTW-2.1.5 library,
   contained there in file transpose_mpi.c
*/


static void exchange_elements(TRANSPOSE_EL_TYPE * buf1, TRANSPOSE_EL_TYPE * buf2, int n)
{
  int i;
  TRANSPOSE_EL_TYPE t0, t1, t2, t3;

  for(i = 0; i < (n & 3); ++i)
    {
      t0 = buf1[i];
      buf1[i] = buf2[i];
      buf2[i] = t0;
    }
  for(; i < n; i += 4)
    {
      t0 = buf1[i];
      t1 = buf1[i + 1];
      t2 = buf1[i + 2];
      t3 = buf1[i + 3];
      buf1[i] = buf2[i];
      buf1[i + 1] = buf2[i + 1];
      buf1[i + 2] = buf2[i + 2];
      buf1[i + 3] = buf2[i + 3];
      buf2[i] = t0;
      buf2[i + 1] = t1;
      buf2[i + 2] = t2;
      buf2[i + 3] = t3;
    }
}


static void do_permutation(TRANSPOSE_EL_TYPE * data,
			   int *perm_block_dest, int num_perm_blocks, int perm_block_size)
{
  int start_block;
  int cur_block;
  int new_block;

  /* Perform the permutation in the perm_block_dest array, following
     the cycles around and *changing* the perm_block_dest array
     to reflect the permutations that have already been performed.
     At the end of this routine, we change the perm_block_dest
     array back to its original state. (ugh) */

  for(start_block = 0; start_block < num_perm_blocks; start_block++)
    {
      cur_block = start_block;
      new_block = perm_block_dest[start_block];

      while(new_block > 0 && new_block < num_perm_blocks && new_block != start_block)
	{
	  exchange_elements(data + perm_block_size * start_block,
			    data + perm_block_size * new_block, perm_block_size);
	  perm_block_dest[cur_block] = -1 - new_block;
	  cur_block = new_block;
	  new_block = perm_block_dest[cur_block];
	}

      if(new_block == start_block)
	perm_block_dest[cur_block] = -1 - new_block;
    }

  /* reset the permutation array (ugh): */
  for(start_block = 0; start_block < num_perm_blocks; ++start_block)
    perm_block_dest[start_block] = -1 - perm_block_dest[start_block];
}


static void local_transpose_copy(TRANSPOSE_EL_TYPE * src,
				 TRANSPOSE_EL_TYPE * dest, int el_size, int nx, int ny)
{
  int x, y;

  if(el_size == 1)
    {
#pragma omp parallel for private (y)
      for(x = 0; x < nx; ++x)
	for(y = 0; y < ny; ++y)
	  dest[y * nx + x] = src[x * ny + y];
    }
  else if(el_size == 2)
    {
#pragma omp parallel for private (y)
      for(x = 0; x < nx; ++x)
	for(y = 0; y < ny; ++y)
	  {
	    dest[y * (2 * nx) + 2 * x] = src[x * (2 * ny) + 2 * y];
	    dest[y * (2 * nx) + 2 * x + 1] = src[x * (2 * ny) + 2 * y + 1];
	  }
    }
  else
    {
#pragma omp parallel for private (y)
      for(x = 0; x < nx; ++x)
	for(y = 0; y < ny; ++y)
	  memcpy(&dest[y * (el_size * nx) + (el_size * x)],
		 &src[x * (el_size * ny) + (el_size * y)], el_size * sizeof(TRANSPOSE_EL_TYPE));
    }

}



/* Out-of-place version of transpose_mpi (or rather, in place using
   a scratch array): */

static void transpose_mpi_out_of_place(transpose_mpi_plan p, int el_size,
				       TRANSPOSE_EL_TYPE * local_data, TRANSPOSE_EL_TYPE * work)
{
  local_transpose_copy(local_data, work, el_size, p->local_nx, p->ny);

  if(p->all_blocks_equal)
    {
      if(p->send_block_size != p->recv_block_size)
	terminate("unequal sizes");

      /* we replace the original MPI_Alltoall with a hypercube */
      /*
         MPI_Alltoall(work, p->send_block_size * el_size, p->el_type,
         local_data, p->recv_block_size * el_size, p->el_type, p->comm);
       */

      int ngrp, recvTask;
      MPI_Status status;
      MPI_Aint ext;
      MPI_Type_extent(p->el_type, &ext);

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      MPI_Sendrecv(work + recvTask * p->send_block_size * el_size,
			   p->send_block_size * el_size, p->el_type,
			   recvTask, TAG_GRAV_A,
			   local_data + recvTask * p->recv_block_size * el_size,
			   p->recv_block_size * el_size, p->el_type,
			   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
	    }
	}

      memcpy(local_data + ThisTask * p->recv_block_size * el_size,
	     work + ThisTask * p->send_block_size * el_size, p->recv_block_size * el_size * ext);
    }
  else
    {
      terminate("this routine only works if the FFT mesh size can be divided by the CPU number");
    }

  do_permutation(local_data, p->perm_block_dest, p->num_perm_blocks, p->perm_block_size * el_size);
}


void my_transpose_mpi(transpose_mpi_plan p, int el_size,
		      TRANSPOSE_EL_TYPE * local_data, TRANSPOSE_EL_TYPE * work)
{
  /* if local_data and work are both NULL, we have no way of knowing
     whether we should use in-place or out-of-place transpose routine;
     if we guess wrong, MPI_Alltoall will block.  We prevent this
     by making sure that transpose_mpi_get_local_storage_size returns
     at least 1. */

  if(!local_data && !work)
    terminate("need work space for this routine");

  if(work)
    {
      transpose_mpi_out_of_place(p, el_size, local_data, work);
    }
  else if(p->local_nx > 0 || p->local_ny > 0)
    {
      terminate("need work space for this routine");
    }
}

#endif
#endif
