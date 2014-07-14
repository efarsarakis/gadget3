
#ifndef FFTW3			/* the following is only for FFTW2 */

#ifdef  FFTW2_THREADS

/* the routines in this file are rewrites of stuff in the FFTW-2.1.5 library,
   contained there in file rfftwnd_mpi.c
*/


#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>


#ifdef NOTYPEPREFIX_FFTW
#include      <rfftw_mpi.h>
#include      <rfftw_threads.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include <drfftw_mpi.h>
#include <drfftw_threads.h>
#else
#include <srfftw_mpi.h>
#include <srfftw_threads.h>
#endif
#endif

#include "../allvars.h"
#include "../proto.h"


void my_transpose_mpi(transpose_mpi_plan p, int el_size,
		      TRANSPOSE_EL_TYPE * local_data, TRANSPOSE_EL_TYPE * work);



struct data_for_thread
{
  fftwnd_plan p;
  fftw_real *loc_data;
  fftw_complex *loc_complex_data;
  int howmany, istride, idist;
  int startindex, lastindex, n_fields, nx;
};


/******************************************************************************************/

static void *aux_threads_rfftwnd_real_to_complex(void *p)
{
  struct data_for_thread *thread_data = (struct data_for_thread *) p;

  rfftwnd_real_to_complex(thread_data->p,
			  thread_data->howmany,
			  thread_data->loc_data, thread_data->istride, thread_data->idist, NULL, 0, 0);
  return NULL;
}



static void threads_rfftwnd_real_to_complex(fftwnd_plan p, int howmany, fftw_real * in, int istride,
					    int idist)
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
      thread_data[j].startindex = j * (howmany / NUM_THREADS);
      thread_data[j].lastindex = (j + 1) * (howmany / NUM_THREADS) - 1;

      thread_data[j].loc_data = in + idist * thread_data[j].startindex;
      thread_data[j].howmany = thread_data[j].lastindex - thread_data[j].startindex + 1;

      thread_data[j].istride = istride;
      thread_data[j].idist = idist;
      thread_data[j].p = p;

      pthread_create(&mythreads[j], &attr, aux_threads_rfftwnd_real_to_complex, &thread_data[j]);
    }

  main_thread.startindex = (NUM_THREADS - 1) * (howmany / NUM_THREADS);
#else
  main_thread.startindex = 0;
#endif

  main_thread.lastindex = howmany - 1;

  main_thread.loc_data = in + idist * main_thread.startindex;
  main_thread.howmany = main_thread.lastindex - main_thread.startindex + 1;
  main_thread.istride = istride;
  main_thread.idist = idist;
  main_thread.p = p;

  aux_threads_rfftwnd_real_to_complex(&main_thread);

#ifdef NUM_THREADS
  for(j = 0; j < NUM_THREADS - 1; j++)
    pthread_join(mythreads[j], NULL);

  pthread_attr_destroy(&attr);
#endif
}




/******************************************************************************************/
static void *aux_threads_rfftwnd_complex_to_real(void *p)
{
  struct data_for_thread *thread_data = (struct data_for_thread *) p;

  rfftwnd_complex_to_real(thread_data->p,
			  thread_data->howmany,
			  thread_data->loc_complex_data, thread_data->istride, thread_data->idist, NULL, 0,
			  0);
  return NULL;
}



static void threads_rfftwnd_complex_to_real(fftwnd_plan p, int howmany, fftw_complex * in, int istride,
					    int idist)
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
      thread_data[j].startindex = j * (howmany / NUM_THREADS);
      thread_data[j].lastindex = (j + 1) * (howmany / NUM_THREADS) - 1;

      thread_data[j].loc_complex_data = in + idist * thread_data[j].startindex;
      thread_data[j].howmany = thread_data[j].lastindex - thread_data[j].startindex + 1;

      thread_data[j].istride = istride;
      thread_data[j].idist = idist;
      thread_data[j].p = p;

      pthread_create(&mythreads[j], &attr, aux_threads_rfftwnd_complex_to_real, &thread_data[j]);
    }

  main_thread.startindex = (NUM_THREADS - 1) * (howmany / NUM_THREADS);
#else
  main_thread.startindex = 0;
#endif

  main_thread.lastindex = howmany - 1;

  main_thread.loc_complex_data = in + idist * main_thread.startindex;
  main_thread.howmany = main_thread.lastindex - main_thread.startindex + 1;
  main_thread.istride = istride;
  main_thread.idist = idist;
  main_thread.p = p;

  aux_threads_rfftwnd_complex_to_real(&main_thread);

#ifdef NUM_THREADS
  for(j = 0; j < NUM_THREADS - 1; j++)
    pthread_join(mythreads[j], NULL);

  pthread_attr_destroy(&attr);
#endif
}




/******************************************************************************************/

static void *aux_threads_fftw(void *p)
{
  struct data_for_thread *thread_data = (struct data_for_thread *) p;
  int fft_iter;

  for(fft_iter = thread_data->startindex; fft_iter <= thread_data->lastindex; fft_iter++)
    fftw(thread_data->p,
	 thread_data->n_fields,
	 ((fftw_complex *) thread_data->loc_data) + (thread_data->nx * thread_data->n_fields) * fft_iter,
	 thread_data->n_fields, 1, NULL, 0, 0);

  return NULL;
}



static void threads_fftw(fftwnd_plan p, int howmany, fftw_real * in, int n_fields, int nx)
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
      thread_data[j].startindex = j * (howmany / NUM_THREADS);
      thread_data[j].lastindex = (j + 1) * (howmany / NUM_THREADS) - 1;

      thread_data[j].loc_data = in;
      thread_data[j].howmany = thread_data[j].lastindex - thread_data[j].startindex + 1;

      thread_data[j].n_fields = n_fields;
      thread_data[j].nx = nx;
      thread_data[j].p = p;

      pthread_create(&mythreads[j], &attr, aux_threads_fftw, &thread_data[j]);
    }

  main_thread.startindex = (NUM_THREADS - 1) * (howmany / NUM_THREADS);
#else
  main_thread.startindex = 0;
#endif

  main_thread.lastindex = howmany - 1;

  main_thread.loc_data = in;
  main_thread.howmany = main_thread.lastindex - main_thread.startindex + 1;
  main_thread.n_fields = n_fields;
  main_thread.nx = nx;
  main_thread.p = p;

  aux_threads_fftw(&main_thread);

#ifdef NUM_THREADS
  for(j = 0; j < NUM_THREADS - 1; j++)
    pthread_join(mythreads[j], NULL);

  pthread_attr_destroy(&attr);
#endif

}



/******************** Computing the Transform *******************/


static void first_dim_aux(rfftwnd_mpi_plan p, int n_fields, fftw_real * local_data)
{
  int local_ny = p->p_transpose->local_ny;
  int nx = p->p_fft_x->n;
  fftw_complex *work_1d = p->work ? p->work : p->p_fft->work;

  if(work_1d)
    terminate("this patched threaded version is only going to work for in-place");

  n_fields *= p->p_fft->n_after[0];	/* dimensions after y 
					   no longer need be considered
					   separately from n_fields */
  if(n_fields > 1)
    {
      fftw_plan p_fft_x = p->p_fft_x;

      threads_fftw(p_fft_x, local_ny, local_data, n_fields, nx);

      /* the above replaced the following loop over fftw() */
      /*
         int fft_iter;
         for (fft_iter = 0; fft_iter < local_ny; ++fft_iter)
         fftw(p_fft_x,  n_fields,
         ((fftw_complex *) local_data) + (nx * n_fields) * fft_iter, 
         n_fields, 1, NULL, 1, 0);
       */
    }
  else
    {
      terminate("not implemented\n");
    }
}


static void other_dims_aux(rfftwnd_mpi_plan p, int n_fields, fftw_real * local_data)
{
  int local_nx = p->p_transpose->local_nx;
  int n_after_x = p->p_fft->n[0] * p->p_fft->n_after[0];

  if(n_fields > 1)
    {
      terminate("n_fields > 1 not implemented in this routine\n");
    }
  else
    {
      if(p->p_fft->dir == FFTW_REAL_TO_COMPLEX)
	{
	  /*
	     rfftwnd_real_to_complex(p->p_fft, local_nx, local_data, 1, 2*n_after_x,  NULL, 0, 0);
	   */

	  threads_rfftwnd_real_to_complex(p->p_fft, local_nx, local_data, 1, 2 * n_after_x);
	}
      else
	{
	  /*
	     rfftwnd_complex_to_real(p->p_fft, local_nx, (fftw_complex *) local_data, 1, n_after_x, NULL, 0, 0);
	   */

	  threads_rfftwnd_complex_to_real(p->p_fft, local_nx, (fftw_complex *) local_data, 1, n_after_x);
	}

    }
}



void rfftwnd_mpi_threads(rfftwnd_mpi_plan p,
			 int n_fields, fftw_real * local_data, fftw_real * work,
			 fftwnd_mpi_output_order output_order)
{
  int el_size = (sizeof(fftw_complex) / sizeof(TRANSPOSE_EL_TYPE)) * n_fields * p->p_fft->n_after[0];

  if(n_fields <= 0)
    return;


  if(p->p_fft->dir == FFTW_REAL_TO_COMPLEX)
    {
      /* First, transform dimensions after the first, which are
         local to this process: */
      other_dims_aux(p, n_fields, local_data);

      /* Second, transpose the first dimension with the second dimension
         to bring the x dimension local to this process: */

      my_transpose_mpi(p->p_transpose, el_size, (TRANSPOSE_EL_TYPE *) local_data, (TRANSPOSE_EL_TYPE *) work);


      /* Third, transform the x dimension, which is now 
         local and contiguous: */

      first_dim_aux(p, n_fields, local_data);

      /* transpose back, if desired: */
      if(output_order == FFTW_NORMAL_ORDER)
	my_transpose_mpi(p->p_transpose_inv, el_size,
			 (TRANSPOSE_EL_TYPE *) local_data, (TRANSPOSE_EL_TYPE *) work);
    }
  else
    {
      /* we have to do the steps in reverse order for c2r transform: */

      /* NOTE: we assume that the same output_order is used for both
         the forward and backward transforms: */

      /* First, if necessary, transpose to get x dimension local: */
      if(output_order == FFTW_NORMAL_ORDER)
	my_transpose_mpi(p->p_transpose, el_size, (TRANSPOSE_EL_TYPE *) local_data,
			 (TRANSPOSE_EL_TYPE *) work);

      /* Second, transform the x dimension, which is now 
         local and contiguous: */
      first_dim_aux(p, n_fields, local_data);

      /* Third, transpose the first dimension with the second dimension
         to bring the others dimensions local to this process: */
      my_transpose_mpi(p->p_transpose_inv, el_size,
		       (TRANSPOSE_EL_TYPE *) local_data, (TRANSPOSE_EL_TYPE *) work);

      /* last, transform dimensions after the first, which are
         local to this process: */
      other_dims_aux(p, n_fields, local_data);
    }
}

#endif
#endif
