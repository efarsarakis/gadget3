#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef SORT_FROM_L3

#define TRANSFER_SIZE_LIMIT  1000000000

#define TAG_TRANSFER  100

/* #define CHECK_LOCAL_RANK
 */
static void msort_serial_with_tmp(sort_type * base, size_t n, size_t s,
                                  int (*compar) (const void *, const void *), sort_type * t, int level,
                                  int maxlevel);
void parallel_sort_test_order(char *base, size_t nmemb, size_t size,
                              int (*compar) (const void *, const void *));


static void get_local_rank(char *element,
			   size_t tie_braking_rank,
			   char *base,
			   size_t nmemb, size_t size, size_t noffs_thistask,
			   long long left, long long right,
			   size_t * loc, int (*compar) (const void *, const void *));

#ifdef CHECK_LOCAL_RANK
static void check_local_rank(char *element,
			     size_t tie_braking_rank,
			     char *base,
			     size_t nmemb, size_t size,
			     size_t noffs_thistask,
			     long long left, long long right,
			     size_t loc, int (*compar) (const void *, const void *));
#endif

void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *))
{
  int i, j, max_task, ranks_not_found, Local_ThisTask, Local_NTask, Local_PTask, Color;
  MPI_Comm MPI_CommLocal;

  /* we create a communicator that contains just those tasks with nmemb > 0. This makes 
   *  it easier to deal with CPUs that do not hold any data.
   */

  if(nmemb)
  {
    Color = 1;
    serial_sort((char*)base, nmemb, size, compar);
  }
  else
    Color = 0;

  MPI_Comm_split(MPI_COMM_WORLD, Color, 0, &MPI_CommLocal);
  MPI_Comm_rank(MPI_CommLocal, &Local_ThisTask);
  MPI_Comm_size(MPI_CommLocal, &Local_NTask);

  if(Local_NTask > 1 && Color == 1)
    {
      for(Local_PTask = 0; Local_NTask > (1 << Local_PTask); Local_PTask++);

      size_t *nlist = (size_t *) mymalloc("nlist", Local_NTask * sizeof(size_t));
      size_t *noffs = (size_t *) mymalloc("noffs", Local_NTask * sizeof(size_t));

      MPI_Allgather(&nmemb, sizeof(size_t), MPI_BYTE, nlist, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      for(i = 1, noffs[0] = 0; i < Local_NTask; i++)
	noffs[i] = noffs[i - 1] + nlist[i - 1];

      char *element_guess = (char *) mymalloc("element_guess", (Local_NTask - 1) * size);
      size_t *element_tie_braking_rank =
	(size_t *) mymalloc("element_tie_braking_rank", (Local_NTask - 1) * sizeof(size_t));
      size_t *desired_glob_rank = (size_t *) mymalloc("desired_glob_rank", (Local_NTask - 1) * sizeof(size_t));
      size_t *current_glob_rank = (size_t *) mymalloc("current_glob_rank", (Local_NTask - 1) * sizeof(size_t));
      size_t *current_loc_rank = (size_t *) mymalloc("current_loc_rank", (Local_NTask - 1) * sizeof(size_t));
      long long *range_left = (long long *) mymalloc("range_left", (Local_NTask - 1) * sizeof(long long));
      long long *range_right = (long long *) mymalloc("range_right", (Local_NTask - 1) * sizeof(long long));
      int *max_loc = (int *) mymalloc("max_loc", (Local_NTask - 1) * sizeof(int));
      size_t max_nmemb;

      /* find the largest nmemb value and initialize the guesses with the median element found there */
      for(i = 0, max_task = max_nmemb = 0; i < Local_NTask; i++)
	if(max_nmemb < nlist[i])
	  {
	    max_nmemb = nlist[i];
	    max_task = i;
	  }

      if(Local_ThisTask == max_task)
	memcpy(element_guess, (char *) base + size * (max_nmemb / 2), size);

      MPI_Bcast(element_guess, size, MPI_BYTE, max_task, MPI_CommLocal);

      for(i = 1; i < Local_NTask - 1; i++)
	memcpy(element_guess + i * size, element_guess, size);

      for(i = 0; i < Local_NTask - 1; i++)
	{
	  desired_glob_rank[i] = noffs[i + 1];
	  element_tie_braking_rank[i] = (max_nmemb / 2) + noffs[max_task];
	  max_loc[i] = max_task;

	  range_left[i] = 0;	/* first element that it can be */
	  range_right[i] = nmemb;	/* first element that it can not be */

	  current_glob_rank[i] = 0;
	}

      int iter = 0;

      do
	{
	  for(i = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])
		{
		  get_local_rank(element_guess + i * size, element_tie_braking_rank[i],
				 (char *) base, nmemb, size, noffs[Local_ThisTask],
				 range_left[i], range_right[i], &current_loc_rank[i], compar);

#ifdef CHECK_LOCAL_RANK
		  check_local_rank(element_guess + i * size, element_tie_braking_rank[i],
				   (char *) base, nmemb, size, noffs[Local_ThisTask],
				   range_left[i], range_right[i], current_loc_rank[i], compar);
#endif
		}
	    }

	  /* now compute the global ranks by summing the local ranks */
	  size_t *list = (size_t *) mymalloc("list", Local_NTask * sizeof(size_t));
	  for(i = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])
		{
		  MPI_Allgather(&current_loc_rank[i], sizeof(size_t), MPI_BYTE, list, sizeof(size_t),
				MPI_BYTE, MPI_CommLocal);

		  for(j = 0, current_glob_rank[i] = 0; j < Local_NTask; j++)
		    current_glob_rank[i] += list[j];
		}
	    }
	  myfree(list);

	  for(i = 0, ranks_not_found = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])	/* here we're not yet done */
		{
		  ranks_not_found++;

		  if(current_glob_rank[i] < desired_glob_rank[i])
		    {
		      range_left[i] = current_loc_rank[i];

		      if(Local_ThisTask == max_loc[i])	// && range_right[i] == range_left[i]+1)
			{
			  range_left[i]++;
			}
		    }

		  if(current_glob_rank[i] > desired_glob_rank[i])
		    range_right[i] = current_loc_rank[i];
		}
	    }

	  /* now we need to determine new element guesses */

	  long long *range_len_list = (long long *) mymalloc("range_len_list", Local_NTask * sizeof(long long));

	  for(i = 0; i < Local_NTask - 1; i++)
	    {
	      if(current_glob_rank[i] != desired_glob_rank[i])	/* here we're not yet done */
		{
		  long long max_range_len, range_len = range_right[i] - range_left[i] + 1;
		  MPI_Allgather(&range_len, sizeof(long long), MPI_BYTE, range_len_list, sizeof(long long),
				MPI_BYTE, MPI_CommLocal);

		  /* find the largest one, and the cpu that holds it */
		  for(j = 0, max_loc[i] = max_range_len = 0; j < Local_NTask; j++)
		    if(max_range_len < range_len_list[j])
		      {
			max_range_len = range_len_list[j];
			max_loc[i] = j;
		      }

		  if(Local_ThisTask == max_loc[i])
		    {
		      long long mid = (range_left[i] + range_right[i]) / 2;

		      memcpy(element_guess + i * size, (char *) base + mid * size, size);
		      element_tie_braking_rank[i] = mid + noffs[Local_ThisTask];
		    }

		  MPI_Bcast(element_guess + i * size, size, MPI_BYTE, max_loc[i], MPI_CommLocal);
		  MPI_Bcast(&element_tie_braking_rank[i], sizeof(size_t), MPI_BYTE, max_loc[i],
			    MPI_CommLocal);
		}
	    }

	  myfree(range_len_list);

	  iter++;

	  if(iter > 400 + NTask && Local_ThisTask == 0)
	    {
	      printf("iter=%d: ranks_not_found=%d  Local_NTask=%d\n", iter, ranks_not_found, Local_NTask);
	      fflush(stdout);
	      if(iter > 500 + NTask)
		terminate("can't find the split points. That's odd");
	    }
	}
      while(ranks_not_found);


      /* At this point we have found all the elements corresponding to the desired split points */
      /* we can now go ahead and determine how many elements of the local CPU have to go to each other CPU */

      size_t *Send_count = (size_t *) mymalloc("Send_count", Local_NTask * sizeof(size_t));
      size_t *Recv_count = (size_t *) mymalloc("Recv_count", Local_NTask * sizeof(size_t));
      size_t *Send_offset = (size_t *) mymalloc("Send_offset", Local_NTask * sizeof(size_t));
      size_t *Recv_offset = (size_t *) mymalloc("Recv_offset", Local_NTask * sizeof(size_t));

      for(i = 0; i < Local_NTask; i++)
	Send_count[i] = 0;

      int target = 0;

      for(i = 0; i < nmemb; i++)
	{
	  while(target < Local_NTask - 1)
	    {
	      int cmp = compar((char *) base + i * size, element_guess + target * size);
	      if(cmp == 0)
		{
		  if(i + noffs[Local_ThisTask] < element_tie_braking_rank[target])
		    cmp = -1;
		  else if(i + noffs[Local_ThisTask] > element_tie_braking_rank[target])
		    cmp = +1;
		}
	      if(cmp >= 0)
		target++;
	      else
		break;
	    }
	  Send_count[target]++;
	}

      MPI_Alltoall(Send_count, sizeof(size_t), MPI_BYTE, Recv_count, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      size_t nimport;

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < Local_NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      if(nimport != nmemb)
	terminate("nimport != nmemb");

      char *basetmp = (char *) mymalloc("basetmp", nmemb * size);

      int ngrp, recvTask;

      /* exchange the data */
      for(ngrp = 0; ngrp < (1 << Local_PTask); ngrp++)
	{
	  recvTask = Local_ThisTask ^ ngrp;

	  if(recvTask < Local_NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
#ifndef MPISENDRECV_SIZELIMIT
		  if(Send_count[recvTask] > TRANSFER_SIZE_LIMIT || Recv_count[recvTask] > TRANSFER_SIZE_LIMIT)
		    terminate("we are above TRANSFER_SIZE_LIMIT");
#endif
		  MPI_Sendrecv((char *) base + Send_offset[recvTask] * size,
			       Send_count[recvTask] * size, MPI_BYTE,
			       recvTask, TAG_TRANSFER,
			       basetmp + Recv_offset[recvTask] * size,
			       Recv_count[recvTask] * size, MPI_BYTE,
			       recvTask, TAG_TRANSFER, MPI_CommLocal, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* instead of doing another serial sort at the end, we could also exploit that the different incoming pieces
         are already individually sorted. This can be used in a variant of merge sort to do the final step a bit quicker */

      memcpy(base, basetmp, nmemb * size);
      myfree(basetmp);
      serial_sort((char *)base, nmemb, size, compar);


      myfree(Recv_offset);
      myfree(Send_offset);
      myfree(Recv_count);
      myfree(Send_count);

      myfree(max_loc);
      myfree(range_right);
      myfree(range_left);
      myfree(current_loc_rank);
      myfree(current_glob_rank);
      myfree(desired_glob_rank);
      myfree(element_tie_braking_rank);
      myfree(element_guess);
      myfree(noffs);
      myfree(nlist);
    }

  MPI_Comm_free(&MPI_CommLocal);

}



static void get_local_rank(char *element,	/* element of which we want the rank */
			   size_t tie_braking_rank,	/* the inital global rank of this element (needed for braking ties) */
			   char *base,	/* base address of local data */
			   size_t nmemb, size_t size,	/* number and size of local data */
			   size_t noffs_thistask,	/* cumulative length of data on lower tasks */
			   long long left, long long right,	/* range of elements on local task that may hold the element */
			   size_t * loc,	/* output: local rank of the element */
			   int (*compar) (const void *, const void *))	/* user-specified  comparison function */
{
  if(right < left)
    terminate("right < left");

  if(left == 0 && right == nmemb + 1)
    {
      if(compar(base + (nmemb - 1) * size, element) < 0)
	{
	  *loc = nmemb;
	  return;
	}
      else if(compar(base, element) > 0)
	{
	  *loc = 0;
	  return;
	}
    }


  if(right == left)		/* looks like we already converged to the proper rank */
    {
      *loc = left;
    }
  else
    {
      if(compar(base + (right - 1) * size, element) < 0)	/* the last element is smaller, hence all elements are on the left */
	*loc = (right - 1) + 1;
      else if(compar(base + left * size, element) > 0)	/* the first element is already larger, hence no element is on the left */
	*loc = left;
      else
	{
	  while(right > left)
	    {
	      long long mid = ((right - 1) + left) / 2;

	      int cmp = compar(base + mid * size, element);
	      if(cmp == 0)
		{
		  if(mid + noffs_thistask < tie_braking_rank)
		    cmp = -1;
		  else if(mid + noffs_thistask > tie_braking_rank)
		    cmp = +1;
		}

	      if(cmp == 0)	/* element has exactly been found */
		{
		  *loc = mid;
		  break;
		}

	      if((right - 1) == left)	/* elements is not on this CPU */
		{
		  if(cmp < 0)
		    *loc = mid + 1;
		  else
		    *loc = mid;
		  break;
		}

	      if(cmp < 0)
		{
		  left = mid + 1;
		}
	      else
		{
		  if((right - 1) == left + 1)
		    {
		      if(mid != left)
			{
			  printf("-->left=%lld  right=%lld\n", left, right);
			  terminate("can't be");
			}
		      *loc = left;
		      break;
		    }

		  right = mid;
		}
	    }
	}
    }
}


#ifdef CHECK_LOCAL_RANK
static void check_local_rank(char *element,	/* element of which we want the rank */
			     size_t tie_braking_rank,	/* the inital global rank of this element (needed for braking ties) */
			     char *base,	/* base address of local data */
			     size_t nmemb, size_t size,	/* number and size of local data */
			     size_t noffs_thistask,	/* cumulative length of data on lower tasks */
			     long long left, long long right,	/* range of elements on local task that may hold the element */
			     size_t loc, int (*compar) (const void *, const void *))	/* user-specified  comparison function */
{

  int i;
  long long count = 0;

  for(i = 0; i < nmemb; i++)
    {
      int cmp = compar(base + i * size, element);

      if(cmp == 0)
	{
	  if(noffs_thistask + i < tie_braking_rank)
	    cmp = -1;
	}

      if(cmp < 0)
	count++;
    }

  if(count != loc)
    {
      printf("Task=%d: loc=%lld count=%lld  left=%lld right=%lld  nmemb=%lld\n",
	     ThisTask, (long long) loc, count, left, right, (long long) nmemb);
      terminate("inconsisteny");
    }
}
#endif



void parallel_sort_test_order(char *base, size_t nmemb, size_t size,
			      int (*compar) (const void *, const void *))
{
  int i, recv, send;
  size_t *nlist;

  nlist = (size_t *) mymalloc("nlist", NTask * sizeof(size_t));

  MPI_Allgather(&nmemb, sizeof(size_t), MPI_BYTE, nlist, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);


  for(i = 0, recv = -1; i < ThisTask && nmemb > 0; i++)
    if(nlist[i] > 0)
      recv = i;

  for(i = ThisTask + 1, send = -1; nmemb > 0 && i < NTask; i++)
    if(nlist[i] > 0)
      {
	send = i;
	break;
      }

  char *element = mymalloc("element", size);

  MPI_Request requests[2];
  int nreq = 0;

  if(send >= 0)
    MPI_Isend(base + (nmemb - 1) * size, size, MPI_BYTE, send, TAG_TRANSFER, MPI_COMM_WORLD,
	      &requests[nreq++]);

  if(recv >= 0)
    MPI_Irecv(element, size, MPI_BYTE, recv, TAG_TRANSFER, MPI_COMM_WORLD, &requests[nreq++]);

  MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

  if(recv >= 0)
    {
      for(i = 0; i < nmemb; i++)
	{
	  if(compar(element, base + i * size) > 0)
	    terminate("wrong order");
	}
    }

  myfree(element);
  myfree(nlist);
}

#ifdef OPENMP
void serial_sort_omp(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *))
{
  size_t n[OPENMP], start[OPENMP];
  size_t *ln, *lstart;
  sort_type *tmp;
  int tid, i;
  int lev, nthreads, maxlevel;
  sort_type *b2, *local_tmp;

  maxlevel = floorLog2(OPENMP);
  if((size % sizesort) != 0)
    terminate("size of structure must be a multiple of sizeof(sort_type)");

  ln = mymalloc("ln", sizeof(size_t) * OPENMP * (maxlevel + 1));
  lstart = mymalloc("lstart", sizeof(size_t) * OPENMP * (maxlevel + 1));
  tmp = (sort_type *) mymalloc("tmp", size * nmemb);

  for(lev = maxlevel; lev >= 0; lev--)
    {
      nthreads = 1 << (maxlevel - lev);
      for(i = 0; i < nthreads; i += 2)
	{
	  if(lev == maxlevel)
	    {
	      ln[lev * OPENMP + i] = nmemb;
	      lstart[lev * OPENMP + i] = 0;
	    }
	  else
	    {
	      ln[lev * OPENMP + i] = ln[(lev + 1) * OPENMP + i / 2] / 2;
	      lstart[lev * OPENMP + i] = lstart[(lev + 1) * OPENMP + i / 2];
	      ln[lev * OPENMP + i + 1] = ln[(lev + 1) * OPENMP + i / 2] - ln[lev * OPENMP + i];
	      lstart[lev * OPENMP + i + 1] = lstart[(lev + 1) * OPENMP + i / 2] + ln[lev * OPENMP + i];
	    }
	}
    }

  for(lev = 0; lev <= maxlevel; lev++)
    {
      nthreads = 1 << (maxlevel - lev);

      if(ThisTask == 0)
	printf("Task:%d using %d threads in level %d\n", ThisTask, nthreads, lev);

#pragma omp parallel for private(b2, local_tmp)
      for(tid = nthreads - 1; tid >= 0; tid--)
	{
	  local_tmp = tmp + lstart[lev * OPENMP + tid] * size / sizesort;
	  b2 = (sort_type *) base + lstart[lev * OPENMP + tid] * size / sizesort;

	  msort_serial_with_tmp(b2, ln[lev * OPENMP + tid], size, compar, local_tmp, 0,
				(lev == 0) ? 10000 : 1);
	}

    }

  myfree(tmp);
  myfree(lstart);
  myfree(ln);
}

int floorLog2(unsigned int n)
{
  int pos = 0;
  if(n >= 1 << 16)
    {
      n >>= 16;
      pos += 16;
    }
  if(n >= 1 << 8)
    {
      n >>= 8;
      pos += 8;
    }
  if(n >= 1 << 4)
    {
      n >>= 4;
      pos += 4;
    }
  if(n >= 1 << 2)
    {
      n >>= 2;
      pos += 2;
    }
  if(n >= 1 << 1)
    {
      pos += 1;
    }
  return ((n == 0) ? (-1) : pos);
}

#endif

void serial_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *))
{
  const size_t storage = nmemb * size;

  if((size % sizesort) != 0)
    terminate("size of structure must be a multiple of sizeof(sort_type)");
  sort_type *tmp = (sort_type *) mymalloc("int *tmp", storage);

  msort_serial_with_tmp((sort_type *) base, nmemb, size, compar, tmp, 0, 10000);

  myfree(tmp);
}


void msort_serial_with_tmp(sort_type * base, size_t n, size_t s,
				  int (*compar) (const void *, const void *), sort_type * t, int level,
				  int maxlevel)
{
  sort_type *tmp;
  sort_type *b1, *b2;
  size_t n1, n2;

  if(n <= 1 || level >= maxlevel)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = base;
  b2 = base + n1 * s / sizesort;

  msort_serial_with_tmp(b1, n1, s, compar, t, level + 1, maxlevel);
  msort_serial_with_tmp(b2, n2, s, compar, t, level + 1, maxlevel);

  tmp = t;
  size_t ss, st;
  st = s / sizesort;
  while(n1 > 0 && n2 > 0)
    {
      ss = st;
      if(compar(b1, b2) < 0)
	{
	  --n1;
	  while(ss--)
	    *tmp++ = *b1++;

	}
      else
	{
	  --n2;
	  while(ss--)
	    *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * s);
  memcpy(base, t, (n - n2) * s);
}

void qs(sort_type * base, size_t n, size_t size, int (*compar) (const void *, const void *))
{
  register int i, j, k;		//, ss;
  register int ss, st;
  register sort_type *pivot, t;

  if(n <= 1)
    return;

//  printf("ad, n:%d\n", n);
//  pivot = malloc(size);
//  t = malloc(size);
//  memcpy(pivot, base, size);

  st = size / sizesort;
  pivot = (base + 2 * st);



  i = 0;
  j = n;
  for(;;)
    {
      while(compar(base + (++i) * st, pivot) < 0 && i < n)
	{
	}
      while(compar(base + (--j) * st, pivot) > 0 && j >= 0)
	{
	}
      if(i >= j)
	break;

      ss = st;
      k = 0;
      while(ss--)
	{
	  t = *(base + i * st + k);
	  *(base + i * st + k) = *(base + j * st + k);
	  *(base + j * st + k) = t;
	  k++;
	}

//    memcpy(t, base + i*st ,size);
//    memcpy(base + i*st, base + j*st, size);
//    memcpy(base + j*st, t ,size);
    }


  ss = st;
  k = 0;
  while(ss--)
    {
      t = *(base + (i - 1) * st + k);
      *(base + (i - 1) * st + k) = *(base + k);
      *(base + k) = t;
      k++;
    }

  //printf("n %d\n", n);
  // for(k = 0; k < n; k++)
  //   printf("%d %lld\n", k, (long long) (*(base + k * st)));
  // printf(" --------------------------------\n");

//   memcpy(t, base + (i-1)*st, size);
//   memcpy(base + (i-1)*st, base, size);
//   memcpy(base, t, size);

//#pragma omp task
  qs(base, i - 1, size, compar);
  qs(base + i * st, n - i, size, compar);
//#pragma omp taskwait

}

typedef long WORD;
#define W sizeof(WORD)		/* must be a power of 2 */
#define SWAPINIT(a, es) swaptype =                \
    ((a-(char*)0) | (es)) % W ? 2 : es > W ? 1 : 0

#define exch(a, b, t) (t = a, a = b, b = t);
#define swap(a, b)                                \
   swaptype != 0 ? swapfunc(a, b, es, swaptype) : \
   (void)exch(*(WORD*)(a), *(WORD*)(b), t)

#define vecswap(a, b, n) if (n > 0) swapfunc(a, b, n, swaptype)

#include <stddef.h>
static void swapfunc(char *a, char *b, size_t n, int swaptype)
{
  if(swaptype <= 1)
    {
      WORD t;
      for(; n > 0; a += W, b += W, n -= W)
	exch(*(WORD *) a, *(WORD *) b, t);
    }
  else
    {
      char t;
      for(; n > 0; a += 1, b += 1, n -= 1)
	exch(*a, *b, t);
    }
}

#define PVINIT(pv, pm)                          \
   if (swaptype != 0) { pv = a; swap(pv, pm); } \
   else { pv = (char*)&v; v = *(WORD*)pm; }

#define min(a, b) ((a) < (b) ? (a) : (b))

static char *med3(char *a, char *b, char *c, int (*cmp) ())
{
  return cmp(a, b) < 0 ?
    (cmp(b, c) < 0 ? b : cmp(a, c) < 0 ? c : a) : (cmp(b, c) > 0 ? b : cmp(a, c) > 0 ? c : a);
}

void quicksort(void *av, size_t n, size_t es, int (*cmp) (const void *, const void *))
{
  char *a, *pa, *pb, *pc, *pd, *pl, *pm, *pn, *pv;
  int r, swaptype;
  WORD t, v;
  size_t s;

  a = (char *) av;


  SWAPINIT(a, es);
  if(n < 7)
    {				/* Insertion sort on smallest arrays */
      for(pm = a + es; pm < a + n * es; pm += es)
	for(pl = pm; pl > a && cmp(pl - es, pl) > 0; pl -= es)
	  swap(pl, pl - es);
      return;
    }
  pm = a + (n / 2) * es;	/* Small arrays, middle element */
  if(n > 7)
    {
      pl = a;
      pn = a + (n - 1) * es;
      if(n > 40)
	{			/* Big arrays, pseudomedian of 9 */
	  s = (n / 8) * es;
	  pl = med3(pl, pl + s, pl + 2 * s, cmp);
	  pm = med3(pm - s, pm, pm + s, cmp);
	  pn = med3(pn - 2 * s, pn - s, pn, cmp);
	}
      pm = med3(pl, pm, pn, cmp);	/* Mid-size, med of 3 */
    }
  PVINIT(pv, pm);		/* pv points to partition value */
  pa = pb = a;
  pc = pd = a + (n - 1) * es;
  for(;;)
    {
      while(pb <= pc && (r = cmp(pb, pv)) <= 0)
	{
	  if(r == 0)
	    {
	      swap(pa, pb);
	      pa += es;
	    }
	  pb += es;
	}
      while(pc >= pb && (r = cmp(pc, pv)) >= 0)
	{
	  if(r == 0)
	    {
	      swap(pc, pd);
	      pd -= es;
	    }
	  pc -= es;
	}
      if(pb > pc)
	break;
      swap(pb, pc);
      pb += es;
      pc -= es;
    }
  pn = a + n * es;
  s = min(pa - a, pb - pa);
  vecswap(a, pb - s, s);
  s = min(pd - pc, pn - pd - es);
  vecswap(pb, pn - s, s);


#pragma omp task
  if((s = pb - pa) > es)
    quicksort(a, s / es, es, cmp);
  if((s = pd - pc) > es)
    quicksort(pn - s, s / es, es, cmp);
#pragma omp taskwait
}


#ifdef OMP_SORT
void omp_qsort(void *base, size_t nmemb, size_t size, int (*compar) ())
{

  if(nmemb > OMP_SORT)
    {
#ifdef VERBOSE      
      if(ThisTask == 0)
	printf("Parallel quicksort\n");
#endif

#pragma omp parallel
#pragma omp single
      {
	quicksort(base, nmemb, size, compar);
      }
    }
  else
    {
      qsort(base, nmemb, size, compar);
    }
}
#endif

#endif
