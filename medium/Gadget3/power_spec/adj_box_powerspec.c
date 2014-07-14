#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>

#include "../allvars.h"
#include "../proto.h"
#include "../forcetree.h"

#if defined(PERIODIC) && defined(ADJ_BOX_POWERSPEC)

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif


#define  FOURIER_GRID  (All.FourierGrid)
#define  FOURIER_GRID2 (2*(FOURIER_GRID/2 + 1))


typedef long long large_array_offset;

static rfftwnd_mpi_plan fft_forward_plan;
static int slabstart_x, nslab_x, slabstart_y, nslab_y;

static int fftsize, maxfftsize;
static fftw_real *velfield[3], *workspace;
static fftw_real *curlfield;
static fftw_real *vorticityfield[3];
static fftw_complex *fft_of_field;

static float *powerspec_nearest_distance, *powerspec_nearest_hsml;


static struct data_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *DataIn, *DataGet;

static struct data_out
{
  MyFloat Distance;
  MyDouble Vel[3];
  MyDouble Curl;
  MyDouble Vorticity[3];
}
 *DataResult, *DataOut;



#define BINS_PS  2000		/* number of bins for power spectrum computation */

static long long CountModes[BINS_PS];
static double SumPower[BINS_PS];
static double Power[BINS_PS];
static double Kbin[BINS_PS];
static double K0, K1;
static double binfac;
static double vel_disp[3];

void adj_box_compute_power_spectrum(fftw_real *field, int count);

void adj_box_powerspec(void)
{
  char fname[500];
  int i;

  mpi_printf("Starting velocity and velocity curl powerspec computation\n");

  double tstart = second();

  /* Set up the FFTW plan  */
  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, FOURIER_GRID, FOURIER_GRID, FOURIER_GRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */
  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);
  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* allocate the memory to hold the FFT fields */
  velfield[0] = (fftw_real *) mymalloc("velgrid[0]", maxfftsize * sizeof(fftw_real));
  velfield[1] = (fftw_real *) mymalloc("velgrid[1]", maxfftsize * sizeof(fftw_real));
  velfield[2] = (fftw_real *) mymalloc("velgrid[2]", maxfftsize * sizeof(fftw_real));
  curlfield = (fftw_real *) mymalloc("curlgrid", maxfftsize * sizeof(fftw_real));
  vorticityfield[0] = (fftw_real *) mymalloc("vorticitygrid[0]", maxfftsize * sizeof(fftw_real));
  vorticityfield[1] = (fftw_real *) mymalloc("vorticitygrid[1]", maxfftsize * sizeof(fftw_real));
  vorticityfield[2] = (fftw_real *) mymalloc("vortycitygrid[2]", maxfftsize * sizeof(fftw_real));

  workspace = (fftw_real *) mymalloc("workspace", maxfftsize * sizeof(fftw_real));

  memset(velfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(velfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(velfield[2], 0, maxfftsize * sizeof(fftw_real));
  memset(curlfield, 0, maxfftsize * sizeof(fftw_real));
  memset(vorticityfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(vorticityfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(vorticityfield[2], 0, maxfftsize * sizeof(fftw_real));

  adj_box_powerspec_obtain_fields();
  adj_box_powerspec_calc_dispersion();


  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }


  adj_box_compute_power_spectrum(velfield[0], 1);
  adj_box_compute_power_spectrum(velfield[1], 0);
  adj_box_compute_power_spectrum(velfield[2], 0);
  adj_box_powerspec_gather_results();

  sprintf(fname, "%s/vel_powerspec_%03d.txt", All.OutputDir, RestartSnapNum);
  adj_box_powerspec_save(fname);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  adj_box_compute_power_spectrum(curlfield, 1);
  adj_box_powerspec_gather_results();

  sprintf(fname, "%s/curl_powerspec_%03d.txt", All.OutputDir, RestartSnapNum);
  adj_box_powerspec_save(fname);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  adj_box_compute_power_spectrum(vorticityfield[0], 1);
  adj_box_compute_power_spectrum(vorticityfield[1], 0);
  adj_box_compute_power_spectrum(vorticityfield[2], 0);
  adj_box_powerspec_gather_results();

  sprintf(fname, "%s/vorticity_powerspec_%03d.txt", All.OutputDir, RestartSnapNum);
  adj_box_powerspec_save(fname);

  myfree(workspace);
  myfree(vorticityfield[2]);
  myfree(vorticityfield[1]);
  myfree(vorticityfield[0]);
  myfree(curlfield);
  myfree(velfield[2]);
  myfree(velfield[1]);
  myfree(velfield[0]);

  rfftwnd_mpi_destroy_plan(fft_forward_plan);

  double tend = second();

  if(ThisTask == 0)
    {
      printf("end gas power spectra. took %g seconds\n", timediff(tstart, tend));
      fflush(stdout);
    }
}

void adj_box_compute_power_spectrum(fftw_real *field, int count)
{
  double k2, kx, ky, kz;
  int x, y, z, zz, ip;

  double BoxSize = All.BoxWidth;

  K0 = 2 * M_PI / BoxSize;	/* minimum k */
  K1 = K0 * FOURIER_GRID / 2;		/* maximum k */
  binfac = BINS_PS / (log(K1) - log(K0));

      /* Do the FFT of the field  */

      rfftwnd_mpi(fft_forward_plan, 1, field, workspace, FFTW_TRANSPOSED_ORDER);

      fft_of_field = (fftw_complex *) field;

      double posum = 0, posum_all = 0;

      for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
	for(x = 0; x < FOURIER_GRID; x++)
	  for(z = 0; z < FOURIER_GRID; z++)
	    {
	      zz = z;
	      if(z >= FOURIER_GRID / 2 + 1)
		zz = FOURIER_GRID - z;

	      if(x > FOURIER_GRID / 2)
		kx = x - FOURIER_GRID;
	      else
		kx = x;
	      if(y > FOURIER_GRID / 2)
		ky = y - FOURIER_GRID;
	      else
		ky = y;
	      if(z > FOURIER_GRID / 2)
		kz = z - FOURIER_GRID;
	      else
		kz = z;

	      k2 = kx * kx + ky * ky + kz * kz;

	      ip = FOURIER_GRID * (FOURIER_GRID / 2 + 1) * (y - slabstart_y) + (FOURIER_GRID / 2 + 1) * x + zz;

	      double po = (fft_of_field[ip].re * fft_of_field[ip].re
			   + fft_of_field[ip].im * fft_of_field[ip].im) / pow(FOURIER_GRID, 6);

	      posum += po;

	      if(k2 > 0)
		{
		  if(k2 < (FOURIER_GRID / 2.0) * (FOURIER_GRID / 2.0))
		    {
		      double k = sqrt(k2) * 2 * M_PI / BoxSize;

		      if(k >= K0 && k < K1)
			{
			  int bin = log(k / K0) * binfac;

			  SumPower[bin] += po;

			  if(count)
			    CountModes[bin] += 1;
			}
		    }
		}
	    }

      MPI_Allreduce(&posum, &posum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      /*if(ThisTask == 0)
	printf("DIM=%d: vel_disp=%g   posum = %g\n", dim, vel_disp[dim], posum_all);*/
}

void adj_box_powerspec_gather_results()
{
  /* Now compute the power spectrum */

  long long int *countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  double *powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));
  int i, n;

  MPI_Allgather(CountModes, BINS_PS * sizeof(long long), MPI_BYTE,
		countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      CountModes[i] = 0;
      for(n = 0; n < NTask; n++)
	CountModes[i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPower, BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      for(n = 0; n < NTask; n++)
	SumPower[i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);

  for(i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[i] > 0)
	Power[i] = SumPower[i] / CountModes[i];
      else
	Power[i] = 0;
    }
}

void adj_box_powerspec_save(const char *fname)
{
  FILE *fd;
  char buf[500];
  int i;

  if(ThisTask == 0)
    {

      if(!(fd = fopen(fname, "w")))
	{
	  sprintf(buf, "can't open file `%s`\n", fname);
	  terminate(buf);
	}

      fprintf(fd, "%g\n", All.Time);
      i = FOURIER_GRID;
      fprintf(fd, "%d\n", i);
      i = BINS_PS;
      fprintf(fd, "%d\n", i);

      fprintf(fd, "%g\n", vel_disp[0]);
      fprintf(fd, "%g\n", vel_disp[1]);
      fprintf(fd, "%g\n", vel_disp[2]);

      for(i = 0; i < BINS_PS; i++)
	{
	  fprintf(fd, "%g %g %g %g\n", Kbin[i], Power[i], (double) CountModes[i], SumPower[i]);
	}

      fclose(fd);
    }
}



/* this function determines the velocity fields by using the nearest cell's values 
 */
double adj_box_powerspec_obtain_fields(void)
{
  int j, dummy;
  long long ntot, npleft;
  int ndone, ndone_flag, ngrp, sendTask, recvTask, place, nexport, nimport, iter;

  double tstart = second();

  double BoxSize = All.BoxWidth;
  double BoxEdge[3] = {All.BoxCenter_x - 0.5 * BoxSize, All.BoxCenter_y - 0.5 * BoxSize, All.BoxCenter_z - 0.5 * BoxSize};

  if(ThisTask == 0)
    {
      printf("Start finding nearest gas-particle for mesh-cell centers (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  large_array_offset i, n, Ncount = ((large_array_offset) nslab_x) * (FOURIER_GRID * FOURIER_GRID);	/* number of grid points on the local slab */

  powerspec_nearest_distance =
    (float *) mymalloc("powerspec_vel_nearest_distance", sizeof(float) * Ncount);
  powerspec_nearest_hsml = (float *) mymalloc("powerspec_vel_nearest_hsml", sizeof(float) * Ncount);

  for(n = 0; n < Ncount; n++)
    {
      powerspec_nearest_distance[n] = 1.0e30;
      powerspec_nearest_hsml[n] = BoxSize / pow(All.TotN_gas, 1.0 / 3);
    }

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", Ncount * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct data_in) + sizeof(struct data_out) +
					     sizemax(sizeof(struct data_in), sizeof(struct data_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  iter = 0;
  /* we will repeat the whole thing for those points where we didn't find enough neighbours */
  do
    {
      i = 0;			/* begin with this index */

      do
	{
	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */
	  for(nexport = 0; i < Ncount; i++)
	    {
	      if(powerspec_nearest_distance[i] > 1.0e29)
		{
		  if(adj_box_powerspec_find_nearest_evaluate(i, 0, &nexport, Send_count) < 0)
		    break;
		}
	    }

	  mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

	  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  DataGet = (struct data_in *) mymalloc("DataGet", nimport * sizeof(struct data_in));
	  DataIn = (struct data_in *) mymalloc("DataIn", nexport * sizeof(struct data_in));

	  if(ThisTask == 0)
	    {
	      printf("still finding nearest... (presently allocated=%g MB)\n",
		     AllocatedBytes / (1024.0 * 1024.0));
	      fflush(stdout);
	    }

	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      int xx = place / (FOURIER_GRID * FOURIER_GRID);
	      int yy = (place - xx * FOURIER_GRID * FOURIER_GRID) / FOURIER_GRID;
	      int zz = (place - xx * FOURIER_GRID * FOURIER_GRID - yy * FOURIER_GRID);
	      xx += slabstart_x;

	      double x = (xx + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[0];
	      double y = (yy + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[1];
	      double z = (zz + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[2];

	      DataIn[j].Pos[0] = x;
	      DataIn[j].Pos[1] = y;
	      DataIn[j].Pos[2] = z;
	      DataIn[j].Hsml = powerspec_nearest_hsml[place];

	      memcpy(DataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	    }

	  /* exchange particle data */
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct data_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct data_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(DataIn);
	  DataResult = (struct data_out *) mymalloc("DataResult", nimport * sizeof(struct data_out));
	  DataOut = (struct data_out *) mymalloc("DataOut", nexport * sizeof(struct data_out));

	  for(j = 0; j < nimport; j++)
	    adj_box_powerspec_find_nearest_evaluate(j, 1, &dummy, &dummy);

	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct data_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct data_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }

	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      if(DataOut[j].Distance < powerspec_nearest_distance[place])
		{
		  powerspec_nearest_distance[place] = DataOut[j].Distance;

		  int ii = place / (FOURIER_GRID * FOURIER_GRID);
		  int jj = (place - ii * FOURIER_GRID * FOURIER_GRID) / FOURIER_GRID;
		  int kk = (place - ii * FOURIER_GRID * FOURIER_GRID - jj * FOURIER_GRID);
		  int ip = FOURIER_GRID2 * (FOURIER_GRID * ii + jj) + kk;

		  velfield[0][ip] = DataOut[j].Vel[0];
		  velfield[1][ip] = DataOut[j].Vel[1];
		  velfield[2][ip] = DataOut[j].Vel[2];
		  curlfield[ip] = DataOut[j].Curl;
		  vorticityfield[0][ip] = DataOut[j].Vorticity[0];
		  vorticityfield[1][ip] = DataOut[j].Vorticity[1];
		  vorticityfield[2][ip] = DataOut[j].Vorticity[2];
		}
	    }

	  if(i >= Ncount)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	  myfree(DataOut);
	  myfree(DataResult);
	  myfree(DataGet);
	}
      while(ndone < NTask);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < Ncount; i++)
	{
	  if(powerspec_nearest_distance[i] > 1.0e29)
	    {
	      /* need to redo this particle */
	      npleft++;
	      powerspec_nearest_hsml[i] *= 2.0;
	      if(iter >= MAXITER - 10)
		{
		  int xx = i / (FOURIER_GRID * FOURIER_GRID);
		  int yy = (i - xx * FOURIER_GRID * FOURIER_GRID) / FOURIER_GRID;
		  int zz = (i - xx * FOURIER_GRID * FOURIER_GRID - yy * FOURIER_GRID);
		  xx += slabstart_x;

		  double x = (xx + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[0];
		  double y = (yy + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[1];
		  double z = (zz + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[2];

		  printf("i=%d task=%d Hsml=%g  pos=(%g|%g|%g)\n",
			 (int) i, ThisTask, powerspec_nearest_hsml[i], x, y, z);
		  fflush(stdout);
		}
	    }
	  else
	    {
	      powerspec_nearest_distance[i] = 0;	/* we not continue to search for this particle */
	    }
	}

      sumup_longs(1, &npleft, &ntot);
      if(ntot > 0)
	{
	  iter++;
	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("powespec_vel nearest iteration %d: need to repeat for %lld particles.\n", iter, ntot);
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    terminate("failed to converge");
	}
    }
  while(ntot > 0);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  myfree(powerspec_nearest_hsml);
  myfree(powerspec_nearest_distance);

  mpi_printf("done finding velocity field\n");

  double tend = second();
  return timediff(tstart, tend);
}


void adj_box_powerspec_calc_dispersion(void)
{
  int dim, i, j, k;

  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i = 0; i < nslab_x; i++)
	for(j = 0; j < FOURIER_GRID; j++)
	  for(k = 0; k < FOURIER_GRID; k++)
	    {
	      int ip = FOURIER_GRID2 * (FOURIER_GRID * i + j) + k;

	      vsum += velfield[dim][ip];
	    }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(FOURIER_GRID, 3);

      for(i = 0; i < nslab_x; i++)
	for(j = 0; j < FOURIER_GRID; j++)
	  for(k = 0; k < FOURIER_GRID; k++)
	    {
	      int ip = FOURIER_GRID2 * (FOURIER_GRID * i + j) + k;

	      velfield[dim][ip] -= vmean;
	    }

      for(i = 0; i < nslab_x; i++)
	for(j = 0; j < FOURIER_GRID; j++)
	  for(k = 0; k < FOURIER_GRID; k++)
	    {
	      int ip = FOURIER_GRID2 * (FOURIER_GRID * i + j) + k;

	      vdisp += velfield[dim][ip] * velfield[dim][ip];
	    }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      vel_disp[dim] = vdisp_all / pow(FOURIER_GRID, 3);
    }
}



int adj_box_powerspec_find_nearest_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, index, listindex = 0;
  int startnode, numngb_inbox;
  double h, r2max;
  double dx, dy, dz, r2;
  MyDouble pos[3];

  if(mode == 0)
    {
      int xx = target / (FOURIER_GRID * FOURIER_GRID);
      int yy = (target - xx * FOURIER_GRID * FOURIER_GRID) / FOURIER_GRID;
      int zz = (target - xx * FOURIER_GRID * FOURIER_GRID - yy * FOURIER_GRID);
      xx += slabstart_x;

      double BoxSize = All.BoxWidth;
      double BoxEdge[3] = {All.BoxCenter_x - 0.5 * BoxSize, All.BoxCenter_y - 0.5 * BoxSize, All.BoxCenter_z - 0.5 * BoxSize};

      double x = (xx + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[0];
      double y = (yy + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[1];
      double z = (zz + 0.5) / FOURIER_GRID * BoxSize + BoxEdge[2];

      pos[0] = x;
      pos[1] = y;
      pos[2] = z;
      h = powerspec_nearest_hsml[target];
    }
  else
    {
      pos[0] = DataGet[target].Pos[0];
      pos[1] = DataGet[target].Pos[1];
      pos[2] = DataGet[target].Pos[2];
      h = DataGet[target].Hsml;
    }

  index = -1;
  r2max = 1.0e60;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = adj_box_powerspec_treefind(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];
	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

	      /*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;

	      r2 = dx * dx + dy * dy + dz * dz;
	      if(r2 < r2max && r2 < h * h)
		{
		  index = j;
		  r2max = r2;
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
      if(index >= 0)
	{
	  if(index >= N_gas)
	    terminate("index >= N_gas");

	  powerspec_nearest_distance[target] = sqrt(r2max);

	  int i = target / (FOURIER_GRID * FOURIER_GRID);
	  int j = (target - i * FOURIER_GRID * FOURIER_GRID) / FOURIER_GRID;
	  int k = (target - i * FOURIER_GRID * FOURIER_GRID - j * FOURIER_GRID);
	  int ip = FOURIER_GRID2 * (FOURIER_GRID * i + j) + k;

	  velfield[0][ip] = P[index].Vel[0];
	  velfield[1][ip] = P[index].Vel[1];
	  velfield[2][ip] = P[index].Vel[2];

	  curlfield[ip] = SphP[index].r.CurlVel;

	  vorticityfield[0][ip] = SphP[index].Vorticity[0];
	  vorticityfield[1][ip] = SphP[index].Vorticity[1];
	  vorticityfield[2][ip] = SphP[index].Vorticity[2];
	}
    }
  else
    {
      if(index >= 0)
	{
	  if(index >= N_gas)
	    terminate("index >= N_gas");

	  DataResult[target].Distance = sqrt(r2max);
	  DataResult[target].Vel[0] = P[index].Vel[0];
	  DataResult[target].Vel[1] = P[index].Vel[1];
	  DataResult[target].Vel[2] = P[index].Vel[2];

	  DataResult[target].Curl = SphP[index].r.CurlVel;

	  DataResult[target].Vorticity[0] = SphP[index].Vorticity[0];
	  DataResult[target].Vorticity[1] = SphP[index].Vorticity[1];
	  DataResult[target].Vorticity[2] = SphP[index].Vorticity[2];
	}
      else
	DataResult[target].Distance = 2.0e30;
    }
  return 0;
}




int adj_box_powerspec_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			   int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  MyDouble dx, dy, dz, dist, r2;

#define FACT2 0.86602540
#ifdef PERIODIC
  MyDouble xtmp;
#endif
  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type != 0)
	    continue;

	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(mode == 0)
		{
		  if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= All.BunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13005);	/* in this case, the buffer is too small to process even a single particle */
			  for(task = 0; task < NTask; task++)
			    nsend_local[task] = 0;
			  for(no = 0; no < nexport_save; no++)
			    nsend_local[DataIndexTable[no].Task]++;
			  return -1;
			}
		      Exportnodecount[task] = 0;
		      Exportindex[task] = *nexport;
		      DataIndexTable[*nexport].Task = task;
		      DataIndexTable[*nexport].Index = target;
		      DataIndexTable[*nexport].IndexGet = *nexport;
		      *nexport = *nexport + 1;
		      nsend_local[task]++;
		    }

		  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		    DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      if(mode == -1)
		{
		  *nexport = 1;
		}

	      no = Nextnode[no - MaxNodes];
	      continue;

	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}




#endif
