#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef MPISENDRECV_CHECKSUM

#undef MPI_Sendrecv


int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		       int dest, int sendtag, void *recvbufreal, int recvcount,
		       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status)
{
  int checksumtag = 1000, errtag = 2000;
  int i, iter = 0, err_flag, err_flag_imported, size_sendtype, size_recvtype, Local_ThisTask, Local_NTask;
  long long sendCheckSum, recvCheckSum, importedCheckSum;
  unsigned char *p, *buf, *recvbuf;

  if(dest != source)
    endrun(3);

  MPI_Comm_rank(comm, &Local_ThisTask);
  MPI_Comm_size(comm, &Local_NTask);

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  if(dest == Local_ThisTask)
    {
      memcpy(recvbufreal, sendbuf, recvcount * size_recvtype);
      return 0;
    }


  if(!(buf = mymalloc("buf", recvcount * size_recvtype + 1024)))
    endrun(6);

  for(i = 0, p = buf; i < recvcount * size_recvtype + 1024; i++)
    *p++ = 255;

  recvbuf = buf + 512;

  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
	       recvbuf, recvcount, recvtype, source, recvtag, comm, status);

  for(i = 0, p = buf; i < 512; i++, p++)
    {
      if(*p != 255)
	{
	  printf
	    ("MPI-ERROR: Task=%d/%s: Recv occured before recv buffer. message-size=%d from %d, i=%d c=%d\n",
	     Local_ThisTask, getenv("HOST"), recvcount, dest, i, *p);
	  fflush(stdout);
	  endrun(6);
	}
    }

  for(i = 0, p = recvbuf + recvcount * size_recvtype; i < 512; i++, p++)
    {
      if(*p != 255)
	{
	  printf
	    ("MPI-ERROR: Task=%d/%s: Recv occured after recv buffer. message-size=%d from %d, i=%d c=%d\n",
	     Local_ThisTask, getenv("HOST"), recvcount, dest, i, *p);
	  fflush(stdout);
	  endrun(6);
	}
    }


  for(i = 0, p = sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
    sendCheckSum += *p;

  importedCheckSum = 0;

  if(dest > Local_ThisTask)
    {
      if(sendcount > 0)
	MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
      if(recvcount > 0)
	MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, comm, status);
    }
  else
    {
      if(recvcount > 0)
	MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, comm, status);
      if(sendcount > 0)
	MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
    }

  checksumtag++;

  for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount * size_recvtype; i++, p++)
    recvCheckSum += *p;


  err_flag = err_flag_imported = 0;

  if(recvCheckSum != importedCheckSum)
    {
      printf
	("MPI-ERROR: Receive error on task=%d/%s from task=%d, message size=%d, sendcount=%d checksums= %d %d  %d %d. Try to fix it...\n",
	 Local_ThisTask, getenv("HOST"), source, recvcount, sendcount, (int) (recvCheckSum >> 32),
	 (int) recvCheckSum, (int) (importedCheckSum >> 32), (int) importedCheckSum);
      fflush(stdout);

      err_flag = 1;
    }

  if(dest > Local_ThisTask)
    {
      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
    }
  else
    {
      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
    }
  errtag++;

  if(err_flag > 0 || err_flag_imported > 0)
    {
      printf("Task=%d is on %s, wants to send %d and has checksum=%d %d of send data\n",
	     Local_ThisTask, getenv("HOST"), sendcount, (int) (sendCheckSum >> 32), (int) sendCheckSum);
      fflush(stdout);

      do
	{
	  sendtag++;
	  recvtag++;

	  for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
	    *p = 0;

	  if((iter & 1) == 0)
	    {
	      if(dest > Local_ThisTask)
		{
		  if(sendcount > 0)
		    MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, comm);
		  if(recvcount > 0)
		    MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, comm, status);
		}
	      else
		{
		  if(recvcount > 0)
		    MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, comm, status);
		  if(sendcount > 0)
		    MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, comm);
		}
	    }
	  else
	    {
	      if(iter > 5)
		{
		  printf("we're trying to send each byte now on task=%d (iter=%d)\n", Local_ThisTask, iter);
		  fflush(stdout);
		  if(dest > Local_ThisTask)
		    {
		      for(i = 0, p = sendbuf; i < sendcount * size_sendtype; i++, p++)
			MPI_Ssend(p, 1, MPI_BYTE, dest, i, comm);
		      for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
			MPI_Recv(p, 1, MPI_BYTE, dest, i, comm, status);
		    }
		  else
		    {
		      for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
			MPI_Recv(p, 1, MPI_BYTE, dest, i, comm, status);
		      for(i = 0, p = sendbuf; i < sendcount * size_sendtype; i++, p++)
			MPI_Ssend(p, 1, MPI_BYTE, dest, i, comm);
		    }
		}
	      else
		{
		  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
			       recvbuf, recvcount, recvtype, source, recvtag, comm, status);
		}
	    }

	  importedCheckSum = 0;

	  for(i = 0, p = sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
	    sendCheckSum += *p;

	  printf("Task=%d gas send_checksum=%d %d\n", Local_ThisTask, (int) (sendCheckSum >> 32),
		 (int) sendCheckSum);
	  fflush(stdout);

	  if(dest > Local_ThisTask)
	    {
	      if(sendcount > 0)
		MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
	      if(recvcount > 0)
		MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag,
			 comm, status);
	    }
	  else
	    {
	      if(recvcount > 0)
		MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag,
			 comm, status);
	      if(sendcount > 0)
		MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
	    }

	  for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount; i++, p++)
	    recvCheckSum += *p;

	  err_flag = err_flag_imported = 0;

	  if(recvCheckSum != importedCheckSum)
	    {
	      printf
		("MPI-ERROR: Again (iter=%d) a receive error on task=%d/%s from task=%d, message size=%d, checksums= %d %d  %d %d. Try to fix it...\n",
		 iter, Local_ThisTask, getenv("HOST"), source, recvcount, (int) (recvCheckSum >> 32),
		 (int) recvCheckSum, (int) (importedCheckSum >> 32), (int) importedCheckSum);
	      fflush(stdout);
	      err_flag = 1;
	    }

	  if(dest > Local_ThisTask)
	    {
	      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
	      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
	    }
	  else
	    {
	      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
	      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
	    }

	  if(err_flag == 0 && err_flag_imported == 0)
	    break;

	  errtag++;
	  checksumtag++;
	  iter++;
	}
      while(iter < 10);

      if(iter >= 10)
	{
	  char buf[1000];
	  int length;
	  FILE *fd;

	  sprintf(buf, "send_data_%d.dat", Local_ThisTask);
	  fd = fopen(buf, "w");
	  length = sendcount * size_sendtype;
	  fwrite(&length, 1, sizeof(int), fd);
	  fwrite(sendbuf, sendcount, size_sendtype, fd);
	  fclose(fd);

	  sprintf(buf, "recv_data_%d.dat", Local_ThisTask);
	  fd = fopen(buf, "w");
	  length = recvcount * size_recvtype;
	  fwrite(&length, 1, sizeof(int), fd);
	  fwrite(recvbuf, recvcount, size_recvtype, fd);
	  fclose(fd);

	  printf("MPI-ERROR: Even 10 trials proved to be insufficient on task=%d/%s. Stopping\n",
		 Local_ThisTask, getenv("HOST"));
	  fflush(stdout);
	  endrun(10);
	}
    }

  memcpy(recvbufreal, recvbuf, recvcount * size_recvtype);

  myfree(buf);

  return 0;
}

#endif