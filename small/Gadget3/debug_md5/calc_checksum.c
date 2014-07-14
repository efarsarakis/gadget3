
#include "../allvars.h"
#include "../proto.h"

#ifdef CHECKSUM_DEBUG

#include "Md5.h"

void calc_memory_checksum(void *base, size_t bytes)
{
  MD5_CTX sum;
  union 
  {  unsigned char digest[16];   
    int val[4]; 
  } u, uglob;
  

  MD5Init(&sum);
  MD5Update(&sum, base, bytes);
  MD5Final(&sum);

  int i;
  
  for (i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Step=%d  MD5=", All.NumCurrentTiStep);
      for (i = 0; i < 16; i++)
	printf ("%02x", uglob.digest[i]);
      printf("\n");
    }
}


#endif

