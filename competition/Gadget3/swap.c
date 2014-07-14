#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include "utilities.h"
#include "swap.h"

double SwapDouble(double Val)
{
  double nVal;
  int i;
  const char *readFrom = (const char *) &Val;
  char *writeTo = ((char *) &nVal) + sizeof(nVal);

  for(i = 0; i < sizeof(Val); ++i)
    {
      *(--writeTo) = *(readFrom++);
    }
  return nVal;
}

float SwapFloat(float Val)
{
  float nVal;
  int i;
  const char *readFrom = (const char *) &Val;
  char *writeTo = ((char *) &nVal) + sizeof(nVal);

  for(i = 0; i < sizeof(Val); ++i)
    {
      *(--writeTo) = *(readFrom++);
    }
  return nVal;
}

int SwapInt(int Val)
{
  int nVal;
  int i;
  const char *readFrom = (const char *) &Val;
  char *writeTo = ((char *) &nVal) + sizeof(nVal);

  for(i = 0; i < sizeof(Val); ++i)
    {
      *(--writeTo) = *(readFrom++);
    }
  return nVal;
}

int CheckSwap(char *fname, int *swap)
{
  FILE *fd;
  off_t fsize, fpos;
  int blocksize, blockend;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't open file `%s'.\n", fname);
      return -1;
    }

  fseeko(fd, 0, SEEK_END);
  fsize = ftello(fd);

  *swap = 0;
  fpos = 0;
  fseeko(fd, 0, SEEK_SET);
  safe_fread(&blocksize, sizeof(int), 1, fd);
  while(!feof(fd))
    {
      if(fpos + blocksize + 4 > fsize)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4 + blocksize;
      fseeko(fd, fpos, SEEK_SET);
      safe_fread(&blockend, sizeof(int), 1, fd);
      if(blocksize != blockend)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4;
      if(!fread(&blocksize, sizeof(int), 1, fd))
	break;
    }

  if(*swap == 0)
    {
      fclose(fd);
      return 0;
    }

  fpos = 0;
  fseeko(fd, 0, SEEK_SET);
  safe_fread(&blocksize, sizeof(int), 1, fd);
  while(!feof(fd))
    {
      blocksize = SwapInt(blocksize);
      if(fpos + blocksize + 4 > fsize)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4 + blocksize;
      fseeko(fd, fpos, SEEK_SET);
      safe_fread(&blockend, sizeof(int), 1, fd);
      blockend = SwapInt(blockend);
      if(blocksize != blockend)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4;
      if(!fread(&blocksize, sizeof(int), 1, fd))
	break;
    }

  fclose(fd);
  return 0;
}
