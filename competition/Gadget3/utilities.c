#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "utilities.h"
#include "proto.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

void myprintf(const char *format, ...)
{
#ifdef PARALLEL
  printf("Task %03d: ", ThisTask);
#endif
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

char *util_fgets(char *str, int num, FILE * stream, char *file, int line)
{
  char *ret = fgets(str, num, stream);
  if(ret == NULL)
    {
      printf("error: fgets in file %s at line %d\n", file, line);
      endrun(200);
    }

  return ret;
}

size_t util_fread(void *ptr, size_t size, size_t count, FILE * stream, char *file, int line)
{
  size_t result = fread(ptr, size, count, stream);
  if(result != count)
    {
      printf("error: fread in file %s at line %d\n", file, line);
      endrun(201);
    }

  return result;
}
