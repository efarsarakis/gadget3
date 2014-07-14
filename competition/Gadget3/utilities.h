#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdlib.h>
#include <stdio.h>
#include "gadgetconfig.h"

#ifdef PARALLEL
int ThisTask;
int NTask;
#endif

#define safe_fgets( str, num, stream ) util_fgets( str, num, stream, __FILE__, __LINE__ )
#define safe_fread( ptr, size, count, stream ) util_fread( ptr, size, count, stream, __FILE__, __LINE__ )

void myprintf( const char* format, ... );
char * util_fgets( char *str, int num, FILE *stream, char *file, int line );
size_t util_fread( void *ptr, size_t size, size_t count, FILE * stream, char *file, int line );

__inline static int imin( int a, int b ) { return a < b ? a : b; }
__inline static int imax( int a, int b ) { return a > b ? a : b; }
__inline static double dmin( double a, double b ) { return a < b ? a : b; }
__inline static double dmax( double a, double b ) { return a > b ? a : b; }

#define VECT_NORM(N, V) ({double VN; int vcc; for(vcc = 0, VN = 0; vcc < (N); vcc++) VN+=((V)[vcc]*(V)[vcc]); VN;})

#endif /* UTILITIES_H */
