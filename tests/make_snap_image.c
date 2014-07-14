#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define PIXELS  768

const unsigned char red_table[256] = {
  0,   1,   2,   4,   5,   7,   8,  10, 
  11,  13,  14,  15,  17,  18,  20,  21, 
  23,  24,  26,  27,  28,  30,  31,  33, 
  34,  36,  37,  39,  40,  42,  43,  44, 
  46,  47,  49,  50,  52,  53,  55,  56, 
  57,  59,  60,  62,  63,  65,  66,  68, 
  69,  70,  72,  73,  75,  76,  78,  79, 
  81,  82,  84,  85,  86,  88,  89,  91, 
  92,  94,  95,  97,  98,  99, 101, 102, 
  104, 105, 107, 108, 110, 111, 113, 114, 
  115, 117, 118, 120, 121, 123, 124, 126, 
  127, 128, 130, 131, 133, 134, 136, 137, 
  139, 140, 141, 143, 144, 146, 147, 149, 
  150, 152, 153, 155, 156, 157, 159, 160, 
  162, 163, 165, 166, 168, 169, 170, 172, 
  173, 175, 176, 178, 179, 181, 182, 184, 
  185, 186, 188, 189, 191, 192, 194, 195, 
  197, 198, 199, 201, 202, 204, 205, 207, 
  208, 210, 211, 212, 214, 215, 217, 218, 
  220, 221, 223, 224, 226, 227, 228, 230, 
  231, 233, 234, 236, 237, 239, 240, 241, 
  243, 244, 246, 247, 249, 250, 252, 253, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255, 
  255, 255, 255, 255, 255, 255, 255, 255};

const unsigned char green_table[256] = {
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   1,   3,   5,   7,   9,  11,  13, 
  15,  17,  18,  20,  22,  24,  26,  28, 
  30,  32,  34,  35,  37,  39,  41,  43, 
  45,  47,  49,  51,  52,  54,  56,  58, 
  60,  62,  64,  66,  68,  69,  71,  73, 
  75,  77,  79,  81,  83,  85,  86,  88, 
  90,  92,  94,  96,  98, 100, 102, 103, 
  105, 107, 109, 111, 113, 115, 117, 119, 
  120, 122, 124, 126, 128, 130, 132, 134, 
  136, 137, 139, 141, 143, 145, 147, 149, 
  151, 153, 154, 156, 158, 160, 162, 164, 
  166, 168, 170, 171, 173, 175, 177, 179, 
  181, 183, 185, 187, 188, 190, 192, 194, 
  196, 198, 200, 202, 204, 205, 207, 209, 
  211, 213, 215, 217, 219, 221, 222, 224, 
  226, 228, 230, 232, 234, 236, 238, 239, 
  241, 243, 245, 247, 249, 251, 253, 255};

const unsigned char blue_table[256] = {
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   0, 
  0,   0,   0,   0,   0,   0,   0,   3, 
  7,  11,  15,  19,  23,  27,  31,  35, 
  39,  43,  47,  51,  54,  58,  62,  66, 
  70,  74,  78,  82,  86,  90,  94,  98, 
  102, 105, 109, 113, 117, 121, 125, 129, 
  133, 137, 141, 145, 149, 153, 156, 160, 
  164, 168, 172, 176, 180, 184, 188, 192, 
  196, 200, 204, 207, 211, 215, 219, 223, 
  227, 231, 235, 239, 243, 247, 251, 255};



/* Writes a square image in 8-bit/color PPM format.
 */
void write_image(char *fname, int pixels, double *red, double *green, double *blue)
{
  FILE *fd;

  if((fd = fopen(fname, "w")))
    {
      printf("writing PPM image '%s'.\n", fname);
      fflush(stdout);


      int row, col;

      fprintf(fd, "P6\n%d %d\n%d\n", pixels, pixels, 255);

      for(row = 0; row < pixels; row++)
	for(col = 0; col < pixels; col++)
	  {
	    unsigned char rgb[3];
	    rgb[0] = red[row * pixels + col];
	    rgb[1] = green[row * pixels + col];
	    rgb[2] = blue[row * pixels + col];

	    fwrite(rgb, 3, sizeof(char), fd);
	  }
      fclose(fd);
    }
  else
    {
      printf("file %s can't be opened\n", fname);
      exit(1);
    }
}



struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[96];        /* fills to 256 Bytes */
} header1;


static float *Pos;
int    NumPart;


/* this routine loads particle positions of type=1 from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.)
 */
void load_snapshot_positions(char *fname)
{
  FILE *fd;
  int i, j, k, dummy, ntot_withmasses;
  int t, n, off, pc, pc_new, pc_sph;


  if(!(fd = fopen(fname, "r")))
    {
      fprintf(stderr, "can't open file `%s`\n", fname);
      exit(0);
    }

  printf("reading `%s' ...\n", fname);
  fflush(stdout);

  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  
  NumPart = header1.npart[1];
  if(!(Pos = malloc(NumPart * 3 * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  fread(&dummy, sizeof(dummy), 1, fd);
  fread(Pos, sizeof(float), 3 * NumPart, fd);
  fread(&dummy, sizeof(dummy), 1, fd);

  printf("reading done.\n");
  fflush(stdout);
}




int main(int argc, char **argv)
{
  int i, j, n;

  if(argc != 3)
    {
      printf("call with:  <snapshotfile>  <output_imagefile>  as arguments\n");
      exit(1);
    }

  load_snapshot_positions(argv[1]);


  printf("binning...\n");
  fflush(stdout);

  double *rho =  malloc(PIXELS * PIXELS * sizeof(double));

  /* set density field to zero */
  for(i = 0; i < PIXELS; i++)
    for(j = 0; j < PIXELS; j++)
      rho[i * PIXELS + j] = 0; 

  /* grid spacing */
  double h = 250.0 / PIXELS;

  for(n = 0; n < NumPart; n++)
    {
      double x = Pos[n*3  + 0];
      double y = Pos[n*3  + 1];
      double z = Pos[n*3  + 2];

      if(z < 250.0 / 10)
	{
	  /* calculate principal bin */
	  double xx = x / h;
	  int i = (int) xx;
	  double yy = y / h;
	  int j = (int) yy;

	  /* weights for assignment */
	  double u = xx - i;
	  double v = yy - j;

	  /* determine neighbouring bins */
	  int ii = i + 1;
	  if(ii >= PIXELS)
	    ii = 0;
	  int jj = j + 1;
	  if(jj >= PIXELS)
	    jj = 0;

	  /* now do the CIC assignment of a unit mass */
	  rho[i * PIXELS + j] += (1 - u) * (1 - v) / (h * h);
	  rho[ii * PIXELS + j] += (u) * (1 - v) / (h * h);
	  rho[i * PIXELS + jj] += (1 - u) * (v) / (h * h);
	  rho[ii * PIXELS + jj] += (u) * (v) / (h * h);
	}
    }


  /* allocate some storage for the image, and then read it */

  double *red = malloc(PIXELS * PIXELS * sizeof(double));
  double *green = malloc(PIXELS * PIXELS * sizeof(double));
  double *blue = malloc(PIXELS * PIXELS * sizeof(double));

  double ma = 2000.0;
  double mi = 2.0;

  for(n = 0; n < PIXELS * PIXELS; n++)
    {
      double val = rho[n];
      int idx = 0;

      if(val > 0)
	{
	  double v = log(val/mi) / log(ma/mi) * 256.0;

	  if(v < 0)
	    idx = 0;
	  else if (v < 255)
	    idx = v;
	  else
	    idx = 255;
	}

      red[n] = red_table[idx];
      green[n] = green_table[idx];
      blue[n] = blue_table[idx];
    }

  write_image(argv[2], PIXELS, red, green, blue);

  exit(0);
}
