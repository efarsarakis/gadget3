#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=Gadget_new
# 
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype 
#     This file has priority over your shell variable.:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype 
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype

ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif

CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat gadgetconfig.h)


CC       = mpicc        # sets the C-compiler (default)
CXX       = mpiCC       # sets the C++-compiler (default)

FC 	 = mpif90

OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =   mpiCC     # default C++ compiler
OPTIMIZE =   -g 
OPT     +=  -DX86FIX   # only needed for 32-bit intel/amd systems
endif

ifeq (VORONOI,$(findstring VORONOI,$(CONFIGVARS)))
GMPLIB   =  -lgmp
endif

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES =  -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
#  ifeq ($(SYSTYPE),"all-intel")
#    FFTW_LIBNAMES = -lfftw3x_cdft -lfftw3xc_intel -lfftw2x_cdft_DOUBLE -lfftw2xc_double_intel
#  else
    FFTW_LIBNAMES =  -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
#  endif
else
#  ifeq ($(SYSTYPE),"all-intel")
#    FFTW_LIBNAMES = -lfftw3x_cdft -lfftw3xc_intel -lfftw2x_cdft_SINGLE -lfftw2xc_single_intel
#  else
    FFTW_LIBNAMES =  -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
#  endif
endif
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"all-intel")

$(info ****************** Using all-intel config ********************)

CC       = mpiicc
CXX      = mpiicpc
FC       = mpiifort
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =   mpiicpc
endif
GSL_INCL = -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS = -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL= -I/shared/apps/gadget/libraries/fftw2-intelmpi/include/ #-I/home/intel/cluster-studio/mkl/include/ -I/home/intel/cluster-studio/mkl/include/fftw
FFTW_LIBS= -L/shared/apps/gadget/libraries/fftw2-intelmpi/lib/
#FFTW_LIBS=  -L/home/intel/cluster-studio/mkl/lib/intel64/ -lmkl_blacs_lp64 -lmkl_sequential -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64
MPICHLIB =
OPTIMIZE = -O3 -g -Wall #-openmp -lpthread #-liomp5 -mt_mpi -DFFTW_ENABLE_FLOAT -DMKL_SINGLE
endif

ifeq ($(SYSTYPE),"mvapich-pgi")
GSL_INCL = -I/shared/apps/gadget/libraries/gsl-gcc/include/
GSL_LIBS = -L/shared/apps/gadget/libraries/gsl-gcc/lib/
FFTW_INCL= -I/shared/apps/gadget/libraries/fftw2-mvapichPGI/include/
FFTW_LIBS= -L/shared/apps/gadget/libraries/fftw2-mvapichPGI/lib/
MPICHLIB =

OPTIMIZE = -fastsse
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline -Mvect=sse -Munroll=c:1 -Mcache_align
#OPTIMIZE = -O4 -Mipa=fast,inline -Munroll
#OPTIMIZE = -O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mscalarsse -Mcache_align -Mflushz -Mpre
#OPTIMIZE = -fastsse -Mipa=fast,inline
#OPTIMIZE = -fastsse -Mipa=fast,inline -Mvect=sse,short
#OPTIMIZE = -fastsse -Mipa=fast,inline -Munroll
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline,libc,libopt,libinline -Mvect=simd -Mcache_align -Munroll
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline -Mvect=sse -Munroll=c:1 -Mcache_align -g
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline,libc,libopt,libinline -Mvect=simd -Mcache_align -Munroll -g
endif

ifeq ($(SYSTYPE),"pgprof")
GSL_INCL = -I/shared/apps/gadget/libraries/gsl-gcc/include/
GSL_LIBS = -L/shared/apps/gadget/libraries/gsl-gcc/lib/
FFTW_INCL= -I/shared/apps/gadget/libraries/fftw2-mvapichPGI/include/
FFTW_LIBS= -L/shared/apps/gadget/libraries/fftw2-mvapichPGI/lib/
MPICHLIB =

OPTIMIZE = -fastsse -pg
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline -Mvect=sse -Munroll=c:1 -Mcache_align
#OPTIMIZE = -O4 -Mipa=fast,inline -Munroll
#OPTIMIZE = -O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mscalarsse -Mcache_align -Mflushz -Mpre
#OPTIMIZE = -fastsse -Mipa=fast,inline
#OPTIMIZE = -fastsse -Mipa=fast,inline -Mvect=sse,short
#OPTIMIZE = -fastsse -Mipa=fast,inline -Munroll
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline,libc,libopt,libinline -Mvect=simd -Mcache_align -Munroll
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline -Mvect=sse -Munroll=c:1 -Mcache_align -g
#OPTIMIZE = -fastsse -O4 -Mipa=fast,inline,libc,libopt,libinline -Mvect=simd -Mcache_align -Munroll -g
endif



ifeq ($(SYSTYPE),"acc-boston")
GSL_INCL = -I/shared/apps/gadget/libraries/gsl-gcc/include/
GSL_LIBS = -L/shared/apps/gadget/libraries/gsl-gcc/lib/
FFTW_INCL= -I/shared/apps/gadget/libraries/fftw2-mvapichPGI/include/
FFTW_LIBS= -L/shared/apps/gadget/libraries/fftw2-mvapichPGI/lib/
MPICHLIB =

OPTIMIZE = -fastsse -ta=nvidia,kepler -acc -Minfo=accel -Mipa=inline -Msafeptr

endif

ifeq ($(SYSTYPE),"archer")
CC       = cc
CXX      = cpc
FC       = cc


OPTIMIZE = -g
endif

ifeq ($(SYSTYPE),"swan-cray")
CC       = cc
CXX      = cpc
FC       = cc
GSL_INCL = -I/lus/scratch/p02045/lib/mygsl/include/
GSL_LIBS = -L/lus/scratch/p02045/lib/mygsl/lib/

OPTIMIZE = -O1
endif

ifeq ($(SYSTYPE),"swan-pgi")
CC       = cc
CXX      = cpc
FC       = cc
GSL_INCL = -I/lus/scratch/p02045/lib/mygsl/include/
GSL_LIBS = -L/lus/scratch/p02045/lib/mygsl/lib/

OPTIMIZE = -O3
endif

ifeq ($(SYSTYPE),"acc-swan")
CC       = cc
CXX      = cpc
FC       = cc
GSL_INCL = -I/lus/scratch/p02045/lib/mygsl/include/
GSL_LIBS = -L/lus/scratch/p02045/lib/mygsl/lib/

OPTIMIZE = -g -target-accel=nvidia35
endif



ifeq ($(SYSTYPE),"mvapich-intel")
GSL_INCL = -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS = -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL= -I/shared/apps/gadget/libraries/fftw2-mvapichIntel/include/
FFTW_LIBS= -L/shared/apps/gadget/libraries/fftw2-mvapichIntel/lib/
MPICHLIB = 
OPTIMIZE = -O3 -g -Wall 
endif

ifeq ($(SYSTYPE),"openmpi-intel")
GSL_INCL = -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS = -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL= -I/shared/apps/gadget/libraries/fftw2-openmpiIntel/include/
FFTW_LIBS= -L/shared/apps/gadget/libraries/fftw2-openmpiIntel/lib/
MPICHLIB = 
OPTIMIZE = -O3 -g -Wall 
endif

ifeq ($(SYSTYPE),"Ranger_intel")
CC       =  mpicc -g -O2 #-xW -ipo -Wall
CXX      =  mpiCC -g -O2 -xW -ipo -Wall
OPTIMIZE =
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"Curie")  #  module load gsl  fftw2/2.1.5  hdf5/1.8.8
CC       =  mpicc -g -O2  -Wno-unknown-pragmas
OPTIMIZE =
GSL_INCL = -I$(GSL_INC_DIR)
GSL_LIBS = -L$(GSL_LIB_DIR)
FFTW_INCL= -I$(FFTW2_INC_DIR)
FFTW_LIBS= -L$(FFTW2_LIB_DIR)
HDF5INCL = -I$(HDF5_INC_DIR) -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_LIB_DIR) -lhdf5 -lz
MPICHLIB =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"Ranger_pgi")
CC       =  mpicc -g -O2 #-tp barcelona-64
CXX      =  mpiCC -g -O2 -tp barcelona-64
OPTIMIZE =
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
endif



ifeq ($(SYSTYPE),"Octopus")
CC       =  mpicc
CXX      =  mpiCC
OPTIMIZE = -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS= -lfftw3
MPICHLIB = -lmpi
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc
CXX      =  mpiCC 
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC  
endif
ifeq (GDE_BIGFLOAT,$(findstring GDE_BIGFLOAT,$(CONFIGVARS)))
CC       =  mpiCC
endif
#gnu
OPT      +=  -DMPICH_IGNORE_CXX_SEEK  
OPTIMIZE =   -O3 -g -Wall -m64
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -I/n/sw/hdf5-1.8.5_gnu/include
HDF5LIB  =  -L/n/sw/hdf5-1.8.5_gnu/lib -lhdf5
#intel
#OPT      +=  -DMPICH_IGNORE_CXX_SEEK  -DNOTYPEPREFIX_FFTW
#OPTIMIZE =   -parallel -O3 -ipo -funroll-loops -Wall
#GSL_INCL =
#GSL_LIBS =
#FFTW_INCL=
#FFTW_LIBS=
#MPICHLIB =
#HDF5INCL =
#HDF5LIB  =  -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"Darwin")
CC       =  mpicc -g -O2 -Wall -Wno-unknown-pragmas
CXX      =  mpiCC -g -O2 -Wall
OPTIMIZE =  
GSL_INCL = -I/sw/include -I/opt/local/include 
GSL_LIBS = -L/sw/lib     -L/opt/local/lib
FFTW_INCL= -I/sw/include -I/opt/local/include
FFTW_LIBS= -L/sw/lib     -L/opt/local/lib
MPICHLIB = -lmpi 
HDF5INCL = -I/sw/lib     -I/opt/local/include -DH5_USE_16_API 
HDF5LIB  = -L/sw/lib -L/opt/local/lib  -lhdf5 -lz 
endif


# modules for Genius
# module load hdf5/1.8.4p1

ifeq ($(SYSTYPE),"Genius")
CC       =  mpixlc_r
CXX      =  mpixlcxx_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450
GSL_INCL =  -I/u/vrs/BlueGene/Libs/include
GSL_LIBS =  -L/u/vrs/BlueGene/Libs/lib -static
FFTW_INCL=  -I/u/vrs/BlueGene/Libs/include
FFTW_LIBS=  -L/u/vrs/BlueGene/Libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"BGP_toy")
CC       =  mpixlc_r
CXX      =  mpixlcxx_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450
GSL_INCL =  -I/gpfs/scratch/userinternal/cin8263a/lib/gsl-1.9/include/
GSL_LIBS =  -L/gpfs/scratch/userinternal/cin8263a/lib/gsl-1.9/lib/ -static
FFTW_INCL=  -I/gpfs/scratch/userinternal/cin8263a/lib/fftw-2.1.5/include
FFTW_LIBS=  -L/gpfs/scratch/userinternal/cin8263a/lib/fftw-2.1.5/lib
MPICHLIB =
HDF5INCL =  -I/gpfs/scratch/userinternal/cin0753a/hdf5_2/hdf5/hdf5-1.8.5-patch1/hdf5/include -I/gpfs/scratch/userinternal/cin0753a/zlib/zlib-1.2.5/include
HDF5LIB  =  -L/gpfs/scratch/userinternal/cin0753a/hdf5_2/hdf5/hdf5-1.8.5-patch1/hdf5/lib -static -lhdf5 -L/gpfs/scratch/userinternal/cin0753a/zlib/zlib-1.2.5/lib -lz
endif

# modules for Jugene
# module load gsl/1.14
# module load fftw/v2.1.5
# module load hdf5/1.8.4_450

ifeq ($(SYSTYPE),"Jugene")
CC       =  mpixlc_r
CXX      =  mpixlcxx_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450 -qsmp=omp:noauto -qreport
GSL_INCL =  -I$(GSL_DIR)/include
GSL_LIBS =  -L$(GSL_DIR)/lib
FFTW_INCL=  -I$(FFTW_DIR)/include
FFTW_LIBS=  -L$(FFTW_DIR)/lib
MPICHLIB =
HDF5INCL =  $(HDF5_INCLUDE)
HDF5LIB  =  -L/bgsys/local/hdf5/v1.8.4/lib_450 -lhdf5 -lhdf5_hl -L/bgsys/local/zlib/lib/ -lz -L/bgsys/local/szip/lib/ -lsz # $(HDF5_LIB)
ZLIB     =  -L/bgsys/local/zlib/lib/ -lz
SZIP     =  -L/bgsys/local/szip/lib/ -lsz
OPT      +=  -DNOCALLSOFSYSTEM
endif

# modules for Fermi
# module load bgq-xl/1.0
# module load gsl/1.15--bgq-xl--1.0
# module load hdf5/1.8.9_ser--bgq-xl--1.0
# module load zlib/1.2.7--bgq-gnu--4.4.6
# module load szip/2.1--bgq-xl--1.0

ifeq ($(SYSTYPE),"Fermi")
OPT      +=  -DOLD_HDF5
CC       =  mpixlc
CXX      =  mpixlcxx
OPTIMIZE =  -O2 -qreport -qsmp=omp:noauto -qthreaded
GSL_INCL =  -I$(GSL_INCLUDE)
GSL_LIBS =  -L$(GSL_LIB)            
FFTW_INCL=  -I/gpfs/scratch/userexternal/gmurante/local/include
FFTW_LIBS=  -L/gpfs/scratch/userexternal/gmurante/local/lib
#FFTW_INCL=  -I/gpfs/scratch/userexternal/kdolag00/include
#FFTW_LIBS=  -L/gpfs/scratch/userexternal/kdolag00/lib
#FFTW_LIBNAMES = -lfftw3f_mpi -lfftw3f_threads -lfftw3f_omp -lfftw3f             
MPICHLIB =           
HDF5INCL =  -I$(HDF5_HOME)/include
HDF5LIB  =  -L$(HDF5_HOME)/lib -lhdf5 -lhdf5_hl -L$(ZLIB_HOME)/lib/ -lz -L$(SZIP_HOME)/lib/ -lsz # $(HDF5_LIB)
ZLIB     =  -L$(ZLIB_HOME)/lib/ -lz
SZIP     =  -L$(SZIP_HOME)/lib/ -lsz                    
OPT      +=  -DNOCALLSOFSYSTEM
endif            

ifeq ($(SYSTYPE),"Solaris_amd64")
CC       = /opt/SUNWhpc/bin/mpicc
CXX      = /opt/SUNWhpc/bin/mpiCC
OPTIMIZE = -xO4 -xchip=native -xarch=amd64a
#OPTIMIZE =   -g -xarch=amd64a
GSL_INCL = $(GSL64INCL)
GSL_LIBS = $(GSL64LIB)
FFTW_INCL= $(FFTW64INCL)
FFTW_LIBS= $(FFTW64LIB)
FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
R_PATH   = $(R64LIB)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL)  $(GSL_LIBS)
HDF5INCL =
HDF5LIB  =
LIBS     =   -lm  $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) $(R64LIB)
endif


ifeq ($(SYSTYPE),"MBM") 
CC       =   mpicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/sw/tap/include
GSL_LIBS =  -L/sw/tap/lib  -Xlinker -R -Xlinker /sw/tap/lib
FFTW_INCL=  -I/sw/tap/include
FFTW_LIBS=  -L/sw/tap/lib  -Xlinker -R -Xlinker /sw/tap/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif

# modules for Magny
# module add mvapich2/gcc/64/1.4.1-qlc

ifeq ($(SYSTYPE),"Magny") 
CC       =   mpicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/hits/tap/sw/libs/include
GSL_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
FFTW_INCL=  -I/hits/tap/sw/libs/include
FFTW_LIBS=  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
GMP_INCL =  -I/hits/tap/sw/libs/include
GMP_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
MPICHLIB =
HDF5INCL =  -I/hits/tap/sw/libs/include -DH5_USE_16_API
HDF5LIB  =  -L/hits/tap/sw/libs/lib -lhdf5 
OPT      +=  -DNOCALLSOFSYSTEM
OPT      +=  -DIMPOSE_PINNING
endif

ifeq ($(SYSTYPE),"OpenSuse")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=
MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"OpenSuse64")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB = -L/usr/lib64/mpi/gcc/openmpi/lib64 -Xlinker -R -Xlinker /usr/lib64/mpi/gcc/openmpi/lib64 -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"SUNNYVALE")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpicxx
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib64
FFTW_INCL= -I${MCKENZIE_FFTW_INC_PATH}
FFTW_LIBS= -L${MCKENZIE_FFTW_LIB_PATH}
MPICHLIB =  -L/opt/lam-7.1.2-intel/lib -lmpi
HDF5INCL = -I${MCKENZIE_HDF_INC_PATH}
HDF5LIB  = -L${MCKENZIE_HDF_LIB_PATH} -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"McKenzie")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpicxx
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib
FFTW_INCL= -I/opt/fftw/intel_8.1/2.1.5/include
FFTW_LIBS= -L/opt/fftw/intel_8.1/2.1.5/lib #-ldrfftw_mpi
MPICHLIB = -L/opt/lam-7.1.2b24-g77/lib -lmpi
HDF5INCL = -I/opt/hdf5-oscar-1.6.4/include
HDF5LIB  = -L/opt/hdf5-oscar-1.6.4/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 


ifeq ($(SYSTYPE),"Stella")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/home/schaye/libs/include
GSL_LIBS =  -L/home/schaye/libs/lib -static
FFTW_INCL=  -I/home/schaye/libs/include
FFTW_LIBS=  -L/home/schaye/libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif


# modules for OPA  
# module load intel
# module load impi
# module load fftw/2.1.5
# module load gsl/1.14
# module load hdf5-serial

ifeq ($(SYSTYPE),"OPA-Cluster64-Intel")  
CC       =   mpiicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64 -wd981 -wd2259 -wd1572 -ipo4 -funroll-loops -finline -finline-functions -finline-limit=200
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR) -Xlinker -R -Xlinker $(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib -Xlinker -R -Xlinker $(FFTW_HOME)/lib
MPICHLIB =  
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB  = -L$(HDF5_HOME)/lib -Xlinker -R -Xlinker $(HDF5_HOME)/lib -lhdf5 -lz
#OPT      +=  -DOLD_HDF5
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"opus")
CC       =   mpicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/snfs1/cgwork/dgad0043/include
GSL_LIBS =  -L/snfs1/cgwork/dgad0043/lib
FFTW_INCL=  -I/snfs1/cgwork/dgad0043/include
FFTW_LIBS=  -L/snfs1/cgwork/dgad0043/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"Sauron-gcc")
CC       =   mpicc.gcc   # sets the C-compiler
OPTIMIZE =   -O3 -funroll-loops -march=k8 -msse2 -static
GSL_INCL =   -I/usr/local/gsl.gcc/include
GSL_LIBS =   -L/usr/local/gsl.gcc/lib -static -lgsl -lgslcblas
FFTW_INCL=   -I/usr/local/fftw.gcc/include
FFTW_LIBS=   -L/usr/local/fftw.gcc/lib -static -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Sauron")
CC       =  mpicc  -m64 # sets the C-compiler
CXX      =  mpiCC  -m64
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC  -m64
endif
OPTIMIZE =   -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


ifeq ($(SYSTYPE),"OpteronMPA-Intel")
CC       =  /usr/local/appl/openmpi-intel/bin/mpicc  -m64 # sets the C-compiler
CXX      =  /usr/local/appl/openmpi-intel/bin/mpiCC  -m64
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  /usr/local/appl/openmpi-intel/bin/mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"OpteronMPA-Gnu")
CC       =  /usr/local/appl/openmpi-gnu/bin/mpicc  -m64 # sets the C-compiler
CXX      =  /usr/local/appl/openmpi-gnu/bin/mpiCC  -m64
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  /usr/local/appl/openmpi-gnu/bin/mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpiCC
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -g -Wall -fopenmp 
# GSL_INCL =  -I/usr/common/pdsoft/include
# GSL_LIBS =  -L/usr/common/pdsoft/lib
GSL_INCL =  -I/afs/mpa/home/volker/Libs/include
GSL_LIBS =  -L/afs/mpa/home/volker/Libs/lib
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz 
OPT     +=  -DOLD_HDF5
endif


ifeq ($(SYSTYPE),"Mako")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpiCC
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -O3 -march=athlon-mp  -mfpmath=sse
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


# modules for VIP
# module load fftw/2.1.5
# module load gsl/1.14
# module load hdf5

ifeq ($(SYSTYPE),"VIP")
CC       =  mpcc_r -g -qcpluscmt # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
CCX      =  mpCC_r -g
FC       =  mpxlf_r
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPT      += -DVIP
FOPT     += -WF,-DVIP
OPTIMIZE =  -q64 -O2 -qarch=auto -qtune=auto -qinline -qstrict -qipa -qhot -qipa=level=2 -qlanglvl=extc99 -qalias=ansi -qinline=200 # -qfloat=hdflt
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
OPT     +=   -DPOWER6 -DDO_NOT_BRACH_IF -DPAD_STRUCTURES
endif

# modules for SARA
# module load c/ibm/11.1
# module load fftw2/2.1.5
# module load gsl/1.11
# module load hdf5

ifeq ($(SYSTYPE),"SARA")
CC       =  mpcc -g -qcpluscmt # -qsmp -qreport # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
CXX      =  mpCC -g
FC       =  mpxlf_r
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(CONFIGVARS)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPT      += -DVIP
FOPT     += -WF,-DVIP
OPTIMIZE =  -q64 -O5 -qarch=auto -qtune=auto -qinline -qstrict -qhot -qipa=level=2 -qlanglvl=extc99 -qalias=ansi -qinline=200 # -qfloat=hsflt
GSL_INCL = -I$(SARA_GSL_INCLUDE)
GSL_LIBS = -L$(SARA_GSL_LIB)
FFTW_INCL= -I$(SARA_FFTW2_INCLUDE)
FFTW_LIBS= -L$(SARA_FFTW2_LIB)
MPICHLIB =
HDF5INCL = -I$(SARA_HDF5_INCLUDE)
HDF5LIB  = -L$(SARA_HDF5_LIB) -lhdf5 -lz -static
OPT     +=   -DPOWER6 -DDO_NOT_BRACH_IF -DPAD_STRUCTURES
endif


ifeq ($(SYSTYPE),"CINECA_SP6")
CC       =   mpcc # sets the gcc C-compiler
CXX      =   mpCC
#to be used for deugging: -g -qcheck=all -qdbxextra -qheapdebug -qkeepparm -qoptdebug
ifndef DEBUG
OPTIMIZE =  -O5 -q64  -qhot -qstrict -qarch=pwr6 -qtune=pwr6 -qcpluscmt -qlanglvl=extc99 -qalias=ansi -qinline=200 -qipa=level=2 # :inline=noauto:list
else
OPTIMIZE = $(DEBUG) -qarch=pwr6 -qtune=pwr6 -qcpluscmt
endif

GSL_INCL = -I/cineca/prod/libraries/gsl/1.9/xl--10.1/include
GSL_LIBS = -L/cineca/prod/libraries/gsl/1.9/xl--10.1/lib
FFTW_INCL= -I/cineca/prod/libraries/fftw/2.1.5/xl--10.1/include
FFTW_LIBS= -L/cineca/prod/libraries/fftw/2.1.5/xl--10.1/lib
HDF5INCL = -I/sp6/userexternal/sborgani/code/HDF5/include
HDF5LIB  = -L/sp6/userexternal/sborgani/code/HDF5/lib -lhdf5
MPICHLIB =
OPT     +=   -DPOWER6 -DDO_NOT_BRACH_IF -DPAD_STRUCTURES
endif


ifeq ($(SYSTYPE),"Ingeld_LinuxCluster")
CC       =   mpicc  # sets the C-compiler
CXX      =   mpiCC
OPTIMIZE =   -O3 -Wall

GSL_INCL = -I/home/patricia/my-gsl/include
GSL_LIBS = -L/home/patricia/my-gsl/lib -static
FFTW_INCL= -I/home/patricia/my-fftw/include
FFTW_LIBS= -L/home/patricia/my-fftw/lib
endif


ifeq ($(SYSTYPE),"hpcf")
CC       =   mpicc  # sets the C-compiler
CXX      =   mpiCC
OPT     +=  -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPTIMIZE =  -O3
GSL_INCL = -I/home/gadget/Libs/include
GSL_LIBS = -L/home/gadget/Libs/lib
FFTW_INCL= -I/home/gadget/Libs/include
FFTW_LIBS= -L/home/gadget/Libs/lib
endif


ifeq ($(SYSTYPE),"Warp")
CC       =   /opt/intel/impi/3.1/bin64/mpiicc
CXX      =   /opt/intel/impi/3.1/bin64/mpiiCC
OPTIMIZE =  -xW
GSL_INCL = -I/opt/gsl/impi/include/gsl
GSL_LIBS = -L/opt/gsl/impi/lib64 -Wl,-rpath,/opt/gsl/impi/lib64
FFTW_INCL= -I/opt/fftw/impi_intel-3.1/include
FFTW_LIBS= -L/opt/fftw/impi_intel-3.1/lib64 -Wl,-rpath,/opt/fftw/impi_intel-3.1/lib64
MPICHLIB = -L/opt/intel/impi/3.1/lib64 -Wl,-rpath,/opt/intel/impi/3.1/lib64
endif


ifeq ($(SYSTYPE),"Kraken")
OPT     +=   -DNOCALLSOFSYSTEM
CC       =   cc # sets the C-compiler
OPTIMIZE =   -fast
GSL_INCL =  -I/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/include
GSL_LIBS =   -L/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/lib
FFTW_INCL=  -I/opt/fftw/2.1.5/cnos/include
FFTW_LIBS=  -L/opt/fftw/2.1.5/cnos/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
endif

ifeq ($(SYSTYPE),"Kraken_FFTW3")
OPT     +=   -DNOCALLSOFSYSTEM
CC       =   cc # sets the C-compiler
OPTIMIZE =   -fast
GSL_INCL =  -I/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/include
GSL_LIBS =   -L/sw/xt5/gsl/1.11/cnl2.1_pgi7.2.3/lib
FFTW_INCL=  -I/sw/xt/fftw/3.3_alpha/cnl2.2_pgi9.0.3/include
FFTW_LIBS=  -L/sw/xt/fftw/3.3_alpha/cnl2.2_pgi9.0.3/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))    # fftw installed with type prefix?
  FFTW_LIBNAMES = $(FFTW_LIBS) -lfftw3_mpi -lfftw3_threads -lfftw3
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))
  FFTW_LIBNAMES =  -lfftw3_mpi -lfftw3_threads -lfftw3
else
  FFTW_LIBNAMES =  -lfftw3f_mpi -lfftw3f_threads -lfftw3f
endif
endif

endif

ifeq ($(SYSTYPE),"INAF_IRA")
OPT     += -DNOCALLSOFSYSTEM
CC       = mpicc # sets the C-compiler
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/homes/donnert/Libs/include
GSL_LIBS = -L/homes/donnert/Libs/lib
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB =
HDF5INCL =
HDF5LIB  = -lhdf5 
endif


ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
endif


ifeq (CHEMCOOL,$(findstring CHEMCOOL,$(CONFIGVARS)))
FOPT	+=  -WF,-DCHEMCOOL
OPT	+=  -DCHEMISTRYNETWORK=1
FOPT	+=  -WF,-DCHEMISTRYNETWORK=1
endif


ifneq (MODGRAV,$(findstring MODGRAV,$(CONFIGVARS))) # standard gravity object files
GRAVITY_OBJS  =  forcetree.o darkenergy.o pm_periodic.o accel.o
else # replacements used in the modified gravity model 
GRAVITY_OBJS  =  modgrav/forcetree.o modgrav/darkenergy.o modgrav/pm_periodic.o modgrav/accel.o
endif

L3_OBJS = 

ifeq (SORT_FROM_L3,$(findstring SORT_FROM_L3,$(CONFIGVARS))) # add L3 parallel sorting
L3_OBJS  +=  ImportFromL3/parallel_sort_l3.o
endif

ifeq (PM_FROM_L3,$(findstring PM_FROM_L3,$(CONFIGVARS))) # add L3 pm part
L3_OBJS  +=  ImportFromL3/pm_periodic_l3.o ImportFromL3/rfftwnd_mpi_threads.o ImportFromL3/transpose_mpi_threads.o
endif


OPTIONS = $(OPTIMIZE) $(OPT) 

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = GadgetMedium

OBJS  = endrun.o mpi_util.o restart.o subfind_loctree.o fof.o  mymalloc.o  run.o  subfind_nearesttwo.o \
 allocate.o ngb.o second_order.o subfind_potential.o allvars.o forcetree_update.o  parallel_sort.o sfr_eff.o subfind_reshuffle.o \
 begrun.o  global.o parallel_sort_special.o sfr_mhm.o subfind_serial.o checksummed_sendrecv.o gravtree.o peano.o sizelimited_sendrecv.o subfind_so.o \
 compare_partitions.o gravtree_forcetest.o pinning.o subfind.o subfind_vars.o compile_time_info.o hydra.o pm_nonperiodic.o subfind_alternative_collective.o swap.o \
 cooling.o init.o subfind_collective.o system.o io.o potential.o subfind_cont.o timestep.o \
 density.o kicks.o predict.o subfind_density.o twopoint.o domain.o longrange.o read_ic.o subfind_distribute.o utilities.o \
 driftfac.o main.o read_ic_cluster.o subfind_findlinkngb.o


OBJS	+= $(GRAVITY_OBJS)

OBJS	+= $(L3_OBJS)

INCL    += allvars.h      cooling.h      domain.h       forcetree.h    kernel.h       proto.h        swap.h         utilities.h \
	   assert.h       dd.h           fof.h           myqsort.h      subfind.h      tags.h         vector.h  Makefile

ifeq (SIDM,$(findstring SIDM,$(CONFIGVARS)))
OBJS    +=  sidm.o
INCL    +=
endif

ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
OBJS	+=  utilities.o integrate.o network_solver.o network.o 
INCL	+=  utilities.h integrate.h network_solver.h network.h
endif

ifeq (MODGRAV,$(findstring MODGRAV,$(CONFIGVARS)))
OBJS	+=  modgrav/modgrav_fieldsolve.o modgrav/modgrav_background.o modgrav/modgrav_cic.o modgrav/modgrav_effective_mass.o 
INCL	+=  modgrav/modgrav_cic.h modgrav/modgrav_effective_mass.h 
endif

ifeq (WINDTUNNEL,$(findstring WINDTUNNEL,$(CONFIGVARS)))
OBJS	+=  windtunnel/tunnel.o
INCL	+=  
endif

ifeq (CS_MODEL,$(findstring CS_MODEL,$(CONFIGVARS)))
OBJS  += cs_model/cs_metals.o cs_model/cs_sfr.o cs_model/cs_enrichment.o cs_model/cs_findhotngbs.o cs_model/cs_update_weights.o cs_model/cs_ngbfind.o
INCL  += cs_model/cs_metals.h
endif

ifeq (VS_TURB,$(findstring VS_TURB,$(CONFIGVARS)))
OBJS	+=  turb/turb_driving.o turb/turb_powerspectra.o
INCL	+=  
endif

ifeq (AB_TURB,$(findstring AB_TURB,$(CONFIGVARS)))
OBJS	+= turb/ab_turb.o  turb/turb_driving.o turb/turb_powerspectra.o
endif

ifeq (TURB_DRIVING,$(findstring TURB_DRIVING,$(CONFIGVARS)))
OBJS	+= fb_driveturb/turb_driving.o
endif

ifeq (LT_METAL_COOLING,$(findstring LT_METAL_COOLING,$(CONFIGVARS)))
OBJS  += lt_sfr/lt_wal_cooling.o lt_sfr/lt_sn.o lt_sfr/lt_io.o lt_sfr/lt_utils.o ImfLib/libimf.o ImfLib/libimf_vars.o 
INCL  += lt_sfr/lt_sn.h lt_sfr/lt.h lt_sfr/lt_error_codes.h lt_sfr/lt_wal_cooling.h 
ifeq (LT_SEv_INFO,$(findstring LT_SEv_INFO,$(CONFIGVARS)))
OBJS  += lt_sfr/lt_sfr.o
else
ifeq (LT_ZAGE,$(findstring LT_ZAGE,$(CONFIGVARS)))
OBJS  += lt_sfr/lt_sfr.o
else
OBJS  += lt_sfr/lt_sfr_light.o
endif
endif
endif

ifeq (BP_REAL_CRs,$(findstring BP_REAL_CRs,$(CONFIGVARS))) # add bp cr part
OBJS += bp_cosmic_rays/bp_cosmic_rays.o
INCL += bp_cosmic_rays/bp_cosmic_rays.h
endif

ifeq (UM_CHEMISTRY,$(findstring UM_CHEMISTRY,$(CONFIGVARS)))
OBJS += um_chemistry/um_chemistry_noneq.o
endif

ifeq (UM_METAL_COOLING,$(findstring UM_METAL_COOLING,$(CONFIGVARS)))
OBJS += um_chemistry/um_metal_cooling.o 
endif

ifeq (ADAPTGRAVSOFT,$(findstring ADAPTGRAVSOFT,$(CONFIGVARS)))
INCL += AdaptGravSoft/ads_ngb.o AdaptGravSoft/ads_density.o
endif


ifeq (ADJ_BOX_POWERSPEC,$(findstring ADJ_BOX_POWERSPEC,$(CONFIGVARS)))
OBJS += power_spec/adj_box_powerspec.o 
INCL += power_spec/adj_box_powerspec_proto.h 
endif



ifeq (CHEMCOOL,$(findstring CHEMCOOL,$(CONFIGVARS)))

OBJS  += chemcool.o raytrace.o \

INCL  += chemcool_consts.h \

FOBJS = calc_photo.o calc_temp.o cheminmo.o compute_heating.o compute_md_kappa.o const_rates.o cool_func.o \
	cool_util.o coolinmo.o dvode.o evolve_abundances.o jac.o photoinit_lowZ.o \
	rate_eq_primordial.o spline.o validate_rates.o \

FINCL = cool.h fs_data.h mol_data.h non_eq.h shield_data.h Makefile \

else
FC    = $(CC)
endif

ifeq (EAGLE,$(findstring EAGLE,$(CONFIGVARS)))
OBJS  += eagle/eagle_cooling.o 
INCL  += eagle/eagle_cooling.h 
endif


CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL)

ifeq (VIP,$(findstring VIP,$(CONFIGVARS)))
FFLAGS = $(FOPTIONS)
else
FFLAGS = $(OPTIONS)
endif


ifeq (ALTERNATIVE_PSORT,$(findstring ALTERNATIVE_PSORT,$(CONFIGVARS)))
OBJS  += fof_alt_psort.o psort-1.0/error_handling.o
CXXFLAGS = $(CFLAGS)
FC    = $(CXX)
endif

FFTW = $(FFTW_LIBS)  $(FFTW_LIBNAMES) 


#LIBS   = -lm -openmp $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl $(FFTW) -L/home/intel/cluster-studio/mkl/lib/intel64/ -lmkl_core -lmkl_pgi_thread -lmkl_blacs_intelmpi_lp64
LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW)

ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS))) 
LIBS   +=  -lpthread
endif

ifeq (VORONOI,$(findstring VORONOI,$(CONFIGVARS)))
LIBS   +=  $(GMP_LIBS) $(GMPLIB)
endif

ifeq (GDE_BIGFLOAT,$(findstring GDE_BIGFLOAT,$(CONFIGVARS)))
LIBS   +=  -L/n/home12/mvogelsberger/opt/lib -lhpa
CFLAGS +=  -I/n/home12/mvogelsberger/opt/include 
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c


$(FOBJS): $(FINCL)

complile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) gadgetconfig.h *.oo


