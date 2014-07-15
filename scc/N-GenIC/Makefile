EXEC   = N-GenIC 

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile

#OPT   +=  -DDOUBLEPRECISION  # Set this to calulate and save the particles in double precision

#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type

#OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)

#OPT   +=  -DNEUTRINOS  # this will produce a second component as slight neutrinos (needs to be in initial glass)
#OPT   +=  -DNEUTRINO_PAIRS  # this will produce an additional partner for every neutrino with opposite thermal velocities


#SYSTYPE="Curie"
#SYSTYPE="Ranger"
#SYSTYPE="MBM"
#SYSTYPE="Magny"
#SYSTYPE="odyssey"
#SYSTYPE="Stella"
#SYSTYPE="VIP"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="Solaris"

#SYSTYPE="all-intel"
#SYSTYPE="archer"
SYSTYPE="mvapich-pgi"
#SYSTYPE="mvapich-intel"
#SYSTYPE="openmpi-intel"

FFTW_INCL = -I/shared/apps/gadget/libraries/fftw2/include/
FFTW_LIBS = -L/shared/apps/gadget/libraries/fftw2/lib/

CC       =   mpicc        # sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)
MPICHLIB =   -lmpich

ifeq ($(SYSTYPE),"all-intel")

$(info ****************** Using intel mpi config ********************)

CC       =   mpiicc
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/shared/apps/gadget/libraries/gsl-intel/include/ 
GSL_LIBS =  -L/shared/apps/gadget/libraries/gsl-intel/lib/
#FFTW_INCL=  -I/home/intel/cluster-studio/mkl/include/fftw/
#FFTW_LIBS=  -L/home/intel/cluster-studio/mkl/lib/intel64/ -lmkl_blacs_lp64 -lmkl_sequential -lmkl_core -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm
FFTW_INCL=  -I/shared/apps/gadget/libraries/fftw2-intelmpi/include/ #-openmp -I/home/intel/cluster-studio/mkl/include/
FFTW_LIBS=  -L/shared/apps/gadget/libraries/fftw2-intelmpi/lib/
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 
endif

ifeq ($(SYSTYPE),"archer")

$(info ****************** Using archer config ********************)

CC       =   cc
OPTIMIZE =  -O3
GSL_INCL =  -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS =  -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL=  -I/shared/apps/gadget/libraries/fftw2-mvapichPGI/include/
FFTW_LIBS=  -L/shared/apps/gadget/libraries/fftw2-mvapichPGI/lib/
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5
endif



ifeq ($(SYSTYPE),"mvapich-pgi")

$(info ****************** Using mvapich-pgi config ********************)

CC       =   mpicc
OPTIMIZE =  -O3 -g -m64
GSL_INCL =  -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS =  -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL=  -I/shared/apps/gadget/libraries/fftw2-mvapichPGI/include/
FFTW_LIBS=  -L/shared/apps/gadget/libraries/fftw2-mvapichPGI/lib/
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5
endif

ifeq ($(SYSTYPE),"mvapich-intel")

$(info ****************** Using mvapich-intel config ********************)

CC 	= mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS =  -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL=  -I/shared/apps/gadget/libraries/fftw2-mvapichIntel/include/
FFTW_LIBS=  -L/shared/apps/gadget/libraries/fftw2-mvapichIntel/lib/
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"openmpi-intel")

$(info ****************** Using openmpi-intel config ********************)

CC      = mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/shared/apps/gadget/libraries/gsl-intel/include/
GSL_LIBS =  -L/shared/apps/gadget/libraries/gsl-intel/lib/
FFTW_INCL=  -I/shared/apps/gadget/libraries/fftw2-openmpiIntel/include/
FFTW_LIBS=  -L/shared/apps/gadget/libraries/fftw2-openmpiIntel/lib/
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"OPA-Cluster64-Intel")  # module load mvapich2-1.2-sdr-intel/11.0
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib  -Xlinker -R -Xlinker /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"OPA-Cluster64-Gnu")  # module load mvapich2-1.2-sdr-gnu/4.1.2
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib  -Xlinker -R -Xlinker /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"VIP")
CC       =  mpcc_r -g -qcpluscmt # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(OPT)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE n\
ot needed)
endif
OPTIMIZE =  -q64 # -O3 -qarch=auto -qtune=auto -qinline # -qhot -qipa=level=2 -qfloat=hdflt
GSL_INCL = -I/afs/rzg/u/kud/include
GSL_LIBS = -L/afs/rzg/u/kud/lib -bstackpsize:64k -bdatapsize:64k -btextpsize:64k
FFTW_INCL= -I/afs/rzg/u/kud/include
FFTW_LIBS= -L/afs/rzg/u/kud/lib  -q64 # -qipa
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT     +=  # -DPOWER6 -DALT_QSORT -DPAD_STRUCTURES
endif


ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =   mpicci   
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/lib                -static
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_linux/include
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_linux/lib
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =   /afs/ipp-garching.mpg.de/u/vrs/bin/mpiccf   
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/lib
FFTW_INCL= -I/afs/ipp-garching.mpg.de/u/vrs/fftw_linux_gcc3.2/include
FFTW_LIBS= -L/afs/ipp-garching.mpg.de/u/vrs/fftw_linux_gcc3.2/lib  
endif

ifeq ($(SYSTYPE),"Solaris")
CC       =   mpcc   # sets the C-compiler
OPTIMIZE =   -i -fast -xvector -xarch=v9b -xchip=ultra3 -xcache=64/32/4:8192/512/1 -I/opt/local/include

GSL_INCL = -I/opt/local/include/gsl
GSL_LIBS = -L/opt/local/lib/sparcv9               
FFTW_INCL= -I/opt/local/include
FFTW_LIBS= -L/opt/local/lib/sparcv9
endif

OPTIONS =  $(OPT)


ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
endif



LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

ifeq ($(SYSTYPE),"Solaris")
LIBS   =   -R/opt/local/lib/sparcv9 -lm  -lmpi   $(GSL_LIBS) -lgsl -lgslcblas  $(FFTW_LIB)
endif



CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



