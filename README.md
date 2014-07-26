gadget3
=======

=============================================================
=============================================================
Gadget3 OpenACC version
=============================================================

##6/3/2014 01.53
Successfully compiled and ran GADGET3. Downloaded and manually compiled and installed FFTW2 for both single and double precision. GSL was installed from yum. 
FFTW2 config:

	$ ./configure --enable-type-prefix --enable-mpi
	$ make
	$ sudo make install
	$ make clean
	$ ./configure --enable-float --enable-type-prefix --enable-mpi
	$ make
	$ sudo make install

N-GenIC compilation:

1. Change SYSTYPE
2. $ make

N-GenIC running:

	$ mpirun -np 8 ./N-GenIC ics_small.param
	$ mpirun -f mfile4nodes -np 100 -ppn 25 -bind-to socket -map-by hwthread ./N-GenIC ics_medium.param 

Note that folders ../ICs and ../output must exist.

GADGET3 compilation:

Systypes have been created which in the Makefile sets the appropriate library paths. 

	$ make CONFIG=Config-Small.sh  EXEC=Gadget3-small-double

GADGET3 running:

	$ mpirun -f mfile5nodes -np 100 -ppn 20 -bind-to socket -map-by hwthread ./Gadget3-small-double param-small.txt
	
	or...?
	
	$ mpirun -mca btl openib.self --allow-run-as-root -x LD_LIBRARY_PATH -hostfile mfile -np 32 --map-by node ./GadgetSmall param-small.txt
	
	
	
##Update

30/05/2014:

Gadget has been ported to the final cluster configuration. It shows good speedup, with the small test case performing at ~13sec on 16 cores / 1 node, ~7sec on 32 cores / 2 nodes, ~5sec on 48 cores / 3 nodes.


