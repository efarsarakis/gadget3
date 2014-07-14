#!/bin/bash
#MSUB -@ volker.springel@h-its.org:begin,end
#MSUB -r GB        # Request name
#MSUB -n 4096     # Number of tasks to use
#MSUB -T 1800      # max elapsed time limit in seconds
#MSUB -o log-%I.o  # Standard output. %I is the job id
#MSUB -e log-%I.e  #  Error output.    %I is the job id
#MSUB -A pa0918 
#MSUB -q standard

set -x
cd ${BRIDGE_MSUB_PWD}
echo ${BRIDGE_MSUB_PWD}

#ccc_mprun  ./Gadget3-Small   param-small.txt
#ccc_mprun  ./Gadget3-Medium  param-medium.txt
ccc_mprun  ./Gadget3-Large   param-large.txt
