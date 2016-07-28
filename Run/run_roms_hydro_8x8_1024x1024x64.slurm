#!/bin/bash
#MSUB -r Croco_nhmg_seamount    # Request na me
#MSUB -n 64               # Number of tasks to use
#MSUB -T 180            # Elapsed time limit in seconds
#MSUB -o nh_%I.txt       # Standard output.%I is the job id
#MSUB -e nh_%I.txt       # Error output.%I is the job id
#MSUB -q standard
#MSUB -A gen9178         # IMPORTANT to specifuy here the correct project ID

# Project ID

set -x

module unload intel

module load intel/16.0.2.181 bullxmpi netcdf/4.2

########################################
###         Number of processes         ##
### HAVE TO BE EQUAL TO -PBS -l select- ##
###   at the begining of this script    ##
##########################################
export NB_PROC=64

####################################
##  Define files and run parameters #
#####################################
export SRC_DIR=/ccc/cont003/home/gen9178/griman/CROCO-NH/Run
export CODFILE=./roms_hydro
export MODEL=roms.in.seamount
export RUNEXP=Seamount_hydro
export EXEC_DIR=/ccc/work/cont003/gen9178/griman/CROCO-NH

###########################
## Move to the Scratchdir #
###########################
cd $EXEC_DIR/$RUNEXP

###################################
## Get the executable code (roms) #
###################################
echo "Getting $CODFILE from $SRC_DIR"
cp -f $SRC_DIR/$CODFILE $EXEC_DIR/$RUNEXP/.
### to be an unix executable ##
chmod 755 $CODFILE

######################################
## Get the input arguments (roms.in) #
######################################
echo "Getting $MODEL from $SRC_DIR"
cp -f $SRC_DIR/$MODEL $EXEC_DIR/$RUNEXP/.

######################################################
############  HERE ROMS IS STARTING TO RUN ###########
############            ON CAPARMOR        ###########
######################################################
ls -rtl
 date
   time ccc_mprun -n ${NB_PROC} ${CODFILE}  ${MODEL} > ${MODEL}.out
 date
ls -rtl
