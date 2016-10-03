#!/bin/bash

#############################################################################
#
# USAGE 
#       bash create_input_command_batcher
#              
# DESCRIPTION
#        This bash script is used to generate the input file which will be
#		 passed to Batcher.
#		 The user needs to modify the Executable path to PROPKA (the EXECUTABLE
#		 variable, the path to the working directory (DIRECTORY variable) and the  
#		 last PDB file number in the for loop (here 39999).
#
#############################################################################


# Name of the input command file to be generated
INPUT_COMMAND_FILE="batcher_input_command"
shift
# Get the executable name 
EXECUTABLE=/hpc/home/ejl62/propka-3.1-master/propka31
shift

DIRECTORY=/hpc/home/ejl62/nme/nme_2/analysis_apo_1_2/pka
shift

# Create a command line per input parameter "inputi"
for var in {0..39999..1};
do
    echo "$EXECUTABLE $DIRECTORY/pdbfiles/New_$var.pdb > $var.log " >> $INPUT_COMMAND_FILE
done





