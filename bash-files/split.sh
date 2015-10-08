#!/usr/bin/env bash

## This script splits a given OGLE data file into smaller files for parallel
## processing on the NeSI cluster.

#------------------------------------------------------------------------------
# Arguments: <vartype>  --specify basename (no suffix) of variable type when
#                         calling script
#------------------------------------------------------------------------------

CORES=16  # Default value for running on SandyBridge architecture is 16 cores; change as appropriate

for i in anocepheid cepheid cepheid_II Del_Sct DPV RR_Lyra LPV; do

    VARTYPE=$i
    echo "vartype: $VARTYPE"

    FILENAME=$VARTYPE.csv

    # Read header and write to file
    sed -n 1,7p $FILENAME > $VARTYPE-header.txt

    # Read data and write to file
    sed -n '8,$ p' $FILENAME > $VARTYPE-data.txt

    # Calculate no. of lines per file
    N_LINES=$[($(wc -l < $VARTYPE-data.txt)+15)/$CORES]
    echo "lines/file: $N_LINES"

    # Split data file and write to temporary files
    cat $VARTYPE-data.txt | awk -v lines=$N_LINES -v fmt="$VARTYPE-temp-%d.txt" '{print>sprintf(fmt,1+int((NR-1)/lines))}'

    # Concatenate header onto each of the split files
    for i in {1..16}; do
        cat $VARTYPE-header.txt $VARTYPE-temp-$i.txt > $VARTYPE-$i.txt
    done

    # Remove temporary files
    rm $VARTYPE-temp-*.txt

done
