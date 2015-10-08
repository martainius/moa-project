#!/usr/bin/env bash

## This is a short script to concatenate each of the MOA subfield data files
## into a single arff test file for use with the RF classifier.
##
## [Ignore error message "cat: gb*-R-*.arff: No such file or directory". 
## These are just subfields where the photometry fit was bad and were excluded 
## from the test data.]

# Loop over fields
for i in {1..22}; do
    # Loop over CCD chip no.
    for j in {1..10}; do
        # If first iteration write new file, else append data
        if [ $i = 1 ] && [ $j = 1 ]; then
            echo "First loop"
            cat gb$i-R-$j.arff > data.arff
        else
            cat gb$i-R-$j.arff >> data.arff
        fi
    done
done

# Write headers onto data and save final test files
cat header.arff data.arff > testing.arff
cat header-sub.arff data.arff > testing-sub.arff

echo "Done!"
