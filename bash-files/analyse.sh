#!/usr/bin/env bash

# This is a short script which reads a given class file and writes each
# variable into a separate file, allowing more in-depth analysis of the results.i

filename=$1

for i in "1:cepheid" "2:cepheid_" "3:anocephe" "4:LPV" "5:DPV" "6:RR_Lyra" "7:Del_Sct" "8:EB"; do

    search_string=$i

    IFS=':' read -a split_string <<< "$search_string"
    echo ${split_string[0]}
    echo ${split_string[1]} 
    
    export vartype=${split_string[1]}
    export search_string=$i

    echo "Searching for ${search_string} in ${filename}"

    awk '$3==ENVIRON["search_string"]{print}' ${filename} > ./gen/${vartype}.dat

done
