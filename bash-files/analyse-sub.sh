#!/usr/bin/env bash

# This is a short script which reads a given class file and writes each
# variable into a separate file, allowing more in-depth analysis of the results.i

filename=$1

for i in "1:cepheid" "2:anocephe" "3:BLHer" "4:RVTau" "5:WVir" "6:pWvir" "7:RRab" "8:RRc" "9:RRd" "10:RRe" "11:OSARG" "12:Mira" "13:SRV" "14:DPV" "15:Del_Sct" "16:EB"; do

    search_string=$i

    IFS=':' read -a split_string <<< "$search_string"
    echo ${split_string[0]}
    echo ${split_string[1]} 
    
    export vartype=${split_string[1]}
    export search_string=$i

    echo "Searching for ${search_string} in ${filename}"

    awk '$3==ENVIRON["search_string"]{print}' ${filename} > ./sub/${vartype}.dat

done
