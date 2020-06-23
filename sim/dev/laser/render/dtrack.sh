#!/bin/bash
N=$1
for filename in ./defect*.dat; do
    #echo $filename
    time=$(echo $filename |grep -o "\([0-9]*\)")
    #echo $time'   '$N'    '$filename
    ./defectT.o $filename $N $time 
done


