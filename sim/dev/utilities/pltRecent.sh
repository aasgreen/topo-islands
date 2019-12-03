#!/bin/bash

dName=`ls -d dataFolder/run* | sort -n -r |head -1`
echo $dName
cd $dName/data*
#bRE='beta-([0-9].?[0-9]*)'
#name=`pwd`
#if [[ $name =~ $bRE ]]
#then
#    beta="${BASH_REMATCH[1]}"
#    echo "$beta"
#    cmd=`awk '{print $2/$3}'<<< " plot 1.0 $beta"`
#    sed -i "3i\"meanE.dat\" u 1:($cmd) " ePlot.gnu
#    more ePlot.gnu
#fi
gnuplot ePlot.gnu

