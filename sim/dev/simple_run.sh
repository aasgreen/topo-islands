#!/bin/bash
mkdir dataFolder/
cp simulations/* dataFolder/
cd dataFolder
gfortran -O3 -o defect.o cb.f90
gfortran -O3 -o defectT.o aveDefect.f90

#mkdir data
#arguments g, beta, N, endT 
python param.py #uncomment to have python generate list of parameters
file="./params.txt"
decross='0'
runName=$(echo $(date +%Y%m%d_%H%M%S))
mkdir run"$runName"
cd run"$runName"
mkdir movies
cp ../params.txt ./
while read -r params
do
    k=`echo $params| awk '{print $1}'`
    beta=`echo $params| awk '{print $2}'`
    mu=`echo $params| awk '{print $3}'`
    N=`echo $params| awk '{print $4}'`
    seed=`echo $params| awk '{print $6}'`
 #   mkdir ./k-"$k"
 #   cd ./k-"$k"
 #   mkdir ./beta-"$beta"
 #   cd beta-"$beta"
 #   mkdir ./mu-"$mu"
 #   cd ./mu-"$mu"
    echo $params > param.txt
    dName=`echo "./data-k-$k-beta-$beta-mu-$mu-seed-$seed"` 
    mkdir $dName
    #mkdir ./g-"$g"-beta-"$beta"-mu-"$mu"-data 
    cd $dName

    cp ../../defect.o ./
#1    cp ../../ePlot.gnu ./
#1    cp ../../imgGen.py ./
    ./defect.o $params | tee -a meanE.dat
    echo 'start video encoding program'
    cp ../../dtrack.sh ./
    cp ../../defectT.o ./
#1    python ./imgGen.py $decross
    #./dtrack.sh $N >> output.txt
   # mv ./output.txt ../k-"$k"_b-"$beta"_mu-"$mu"-defectTracks.txt
    cd ../
    pwd
#1    cp ../defectMovie.py ./
#1    python defectMovie.py $dName
#1    cp defect.mp4 ./movies/k-"$k"_b-"$beta"_mu-"$mu"-seed-"$seed"defect.mp4
#1    cp dT.mp4 ./movies/k-"$k"_b-"$beta"_mu-"$mu"-tracking.mp4
#1    mv defect.mp4 k-"$k"_b-"$beta"_mu-"$mu"-seed-"$seed"defect.mp4
done <$file
#./defect.o 10, .4, 64, 1001
#mv *.dat data/

