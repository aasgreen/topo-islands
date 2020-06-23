#!/bin/bash
mkdir dataFolder/
cd dataFolder/
runName=$(echo $(date +%Y%m%d_%H%M%S))
mkdir run"$runName"
cd run"$runName"
mkdir movies
cp -r ../../simulations ./
cp -r ../../render ./
cp -r ../../utilities ./

#cp ../params.txt ./
#cp ../*.f90 ./

#cp simulations/* dataFolder/
#cp render/* dataFolder/
#cp utilities/* dataFolder/

#cd dataFolder

gfortran -O3 -o laserxy.o simulations/laserxy.f90 -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5
gfortran -O3 -o laserxyT.o simulations/aveDefect.f90

#mkdir data
#arguments g, beta, N, endT 
python simulations/param.py #uncomment to have python generate list of parameters
file="./params.txt"
decross='0'
while read -r params
do
    k=`echo $params| awk '{print $1}'`
    beta=`echo $params| awk '{print $2}'`
    mu=`echo $params| awk '{print $3}'`
    N=`echo $params| awk '{print $4}'`
    endT=`echo $params| awk '{print $5}'`
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
    
    cp ../laserxy.o ./
    cp ../utilities/ePlot.gnu ./
    cp ../simulations/imgGen.py ./
    ./laserxy.o $params | tee -a meanE.dat
    echo 'start video encoding program'
    cp ../render/dtrack.sh ./
    cp ../laserxyT.o ./
    cp ../utilities/defectLoc.o ./
    totalNTimes=`echo $(wc datanames.dat | awk '{print $1}')`
    ./defectLoc.o dsetf.h5 datanames.dat $N $totalNTimes
    python ./imgGen.py $decross
    #./dtrack.sh $N >> output.txt
   # mv ./output.txt ../k-"$k"_b-"$beta"_mu-"$mu"-defectTracks.txt
    cd ../
    pwd
    cp render/defectMovie.py ./
    python defectMovie.py $dName
    cp defect.mp4 ./movies/k-"$k"_b-"$beta"_mu-"$mu"-seed-"$seed"defect.mp4
    cp dT.mp4 ./movies/k-"$k"_b-"$beta"_mu-"$mu"-tracking.mp4
    mv defect.mp4 k-"$k"_b-"$beta"_mu-"$mu"-seed-"$seed"defect.mp4
    vlc.exe k-"$k"_b-"$beta"_mu-"$mu"-seed-"$seed"defect.mp4
done <$file

