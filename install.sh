#!/bin/bash

NUM_CORES=1

while [[ $# -gt 0 ]];
do
key="$1"

case $key in
    -n)
    NUM_CORES="$2"
    shift # past argument
    ;;
    -m)
    NUM_CORES="0"
    shift # past argument
    ;;
    *)
        # unknown option
    ;;
esac
shift # past argument or value
done

cd src

make clean && make 

cd ..
if [ ! -d "data" ]; then
    mkdir "data"
    mkdir "data/jmean"
    mkdir "data/im"
    mkdir "data/deposit"
fi

if [ ! -d "build" ]; then
   mkdir "build"
fi
cd build
ndirec="$(pwd)"
cd ..
if [ ! -d "bin" ]; then
   mkdir "bin"
fi
cd bin
bdirc="$(pwd)"
cd ..
cd src


for i in *; do
   if [ "${i}" != "${i%.mod}" ];then
      cp "${i}" "$ndirec"
   fi
   if [ "${i}" != "${i%.o}" ];then
      mv "${i}" "$ndirec"
   fi
done


if [ "$NUM_CORES" = "0" ]; then #just make code
    exit 0
fi

mv mcgrid "$bdirc" && echo " "&& echo "*****Install complete*****" && echo " "

clear
cd ../bin

if [ "$NUM_CORES" = "1" ]; then
    ./mcgrid
else
    mpirun -n "$NUM_CORES" ./mcgrid
fi