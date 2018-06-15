#!/bin/bash


if [[ -e hpccoutf.txt ]];
then
    rm hpccoutf.txt
    printf "removing Output File \n"
    echo hello
fi
if [[ -e hpccoutTP.txt ]];
then   
    rm hpccoutTP.txt hpccoutS.txt hpccoutEP.txt
    printf "removing Unnecessary Files \n"
fi

make arch=intel64_icc
sed -n 37,38p hpccinf.txt 

touch hpccoutTP.txt hpccoutS.txt hpccoutEP.txt
cd RFiles
touch RinputTP8threads.txt RinputEP8threads.txt RinputS8threads.txt
cd ..

printf  "TOTALLY PARALLEL RUNNING\n\n"
OMP_NUM_THREADS=8 mpirun -np 1 ./hpcc
cp -r hpccoutf.txt hpccoutTP.txt
sed -n  84,200p hpccoutTP.txt
cp -r RFiles/Rinput.txt RFiles/RinputTP8threads.txt
 

rm hpccoutf.txt
printf "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n"
printf "\nSERIAL RUNNING\n"
OMP_NUM_THREADS=1 mpirun -np 1 ./hpcc
cp -r hpccoutf.txt hpccoutS.txt
sed -n 84,200p hpccoutS.txt
cp -r RFiles/Rinput.txt RFiles/RinputS8threads.txt

#SEG FAULT ISSUE, fix coming soon#
# echo EMBARASIINGLY PARALLEL RUNNING
# OMP_NUM_THREADS=1 mpirun -np 4 ./hpcc
# cp -r hpccoutf.txt hpccoutEP.txt
# cp -r RFiles/Rinput.txt RFiles/RinputEP.txt
printf "ALL DONE"


