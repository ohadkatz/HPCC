#!/bin/bash


if [[ -e  hpccoutf.txt ]] ;
then
    printf "removing Output File\n"
    rm hpccoutf.txt
fi

if [[ -e hpccoutTP.txt ]] ;
then 
    
    printf "removing Unnecessary Files\n"
    rm hpccoutTP.txt hpccoutS.txt hpccoutEP.txt
fi


make arch=intel64_icc
sed -n 37,38p hpccinf.txt 

touch hpccoutTP.txt hpccoutS.txt hpccoutEP.txt
cd RFiles
touch RinputTP.txt RinputEP.txt RinputS.txt
cd ..

printf  "TOTALLY PARALLEL RUNNING\n\n"
OMP_NUM_THREADS=4 mpirun -np 1 ./hpcc
cp -r hpccoutf.txt hpccoutTP.txt
sed -n  84,200p hpccoutTP.txt
cp -r RFiles/Rinput.txt RFiles/RinputTP.txt
 

rm hpccoutf.txt

printf "SERIAL RUNNING\n"
OMP_NUM_THREADS=1 mpirun -np 1 ./hpcc
cp -r hpccoutf.txt hpccoutS.txt
sed -n 84,200p hpccoutS.txt
cp -r RFiles/Rinput.txt RFiles/RinputS.txt

#SEG FAULT ISSUE, fix coming soon#
# echo EMBARASIINGLY PARALLEL RUNNING
# OMP_NUM_THREADS=1 mpirun -np 4 ./hpcc
# cp -r hpccoutf.txt hpccoutEP.txt
# cp -r RFiles/Rinput.txt RFiles/RinputEP.txt

printf "--ALL--DONE--"
sleep 5 


