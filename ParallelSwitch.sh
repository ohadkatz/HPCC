#!/bin/bash

make arch=intel64_icc
if [[ -e  hpccoutf.txt ]] ;
then
    echo removing Output File
    rm hpccoutf.txt
fi

if [[ -e hpccoutTP.txt ]] ;
then 
    echo 
    echo removing Unnecessary Files
    rm hpccoutTP.txt hpccoutS.txt hpccoutEP.txt
fi

sed -n 37,38p hpccinf.txt 

touch hpccoutTP.txt hpccoutS.txt hpccoutEP.txt
cd RFiles
touch RinputTP.txt RinputEP.txt RinputS.txt
cd ..

echo  TOTALLY PARALLEL RUNNING
OMP_NUM_THREADS=4 mpirun -np 1 ./hpcc

cp -r hpccoutf.txt hpccoutTP.txt
cp -r RFiles/Rinput.txt RFiles/RinputTP.txt

rm hpccoutf.txt
echo 

echo SERIAL RUNNING
OMP_NUM_THREADS=1 mpirun -np 1 ./hpcc
cp -r hpccoutf.txt hpccoutS.txt
cp -r RFiles/Rinput.txt RFiles/RinputS.txt
echo 


#SEG FAULT ISSUE, fix coming soon#
# echo EMBARASIINGLY PARALLEL RUNNING
# OMP_NUM_THREADS=1 mpirun -np 4 ./hpcc
# cp -r hpccoutf.txt hpccoutEP.txt
# cp -r RFiles/Rinput.txt RFiles/RinputEP.txt


echo ALL DONE
sleep 5 


