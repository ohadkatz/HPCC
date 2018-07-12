#!/usr/bin/env bash

OutputDir=~/HPCC/outputs

printf '\e[5t'
cat ".DGEMMtesting.txt"
Runs=(2 4 8)

function Reset(){
    if [ -e hpccoutf.txt ]; then
        > hpccoutf.txt
    fi
    > StarResults.txt
    > ParallelResults.txt 
}

function InputReset(){
    sed -i '11s/.*/1/' hpccinf.txt
    sed -i '11s/$/            Ps/' hpccinf.txt
    sed -i "12s/.*/1/" hpccinf.txt
    sed -i '12s/$/            Qs/' hpccinf.txt
}

function Switch(){
    echo $1
    if [ $1 == 2 ]; then
        sed -i "11s/.*/2/" hpccinf.txt
        sed -i '11s/$/            Ps/' hpccinf.txt
        sed -i "12s/.*/1/" hpccinf.txt
        sed -i '12s/$/            Qs/' hpccinf.txt
    fi
    if [ $1 == 4 ]; then
        sed -i "11s/.*/2/" hpccinf.txt
        sed -i '11s/$/            Ps/' hpccinf.txt
        sed -i "12s/.*/2/" hpccinf.txt
        sed -i '12s/$/            Qs/' hpccinf.txt
    fi
    if [ $1 == 8 ]; then
        sed -i "11s/.*/4/" hpccinf.txt
        sed -i '11s/$/            Ps/' hpccinf.txt
        sed -i "12s/.*/2/" hpccinf.txt
        sed -i '12s/$/            Qs/' hpccinf.txt
    fi
}

InputReset

printf "===================== PART 2 =========================\n"
#=======================Single Node: Serial=================================#
printf "Serial Run on a Single Node! \n"
OMP_NUM_THREADS=1 mpirun -np 1 ./hpcc
#TURN OFF SCALAPACK??
cp hpccoutf.txt $OutputDir/Serial/Serialout.txt
cp StarResults.txt $OutputDir/Serial/SerialResults.txt
Reset
printf "All set with serial! \n"

#====================Single Node: True Parallel=============================#

printf "True Parallel Run on Single Node! \n"
for Thread in "${Runs[@]}"; do 
    OMP_NUM_THREADS=$Thread mpirun -np 1 ./hpcc
    cp hpccoutf.txt $OutputDir/TotParallel/TPout$Thread.txt
    cp StarResults.txt $OutputDir/TotParallel/TPresults$Thread.txt
    Reset
done
printf "Done with True Parallel Run \n"
printf "True Parallel Run (scaLAPACK) on Single Node! \n"

for P in "${Runs[@]}"; do 
    Switch $P
    OMP_NUM_THREADS=1 mpirun -np $P ./hpcc
    cp hpccoutf.txt $OutputDir/TotParallel/ScaLAPACK/ScaLAPACKTPSPout$P.txt
    cp ParallelResults.txt $OutputDir/TotParallel/ScaLAPACK/TPSPresults$P.txt
    Reset
done
printf "All set with True Parallel (scaLAPACK) Run\n"

#================Single Node: Embarassingly Parallel========================#
InputReset
printf "Embarassingly Parallel Run on Single Node!\n"
for Q in "${Runs[@]}"; do 
    Switch $Q
    OMP_NUM_THREADS=1 mpirun -np $Q ./hpcc
    cp hpccoutf.txt $OutputDir/EmbParallel/EPout$Q.txt
    cp StarResults.txt $OutputDir/EmbParallel/EPresults$Q.txt
    Reset
done

InputReset
printf "Embarassingly Parallel All Done!\n"
printf "All tests complete! Goodbye.\n"
printf "To find outputs please go to directory marked outputs! Each portion will have its own Directory \n"

