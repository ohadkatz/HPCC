#!/usr/bin/env bash

OutputDir=~/HPCC/outputs

printf '\e[5t'
cat ".DGEMMtesting.txt"
Runs=(2 4 8)

function InputReset(){
    sed -i '11s/.*/1/' hpccinf.txt
    sed -i '11s/$/            Ps/' hpccinf.txt
    sed -i "12s/.*/1/" hpccinf.txt
    sed -i '12s/$/            Qs/' hpccinf.txt
    > Results.txt
    > FinalResults.txt
    
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
Start="0"
End="0"
InputReset
printf "===================== PART 1 =========================\n"
#=======================Single Node: Serial=================================#
printf "Serial Run on a Single Node! \n"
Start=$(wc -l < Results.txt)
time OMP_NUM_THREADS=1  mpirun -np 1 ./hpcc
End=$(wc -l < Results.txt)
python IOchange.py 1 1 "Serial" $Start $End
echo "==============================SERIAL============================" >> hpccoutf.txt
printf "All set with serial! \n"

#====================Single Node: True Parallel=============================#

printf "True Parallel Run on Single Node! \n"
for Thread in "${Runs[@]}"; do 
    Start=$(wc -l < Results.txt)
    time OMP_NUM_THREADS=$Thread  mpirun -np 1 ./hpcc
    End=$(wc -l < Results.txt)
    python IOchange.py 1 $Thread "True Parallel" $Start $End 
    echo "==========TRUE=======Parallel=========" >> hpccoutf.txt
done
printf "Done with True Parallel Run \n"
printf "True Parallel Run (scaLAPACK) on Single Node! \n"

for P in "${Runs[@]}"; do 
    Switch $P
    Start=$(wc -l < Results.txt)
    time OMP_NUM_THREADS=1  mpirun -np $P ./hpcc
    End=$(wc -l < Results.txt)
    python IOchange.py $P 1 "True Parallel(SCALAPACK)" $Start $End 
    echo "==========True=======Parallel=========SCALAPACK=========" >> hpccoutf.txt
done
printf "All set with True Parallel (scaLAPACK) Run\n"

#================Single Node: Embarassingly Parallel========================#

printf "Embarassingly Parallel Run on Single Node!\n"
for Q in "${Runs[@]}"; do 
    Switch $Q
    Start=$(wc -l < Results.txt)
    time OMP_NUM_THREADS=1 mpirun -np $Q ./hpcc
    End=$(wc -l < Results.txt)
    python IOchange.py $Q 1 "Embarassingly Parallel" $Start $End 
    echo "==========Embarassingly=======Parallel=========" >> hpccoutf.txt
done
printf "Embarassingly Parallel All Done!\n"
printf "All tests complete! Goodbye.\n"
printf "To find outputs please go to directory marked outputs! Each portion will have its own Directory \n"

cp FinalResults.txt $OutputDir
