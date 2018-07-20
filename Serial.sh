#!/usr/bin/env bash
#SBATCH -J hpccSerial
#SBATCH -o SerialResult.o  
#SBATCH -e SerialErr.e  
#SBATCH -p skx-normal  
#SBATCH -N 1   
#SBATCH -n 96  
#SBATCH -t 24:00:00   
#SBATCH --mail-user=ohadkatz@buffalo.edu  
#SBATCH --mail-type=all  
module list 
pwd
date


HPCCdir=~/HPCC
OutputDir=~/HPCC/SingleNode/outputs
cat $HPCCdir/.DGEMMtesting.txt
Runs=(2 4 8 16 24 48 96)

function InputReset(){
    sed -i '11s/.*/1/' $HPCCdir/hpccinf.txt
    sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
    sed -i "12s/.*/1/" $HPCCdir/hpccinf.txt
    sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    > $HPCCdir/Results.txt
    > $HPCCdir/FinalResults.txt
    > $HPCCdir/hpccoutf.txt
}

function Switch(){
    echo $1
    if [ $1 == 2 ]; then
        sed -i "11s/.*/1/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/2/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi
   
    if [ $1 == 4 ]; then
        sed -i "11s/.*/2/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/2/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi

    if [ $1 == 8 ]; then
        sed -i "11s/.*/2/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/4/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi
   
    if [ $1 == 16 ]; then
        sed -i "11s/.*/4/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/4/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi

    if [ $1 == 24 ]; then
        sed -i "11s/.*/4/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/6/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi

    if [ $1 == 48 ]; then
        sed -i "11s/.*/6/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/8/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi

    if [ $1 == 96 ]; then
        sed -i "11s/.*/8/" $HPCCdir/hpccinf.txt
        sed -i '11s/$/            Ps/' $HPCCdir/hpccinf.txt
        sed -i "12s/.*/12/" $HPCCdir/hpccinf.txt
        sed -i '12s/$/            Qs/' $HPCCdir/hpccinf.txt
    fi

}

Start="0"
End="0"
InputReset
printf "===================== PART 1 =========================\n"
#=======================Single Node: Serial=================================#
printf "Serial Run on a Single Node! \n"
Starit=$(wc -l < $HPCCdir/Results.txt)
export IBRUN_TASKS_PER_NODE=1 
$HPCCdir/hpcc
End=$(wc -l < $HPCCdir/Results.txt)
cp -r $HPCCdir/hpccoutf.txt $OutputDir/
python $HPCCdir/IOchange.py 1 1 "Serial" $Start $End
echo "=============================SERIAL============================" >> $HPCCdir/hpccoutf.txt
printf "All set with serial! \n"

#====================Single Node: True Parallel=============================#
printf "True Parallel Run on Single Node! \n"
for Thread in "${Runs[@]}"; do 
   Switch $Thread
   Start=$(wc -l < $HPCCdir/Results.txt)
   export OMP_NUM_THREADS=48
   ./hpcc 
   End=$(wc -l < $HPCCdir/Results.txt)
   cp -r $HPCCdir/hpccoutf.txt $OutputDir/
   python $HPCCdir/IOchange.py 1 $Thread "True Parallel" $Start $End 
   echo "==========TRUE=======Parallel=========" >> $HPCCdir/hpccoutf.txt
done
printf "Done with True Parallel Run \n"
printf "True Parallel Run (scaLAPACK) on Single Node! \n"

for P in "${Runs[@]}"; do 
   Switch $P
   Start=$(wc -l < $HPCCdir/Results.txt)
   export OMP_NUM_THREADS=1
   export IBRUN_TASKS_PER_NODE=$P
   time ibrun $HPCCdir/hpcc
   End=$(wc -l < $HPCCdir/Results.txt)
   cp -r $HPCCdir/hpccoutf.txt $OutputDir/
   python $HPCCdir/IOchange.py $P 1 "True Parallel(SCALAPACK)" $Start $End 
   echo "==========True=======Parallel=========SCALAPACK=========" >> $HPCCdir/hpccoutf.txt
done
printf "All set with True Parallel (scaLAPACK) Run\n"

#================Single Node: Embarassingly Parallel========================#
printf "Embarassingly Parallel Run on Single Node!\n"
for Q in "${Runs[@]}"; do 
   Switch $Q
   Start=$(wc -l < $HPCCdir/Results.txt)
   export OMP_NUM_THREADS=1
   export IBRUN_TASKS_PER_NODE=$Q
   time ibrun $HPCCdir/hpcc
   End=$(wc -l < $HPCCdir/Results.txt)
   cp -r $HPCCdir/hpccoutf.txt $OutputDir/
   python $HPCCdir/IOchange.py $Q 1 "Embarassingly Parallel" $Start $End 
   echo "==========Embarassingly=======Parallel=========" >> $HPCCdir/hpccoutf.txt
done
printf "Embarassingly Parallel All Done!\n"
printf "All Single Node tests complete!\n"

cp -r $HPCCdir/FinalResults.txt $OutputDir/
