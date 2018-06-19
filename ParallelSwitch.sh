#!/bin/bash

#  
#  -- High Performance Computing Parallel Analysis              
#                               
#     Ohad Katz 
#     SUMMER 2018                                           
#     University at Buffalo                                
#     Center for Computational Research                                                    
#
cat .DGEMMtesting.txt
printf "Beginning Tests \n"

make arch=intel64_icc
#
# PART 1: MULTI NODES
#

# 
# True Parallel
#
# printf("True Parallel: \n")
# OMP_NUM_THREADS= 1 mpirun -np 2 ./hpcc
# cp -r results.txt resultsTPmulti.txt
# cp -r hpccoutf.txt hpccoutfTPmulti.txt

# #
# # Embarassingly Parallel
# #
# printf("Embarassingly Parallel: \n")
# OMP_NUM_THREADS=1 mpirun -np 2 ./hpcc
# cp -r results.txt resultsEPmulti.txt
# cp -r hpccoutf.txt hpccoutEPmulti.txt


# #
# # PART 2 Single
# #

# #
# # Serial
# #
# printf ("Serial: \n")
# OMP_NUM_THREADS=1 mpirun -np 1 ./hpcc
# cp -r results.txt resultsS.txt
# cp -r hpccoutf.txt hpccoutfS.txt

# #
# #True Parallel
# #
# printf("True Parallel(single): \n")
# OMP_NUM_THREADS=4 mpirun -np 1 ./hpcc
# cp -r results.txt resultsTPsingle.txt
# cp -r hpccoutf.txt hpccoutfTPsingle.txt

# #
# #Embarassingly Parallel
# #
# printf("Embarassingly Parallel(single): \n")
# mpirun -np 2 ./hpcc



