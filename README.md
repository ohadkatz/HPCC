# HPCC
Extended HPC Challenge Benchmark Tests which allow user inputted range and parallel Marix Multiplication tools.

## Installing Interfaces
make sure intel/impi/VERSIONNUMBER/bin64 is added to your PATH

make sure intel/compilers_and_libraries_VERSION/linux/compiler/lib/intel64 added to LD_LIBRARY_PATH

To check these work run which mpicc, should see path to bin64

Install interface to your mkl/interfaces path(GCC & mpicc needed) using the following line:

        make libintel64 PRECISION=MKL_DOUBLE interface=lp64 compiler=gnu


            
## How To Run Tests    
Fork the Repo or Clone and navigate to HPCC/hpl folder

edit the Make.intel64_iccSP:

    1.) Change MPI path to where your intel path is
    2.) Change MKL path to where your intel path is

go to HPCC directory

run: 
  
    make arch=intel64_iccSP
    
edit hpccinf.txt to your preferred range values


To run benchmark you have a couple options (Make sure to check your maximum allowed processes in your cpu info):
    
    mpicc -np (NUMBER OF PROCESSES) ./hpcc
    
    OMP_NUM_THREADS=(NUMBER OF THREADS) mpicc -np (NUMBER OF PROCESSES) ./hpcc
    
    ./hpcc (Default)
