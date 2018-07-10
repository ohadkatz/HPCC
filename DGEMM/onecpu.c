/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */

#include <hpcc.h>

int
HPCC_StarDGEMM(HPCC_Params *params) {
  int myRank, commSize;
  double localGflops, minGflops, maxGflops, avgGflops;
  int n;
  int rv, errCount, failure, failureAll;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  localGflops = minGflops = maxGflops = avgGflops = 0.0;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  rv = HPCC_TestDGEMM( params, 0 == myRank, &localGflops, &n, &failure, 0, 1 );
  //params->DGEMM_N = n;

  MPI_Reduce( &rv, &errCount, 1, MPI_INT, MPI_SUM, 0, comm );
  MPI_Allreduce( &failure, &failureAll, 1, MPI_INT, MPI_MAX, comm );
  if (failureAll) params->Failure = 1;

  MPI_Reduce( &localGflops, &minGflops, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &localGflops, &avgGflops, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &localGflops, &maxGflops, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
  avgGflops /= (double)commSize;

  MPI_Bcast( &avgGflops, 1, MPI_DOUBLE, 0, comm ); params->StarDGEMMGflops = avgGflops;

  BEGIN_IO( myRank, params->outFname, outputFile);
  fprintf( outputFile, "Node(s) with error %d\n", errCount );
  fprintf( outputFile, "Minimum Gflop/s %.6f\n", minGflops );
  fprintf( outputFile, "Average Gflop/s %.6f\n", avgGflops );
  fprintf( outputFile, "Maximum Gflop/s %.6f\n", maxGflops );
  END_IO( myRank, outputFile );

  return 0;
}

int
HPCC_ParallelDGEMM(HPCC_Params *params) {
  int myRank, commSize;
  double localGflops, minGflops, maxGflops, avgGflops;
  int n;
  int rv, errCount, failure, failureAll;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  localGflops = minGflops = maxGflops = avgGflops = 0.0;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  rv = HPCC_TestDGEMM( params, 0 == myRank, &localGflops, &n, &failure, 0, 0 );
  //params->DGEMM_N = n;

  MPI_Reduce( &rv, &errCount, 1, MPI_INT, MPI_SUM, 0, comm );
  MPI_Allreduce( &failure, &failureAll, 1, MPI_INT, MPI_MAX, comm );
  if (failureAll) params->Failure = 1;

  MPI_Reduce( &localGflops, &minGflops, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &localGflops, &avgGflops, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &localGflops, &maxGflops, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
  avgGflops /= (double)commSize;

  MPI_Bcast( &avgGflops, 1, MPI_DOUBLE, 0, comm ); params->StarDGEMMGflops = avgGflops;

  BEGIN_IO( myRank, params->outFname, outputFile);
  fprintf( outputFile, "Node(s) with error %d\n", errCount );
  fprintf( outputFile, "Minimum Gflop/s %.6f\n", minGflops );
  fprintf( outputFile, "Average Gflop/s %.6f\n", avgGflops );
  fprintf( outputFile, "Maximum Gflop/s %.6f\n", maxGflops );
  END_IO( myRank, outputFile );

  return 0;
}
