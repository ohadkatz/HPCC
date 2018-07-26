/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* mpifft.c
 */

#include <hpcc.h>

#include "hpccfft.h"
#include "wrapmpifftw.h"

double *HPCC_fft_timings_forward, *HPCC_fft_timings_backward;

static void
MPIFFT0(HPCC_Params *params, int doIO, FILE *outFile, FILE *Rfile, MPI_Comm comm, int locN,
        double *UGflops, s64Int_t *Un, double *UmaxErr, int *Ufailure, int VecSize, int Repetitions) {
  int commRank, commSize, failure, flags;
  s64Int_t i, n;
  s64Int_t locn, loc0, alocn, aloc0, tls;
  double maxErr, tmp1, tmp2, tmp3, t0, t1, t2, t3, Gflops;
  double deps;
  fftw_complex *inout, *work;
  fftw_mpi_plan p;
  hpcc_fftw_mpi_plan ip;
  int sAbort, rAbort;
#ifdef USING_FFTW
  int ilocn, iloc0, ialocn, ialoc0, itls;
#endif
  
  failure = 1;
  Gflops = -1.0;
  deps = HPL_dlamch( HPL_MACH_EPS );
  maxErr = 1.0 / deps;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &commRank );

  n = locN;
  int inner_rep= (VecSize > 2000) ? 1 : Repetitions;
  inner_rep=1;
  /* number of processes have been factored out - need to put it back in */
  //n *= commSize;

  //n *= commSize; /* global vector size */

#ifdef USING_FFTW
  /* FFTW ver. 2 only supports vector sizes that fit in 'int' */
  if (n > (1<<30)-1+(1<<30)) {
#ifdef HPCC_FFTW_CHECK32
    goto no_plan;
#else
  if (doIO) {
    fprintf( outFile, "Warning: problem size too large: %ld*%d*%d\n", (long)(n / commSize / commSize), commSize, commSize );
  }
#endif
  }
#endif

#ifdef HPCC_FFTW_ESTIMATE
  flags = FFTW_ESTIMATE;
#else
  flags = FFTW_MEASURE;
#endif
  t1 = -MPI_Wtime();
  for(int i = 0; i < Repetitions; i++){

    p = fftw_mpi_create_plan( comm, n, FFTW_FORWARD, flags );

    if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Tuning", i+1 ,VecSize, t1);
  }
  t1 += MPI_Wtime();
  if (! p) goto no_plan;

#ifdef USING_FFTW
  fftw_mpi_local_sizes( p, &ilocn, &iloc0, &ialocn, &ialoc0, &itls );
  locn = ilocn;
  loc0 = iloc0;
  alocn = ialocn;
  aloc0 = ialoc0;
  tls = itls;
#else
  fftw_mpi_local_sizes( p, &locn, &loc0, &alocn, &aloc0, &tls );
#endif

  inout = (fftw_complex *)HPCC_fftw_malloc( tls * (sizeof *inout) );
  work  = (fftw_complex *)HPCC_fftw_malloc( tls * (sizeof *work) );

  sAbort = 0;
  if (! inout || ! work) sAbort = 1;
  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, comm );
  if (rAbort > 0) {
    fftw_mpi_destroy_plan( p );
    goto comp_end;
  }

  /* Make sure that `inout' and `work' are initialized in parallel if using
     Open MP: this will ensure better placement of pages if first-touch policy
     is used by a distrubuted shared memory machine. */
#ifdef _OPENMP
#pragma omp parallel for
  for (i = 0; i < tls; ++i) {
    c_re( inout[i] ) = c_re( work[i] ) = 0.0;
    c_re( inout[i] ) = c_im( work[i] ) = 0.0;
  }
#endif
  t0 = -MPI_Wtime();
  for(int i = 0 ; i < Repetitions ; i++){
    HPCC_bcnrand( 2 * tls, 53 * commRank * 2 * tls, inout );
    
    if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Gen. Time",i+1,VecSize, t0);
  }
  t0 += MPI_Wtime();

  t2 = -MPI_Wtime();
  for(int i = 0 ; i < Repetitions ; i++){
  
 
    fftw_mpi( p, 1, inout, work );
    if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Computing", i+1 ,VecSize, t2);
  }
  t2 += MPI_Wtime();
  //fftw_mpi_destroy_plan( p );

  ip = HPCC_fftw_mpi_create_plan( comm, n, FFTW_BACKWARD, FFTW_ESTIMATE );

  if (ip) {
    t3 = -MPI_Wtime();
    for(int i = 0 ; i < Repetitions ; i++){
      
      HPCC_fftw_mpi( ip, 1, inout, work );
      
      if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Inverse", i+1 ,VecSize, t3);
    }
    t3 += MPI_Wtime();
    HPCC_fftw_mpi_destroy_plan( ip );
  }

  HPCC_bcnrand( 2 * tls, 53 * commRank * 2 * tls, work ); /* regenerate data */

  maxErr = 0.0;
  for (i = 0; i < locn; ++i) {
    tmp1 = c_re( inout[i] ) - c_re( work[i] );
    tmp2 = c_im( inout[i] ) - c_im( work[i] );
    tmp3 = sqrt( tmp1*tmp1 + tmp2*tmp2 );
    maxErr = maxErr >= tmp3 ? maxErr : tmp3;
  }
  MPI_Allreduce( &maxErr, UmaxErr, 1, MPI_DOUBLE, MPI_MAX, comm );
  maxErr = *UmaxErr;
  if (maxErr / log(n) / deps < params->test.thrsh) failure = 0;

  if (t2 > 0.0) Gflops = 1e-9 * (5.0 * n * log(n) / log(2.0)) / t2;
  
  if (doIO) {
    fprintf( outFile, "\nNumber of nodes: %d\n", commSize );
    fprintf( outFile, "Vector size: %9.3d\n", n );
    fprintf( outFile, "Repetitions: %d\n", Repetitions);
    fprintf( outFile, "Generation time: %9.3f\n", t0 );
    fprintf( outFile, "Tuning: %9.3f\n", t1 );
    fprintf( outFile, "Computing: %9.3f\n", t2 );
    fprintf( outFile, "Inverse FFT: %9.3f\n", t3 );
    fprintf( outFile, "max(|x-x0|): %9.3e\n", maxErr );
    fprintf( outFile, "Gflop/s: %9.3f\n", Gflops );
    // fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* GFLOP",Repetitions,VecSize, Gflops);
    // fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Gen. Time",Repetitions,VecSize, t0);
    // fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Tuning",Repetitions,VecSize, t1);
    // fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Computing",Repetitions,VecSize, t2);
    // fprintf(Rfile,"%s,%d,%d,%f\n","MPIFFT* Inverse",Repetitions,VecSize, t3);
    
  }

  comp_end:

  if (work) HPCC_fftw_free( work );
  if (inout) HPCC_fftw_free( inout );

  no_plan:

  *UGflops = Gflops;
  *Un = n;
  *UmaxErr = maxErr;
  *Ufailure = failure;
}

int
HPCC_MPIFFT(HPCC_Params *params) {
  int commRank, commSize;
  int locN, procCnt, isComputing, doIO, failure = 0;
  s64Int_t n;
  double Gflops = -1.0, maxErr = -1.0;
  MPI_Comm comm;
  FILE *outFile;
  FILE *Rfile;
  int Repetitions = 0;
  int VecSize = 0;
  MPI_Comm_size( MPI_COMM_WORLD, &commSize );
  MPI_Comm_rank( MPI_COMM_WORLD, &commRank );

  doIO = commRank == 0 ? 1 : 0;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    Rfile = fopen( params->Results, "a" );
    if (! outFile) outFile = stderr;
    if (! Rfile) Rfile = stderr;
  }

  /*
  There are two vectors of size 'n'/'commSize': inout, work,
  and internal work: 2*'n'/'commSize'; it's 4 vectors then.

  FFTE requires that the global vector size 'n' has to be at least
  as big as square of number of processes. The square is calculated
  in each factor independently. In other words, 'n' has to have
  at least twice as many 2 factors as the process count, twice as many
  3 factors and twice as many 5 factors.
  */
for(int i_vec=0; i_vec<params->FFT_Size;i_vec++){
  VecSize = params->FFT_UserVector[i_vec];
  Repetitions = params->FFT_repetitions[i_vec];
#ifdef HPCC_FFT_235
  locN = 0; procCnt = commSize + 1;
  do {
    int f[3];

    procCnt--;

    for ( ; procCnt > 1 && HPCC_factor235( procCnt, f ); procCnt--)
      ; /* EMPTY */

    /* Make sure the local vector size is greater than 0 */
    locN = params->FFT_UserVector[i_vec];
    for ( ; locN >= 1 && HPCC_factor235( locN, f ); locN--)
      ; /* EMPTY */
  } while (locN < 1);
#else
  /* Find power of two that is smaller or equal to number of processes */
  for (procCnt = 1; procCnt <= (commSize >> 1); procCnt <<= 1)
    ; /* EMPTY */

  /* Make sure the local vector size is greater than 0 */
  while (1) {
    locN = params->FFT_UserVector[i_vec];
    if (locN) break;
    procCnt >>= 1;
  }
#endif

  isComputing = commRank < procCnt ? 1 : 0;

  HPCC_fft_timings_forward = params->MPIFFTtimingsForward;
  HPCC_fft_timings_backward = params->MPIFFTtimingsBackward;

  if (commSize == procCnt)
    comm = MPI_COMM_WORLD;
  else
    MPI_Comm_split( MPI_COMM_WORLD, isComputing ? 0 : MPI_UNDEFINED, commRank, &comm );

  if (isComputing)
    MPIFFT0( params, doIO, outFile, Rfile, comm, locN, &Gflops, &n, &maxErr, &failure, VecSize, Repetitions );

  if (commSize != procCnt && isComputing && comm != MPI_COMM_NULL)
    MPI_Comm_free( &comm );

  params->MPIFFT_N = n;
  params->MPIFFT_Procs = procCnt;
  params->MPIFFT_maxErr = maxErr;

  MPI_Bcast( &Gflops, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  params->MPIFFTGflops = Gflops;

  params->FFTEnblk = FFTE_NBLK;
  params->FFTEnp = FFTE_NP;
  params->FFTEl2size = FFTE_L2SIZE;

  if (failure)
    params->Failure = 1;
}
  if (doIO) if (outFile != stderr) fclose( outFile );
  if (doIO) if (Rfile != stderr) fclose( Rfile );
  return 0;
}
