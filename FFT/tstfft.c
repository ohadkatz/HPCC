/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* tstfft.c
 */

#include <hpcc.h>

#include "hpccfft.h"

static int
TestFFT1(HPCC_Params *params, int doIO, FILE *outFile, FILE *Rfile, double *UGflops, int *Un, int *Ufailure, int VecSize, int Repetitions) {
  fftw_complex *in, *out;
  fftw_plan p;
  hpcc_fftw_plan ip;
  double Gflops = -1.0;
  double maxErr, tmp1, tmp2, tmp3, t0, t1, t2, t3;
  int i, n, flags, failure = 1;
  double deps = HPL_dlamch( HPL_MACH_EPS );
  int inner_rep= (VecSize > 2000) ? 1 : Repetitions;
  inner_rep=1;
  //if(VecSize<5000) inner_rep=Repetitions;
#ifdef HPCC_FFT_235
  int f[3];

  /* Need 2 vectors for input and output and 1 vector of scratch spaces */
  n = VecSize;

  /* Adjust local size for factors */
  for ( ; HPCC_factor235( n, f ); n--)
    ; /* EMPTY */
#else
  /* Need 2 vectors and vectors' sizes as power of 2 */
  n = VecSize;
#endif

  /* need to use fftw_malloc() so that the returned pointers will be aligned properly for SSE
     instructions on Intel/AMD systems */
  in  = (fftw_complex *)HPCC_fftw_malloc( (sizeof *in)  * n );
  out = (fftw_complex *)HPCC_fftw_malloc( (sizeof *out) * n );

  if (! in || ! out) goto comp_end;

  /* Make sure that `inout' and `work' are initialized in parallel if using
     Open MP: this will ensure better placement of pages if first-touch policy
     is used by a distrubuted shared memory machine. */
#ifdef _OPENMP
#pragma omp parallel for
  
  for (i = 0; i < n; ++i) {
    c_re( in[i] ) = c_re( out[i] ) = 0.0;
    c_re( in[i] ) = c_im( out[i] ) = 0.0;
  }
#endif
  t0 = -MPI_Wtime();
  for (int i=0 ; i < Repetitions; i++){
    HPCC_bcnrand( 2*n, 0, in );
    if (doIO) fprintf(Rfile,"%s,%d,%d,%f\n","FFT* Gen. Time",i+1,VecSize, t0+MPI_Wtime());
  }
 t0 += MPI_Wtime();

#ifdef HPCC_FFTW_ESTIMATE
  flags = FFTW_ESTIMATE;
#else
  flags = FFTW_MEASURE;
#endif
  t1 = -MPI_Wtime();
  for (int i=0 ; i < Repetitions; i++){
    p = fftw_create_plan( n, FFTW_FORWARD, flags );
    if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","FFT* Tuning",i+1,VecSize, t1+MPI_Wtime());
  }
  t1 += MPI_Wtime();

  if (! p) goto comp_end;

  t2 = -MPI_Wtime();
  for (int i=0 ; i < Repetitions; i++){
    fftw_one( p, in, out );
    if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","FFT* Computing",i+1,VecSize, t2+MPI_Wtime());
    if (t2+MPI_Wtime() > 0.0) Gflops = 1e-9 * (5.0 * n * log(n) / log(2.0)) / (t2+MPI_Wtime());
    if(doIO) fprintf(Rfile,"%s,%d,%d,%f\n","FFT* Gflops",Repetitions,VecSize, Gflops);

  }

  t2 += MPI_Wtime();

  //  fftw_destroy_plan(p);
  ip = HPCC_fftw_create_plan( n, FFTW_BACKWARD, FFTW_ESTIMATE );
  if (ip) {

  t3 = -MPI_Wtime();
  for (int i=0 ; i < Repetitions; i++){
    
    HPCC_fftw_one( ip, out, in );
    
    if (doIO) fprintf(Rfile,"%s,%d,%d,%f\n","FFT* Inverse", i+1, VecSize, t3+MPI_Wtime());
  }
  t3 += MPI_Wtime();
  HPCC_fftw_destroy_plan( ip );
  }

  HPCC_bcnrand( 2*(s64Int)n, 0, out ); /* regenerate data */

  maxErr = 0.0;
  for (i = 0; i < n; i++) {
    tmp1 = c_re( in[i] ) - c_re( out[i] );
    tmp2 = c_im( in[i] ) - c_im( out[i] );
    tmp3 = sqrt( tmp1*tmp1 + tmp2*tmp2 );
    maxErr = maxErr >= tmp3 ? maxErr : tmp3;
  }

  if (maxErr / log(n) / deps < params->test.thrsh) failure = 0;
  if (t2 > 0.0) Gflops = 1e-9 * (5.0 * n * log(n) / log(2.0)) / t2;
 
  if (doIO) {
    fprintf( outFile, "\nVector size: %d\n", n );
    fprintf( outFile, "Repetitions: %d\n", Repetitions );
    fprintf( outFile, "Generation time: %9.3f\n", t0 );
    fprintf( outFile, "Tuning: %9.3f\n", t1 );
    fprintf( outFile, "Computing: %9.3f\n", t2 );
    fprintf( outFile, "Inverse FFT: %9.3f\n", t3 );
    fprintf( outFile, "max(|x-x0|): %9.3e\n", maxErr );
    fprintf( outFile, "Gflop/s: %9.3f\n", Gflops );
  }

  comp_end:

  if (out) HPCC_fftw_free( out );
  if (in)  HPCC_fftw_free( in );

  *UGflops = Gflops;
  *Un = n;
  *Ufailure = failure;

  return 0;
}

int
HPCC_TestFFT(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure) {
  int rv, n, failure = 1;
  double Gflops;
  FILE *outFile = NULL;
  FILE *Rfile = NULL;
  int UserVector=0;
  int Repetitions=0;
  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    Rfile = fopen( params->Results , "a" );
    if (! outFile) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
    if(! Rfile) {
      Rfile=stderr;
      fprintf(Rfile, "Cannot output file.\n");
      return 1;
    }
  } else {
    outFile = fopen( "/dev/null", "w" ); /* special filename Unix file systems */
    if (! outFile) {
      outFile = fopen( "nul", "w"); /* special filename on Windows, produces no output */
    }
  }
  for(int i_vec=0; i_vec<params->FFT_Size;i_vec++){
    n = 0;
    Gflops = -1.0;
    UserVector= params->FFT_UserVector[i_vec];
    Repetitions= params->FFT_repetitions[i_vec];
    
    rv = TestFFT1( params, doIO, outFile, Rfile, &Gflops, &n, &failure, UserVector, Repetitions);
  }
    if (doIO) {
      fflush( outFile );
      fflush( Rfile );
    }

    if (outFile)
      fclose( outFile );
    if (Rfile)
      fclose( Rfile );
    if (UGflops) *UGflops = Gflops;
    if (Un) *Un = n;
    if (Ufailure) *Ufailure = failure;

    return rv;
}
