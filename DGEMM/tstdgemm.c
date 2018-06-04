/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* tstdgemm.c
 */

#include <hpcc.h>
#include <sys/time.h>
/* Generates random matrix with entries between 0.0 and 1.0 */
static void
dmatgen(int m, int n, double *a, int lda, int seed) {
  int i, j;
  double *a0 = a, rcp = 1.0 / RAND_MAX;

  srand( seed );

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      a0[i] = rcp * rand();

    a0 += lda;
  }
}

static double
dnrm_inf(int m, int n, double *a, int lda) {
  int i, j, k, lnx;
  double mx, *a0;

  int nx = 10;
  double x[10];

  mx = 0.0;

  for (i = 0; i < m; i += nx) {
    lnx = Mmin( nx, m-i );
    for (k = 0; k < lnx; ++k) x[k] = 0.0;

    a0 = a + i;

    for (j = 0; j < n; ++j) {
      for (k = 0; k < lnx; ++k)
        x[k] += fabs( a0[k] );

      a0 += lda;
    }

    for (k = 0; k < lnx; ++k)
      if (mx < x[k]) mx = x[k];
  }

  return mx;
}

HPCC_DGEMM_Calculation(int n, int doIO, double *UGflops, int *Un, int *Ufailure, int *GFLOPS, int iteration){
  int i,j,lda, ldb, ldc, failure = 1;
  double *a=NULL, *b=NULL, *c=NULL, *x=NULL, *y=NULL, *z=NULL, alpha, beta, sres, cnrm, xnrm;
  double Gflops = 0.0, dn, t0, t1;
  long l_n;
  int seed_a, seed_b, seed_c, seed_x;
  if (n < 0) n = -n; /* if 'n' has overflown an integer */
  l_n = n;
  lda = ldb = ldc = n;

  a = HPCC_XMALLOC( double, l_n * l_n );
  b = HPCC_XMALLOC( double, l_n * l_n );
  c = HPCC_XMALLOC( double, l_n * l_n );

  x = HPCC_XMALLOC( double, l_n );
  y = HPCC_XMALLOC( double, l_n );
  z = HPCC_XMALLOC( double, l_n );

  
  seed_a = (int)time( NULL );
  dmatgen( n, n, a, n, seed_a );

  seed_b = (int)time( NULL );
  dmatgen( n, n, b, n, seed_b );

  seed_c = (int)time( NULL );
  dmatgen( n, n, c, n, seed_c );

  seed_x = (int)time( NULL );
  dmatgen( n, 1, x, n, seed_x );

  alpha = a[n / 2];
  beta  = b[n / 2];


  /* MPI_WTIME ISSUE*/
  t0 = MPI_Wtime();
  HPL_dgemm( HplColumnMajor, HplNoTrans, HplNoTrans, n, n, n, alpha, a, n, b, n, beta, c, n );
  t1 = MPI_Wtime();

  t1 -= t0;
  dn = (double)n;
  if (t1 != 0.0 && t1 != -0.0){
    Gflops = 2.0e-9 * dn * dn * dn / t1;
    GFLOPS[iteration]= Gflops;
    
  }
  else
    Gflops = 0.0;

  cnrm = dnrm_inf( n, n, c, n );
  xnrm = dnrm_inf( n, 1, x, n );

  /* y <- c*x */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, 1.0, c, ldc, x, 1, 0.0, y, 1 );

  /* z <- b*x */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, 1.0, b, ldb, x, 1, 0.0, z, 1 );

  /* y <- alpha * a * z - y */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, alpha, a, lda, z, 1, -1.0, y, 1 );

  dmatgen( n, n, c, n, seed_c );

  /* y <- beta * c_orig * x + y */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, beta, c, ldc, x, 1, 1.0, y, 1 );

  sres = dnrm_inf( n, 1, y, n ) / cnrm / xnrm / n / HPL_dlamch( HPL_MACH_EPS );

  if (z) HPCC_free( z );
  if (y) HPCC_free( y );
  if (x) HPCC_free( x );
  if (c) HPCC_free( c );
  if (b) HPCC_free( b );
  if (a) HPCC_free( a );
}

int
HPCC_TestDGEMM(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure) {
  int i,j, n, lda, ldb, ldc, failure = 1;
  double *a=NULL, *b=NULL, *c=NULL, *x=NULL, *y=NULL,  *z=NULL, alpha, beta, sres, cnrm, xnrm;
  double Gflops = 0.0, dn, t0, t1;
  long l_n;
  FILE *outFile;
  int seed_a, seed_b, seed_c, seed_x;
  int GflopArray[params->DGEMM_N]; 
  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
  }
  
  /*
  * AUTHOR= OHAD KATZ
  * 
  * Added functionality such that DGEMM takes in values from the array of matrices
  * and the amount of repitions and calculates matrix multiplication using DGEMM
  * algorithm and reports it out.
  * 
  * After this is done, the code runs DGEMM like normal to assess the differences 
  * between the two methods.
  */
  for(int i_matrix = 0; i_matrix< params->DGEMM_N; i_matrix++){
      int repetitions= params->DGEMM_MatRep[i_matrix];
      struct timeval start_time, final_time;
      
      gettimeofday(&start_time, NULL);
      for (j = 0 ; j < repetitions; j++){
        
        n = params->DGEMM_MatSize[i_matrix]; 
        HPCC_DGEMM_Calculation(n, doIO, UGflops, Un, Ufailure, GflopArray, i);
        if (! a || ! b || ! c || ! x || ! y || ! z) {
          break;
        }
      
      }
      gettimeofday(&final_time, NULL);
      fprintf(outFile, "Time for array of size  %d : %ld Seconds \n",  params->DGEMM_MatSize[i_matrix], ((final_time.tv_sec * 1000000 + final_time.tv_usec) - (start_time.tv_sec * 1000000 + start_time.tv_usec)));
      fprintf(outFile, "# Repetitions: %d \n", params->DGEMM_MatRep[i_matrix]);
    }
  
  if (doIO) fprintf( outFile, "Scaled residual: %g\n", sres );

  if (sres < params->test.thrsh)
    failure = 0;   


	if (z) HPCC_free( z );
	if (y) HPCC_free( y );
	if (x) HPCC_free( x );
	if (c) HPCC_free( c );
	if (b) HPCC_free( b );
	if (a) HPCC_free( a );

	if (doIO) {
	  fflush( outFile );
	  fclose( outFile );
	}

	if (UGflops) *UGflops = Gflops;
	if (Un) *Un = n;
	if (Ufailure) *Ufailure = failure;
  return 0;
}
