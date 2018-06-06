/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* tstdgemm.c
 */

#include <hpcc.h>
#include <sys/time.h>
#include <limits.h>
/* Generates random matrix with entries between 0.0 and 1.0 */


double Max(double arr[], int size){
  double max= arr[0];
  int i;
  for(i=0; i<size;i++){
    if (arr[i]>max ){
      max= arr[i];
    }
  }
  return max;
}

double Min(double arr[], int size){
double min= arr[0];
int i;
  for(i=0; i<size;i++){
    if (arr[i]<min && arr[i]>0){
      min= arr[i];
    }
    if (arr[i]<=0){
      continue;
    }
  }
  return min;
}

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
double
HPCC_DGEMM_Calculation(int n, int doIO, double *UGflops, int *Un, int *Ufailure, double *GFLOPVAL){
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

  t0 = MPI_Wtime();
  HPL_dgemm( HplColumnMajor, HplNoTrans, HplNoTrans, n, n, n, alpha, a, n, b, n, beta, c, n );
  t1 = MPI_Wtime();
  //printf("n = %d\nt1 =%f\nt0= %f\n",n , t1,t0);
  t1 -= t0;
  dn = (double)n;
  if (t1 != 0.0 && t1 != -0.0){
    Gflops = 2.0e-9 * dn * dn * dn / t1;
    *GFLOPVAL = Gflops;
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
  return sres;
}

int
HPCC_TestDGEMM(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure) {
  int i,j, n, failure = 1;
  double sres, cnrm, xnrm,max,min;
  double Gflop = 0.0,start,end;
  FILE *outFile;
  double timer[params->DGEMM_N];
  double maximums[params->DGEMM_N], minimums[params->DGEMM_N], avg[params->DGEMM_N], sresArr[params->DGEMM_N],stddev[params->DGEMM_N], sum[params->DGEMM_N];
  double avgSquare,sumSquare;
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
  * 
  * 
  *  Sample Output:
  * [------------------------------------------------------------------------------------------------]
  * [ Matrix Size |  Repetitions | Total Time  |  Avg GFLOPS  |  Min GFLOP   |   Max GFLOP   |  srec ]
  * [    128      |     1000     |   30 sec    |     47.123   |    57.391    |    14.532     |  .03  ] 
  * [    256      |     500      |   40 sec    |     50.432   |    61.123    |    16.431     |  .021 ]
  * [    512      |     250      | 1 min 30 sec|     48.123   |    53.912    |    23.13      |  .23  ]
  * [    1024     |     100      | 2 min 10 sec|     49.13    |    57.93     |    21.34      |  .14  ]
  * [------------------------------------------------------------------------------------------------]
  */

  /*Iterate through each Matrix Size and repeat operations on it*/
  for(int i_matrix = 0; i_matrix< params->DGEMM_N; i_matrix++){
      int repetitions= params->DGEMM_MatRep[i_matrix];
      double Gflop;
      avg[i_matrix] = 0; 
      sum[i_matrix] = 0;
      max = 0;
      min = INT_MAX;
      start = MPI_Wtime();
      for (int repnum = 0 ; repnum < repetitions; repnum++){\
        /*Set n to fixed size array in Input File*/
        n = params->DGEMM_MatSize[i_matrix]; 
        sres = HPCC_DGEMM_Calculation(n, doIO, UGflops, Un, Ufailure, &Gflop);
        
        
        if (Gflop>max) max=Gflop;
        if (Gflop<min) min=Gflop;
        
        /*Sum up both Gflops and Gflop^2(For std.Deviation)*/
        sum[i_matrix] += Gflop;
        stddev[i_matrix] += Gflop*Gflop;
        
        // if (! a || ! b || ! c || ! x || ! y || ! z) {
        //   break;
        // }
      }
      
      sresArr[i_matrix]= sres;
      /*
      * Average = sum/N
      * Standard Deviation= Mean(Gflops^2)-Mean(Gflops)^2
      */
      /*Calculations for each fixed size matrix*/
      avg[i_matrix] = sum[i_matrix]/repetitions;
      avgSquare= avg[i_matrix]*avg[i_matrix];
      stddev[i_matrix]= sqrt((stddev[i_matrix]/repetitions)-avgSquare);
      
      end = MPI_Wtime();
      start = end-start;
      
      timer[i_matrix] = start;
      fprintf(outFile, "Scaled Residual: %g\n" , sres);
      fprintf(outFile, "\nTime for array of size %d : %f seconds\n",  params->DGEMM_MatSize[i_matrix], timer[i_matrix]);
      fprintf(outFile, "\n# Repetitions: %d\n", params->DGEMM_MatRep[i_matrix]);
      
      maximums[i_matrix]=max;
      minimums[i_matrix]=min;

      fprintf(outFile, "\nMaximum GFLOP/S: %f\n", max);
      fprintf(outFile, "\nMinimum GFLOP/S: %f\n", min);
      fprintf(outFile, "\nAvg GFLOP/S: %f \n", avg[i_matrix]);
      
      fprintf(outFile, "---------------------------------------------------------\n");
    }
  /*OUTPUT TABLE*/
  fprintf(outFile, "|----------------------------------------------------------------------------------------|\n" );
  fprintf(outFile,"|%-10s %-10s %-10s  %-10s  %-10s%-10s %-10s %-10s\n", "Mat.Size", "Repeat Amt.", "Tot.Time(s)", "Avg GFLOP", "Std.Dev",  "Min GFLOP","Max GFLOP", "Scal.Res|");
  fprintf(outFile, "|----------------------------------------------------------------------------------------|\n" );
  for(i = 0 ; i< params->DGEMM_N; i++){
    fprintf(outFile,"%10d %10d  %10.2f %10.2f %10.2f %10.2f %10.2f  %10.2E\n", params->DGEMM_MatSize[i], params->DGEMM_MatRep[i], timer[i], avg[i],stddev[i], minimums[i],maximums[i], sresArr[i]);
  }
  

  if (doIO) fprintf( outFile, "\nScaled residual: %g\n", sres );

  if (sres < params->test.thrsh)
    failure = 0;   


	// if (z) HPCC_free( z );
	// if (y) HPCC_free( y );
	// if (x) HPCC_free( x );
	// if (c) HPCC_free( c );
	// if (b) HPCC_free( b );
	// if (a) HPCC_free( a );

	if (doIO) {
	  fflush( outFile );
	  fclose( outFile );
	}

	if (UGflops) *UGflops = Gflop;
	if (Un) *Un = n;
	if (Ufailure) *Ufailure = failure;
  return 0;
}
