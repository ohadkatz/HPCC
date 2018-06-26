/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* tstdgemm.c
 */

#include <hpcc.h>
#include <sys/time.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>

//#include <mkl_blacs.h>
/* Generates random matrix with entries between 0.0 and 1.0 */


/*
 *Baseline DGEMM(No optimization)
 *AUTHOR= Ohad Katz
 *UB Center for Computational Research
 */
/* C[i][j]=Σa[i][k]*b[k][j] */
static void 
ODGEMM_Calc(double *a, double *b, double *c, int N){
  int i,j,k;
  double sum;
  for (i=0; i<N; i++){

    for (j=0;j<N; j++){
      c[i+j*N] = 0;
    
      for(k=0;k<N;k++){
        c[i+j*N] += a[i+k*N]*b[k+j*N];
      }
    }
  }
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

/*HPCC Parallel DGEMM
* UNIVERSITY AT BUFFALO
* Author = Ohad Katz
* Mentor = Nikolay Simakov
* Center for Computational Research
*/
double
HPCC_PDGEMM_Scala_Calculation(int n, int doIO, double *UGflops, int *Un, int *Ufailure, double *GFLOPVAL){
  int  i,j,lda, ldb, ldc, failure = 1, ONE=1, ZERO=0;
  double *a=NULL, *b=NULL, *c=NULL, *x=NULL, *y=NULL, *z=NULL, alpha, beta, sres, cnrm, xnrm;
  double Gflops = 0.0, dn,oN, t0, t1;
  long l_n;
  int seed_a, seed_b, seed_c, seed_x;
  long int CONTXT, DescA[9], DescB[9], DescC[9],info, nprocs, iam;
  const long l_one= 1, l_zero= 0, l_negone= -1; /*l = long*/
  const double d_one=1, d_zero=0, d_negone=-1;  /*d = double*/
 
  if (n < 0) n = -n; /* if 'n' has overflown an integer */

  l_n = n;
  lda = ldb = ldc = n;
  const long clda= lda, cldb= ldb, cldc= ldc;

  a = HPCC_XMALLOC( double, l_n * l_n );
  b = HPCC_XMALLOC( double, l_n * l_n );
  c = HPCC_XMALLOC( double, l_n * l_n );

  x = HPCC_XMALLOC( double, l_n );
  y = HPCC_XMALLOC( double, l_n );
  z = HPCC_XMALLOC( double, l_n );

  //const long lld_local = Mmax( numroc_( &n, &n, &n, 0, &n ), 1 );
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
  /*=====================PARALLEL INITILIZATION============================*/

  HPL_blacspinfo( &iam, &nprocs );
  HPL_blacsget(&l_negone, &l_zero, &CONTXT);
  HPL_blacsgridinit(&CONTXT, "R", &nprocs, &l_one);
 

  /*=======================Descriptor Array Init===========================*/
  
  HPL_descinit( DescA, &l_n, &l_n, &l_n, &l_n, &l_zero, &l_zero, &CONTXT, &clda, &info );
  HPL_descinit( DescB, &l_n, &l_n, &l_n, &l_n, &l_zero, &l_zero, &CONTXT, &cldb, &info );
  HPL_descinit( DescC, &l_n, &l_n, &l_n, &l_n, &l_zero, &l_zero, &CONTXT, &cldb, &info );
  

  /*==================Broadcast Arrays Over Processess=====================*/
  
  HPL_pdgeadd("N", &l_n, &l_n, &d_one, a, &l_one, &l_one, DescA, &d_zero, a, &l_one, &l_one, DescA);
  HPL_pdgeadd("N", &l_n, &l_n, &d_one, b, &l_one, &l_one, DescB, &d_zero, b, &l_one, &l_one, DescB);


  /*
  *
  * Parallel DGEMM outline :
  * Trans A, Trans B, m , n , k , alpha, a , IA(First Row 1<=IA<=M_A), JA(First Column 1<=JA<=N_A), b, IB, JB, Desc_b, beta, c , IC, JC, Desc_c
  * 
  * Desc_A,Desc_B, Desc_c = 9 sized array: 
  *       DTYPE_= Descriptor Type
  *       CTX_= BLACS context
  *       M_= # Rows in Global Matrix
  *       N_= # Columns in Global Matrix
  *       MB_= Block Size for Row
  *       NB_ = Block Size for Column
  *       RSRC_ = Process of pxq where 1st row is distributed
  *       CSRC_ = Process of pxq where 1st column is distributed
  *       LLD_= Leading Dimension of a
  * 
  */

  /*======================================Main Calculation==========================================*/
  
  HPL_pdgemm("N","N", &l_n, &l_n, &l_n, &alpha, a, &l_one, &l_one, DescA, b, &l_one, &l_one, DescB, &beta, c, &l_one, &l_one, DescC );
 
  printf("%s", "done with pdgemm");
  t1 = MPI_Wtime();
  t1 -= t0;
  dn = (double)n;
  if (t1 != 0.0 && t1 != -0.0){
    Gflops = 2.0e-9 * dn * dn * dn / t1;
  }

  else
    Gflops = 0.0;
  
  *GFLOPVAL = Gflops;
  
  cnrm = dnrm_inf( n, n, c, n );
  xnrm = dnrm_inf( n, 1, x, n );

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
/*PDGEMM FIN*/


double
HPCC_DGEMM_Calculation(int n, int doIO, double *UGflops, int *Un, int *Ufailure, double *GFLOPVAL){
  int i,j,lda, ldb, ldc, failure = 1;
  double *a=NULL, *b=NULL, *c=NULL, *x=NULL, *y=NULL, *z=NULL, alpha, beta, sres, cnrm, xnrm;
  double Gflops = 0.0, dn,oN, t0, t1;
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
  
  // FOR EDUCATIONAL PURPOSES: Unoptimized Matrix-Matrix Multiplication= ODGEMM_Calc(a,b,c,n);
 
  HPL_dgemm( HplColumnMajor, HplNoTrans, HplNoTrans, n, n, n, alpha, a, n, b, n, beta, c, n );
 
  t1 = MPI_Wtime();
  t1 -= t0;
  dn = (double)n;
  if (t1 != 0.0 && t1 != -0.0){
    Gflops = 2.0e-9 * dn * dn * dn / t1;
    // EDUCATIONAL : oGflops= 2.0e-9 * oN *oN *oN/t1;
   
  }
  else
    Gflops = 0.0;
  
  *GFLOPVAL = Gflops;
  
  cnrm = dnrm_inf( n, n, c, n );
  xnrm = dnrm_inf( n, 1, x, n );

  /* y <- c*x */
  //pdgemv_("N", &l_n, &l_n, &alpha, a, &ONE, &ONE, descA, x, &ONE, &ONE, descx, &ONE, &beta, y, &ONE, &ONE, descy, &ONE);
 
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, 1.0, c, ldc, x, 1, 0.0, y, 1 );

  /* z <- b*x */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, 1.0, b, ldb, x, 1, 0.0, z, 1 );
  //pdgemv_("N", &l_n, &l_n, &alpha, A, &ONE, &ONE, descA, x, &ONE, &ONE, descx, &ONE, &beta, y, &ONE, &ONE, descy, &ONE);
 

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
  double Gflop = 0.0,  start,end;
  FILE *outFile;
  FILE *Rfile;
  double timer[params->DGEMM_N];
  double maximums[params->DGEMM_N], minimums[params->DGEMM_N], avg[params->DGEMM_N], sresArr[params->DGEMM_N],stddev[params->DGEMM_N], sum[params->DGEMM_N];
  double avgSquare,stddevAvg, sumSquare;

    
  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    /*Added*/
    Rfile = fopen( params-> results, "a");
    
    if (! outFile ) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
    /*Added*/
    if(! Rfile) {
      Rfile=stderr;
      fprintf( Rfile, "Cannot output file.\n");
      return 1;
    }
  }
  if(doIO) fprintf(Rfile,"%s", "N,RunID,GFLOPS\n");

  /*
  *
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
  * [ Matrix Size |  Repetitions | Total Time  |  Avg GFLOPS  |  Min GFLOP   |   Max GFLOP   |  sres ]
  * [    128      |     1000     |     1 sec   |     47.123   |    57.391    |    14.532     |  .03  ] 
  * [    256      |      500     |    30 sec   |     50.432   |    61.123    |    16.431     |  .021 ]
  * [    512      |      250     |  1:30 sec   |     48.123   |    53.912    |    23.13      |  .23  ]
  * [    1024     |      100     |  2:10 sec   |     49.13    |    57.93     |    21.34      |  .14  ]
  * [------------------------------------------------------------------------------------------------]
  * 
  */

  /*Iterate through each Matrix Size and repeat operations on it*/
  for(int i_matrix = 0; i_matrix< params->DGEMM_N; i_matrix++){
      int repetitions= params->DGEMM_MatRep[i_matrix];
      double Gflop, OGflop;
      max = 0;
      min = INT_MAX;
      start = MPI_Wtime();
      for (int repnum = 0 ; repnum < repetitions; repnum++){
        /*Set n to fixed size array in Input File*/
        n = params->DGEMM_MatSize[i_matrix]; 
        sres = HPCC_PDGEMM_Scala_Calculation(n, doIO, UGflops, Un, Ufailure, &Gflop);
        
        if (doIO) fprintf(Rfile,"%d,%d,%f\n", params->DGEMM_MatSize[i_matrix],repnum+1,Gflop);
        if (Gflop>max) max=Gflop;
        if (Gflop<min) min=Gflop;
        
        /*Sum up both Gflops and Gflop^2(For std.Deviation)*/
        sum[i_matrix] += Gflop;
        stddev[i_matrix] += Gflop*Gflop;
        
      }
      
      sresArr[i_matrix]= sres;
      /*
      * Average = sum/N
      * Standard Deviation= Mean(Gflops^2)-Mean(Gflops)^2
      */

      /*Calculations for each fixed size matrix*/
      avg[i_matrix] = sum[i_matrix]/repetitions;
      avgSquare= avg[i_matrix]*avg[i_matrix];
      stddevAvg= stddev[i_matrix]/repetitions;
      stddev[i_matrix]=(stddevAvg>avgSquare ? sqrt(stddevAvg-avgSquare): 0);
     
      start = MPI_Wtime()-start;
      
      timer[i_matrix] = start;
      maximums[i_matrix]=max;
      minimums[i_matrix]=min;
      if (doIO){      
      fprintf(outFile, "Scaled Residual: %f\n" , sresArr[i_matrix]);
      fprintf(outFile, "Time for array of size %d : %f seconds\n",  params->DGEMM_MatSize[i_matrix], timer[i_matrix]);
      fprintf(outFile, "# Repetitions: %d\n", params->DGEMM_MatRep[i_matrix]);
      

      fprintf(outFile, "Maximum GFLOP/S: %f\n", maximums[i_matrix]);
      fprintf(outFile, "Minimum GFLOP/S: %f\n", minimums[i_matrix]);
      fprintf(outFile, "Avg GFLOP/S: %f \n", avg[i_matrix]);
      }
      //fprintf(outFile, "---------------------------------------------------------\n");
    }

  /*OUTPUT TABLE*/
  if(doIO){
      fprintf(outFile, "|----------------------------------------------------------------------------------------|\n" );
      fprintf(outFile,"%-10s %-10s %-10s  %-10s  %-10s%-10s %-10s %-10s\n", "Mat.Size", "Repeat Amt.", "Tot.Time(s)", "Avg GFLOP", "Std.Dev",  "Min GFLOP","Max GFLOP", "Scal.Res");
      fprintf(outFile, "|----------------------------------------------------------------------------------------|\n" );
      
      for(i = 1 ; i< params->DGEMM_N; i++){
        fprintf(outFile,"%10d %10d  %10.2f %10.2f %10.2f %10.2f %10.2f  %10.2E\n", params->DGEMM_MatSize[i], params->DGEMM_MatRep[i], timer[i], avg[i],stddev[i], minimums[i],maximums[i], sresArr[i]);
      }
  }

  if (doIO) fprintf( outFile, "\nScaled residual: %g\n", sres );

  if (sres < params->test.thrsh)
    failure = 0;   

	if (doIO) {
	  fflush( outFile );
	  fclose( outFile );
    fflush( Rfile );
    fclose( Rfile );
	}

	if (UGflops) *UGflops = Gflop;
	if (Un) *Un = n;
	if (Ufailure) *Ufailure = failure;
  return 0;
}
