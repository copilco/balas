#include <stdio.h>
//#include <gsl/gsl_cblas.h>
#include <Accelerate/Accelerate.h>     
int
main (void)
{
  int lda = 3;
  
  double A[] = { 0.11, 0.12, 0.13,
		0.21, 0.22, 0.23 };
  
  int ldb = 2;
  
  double B[] = { 1011, 1012,
		1021, 1022,
		1031, 1032 };
     
  int ldc = 2;
  
  double C[] = { 0.00, 0.00,
		0.00, 0.00 };
  
  /* Compute C = A B */
  
  cblas_dgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, 2, 2, 3,
	       1.0, A, lda, B, ldb, 0.0, C, ldc);
  
  printf ("[ %g, %g\n", C[0], C[1]);
  printf ("  %g, %g ]\n", C[2], C[3]);
  
  return 0;  
}
