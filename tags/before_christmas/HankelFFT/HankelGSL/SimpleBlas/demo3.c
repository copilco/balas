#include <stdio.h>
#include <gsl/gsl_cblas.h>
     
int
main (void)
{
  int lda = 3;
  
  double A[] = { 0., 0., 1.,
		1., 0., 0. };
  
  int ldb = 1;
       
  double B[] = { 1011,
		1021,
		1031 };
     
  int ldc = 1;
     
  double C[] = { 0.00,
		0.00 };
     
  /* Compute C = A B */
     
  cblas_dgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, 2, 1, 3,
	       1.0, A, lda, B, ldb, 0.0, C, ldc);
     
  printf ("[ %e, %e \n", C[0], C[1]);
  //  printf ("[ %g, %g %g \n", C[0], C[1], C[2]);
  // printf ("  %g, %g ]\n", C[2], C[3]);
     
  return 0;  
}
