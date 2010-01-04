
#include <stdio.h>
#include <Accelerate/Accelerate.h>     
#include <complex>
#include "fftw3.h"
//#include "grid.h"
//#include "wavef3.h"
int
main (void)
{
	
	double A[] = { 1., 1., 1., 
		           0., 2., 0.,
	               1., 0., 0. };
	
	double B[] = {1., 1., 1.};
	
	double C[] = { 0.00, 0.00,	0.0 };
	
 
	int lda = 3;
  
	int ldb = 1;
	
	int ldc = 1;
	
	fftw_complex *aa,*bb,*cc;
	
	aa = (fftw_complex*) fftw_malloc(9 * sizeof(fftw_complex));
//	if(!aa) cout<<"\nError allocating array, please check ...";

	bb = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
//	if(!bb) cout<<"\nError allocating array, please check ...";

	cc = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
//	if(!cc) cout<<"\nError allocating array, please check ...";

	double *raa[9];
	double *rbb[3];
	double *rcc[3];
	
	for(int i=0;i<9;i++)
    {
        aa[i][0]=A[i];
		aa[i][1]=0.;
		
		raa[i]=&aa[i][0];

    	printf(" aa %e %e\n",aa[i][0],*raa[i]);
	}

		
  for(int i=0;i<3;i++)
  {
	  bb[i][0]=B[i];
	  cc[i][0]=C[i];
	  
	  rbb[i]=&bb[i][0];
	  rcc[i]=&cc[i][0];
  }
  

  /* Compute C = A B */
  
	
  cblas_dgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, lda, 1, lda,
	       1.0, A, lda, B, ldb, 0.0, C, ldc);

 
 
	/*
  cblas_dgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, 1, lda,
				 1.0, &raa, lda, &rbb,ldb,0.0, &rcc, ldc);
*/
	
	printf ("[ %g, %g\n",   C[0], C[1]);
	printf ("  %g,  ]\n", C[2]);
	
	printf ("[ %g, %g\n",   cc[0][0], cc[1][0]);
	printf ("  %g,  ]\n", cc[2][0]);
	
	printf ("[ %g, %g\n",  *rcc[0], *rcc[1]);
	printf ("  %g,  ]\n", *rcc[2]);
	
	
	
	
	 cblas_dgemm (CblasRowMajor, 
	 CblasNoTrans, CblasNoTrans, lda, 1, lda,
	 1.0, aa[0], lda, bb[0],ldb,0.0, cc[0], ldc);
	
  printf ("[ %g, %g\n",   C[0], C[1]);
  printf ("  %g,  ]\n", C[2]);

  printf ("[ %g, %g\n",   cc[0][0], cc[1][0]);
  printf ("  %g,  ]\n", cc[2][0]);

	printf ("[ %g, %g\n",  *rcc[0], *rcc[1]);
	printf ("  %g,  ]\n", *rcc[2]);

//	printf ("[ %g, %g\n",   cc[0][1], cc[1][1]);
//	printf ("  %g, ]\n", cc[2][1] );
 

  return 0;  
}
