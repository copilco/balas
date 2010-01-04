#include <stdio.h>
#include <gsl/gsl_cblas.h>
//#include <Accelerate/Accelerate.h>     
#include <complex>
#define complex complex<double>
#include "fftw3.h"
//#include "grid.h"
//#include "wavef3.h"



#include <complex.h>


int main(){
	complex xx;
	complex yy = complex(1,2.718);
	xx = log(yy/3);
	cout << 1+xx;
}

/*

int
main (void)
{
 
  int lda = 3;
  
	double A[] = { 
		1., 1., 1., 
		0., 0., 0.,
	    0., 0., 0. };
	
	//0.11, 0.12, 0.13,
	//	0.21, 0.22, 0.23 };
  
	int ldb = 1;
	
	double B[] = {1., 0., 0.};//{ 1011, 1012,
		//1021, 1022,
		//1031, 1032 };
     
  int ldc = 1;
  
	double C[] = { 0.00, 0.00,
	0.0 };
	
	fftw_complex *aa,*bb,*cc,*dd;
	
	aa = (fftw_complex*) fftw_malloc(9 * sizeof(fftw_complex));
//	if(!aa) cout<<"\nError allocating array, please check ...";

	bb = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
//	if(!bb) cout<<"\nError allocating array, please check ...";

	cc = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
	
	dd = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
//	if(!cc) cout<<"\nError allocating array, please check ...";

	
	for(int i=0;i<9;i++)
    {
        aa[i]=A[i];
		aa[i]=0;//A[i];
		dd[i]=bb[i]*cc[i];
//		printf("a %e %e %e \n",A[i],real(aa[i]),imag(aa[i]));
    }
	printf("\n");
	
	for(int i=0;i<3;i++)
	{
		bb[i]=B[i];
		bb[i]=B[i];
	//	printf("b %e %e %e\n",B[i],bb[i][0],bb[i][1] );		
	}
	printf("\n");
	
	//aa[0][1]=1.;
	//bb[0][1]=0.;
	

  
  // Compute C = A B /
  
  cblas_dgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, lda, 1, lda,
	       1.0, A, lda, B, ldb, 0.0, C, ldc);

	complex a11(1.,0.);
  cblas_zgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, lda, 1, lda,
	       a11, aa, lda, bb,ldb,0.0, cc, ldc);
  
	
	for(int i=0;i<3;i++)
	{
		printf("c %e %e %e\n",C[i],cc[i][0],cc[i][1]);
	}
	printf("\n");
	
	/*
	printf ("[ %g, %g\n",   C[0], C[1]);
  printf ("  %g,  ]\n", C[2]);

  printf ("[ %g, %g\n",   cc[0][0], cc[1][0]);
  printf ("  %g,  ]\n", cc[2][0]);
	
	printf ("[ %g, %g\n",   cc[0][1], cc[1][1]);
	printf ("  %g, ]\n", cc[2][1] );
	//
// printf ("[ %g, %g\n",   cc[0][1], cc[1][1]);
*/


//return 0;  
//}
