#include <stdio.h>
#include <iostream.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_cblas.h>
//#include <gsl/gsl_cblas.h>
#include <Accelerate/Accelerate.h>     
#include <complex>
#include "fftw3.h"
//#include "grid.h"
//#include "wavef3.h"
int
main (void)
{
	
	int lda = 3;
	
	
	__CLPK_doublecomplex * A;
	A = (__CLPK_doublecomplex*) malloc( 9*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex * B;
	B = (__CLPK_doublecomplex*) malloc( 3*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex * C;
	C = (__CLPK_doublecomplex*) malloc( 3*sizeof(__CLPK_doublecomplex) );

	
	fftw_complex *aa,*bb,*cc;
	
	aa = (fftw_complex*) fftw_malloc(9 * sizeof(fftw_complex));
	//	if(!aa) cout<<"\nError allocating array, please check ...";
	
	bb = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
	//	if(!bb) cout<<"\nError allocating array, please check ...";
	
	cc = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
	
	
	for(int i=0;i<9;i++)
	{
		A[i].r=1.345;
		A[i].i=0.;
		
		aa[i][0]=1.345;
		aa[i][1]=0.;
	}
		
	for(int i=0;i<3;i++)
	{
		B[i].r=2.;
		B[i].i=1.;
		
		C[i].r=0.;
		C[i].i=0.;
		
		bb[i][0]=2.;
		bb[i][1]=1.;

		cc[i][0]=0.;
		cc[i][1]=0.;

	
	}
	
	/*
	 A[] = { 
	 1., 0., 0., 
	 0., 0., 0.,
	 0., 0., 0. };
	 
	 //0.11, 0.12, 0.13,
	 //	0.21, 0.22, 0.23 };
	 */
	int ldb = 1;
	
	//double B[] = {1., 1., 1.};//{ 1011, 1012,
	//1021, 1022,
	//1031, 1032 };
	
	int ldc = 1;
	
	//double C[] = { 0.00, 0.00,
	//0.0 };
	
	//	if(!cc) cout<<"\nError allocating array, please check ...";
	
	/*
	for(int i=0;i<9;i++)
    {
        aa[i][0]=A[i];
		aa[i][1]=A[i];
    }
	 */
	
	//aa[0][1]=1.;
	//bb[0][1]=0.;
	
	//for(int i=0;i<4;i++)
	//	cc[i][0]=C[i];
	
	
	/* Compute C = A B */
	
	__CLPK_doublecomplex alpha;
	__CLPK_doublecomplex beta;
	__CLPK_doublecomplex gamma;
	
	alpha.r=1;
	alpha.i=0;
	
	gamma.r=1;
	gamma.i=0;
	
	
	cblas_zgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, 1, lda,
				 &alpha, A, lda, B, ldb, &gamma, C, ldc);
	
	cblas_zgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, 1, lda,
				 &alpha, aa, lda, bb,ldb,&gamma, cc, ldc);

	for (int hh=0;hh<3;hh++)
	{
		cout << C[hh].r << " "<< C[hh].i<< " -->";
		cout << cc[hh][0] << " * "<< cc[hh][1]<< "\n";
	}
		
	//printf ("[ %g, %g\n",   C[0], C[1]);
	//printf ("  %g,  ]\n", C[2]);
	/*
	
	printf ("[ %g, %g\n",   cc[0][0], cc[1][0]);
	printf ("  %g,  ]\n", cc[2][0]);
	
	printf ("[ %g, %g\n",   cc[0][1], cc[1][1]);
	printf ("  %g, ]\n", cc[2][1] );
	
	*/
	return 0;  
}
