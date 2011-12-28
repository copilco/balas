//g++ -framework Accelerate demo11.cpp  -lfftw3 -o aaa


#include <stdio.h>
#include <iostream.h>
#include <Accelerate/Accelerate.h>     
#include <complex>
#include "fftw3.h"
int main (void)
{
	
	int NN=8000;
	
	int lda = NN;
	
	
	__CLPK_doublecomplex * A;
	A = (__CLPK_doublecomplex*) malloc( NN*NN*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex * B;
	B = (__CLPK_doublecomplex*) malloc( NN*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex * C;
	C = (__CLPK_doublecomplex*) malloc( NN*sizeof(__CLPK_doublecomplex) );

	
	fftw_complex *aa,*bb,*cc;
	
	aa = (fftw_complex*) fftw_malloc(NN*NN * sizeof(fftw_complex));
	//	if(!aa) cout<<"\nError allocating array, please check ...";
	
	bb = (fftw_complex*) fftw_malloc(NN * sizeof(fftw_complex));
	//	if(!bb) cout<<"\nError allocating array, please check ...";
	
	cc = (fftw_complex*) fftw_malloc(NN * sizeof(fftw_complex));
	
	
	for(int i=0;i<NN*NN;i++)
	{
		A[i].r=i*0.0054;
		A[i].i=0.;
		
		aa[i][0]=i*0.0054;
		aa[i][1]=0.;
	}
		
	for(int i=0;i<NN;i++)
	{
		B[i].r=i*0.0013;
		B[i].i=0.;
		
		C[i].r=0.;
		C[i].i=0.;
		
		bb[i][0]=i*0.0013;
		bb[i][1]=0.;

		cc[i][0]=0.;
		cc[i][1]=0.;

	
	}
	
	int ldb = 1;
		
	int ldc = 1;
	
	
	
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

	for (int hh=0;hh<NN;hh++)
	{
		cout << hh <<" "<< C[hh].r << " "<< C[hh].i<< " -->";
		cout << cc[hh][0] << " * "<< cc[hh][1]<< "\n";
	}
		
	
	return 0;  
}
