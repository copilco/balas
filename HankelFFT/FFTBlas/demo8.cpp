#include <stdio.h>
#include <iostream.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_cblas.h>
//#include <gsl/gsl_cblas.h>
//#include <Accelerate/Accelerate.h>     
#include <complex>
#include "fftw3.h"
//#include "grid.h"
//#include "wavef3.h"
int
main (void)
{
	
	int lda = 3;
	
	gsl_complex *A;
	A=new gsl_complex[9];
	
	gsl_complex *B;
	B=new gsl_complex[3];
	
	gsl_complex *C;
	C=new gsl_complex[3];
	
	
	
	fftw_complex *aa,*bb,*cc;
	
	aa = (fftw_complex*) fftw_malloc(9 * sizeof(fftw_complex));
	//	if(!aa) cout<<"\nError allocating array, please check ...";
	
	bb = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
	//	if(!bb) cout<<"\nError allocating array, please check ...";
	
	cc = (fftw_complex*) fftw_malloc(3 * sizeof(fftw_complex));
	
	
	for(int i=0;i<9;i++)
	{
		GSL_SET_COMPLEX(&A[i], 1.345, 0);//A[i]=1;
		aa[i][0]=1.345;
		aa[i][1]=0.;
	}
		
	for(int i=0;i<3;i++)
	{
		GSL_SET_COMPLEX(&B[i], 2., 1);
		GSL_SET_COMPLEX(&C[i], 0, 0);
	
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
	
	gsl_complex alpha;
	gsl_complex beta;
	gsl_complex gamma;
	
	GSL_SET_COMPLEX(&alpha, 1, 0);//A[i]=1;
	GSL_SET_COMPLEX(&beta, 1, 0);//A[i]=1;
	GSL_SET_COMPLEX(&gamma, 0, 0);//A[i]=1;
	
	cblas_zgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, 1, lda,
				 &alpha, A, lda, B, ldb, &gamma, C, ldc);
	
	cblas_zgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, 1, lda,
				 &alpha, aa, lda, bb,ldb,&gamma, cc, ldc);

	for (int hh=0;hh<3;hh++)
	{
		cout << GSL_REAL(C[hh]) << " "<< GSL_IMAG(C[hh])<< " -->";
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
