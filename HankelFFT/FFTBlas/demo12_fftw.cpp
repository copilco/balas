//g++ -framework Accelerate demo11.cpp  -lfftw3 -o aaa


#include <stdio.h>
#include <iostream.h>
#include <Accelerate/Accelerate.h>     
#include <complex>
#include "fftw3.h"
int main (void)
{
	
	FILE *out0,*out1;
	out0=fopen("out0.txt","w");
	out1=fopen("out1.txt","w");
	
	int NN=50;
	
	int lda = NN;
	int MM=50;
	
	__CLPK_doublecomplex * A;
	A = (__CLPK_doublecomplex*) malloc( NN*NN*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex * B;
	B = (__CLPK_doublecomplex*) malloc( NN*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex * C;
	C = (__CLPK_doublecomplex*) malloc( NN*sizeof(__CLPK_doublecomplex) );

	
	fftw_complex *aa,*bb,*cc;
	fftw_complex *bb1,*cc1;
	
	aa = (fftw_complex*) fftw_malloc(NN*NN * sizeof(fftw_complex));
	//	if(!aa) cout<<"\nError allocating array, please check ...";
	
	bb = (fftw_complex*) fftw_malloc(NN * sizeof(fftw_complex));
	//	if(!bb) cout<<"\nError allocating array, please check ...";
	
	cc = (fftw_complex*) fftw_malloc(NN * sizeof(fftw_complex));
	
	
	bb1 = (fftw_complex*) fftw_malloc(NN*MM * sizeof(fftw_complex));	
	cc1 = (fftw_complex*) fftw_malloc(NN*MM * sizeof(fftw_complex));
	
	for(int i=0;i<NN*NN;i++)
	{
		A[i].r=i*0.0054;
		A[i].i=0.;
		
		aa[i][0]=0;//i*0.0054;
		aa[i][1]=0.;
	}
	
	
	/*
	for (int j=0; j<MM; j++)
		for(int i=0;i<NN;i++)
		{
			aa[j*NN+i][0]=j;//i*0.0054;
			aa[j*NN+i][1]=0.;
		}
	*/
	
	 for(int i=0;i<NN;i++)
	 {
	 aa[i*NN+i][0]=2;//i*0.0054;
	 aa[i*NN+i][1]=0.;
	 }
	
		
	for(int i=0;i<NN;i++)
	{
		B[i].r=i*0.001;
		B[i].i=0.;
		
		C[i].r=0.;
		C[i].i=0.;
		
		bb[i][0]=i*0.001;
		bb[i][1]=0.;

		cc[i][0]=0.;
		cc[i][1]=0.;

	
	}
	
	
	for(int j=0;j<MM;j++)
		for(int i=0;i<NN;i++)
		{
			
			bb1[j*NN+i][0]=1.;//i*0.+j*1.;
			bb1[j*NN+i][1]=0.;//i*0.0013;
			
			cc1[j*NN+i][0]=0.;
			cc1[j*NN+i][1]=0.;
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
			   
			   
	
	cblas_zgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, MM/2, lda,
				 &alpha, aa, lda, bb1,MM/2,&gamma, cc1, MM/2);
	
	
	/*
	cblas_zgemm (CblasColMajor, 
				 CblasNoTrans, CblasNoTrans, MM/2, NN, MM/2,
				 &alpha, aa, lda, bb1,MM,&gamma, cc1, MM);
	 */
	
	
	
	/*
	for (int hh=0;hh<NN;hh++)
	{
		cout << hh <<" "<< C[hh].r << " "<< C[hh].i<< " -->";
		cout << cc[hh][0] << " * "<< cc[hh][1]<< "\n";
	}
	 */
	
	for( int kk=0;kk<MM;kk++)
	{
		cout <<  "\n";
		for (int hh=0;hh<NN;hh++)
		{
			//cout << hh <<" "<< C[kk*NN+hh].r << " "<< C[kk*NN+hh].i<< " -->";
			cout << cc1[kk*NN+hh][0] << " * "<< cc1[kk*NN+hh][1]<< "\n";
			fprintf(out0,"%e\n", bb1[kk*NN+hh][0] );
			fprintf(out1,"%e\n", cc1[kk*NN+hh][0] );
		}
	}	
	
	return 0;  
}
