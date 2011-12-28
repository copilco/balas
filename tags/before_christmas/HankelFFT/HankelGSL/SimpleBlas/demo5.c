#include <stdio.h>
#include <gsl/gsl_cblas.h>
#include <stdlib.h>
#include <vector.h>

int main (void)
{
	int N=20;
	
	vector<double> aa, bb, cc;
	aa.resize(N*N,0.);
	bb.resize(N,0.);
	cc.resize(N,0.);
	
	double *a, *b, *c;
	a = (double *)calloc( N*N, sizeof( double ) );
	b = (double *)calloc( N, sizeof( double ) );
	c = (double *)calloc( N, sizeof( double ) );
	if( a == NULL || b == NULL || c == NULL ) {
		printf( "\n Can't allocate memory for arrays\n");
		return 0;
	}
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			a[j*N+i]=0.;
	
	for(int i=0;i<N;i++)
    {
		a[i*N+i]=1.;
		b[i]=i;
		c[i]=0.;
		
		aa[i*N+i]=1.;
		bb[i]=i;
		cc[i]=0.;	
    }
    
	
	int lda=N;
	int ldb=1;
	int ldc=1;
	
	cblas_dgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, N, 1, N,
				 1.0, a, lda, b, ldb, 0.0, c, ldc);
	
	
	cblas_dgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, N, 1, N,
				 1.0, aa, lda, bb, ldb, 0.0, cc, ldc);
	
	//  cblas_dgemm(CblasRowMajor, CblasNoTrans, transB, N, 1, N, 1.,
	//	      a, lda, b, ldb, 1, c, ldc);
	
	for(int i=0;i<N;i++)
		printf("%e -> %e\n",c[i],cc[i]);
	
	
}
