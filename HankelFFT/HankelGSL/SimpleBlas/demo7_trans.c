#include <stdio.h>
#include <gsl/gsl_cblas.h>
#include <stdlib.h>

int main (void)
{
	int N=4;
	int MM=3;
	
	
	double *a, *b, *c;
	a = (double *)calloc( N*N, sizeof( double ) );
	b = (double *)calloc( N*MM, sizeof( double ) );
	c = (double *)calloc( N*MM, sizeof( double ) );
	if( a == NULL || b == NULL || c == NULL )
	{
		printf( "\n Can't allocate memory for arrays\n");
		return 0;
	}
	
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			if(j<i)
				a[j*N+i]=0.;
			else
				a[j*N+i]=1.;
		}
	
	for(int j=0;j<N;j++)
		a[j*N+j]=1.;
	
	for(int j=0;j<N;j++)
		for(int i=0;i<MM;i++)
		{
			b[j*MM+i]=1.;
			c[j*MM+i]=0.;
		//	printf("%e ",b[j*N+i]);
		}    
	
	
	printf("\n");
	printf("\n");
	for(int j=0;j<N;j++)
		for(int i=0;i<MM;i++)
		{
			printf("%e ",b[j*MM+i]);
		} 


	printf("\n");
	printf("\n");
	
	int lda=N;
	int ldb=MM;
	int ldc=MM;
	
	
	cblas_dgemm (CblasColMajor, 
				 CblasNoTrans, CblasNoTrans, N, MM, N,
				 1.0, a, lda, b, ldb, 0.0, c, ldc);
	
	
	/*
	cblas_dgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, N, MM, N,
				 1.0, a, lda, b, ldb, 0.0, c, ldc);
   */
	
	//  cblas_dgemm(CblasRowMajor, CblasNoTrans, transB, N, 1, N, 1.,
	//	      a, lda, b, ldb, 1, c, ldc);
	
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)
			printf("%e ",a[j*N+i]);
		printf("\n");
	}
	
	printf("\n");
	printf("\n");
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<MM;i++)
			printf("%e ",b[j*MM+i]);
		printf("\n");
	}
	
	printf("\n");
	printf("\n");
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<MM;i++)
			printf("%e ",c[j*MM+i]);
		printf("\n");
	}
	
	
}
