#include <stdio.h>
#include <gsl/gsl_cblas.h>
#include <stdlib.h>
#include "arrai.h"

int main (void)
{
	int N=20;
	arrai a(N*N,0.);
	arrai b(N,0.);
	arrai c(N,0.);
	
	for(int i=0;i<N;i++)
	{
		a.v[i*N+i]=1.;
		b.v[i]=i;
		c.v[i]=0.;
    }
    
	
	int lda=N;
	int ldb=1;
	int ldc=1;
	
	cblas_dgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, N, 1, N,
				 1.0, a.v, lda, b.v, ldb, 0.0, c.v, ldc);

  //  cblas_dgemm(CblasRowMajor, CblasNoTrans, transB, N, 1, N, 1.,
  //	      a, lda, b, ldb, 1, c, ldc);

  for(int i=0;i<N;i++)
    printf("%e \n",c.v[i]);


}
