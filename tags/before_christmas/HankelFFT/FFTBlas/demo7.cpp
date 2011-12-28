/*
 *  demo7.cpp
 *  
 *
 *  Created by camilo Ruiz Méndez on 07/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

//#include "demo7.h"

#include <complex.h>
#include <alloc.h>
#include <Accelerate/Accelerate.h>     
#include <stdio.h>
#include <math.h>
#define complex complex<double>
#define <gsl/gsl_complex.h>

int main()
{
	complex *A;
	complex *B;
	complex *C;
	
	A=(complex *)malloc((9)*sizeof(complex));
	B=(complex *)malloc((3)*sizeof(complex));
	C=(complex *)malloc((3)*sizeof(complex));
	//if(!sxr) printf("\nerror de colocacion en la asignacion complex");
	
	int lda = 3;
	int ldb = 1;
	int ldc = 1;
	
	//	doublecomplex A[]={1., 1., 1.};
	
	for(int i=0;i<9;i++)
    {
        A[i]=complex(1.,0.);		
    	printf(" A %e\n",A[i]);
	}

	
	for(int i=0;i<3;i++)
    {
        B[i]=complex(1.,0.);	
		C[i]=complex(0.,0.);	
    	printf(" B %e C %e\n",real(A[i]),real(B[i]));
	}
	
	complex alpha=complex(1.,0.);
	cblas_zgemm (CblasRowMajor, 
				 CblasNoTrans, CblasNoTrans, lda, 1, lda,
				 *alpha, A, lda, B, ldb, 0.0, C, ldc);
	
	
	zgemm_(char *transa, char *transb, integer *m, integer *
		   n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, 
		   doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c,
		   integer *ldc)
	
}
