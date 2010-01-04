#include <iostream>
#include <cstdlib>
#include "clapack.h"

int main (int argc, const char * argv[]) {
	
	char eigens = 'V';
	char upper= 'U';
	__CLPK_integer lwork = 12;
	__CLPK_integer LDA = 12;
	__CLPK_integer N = 12;
	__CLPK_integer res;
	
	__CLPK_doublecomplex * H2prime;
	H2prime = (__CLPK_doublecomplex*) malloc( LDA*N*sizeof(__CLPK_doublecomplex) );
	
	__CLPK_doublereal *W;
	W = (__CLPK_doublereal*) malloc( LDA*N*sizeof(__CLPK_doublereal) );
	
	__CLPK_doublecomplex *Work;
	Work = (__CLPK_doublecomplex*) malloc( LDA*N*sizeof(__CLPK_doublecomplex) );
	
	double *rwork;
	rwork = (double*) malloc( LDA*N*sizeof(double) );
	
	
	lwork=-1;
	
	__CLPK_integer n=100;
	
	__CLPK_integer n1=1;
//	n1 = (__CLPK_integer*) malloc( 2*sizeof(__CLPK_integer) );
	//=100;
	
	__CLPK_doublecomplex *lower;
	lower = (__CLPK_doublecomplex*) malloc( (n-1)*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex *middle;
	middle = (__CLPK_doublecomplex*) malloc( n*sizeof(__CLPK_doublecomplex) );

	__CLPK_doublecomplex *upp;
	upp = (__CLPK_doublecomplex*) malloc( (n-1)*sizeof(__CLPK_doublecomplex) );
	
	__CLPK_doublecomplex *b;
	b = (__CLPK_doublecomplex*) malloc( n*sizeof(__CLPK_doublecomplex) );
	
	
	__CLPK_integer info=0;
	
	for(int i=0;i<n;i++)
	{
		lower[1].r=1;
		lower[1].i=0;

		middle[1].r=1;
		middle[1].i=0;

		upp[1].r=1;
		upp[1].i=0;

		b[1].r=1;
		b[1].i=1;

	}
	
	for(int h=1;h<8000000;h++)
	{
		for(int i=0;i<n;i++)
		{
			lower[1].r=1;
			lower[1].i=0;
			
			middle[1].r=1;
			middle[1].i=0;
			
			upp[1].r=1;
			upp[1].i=0;
			
			b[1].r=1;
			b[1].i=1;
			
		}
		zgtsv_(&n,&n1,lower, middle, upp, b,&n, &info);
	}
//	zheev_(&eigens, &upper, &N, H2prime , &LDA, W, Work, &lwork, rwork, &res);
	
	std::cout << "res = " << res << "\n";
	std::cout << "info = " << info << "\n";
}
