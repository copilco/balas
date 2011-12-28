#include <stdio.h>
//#include <gsl/gsl_cblas.h>
#include <Accelerate/Accelerate.h>     
#include "grid.h"
#include "wavef3.h"
int
main (void)
{
 
  int lda = 3;
  
	double A[] = { 1., 0., 0., 0., 1.,0.};
	//0.11, 0.12, 0.13,
	//	0.21, 0.22, 0.23 };
  
  int ldb = 2;
  
	double B[] = {1., 1., 1., 1., 0., 0.};//{ 1011, 1012,
		//1021, 1022,
		//1031, 1032 };
     
  int ldc = 2;
  
  double C[] = { 0.00, 0.00,
		0.00, 0.00 };

  grid g1,g2;
  
  g1.set_grid(6,1,1,1.,1.,1.);
  g2.set_grid(4,1,1,1.,1.,1.);
  
  field aa,bb,cc; 
  aa.put_on_grid(g1);
  bb.put_on_grid(g1);
  cc.put_on_grid(g2);
 
  for(int i=0;i<6;i++)
    {
        aa.w[i][0]=A[i];
		bb.w[i][0]=B[i];
		
		aa.w[i][1]=0.;//A[i];
		bb.w[i][1]=0.;//1.;
    }

	aa.w[0][1]=0.;
	bb.w[0][1]=0.;
	
  for(int i=0;i<4;i++)
    cc.w[i][0]=C[i];

  
  /* Compute C = A B */
  
  cblas_dgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, 2, 2, 3,
	       1.0, A, lda, B, ldb, 0.0, C, ldc);

  cblas_dgemm (CblasRowMajor, 
	       CblasNoTrans, CblasNoTrans, 2, 2, 3,
	       1.0, aa.w[0], lda, bb.w[0],ldb,0.0, cc.w[0], ldc);
  
  printf ("[ %g, %g\n",   C[0], C[1]);
  printf ("  %g, %g ]\n", C[2], C[3]);

  printf ("[ %g, %g\n",   cc.w[0][0], cc.w[1][0]);
  printf ("  %g, %g ]\n", cc.w[2][0], cc.w[3][0]);
	
	printf ("[ %g, %g\n",   cc.w[0][1], cc.w[1][1]);
	printf ("  %g, %g ]\n", cc.w[2][1], cc.w[3][1]);
 

  return 0;  
}
