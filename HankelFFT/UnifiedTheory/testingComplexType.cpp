/*
 *  testingComplexType.cpp
 *  
 *
 *  Created by Camilo Ruiz Méndez on 03/01/10.
 *  Copyright 2010 USAL. All rights reserved.
 *
 */


/*
#include <iostream.h>
//#include <gsl/gsl_cblas.h>
//#include <Accelerate/Accelerate.h> //contains data types used
#include <math.h>
//#include "HankelMatrix.h"
#include <complex.h>
#include <fftw3.h>

int main()
{

	
	
	fftw_complex I;
	
	I[0]=0.;
	I[1]=1.;
	
	fftw_complex y;
	y[0]=5.;
	y[1]=3.;
	
	cout << y;
	fftw_complex x;
	x= y * I;//(3+4*I)
	
	
}
*/

#include <iostream.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
int main()
{
	int N=3;
	int K=4;
	int i, j;
	const double TWOPI = 6.2831853071795864769252867665590057683943388;
	fftw_complex in[N], out[N], twids[(K-1)*(N/K-1)];
	fftw_plan plan;
	fftw_complex I;
	I[0]=0.;
	I[1]=1.;
	/* plan N/K FFTs of size K */
//	plan = fftw_plan_many_dft(1, &K, N/K,
//							  in, NULL, N/K, 1, 
//							  out, NULL, 1, K,
//							  FFTW_FORWARD, FFTW_ESTIMATE);
	
	/* precompute twiddle factors (since we usually want more than one FFT) */
	for (j = 1; j < N/K; ++j)
		for (i = 1; i < K; ++i)
			twids[(j-1)*(K-1) + (i-1)] = cexp((I * FFTW_FORWARD * TWOPI/N) * (i*j));
	
			/*
			...initialize in[N] data....
			
			fftw_execute(plan);
			
			for (j = 1; j < N/K; ++j) {
				out[0] += out[j*K];
				for (i = 1; i < K; ++i)
					out[i] += out[i + j*K] * twids[(j-1)*(K-1) + (i-1)];
					}
	*/
	
	
	fftw_destroy_plan(plan);
}