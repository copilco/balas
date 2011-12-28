/*
 *  testing.cpp
 *  
 *
 *  Created by Camilo Ruiz Méndez on 18/12/09.
 *  Copyright 2009 USAL. All rights reserved.
 *
 */


#include <fftw3.h>
#include <stdio.h>

int main()
{
	FILE *out0;
	out0=fopen("out0.txt","w");
	
	int N=100;
	fftw_complex *in, *out;
	
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N* N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N* N);
	
	fftw_plan p[N];
	
	for(int i=0;i<N;i++)
	{
		p[i] = fftw_plan_dft_1d(N, in+i*N, out+i*N, FFTW_FORWARD, FFTW_ESTIMATE);
	}
	
	
	for(int i=0;i<N;i++)
		fftw_execute(p[i]); /* repeat as needed */
	
		
		
	for (int j=1; j<N; j++)
		for(int i=0;i<N;i++)
		{
			fprintf(out0,"%e\n",in[j*N+i][0]*in[j*N+i][0]+in[j*N+i][1]*in[j*N+i][1]);
		}
	
	
	for(int i=0;i<N;i++)
		fftw_destroy_plan(p[i]);
	fftw_free(in); fftw_free(out);
}