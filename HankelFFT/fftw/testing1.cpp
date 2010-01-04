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
#include <math.h>

int main()
{
	FILE *out0,*out1;
	out0=fopen("out0.txt","w");
	out1=fopen("out1.txt","w");
	
	int N=100;
	fftw_complex *in, *out;
	
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N* N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N* N);
	
	fftw_plan p[N];
	
	for(int i=0;i<N;i++)
	{
		p[i] = fftw_plan_dft_1d(N, in+i*N, out+i*N, FFTW_FORWARD, FFTW_ESTIMATE);
	}
	
	
	double dx=0.1;
	for (int j=0;j<N;j++)
	{
		double y=(-N/2+j)*dx;
		for(int i=0;i<N;i++)
		{
			double y=(-N/2+j)*dx;
			double x=(-N/2+i)*dx;
			in[j*N+i][0]=sin(x);//*exp(-y*y);
			in[j*N+i][1]=0.;
		}
	}
	
	for(int i=0;i<N;i++)
		fftw_execute(p[i]); /* repeat as needed */
	
		
		
	for (int j=0; j<N; j++)
		for(int i=0;i<N;i++)
		{
			fprintf(out0,"%e\n",in[j*N+i][0]*in[j*N+i][0]+in[j*N+i][1]*in[j*N+i][1]);
			fprintf(out1,"%e\n",out[j*N+i][0]*out[j*N+i][0]+out[j*N+i][1]*out[j*N+i][1]);
		}
	
	
	for(int i=0;i<N;i++)
		fftw_destroy_plan(p[i]);
	
	fftw_free(in); 
	fftw_free(out);
}