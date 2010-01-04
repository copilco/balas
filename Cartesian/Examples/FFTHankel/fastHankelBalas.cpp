#include <iostream>
#include <math.h>
#include "grid.h"
#include "wavef3.h"
#include <gsl/gsl_sf_bessel.h>


int main()
{
	
	int M=7;
	int m=5;
	
	double K=10;//11;
	double R=5;//2.25;
	double a=K*R/2/pi;
	
	int N=int(pow(2.,ceil(1 + log2(a*m*log(a*M)))));
	printf("N=%d",N);
	
	a=1/a/m;
	printf("a=%e",a);
	double ro=R*exp(-a*N/2);
	double ko=K*exp(-a*N/2);

	
	int nx=1;
	int ny=N;
	int nz=1;
	
	double dx=1.;
	double dy=1.;
	double dz=1.;
	
	grid g1,g2;
	
	g1.set_grid(N,1,1,1.,1.,1.);
	g2.set_grid(N/2,1,1,1.,1.,1.);
	
	field I,kernel,temp,temp2; 
	I.put_on_grid(g1);
	kernel.put_on_grid(g1);
	temp.put_on_grid(g1);
	temp2.put_on_grid(g1);
	
	field k,r,h,hr,H;
	k.put_on_grid(g2);
	r.put_on_grid(g2);
	h.put_on_grid(g2);
	hr.put_on_grid(g2);
	H.put_on_grid(g2);
	
	for(int i=0;i<N;i++)
	{
		I.w[i][0]=0.;
		I.w[i][1]=0.;
		
		kernel.w[i][0]=0.;
		kernel.w[i][1]=0.;
		
		temp.w[i][0]=0.;
		temp.w[i][1]=0.;
		
		temp2.w[i][0]=0.;
		temp2.w[i][1]=0.;
	}
	
	for(int i=0;i<N/2;i++)
	{
		k.w[i][0]=0.;
		k.w[i][1]=0.;
		
		r.w[i][0]=0.;
		r.w[i][1]=0.;
		
		h.w[i][0]=0.;
		h.w[i][1]=0.;
		
		hr.w[i][0]=0.;
		hr.w[i][1]=0.;
		
		H.w[i][0]=0.;
		H.w[i][1]=0.;
	}
	
	for(int i=0;i<N;i++)
	{
		I.w[i][0]   =exp(a*i);
		kernel.w[i][0]=a*ko*ro*I.w[i][0]*gsl_sf_bessel_J0 (ko*ro*I.w[i][0]);
		//kernel.w[i][0]=gsl_sf_bessel_J0 (ko*ro*I.w[i][0]);
	}
	
	fftw_execute(kernel.p3DB);
	
	for(int i=0;i<=N/2;i++)
	{
	    k.w[i][0]=ko*I.w[i][0];                        // samplings
		r.w[i][0]=ro*I.w[i][0];
		h.w[i][0]=exp(-r.w[i][0]*r.w[i][0] );
		hr.w[i][0]=h.w[i][0]*r.w[i][0];
		temp.w[i][0]=hr.w[i][0];
	}
	
	fftw_execute(temp.p3DF);
	
	for(int i=0;i<N;i++)
	{
		temp2.w[i][0]=temp.w[i][0]*kernel.w[i][0]-temp.w[i][1]*kernel.w[i][1];
		temp2.w[i][1]=temp.w[i][0]*kernel.w[i][1]+temp.w[i][1]*kernel.w[i][0];
	}
	
	
	fftw_execute(temp2.p3DF);
	
	for(int i=0;i<=N/2;i++)
	{
		H.w[i][0]=2.*pi*temp2.w[i][0]/k.w[i][0]/N;
		H.w[i][1]=2.*pi*temp2.w[i][1]/k.w[i][0]/N;
	}
	
	
	FILE *out0,*out1,*out2;
	out0=fopen("out0.txt","w");
	out1=fopen("out1.txt","w");
	out2=fopen("out2.txt","w");
	
	FILE *out3,*out4,*out5;
	out3=fopen("out3.txt","w");
	out4=fopen("out4.txt","w");
	out5=fopen("out5.txt","w");
	
	
	for(int i=0;i<N;i++)
	{
		fprintf(out0,"%e %e\n",kernel.w[i][0]/sqrt(N*N),kernel.w[i][1]/sqrt(N*N) );
		fprintf(out1,"%e %e\n",temp.w[i][0],temp.w[i][1]);
		fprintf(out2,"%e %e\n",temp2.w[i][0],temp2.w[i][1]);
	}
	
	
	
	for(int i=0;i<=N/2;i++)
	{
		fprintf(out3,"%e %e %e\n",r.w[i][0],h.w[i][0],h.w[i][1]);
		fprintf(out4,"%e %e %e\n",k.w[i][0],H.w[i][0],H.w[i][1]);
		
	}
	
	
	
	
	
}
