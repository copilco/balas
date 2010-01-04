#include <iostream>
#include <math.h>
#include "grid.h"
#include "wavef3.h"
#include <gsl/gsl_sf_bessel.h>


int main()
{
	
	int N=250;
	double R=4;
	int ord=0;
	
	printf("N=%d\n",N);
	
	
	FILE *in0;
	in0=fopen("besselCeros.txt","r");
	
	FILE *out0,*out1,*out2,*out3,*out4;
	out0=fopen("out0.txt","w");
	out1=fopen("out1.txt","w");
	out2=fopen("out2.txt","w");
	out3=fopen("out3.txt","w");
	out4=fopen("out4.txt","w");
	
	int zer=5*3001;
	int nn=5;//3001;
	int mm=3001;
	
	vector<double> zeros;
	zeros.resize(zer,0.);
	
	
	for(int j=0;j<mm;j++)
	{
		for(int i=0;i<nn;i++)
		{
			double f;
			fscanf(in0,"%lf",&f);
			zeros[j*nn+i]=f;
		}
		fscanf(in0,"\n");
	}
	
	
	for(int j=0;j<1;j++)
		for(int i=0;i<5;i++)
			printf(" %e ",zeros[j*nn+i]);
	
	
	vector<double> c;
	c.resize(N+1,0.);
	
	
	for(int j=0;j<N+1;j++)
	{
		c[j]=zeros[j*nn+(ord)];
		fprintf(out0,"%e ",c[j]);
	}
	
	double	V = c[N]/(2*pi*R);    //% Maximum frequency
	printf("\n V %e",V);
	
	vector<double> r;
	r.resize(N,0.);
	
	vector<double> v;
	v.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		r[i] = c[i]*R/c[N];   //% Radius vector
		v[i] = c[i]/(2*pi*R);
		fprintf(out1,"%e %e\n",r[i],v[i]);
	}

	vector<double> jnm;
	jnm.resize(N*N,0.);
	vector<double> C;
	C.resize(N*N,0.);
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			jnm[j*N+i]=c[i];
			
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			//jnm[j*N+i]=c[i];
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N])/abs(gsl_sf_bessel_J1(jnm[i*N+j]))/abs(gsl_sf_bessel_J1(jnm[j*N+i]));
			//fprintf(out2,"%e\n",C[j*N+i]);
		}
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			fprintf(out2,"%e\n",C[j*N+i]);
		
	
	vector<double> m1;
	m1.resize(N,0.);
	vector<double> m2;
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=abs(gsl_sf_bessel_J1(c[i]))/R;
		m2[i]=m1[i]*R/V;
	}
	
	/********************/
	//Transformation
	
	vector<double> f;
	f.resize(N,0.);
	vector<double> F;
	F.resize(N,0.);
	vector<double> f2;
	f2.resize(N,0.);
	vector<double> F2;
	F2.resize(N,0.);
	vector<double> Fretrieved;
	Fretrieved.resize(N,0.);
	vector<double> fretrieved;
	fretrieved.resize(N,0.);
	
	
	//Define the function
	for(int i=0;i<N;i++)
	{
		if(r[i]<1.)
			f[i]=1.;
		else
			f[i]=0.;
		
		F[i]=f[i]/m1[i];
	}
	
	//Apply the matrix  Forward Hankel
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			F2[j]+=C[j*N+i]*F[i];
	
	//Apply the matrix  Forward Hankel
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			Fretrieved[j]+=C[j*N+i]*F2[i];
	
	
	
	//Fretrieved = C*F2;           %% Inverse hankel transform
	
	//fretrieved = Fretrieved.*m1
	
	
	for(int i=0;i<N;i++)
	{
		f2[i]=F2[i]*m2[i];
		fprintf(out3,"%e %e %e\n",v[i],F2[i],f2[i]);
		fretrieved[i]=Fretrieved[i]*m1[i];
		fprintf(out4,"%e %e %e %e\n",r[i],f[i],fretrieved[i]);
	}
	
	
	/*
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
	
	field I,kernel,temp,temp2,temp3,temp4; 
	I.put_on_grid(g1);
	kernel.put_on_grid(g1);
	temp.put_on_grid(g1);
	temp2.put_on_grid(g1);
	temp3.put_on_grid(g1);
	temp4.put_on_grid(g1);
	
	field k,r,h,hr,H,Hk,hh;
	k.put_on_grid(g2);
	r.put_on_grid(g2);
	h.put_on_grid(g2);
	hr.put_on_grid(g2);
	H.put_on_grid(g2);
	Hk.put_on_grid(g2);
	hh.put_on_grid(g2);
	
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
		
		temp3.w[i][0]=0.;
		temp3.w[i][1]=0.;
		
		temp4.w[i][0]=0.;
		temp4.w[i][1]=0.;
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
		
		Hk.w[i][0]=0.;
		Hk.w[i][1]=0.;
		
		h.w[i][0]=0.;
		h.w[i][1]=0.;
		
	}
	
	for(int i=0;i<N;i++)
	{
		I.w[i][0]   =exp(a*i);
		kernel.w[i][0]=a*ko*ro*I.w[i][0]*gsl_sf_bessel_J0 (ko*ro*I.w[i][0]);
		//kernel.w[i][0]=gsl_sf_bessel_J0 (ko*ro*I.w[i][0]);
	}
	
	fftw_execute(kernel.p3DB);
	
	for(int i=0;i<N/2;i++)
	{
	    k.w[i][0]=ko*I.w[i][0];                        // samplings
		r.w[i][0]=ro*I.w[i][0];
		h.w[i][0]=2.*r.w[i][0]*exp(-r.w[i][0]*r.w[i][0] );
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
	
	for(int i=0;i<N/2;i++)
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
	
	
	
	for(int i=0;i<N/2;i++)
	{
		fprintf(out3,"%e %e %e\n",r.w[i][0],h.w[i][0],h.w[i][1]);
		fprintf(out4,"%e %e %e\n",k.w[i][0],H.w[i][0],H.w[i][1]);
		
	}
	
	 */
	
	/*
	 
	 //Backwards
	
	
	for(int i=0;i<N/2;i++)
	{
	    //k.w[i][0]=ko*I.w[i][0];                        // samplings
		//r.w[i][0]=ro*I.w[i][0];
		//h.w[i][0]=exp(-r.w[i][0]*r.w[i][0] );
		Hk.w[i][0]=H.w[i][0]*k.w[i][0];
		temp3.w[i][0]=Hk.w[i][0];
	}
	
	fftw_execute(temp3.p3DF);
	
	for(int i=0;i<N;i++)
	{
		temp4.w[i][0]=temp3.w[i][0]*kernel.w[i][0]-temp3.w[i][1]*kernel.w[i][1];
		temp4.w[i][1]=temp3.w[i][0]*kernel.w[i][1]+temp3.w[i][1]*kernel.w[i][0];
	}
	
	fftw_execute(temp4.p3DF);
	
	for(int i=0;i<N/2;i++)
	{
		
		hh.w[i][0]=temp4.w[i][0]/2.*pi/r.w[i][0]/N;
		hh.w[i][1]=temp4.w[i][1]/2.*pi/r.w[i][0]/N;
	}
	
	
	for(int i=0;i<N/2;i++)
	{
		fprintf(out5,"%e %e %e\n",r.w[i][0],hh.w[i][0],hh.w[i][1]);
	}
	
	*/
	
	
	
}
