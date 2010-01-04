/*
 *  HankelMatrix.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 01/03/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
//#include <grid.h>
#include "fftw3.h"
#include <math.h>
#include <gsl/gsl_sf_bessel.h>



class HankelMatrix: public grid
	{
	public:
		
		int N;
		double R;
		int ord;

		int nn;//=5;
		int mm;//=3001;
		int zer;//=nn*mm;
		
		vector<double> zeros;
		vector<double> c;
		double	V;
		vector<double> r;
		vector<double> v;
		vector<double> jnm;
		vector<double> C;
		vector<double> m1;
		vector<double> m2;
		
		FILE *in0;
		
		HankelMatrix ();
		HankelMatrix(int _N,double R);
		HankelMatrix(grid g);
	};

HankelMatrix::HankelMatrix ()
{
	N=100;
	R=1.;
	ord=0.;
	
	nn=5;
	mm=3001;
	zer=nn*mm;
	
	//Read the zeros of the Bessel function
	in0=fopen("besselCeros.txt","r");
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
	
	//Fill the zeros onto the c vector
	c.resize(N+1,0.);
	for(int j=0;j<N+1;j++)
		c[j]=zeros[j*nn+(ord)];
	
	double	V = c[N]/(2*pi*R);    //% Maximum frequency
	
	r.resize(N,0.);	
	v.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		r[i] = c[i]*R/c[N];   //% Radius vector
		v[i] = c[i]/(2*pi*R);
	}
	
	jnm.resize(N*N,0.);
	C.resize(N*N,0.);
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			jnm[j*N+i]=c[i];
	
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N])/abs(gsl_sf_bessel_J1(jnm[i*N+j]))/abs(gsl_sf_bessel_J1(jnm[j*N+i]));
		}
	
	m1.resize(N,0.);
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=abs(gsl_sf_bessel_J1(c[i]))/R;
		m2[i]=m1[i]*R/V;
	}
	
}


HankelMatrix::HankelMatrix(int _N,double _R)
{
	N=_N;
	R=_R;
	ord=0.;
	
	nn=5;
	mm=3001;
	zer=nn*mm;
	
	//Read the zeros of the Bessel function
	in0=fopen("besselCeros.txt","r");
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
	
	//Fill the zeros onto the c vector
	c.resize(N+1,0.);
	for(int j=0;j<N+1;j++)
		c[j]=zeros[j*nn+(ord)];
	
	double	V = c[N]/(2*pi*R);    //% Maximum frequency
	
	r.resize(N,0.);	
	v.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		r[i] = c[i]*R/c[N];   //% Radius vector
		v[i] = c[i]/(2*pi*R);
	}
	
	jnm.resize(N*N,0.);
	C.resize(N*N,0.);
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			jnm[j*N+i]=c[i];
	
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N])/abs(gsl_sf_bessel_J1(jnm[i*N+j]))/abs(gsl_sf_bessel_J1(jnm[j*N+i]));
		}
	
	m1.resize(N,0.);
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=abs(gsl_sf_bessel_J1(c[i]))/R;
		m2[i]=m1[i]*R/V;
	}
	
}


HankelMatrix::HankelMatrix(grid g)
{
	N=g.n1;
	R=g.dx1*(g.n1-1);
	ord=0.;
	
	nn=5;
	mm=3001;
	zer=nn*mm;
	
	//Read the zeros of the Bessel function
	in0=fopen("besselCeros.txt","r");
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
	
	//Fill the zeros onto the c vector
	c.resize(N+1,0.);
	for(int j=0;j<N+1;j++)
		c[j]=zeros[j*nn+(ord)];
	
	double	V = c[N]/(2*pi*R);    //% Maximum frequency
	
	r.resize(N,0.);	
	v.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		r[i] = c[i]*R/c[N];   //% Radius vector
		v[i] = c[i]/(2*pi*R);
	}
	
	jnm.resize(N*N,0.);
	C.resize(N*N,0.);
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			jnm[j*N+i]=c[i];
	
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N])/abs(gsl_sf_bessel_J1(jnm[i*N+j]))/abs(gsl_sf_bessel_J1(jnm[j*N+i]));
		}
	
	m1.resize(N,0.);
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=abs(gsl_sf_bessel_J1(c[i]))/R;
		m2[i]=m1[i]*R/V;
	}
	
}



