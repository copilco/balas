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
//#include "fftw3.h"
#include <math.h>
#include <vector.h>
#include "constant.h"
//#include <gsl/gsl_sf_bessel.h>



class HankelMatrix
	{
	public:
		
		int N;
		double R;
		int ord;

		int nn;//=5;
		int mm;//=3001;
		int zer;//=nn*mm;
		
		double * zeros;
		double * c;
		double	V;
		double * r;
		double * v;
		double * jnm;
		double * C;
		double * m1;
		double * m2;
		
		FILE *in0;
		
		HankelMatrix ();
		HankelMatrix(int _N,double R);
		HankelMatrix(int _N,double R,int aa);
		//HankelMatrix(grid g);
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
	
	zeros = (double *)calloc( zer, sizeof( double ) );
	for (int i=0;i<zer;i++)
		zeros[i]=0.;
	
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
	c = (double *)calloc( N+1, sizeof( double ) );
	for (int i=0;i<N+1;i++)
		c[i]=0.;
	
	
	for(int j=0;j<N+1;j++)
		c[j]=zeros[j*nn+(ord)];
	
	double	V = c[N]/(2*pi*R);    //% Maximum frequency
	
	r = (double *)calloc( N, sizeof( double ) );
	v = (double *)calloc( N, sizeof( double ) );

	for (int i=0;i<N;i++)
	{
		r[i]=0.;
		v[i]=0.;
	}
	
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
			C[j*N+i]=(2./c[N])*j0(jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*j0(jnm[j*N+i]*jnm[i*N+j]/c[N])/fabs(j1(jnm[i*N+j]))/fabs(j1(jnm[j*N+i]));
		}
	
	m1.resize(N,0.);
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=fabs(j1(c[i]))/R;
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
//	in0=fopen("bessel0000.txt","r");
	in0=fopen("besselCeros.txt","r");
	zeros.resize(zer,0.);
	
	
	for(int j=0;j<mm;j++)
	{
		for(int i=0;i<nn;i++)
		{
			double f;
			fscanf(in0,"%lf",&f);
			zeros[j*nn+i]=f;
			//printf("-> %e\n",zeros[j*nn+i]);
		}
		fscanf(in0,"\n");

	}
	
	//Fill the zeros onto the c vector
	c.resize(N+1,0.);
	for(int j=0;j<N+1;j++)
	{
		c[j]=zeros[j*nn+(ord)];
		printf("-> %e\n",c[j]);
	}
		
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
			C[j*N+i]=(2./c[N])*j0(jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*j0(jnm[j*N+i]*jnm[i*N+j]/c[N])/fabs(j1(jnm[i*N+j]))/fabs(j1(jnm[j*N+i]));
		}
	
	m1.resize(N,0.);
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=fabs(j1(c[i]))/R;
		m2[i]=m1[i]*R/V;
		printf("** %e\n",m1[i]);
	}
	
}


HankelMatrix::HankelMatrix(int _N,double _R,int aa)
{
	N=_N;
	R=_R;
	ord=0.;
	
	nn=1;
	mm=10001;
	zer=nn*mm;
	
	//Read the zeros of the Bessel function
	in0=fopen("bessel0000.txt","r");
	zeros.resize(zer,0.);
	
	for(int j=0;j<mm;j++)
	{
		//for(int i=0;i<nn;i++)
		//{
			double f;
			fscanf(in0,"%lf\n",&f);
			zeros[j]=f;
		//}
		//fscanf(in0,"\n");
		//printf("read %d\n",j);
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
			C[j*N+i]=(2./c[N])*j0(jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*j0(jnm[j*N+i]*jnm[i*N+j]/c[N])/fabs(j1(jnm[i*N+j]))/fabs(j1(jnm[j*N+i]));
		}
	
	m1.resize(N,0.);
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=fabs(j1(c[i]))/R;
		m2[i]=m1[i]*R/V;
	}
	
}




	


