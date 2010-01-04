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
#include "arrai.h"
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
		
		arrai zeros;
		arrai c;
		double	V;
		arrai r;
		arrai v;
		arrai jnm;
		arrai C;
		arrai m1;
		arrai m2;
		
		FILE *in0;
		
		HankelMatrix(int _N,double R);	
	};


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
	zeros.init(zer,0.);
	
	
	for(int j=0;j<mm;j++)
	{
		for(int i=0;i<nn;i++)
		{
			double f;
			fscanf(in0,"%lf",&f);
			zeros.v[j*nn+i]=f;
			//printf("-> %e\n",zeros[j*nn+i]);
		}
		fscanf(in0,"\n");

	}
	
	//Fill the zeros onto the c vector
	c.init(N+1,0.);
	for(int j=0;j<N+1;j++)
	{
		c.v[j]=zeros.v[j*nn+(ord)];
		//printf("-> %e\n",c.v[j]);
	}
		
	double	V = c.v[N]/(2*pi*R);    //% Maximum frequency
	
	r.init(N,0.);	
	v.init(N,0.);
	
	for(int i=0;i<N;i++)
	{
		r.v[i] = c.v[i]*R/c.v[N];   //% Radius vector
		v.v[i] = c.v[i]/(2*pi*R);
	}
	
	jnm.init(N*N,0.);
	C.init(N*N,0.);
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			jnm.v[j*N+i]=c.v[i];
	
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			C.v[j*N+i]=(2./c.v[N])*j0(jnm.v[j*N+i]*jnm.v[i*N+j]/c.v[N]);
			C.v[j*N+i]=(2./c.v[N])*j0(jnm.v[j*N+i]*jnm.v[i*N+j]/c.v[N])/fabs(j1(jnm.v[i*N+j]))/fabs(j1(jnm.v[j*N+i]));
		}
	
	m1.init(N,0.);
	m2.init(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1.v[i]=fabs(j1(c.v[i]))/R;
		m2.v[i]=m1.v[i]*R/V;
		//printf("** %e\n",m1.v[i]);
	}
	
}






	


