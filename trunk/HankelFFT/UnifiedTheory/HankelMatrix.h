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
		
		int Nr;
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
		
		HankelMatrix(int _Nr,double R);	
	};


HankelMatrix::HankelMatrix(int _Nr,double _R)
{
	Nr=_Nr;
	R=_R;
	ord=0.; //For this application we will only use the zeros order solutions
	
	nn=5;
	mm=3001;
	zer=nn*mm;
	
	//Read the zeros of the Bessel function
	in0=fopen("besselCeros.txt","r");
	zeros.init(1,zer,0.);
	
	//Readt the zeros in the matrix called zeros, useful for higher order Transformations
	for(int j=0;j<mm;j++)
	{
		for(int i=0;i<nn;i++)
		{
			double f;
			fscanf(in0,"%lf",&f);
			zeros.v[j*nn+i][0]=f;
			//printf("-> %e\n",zeros[j*nn+i]);
		}
		fscanf(in0,"\n");

	}
	
	//Fill the zeros onto the c vector
	c.init(1,Nr+1,0.);
	for(int j=0;j<Nr+1;j++)
	{
		c.v[j][0]=zeros.v[j*nn+(ord)][0];
	}

	double	V = c.v[Nr][0]/(2*pi*R);    //% Maximum frequency
	
	//The axes for the rho and the krho coordinate
	r.init(1,Nr,0.);	
	v.init(1,Nr,0.);
	
	for(int i=0;i<Nr;i++)
	{
		r.v[i][0] = c.v[i][0]*R/c.v[Nr][0];   //% Radius vector
		v.v[i][0] = c.v[i][0]/(2*pi*R);
	}
	
	//The bessel funciotns and the matrix itself
	jnm.init(1,Nr*Nr,0.);
	C.init(1,Nr*Nr,0.);
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nr;i++)
			jnm.v[j*Nr+i][0]=c.v[i][0];
	
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nr;i++)
		{
			C.v[j*Nr+i][0]=(2./c.v[Nr][0])*j0(jnm.v[j*Nr+i][0]*jnm.v[i*Nr+j][0]/c.v[Nr][0]);
			C.v[j*Nr+i][0]=(2./c.v[Nr][0])*j0(jnm.v[j*Nr+i][0]*jnm.v[i*Nr+j][0]/c.v[Nr][0])/fabs(j1(jnm.v[i*Nr+j][0]))/fabs(j1(jnm.v[j*Nr+i][0]));
		}
	
	m1.init(1,Nr,0.);
	m2.init(1,Nr,0.);
	
	for(int i=0;i<Nr;i++)
	{
		m1.v[i][0]=fabs(j1(c.v[i][0]))/R;
		m2.v[i][0]=m1.v[i][0]*R/V;
		//printf("** %e\n",m1.v[i]);
	}
	
}






	


