/*
 *  interp.h
 *  
 *
 *  Created by camilo on 30/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "wave.h"
#include "waveUniform.h"

void interpH2U(wave wHank, waveUniform wU)
{
	// Parameters for interpolate
	double sig;
	complex p,qn,un;
	complex *s=new complex[wHank.Nr];
	complex *u=new complex[wHank.Nr-1];
	
	// Parameters for evaluate the polynom
	int klo,khi,k; 
	double h,b,a;
	
	
	////////////////////////////
	// Create the interpolant //
	////////////////////////////
	
	s[1] =complex(-0.5,-0.5);
	u[1]=(3.0/(wHank.r[2]-wHank.r[1]))*((wHank.phiHank[2]-wHank.phiHank[1])/(wHank.r[2]-wHank.r[1]));
	
	for (int i=2;i<=wHank.Nr-1;i++)
	{	
		sig=(wHank.r[i]-wHank.r[i-1])/(wHank.r[i+1]-wHank.r[i-1]);
		p=sig*s[i-1]+2.0;
		s[i]=(sig-1.0)/p;
		u[i]=(wHank.phiHank[i+1]-wHank.phiHank[i])/(wHank.r[i+1]-wHank.r[i]) - (wHank.phiHank[i]-wHank.phiHank[i-1])/(wHank.r[i]-wHank.r[i-1]);
		u[i]=(6.0*u[i]/(wHank.r[i+1]-wHank.r[i-1])-sig*u[i-1])/p;
	}
	
	qn=complex(0.5,0.5);
	un=(3.0/(wHank.r[wHank.Nr]-wHank.r[wHank.Nr-1]))*(-(wHank.phiHank[wHank.Nr]-wHank.phiHank[wHank.Nr-1])/(wHank.r[wHank.Nr]-wHank.r[wHank.Nr-1]));
	
	s[wHank.Nr]=(un-qn*u[wHank.Nr-1])/(qn*s[wHank.Nr-1]+1.0);
	
	for (int k=wHank.Nr-1;k>=1;k--)	
	{
		s[k]=s[k]*s[k+1]+u[k];
	}
	
	////////////////////////////
	// Evaluate the polynom   //
	////////////////////////////
	
	for(int i=0;i<wU.Nr;i++)
	{
		klo=1;
		khi=wHank.Nr; 
		
		while (khi-klo > 1)
		{
			k=(khi+klo) >> 1;
			if (wHank.r[k] > wU.r[i])
				khi=k;
			else klo=k;
		}
		
		h=wHank.r[khi]-wHank.r[klo];
		if (h == 0.0) 
			nrerror ("Bad x array input");
		
		a=(wHank.r[khi]-wU.r[i])/h;
		b=(wU.r[i]-wHank.r[klo])/h;
		wU.phi[i]=a*wHank.phiHank[klo]+b*wHank.phiHank[khi]+((a*a*a-a)*s[klo]+(b*b*b-b)*s[khi])*(h*h)/6.0;
		
	}
	
	
	delete[] s;
	
}



void interpU2H(waveUniform *wU, wave *wHank)
{
	// Parameters for interpolate
	double sig;
	complex p,qn,un;
	complex *s=new complex[wU->Nr];
	complex *u=new complex[wU->Nr-1];
	
	// Parameters for evaluate the polynom
	int klo,khi,k; 
	double h,b,a;
	
	
	////////////////////////////
	// Create the interpolant //
	////////////////////////////
	
	s[0] =complex(-0.5,-0.5);
	u[0]=(3.0/(wU->r[1]-wU->r[0]))*((wU->phi[1]-wU->phi[0])/(wU->r[1]-wU->r[0]));
	
	for (int i=1;i<=wU->Nr-2;i++)
	{	
		sig=(wU->r[i]-wU->r[i-1])/(wU->r[i+1]-wU->r[i-1]);
		p=sig*s[i-1]+2.0;
		s[i]=(sig-1.0)/p;
		u[i]=(wU->phi[i+1]-wU->phi[i])/(wU->r[i+1]-wU->r[i]) - (wU->phi[i]-wU->phi[i-1])/(wU->r[i]-wU->r[i-1]);
		u[i]=(6.0*u[i]/(wU->r[i+1]-wU->r[i-1])-sig*u[i-1])/p;
	}
	
	qn=complex(0.5,0.5);
	un=(3.0/(wU->r[wU->Nr-1]-wU->r[wU->Nr-2]))*(-(wU->phi[wU->Nr-1]-wU->phi[wU->Nr-2])/(wU->r[wU->Nr-1]-wU->r[wU->Nr-2]));
	
	s[wU->Nr-1]=(un-qn*u[wU->Nr-2])/(qn*s[wU->Nr-2]+1.0);
	
	for (int k=wU->Nr-2;k>=0;k--)	
	{
		s[k]=s[k]*s[k+1]+u[k];
	}
	
	////////////////////////////
	// Evaluate the polynom   //
	////////////////////////////
	
	for(int i=0;i<wU->Nr;i++)
	{
		klo=1;
		khi=wU->Nr; 
		
		while (khi-klo > 1)
		{
			k=(khi+klo) >> 1;
			if (wU->r[k] > wHank->r[i])
				khi=k;
			else klo=k;
		}
		
		h=wU->r[khi]-wU->r[klo];
		if (h == 0.0) 
			nrerror ("Bad x array input");
		
		a=(wU->r[khi]-wHank->r[i])/h;
		b=(wHank->r[i]-wU->r[klo])/h;
		wHank->phiHank[i]=a*wU->phi[klo]+b*wU->phi[khi]+((a*a*a-a)*s[klo]+(b*b*b-b)*s[khi])*(h*h)/6.0;
		
	}
	
	
	delete[] s;
	
}

