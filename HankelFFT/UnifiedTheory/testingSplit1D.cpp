/*
 *  testingSplit.cpp
 *  
 *
 *  Created by Camilo Ruiz Méndez on 03/01/10.
 *  Copyright 2010 USAL. All rights reserved.
 *
 */


#include <iostream>
//#include <gsl/gsl_cblas.h>
#include <Accelerate/Accelerate.h> //contains data types used
#include <math.h>
#include "HankelMatrix.h"

#include <complex>
#define complex complex<double>

int main()
{
	
	
	int Nt=600;//50;
		
	double dt=0.1;
	double nisq=pi/dt;
	double dw=dospi/Nt/dt;
	double wfactor=dt/sqrt(dospi);
	
	double *t;
	t = (double*) malloc(sizeof(double) * Nt);
	
	double *w;
	w = (double*) malloc(sizeof(double) * Nt);
	
	for(int i=0;i<Nt;i++)
        t[i]=(-Nt/2.+i)*dt;
	
	for(int i=0;i<Nt/2;i++)
        w[i]=i*dw;
	
	for(int i=Nt/2;i<Nt;i++)
        w[i]=-nisq+(i-Nt/2)*dw;
	
	
	FILE *in0,*in1;
	in0=fopen("besselCeros.txt","r");
	in1=fopen("InputFunction.txt","r");
	
	FILE *out0,*out1,*out2,*out3,*out4,*out5;
	FILE *sout0;
	
	sout0=fopen("spout0.txt","w");
	out0=fopen("pout0.txt","w");
	out1=fopen("pout1.txt","w");
	out2=fopen("pout2.txt","w");
	out3=fopen("pout3.txt","w");
	out4=fopen("pout4.txt","w");
	out5=fopen("pout5.txt","w");
	
	
	arrai F;
	F.init(1,Nt,0.);
	
	fftw_plan pf;
		pf = fftw_plan_dft_1d(Nt, F.v, F.v, FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_plan pb;
		pb = fftw_plan_dft_1d(Nt, F.v, F.v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	
	double K=6.4387e+04;
	double z_max = 25;		//%	Maximum propagation distance
	int Nz = 400;//int(z_max/dz); 
	double dz = z_max/(Nz-1);
	
	for(int i=0;i<Nt;i++)
	{			
		F.v[i][0]=exp(-t[i]*t[i]/40.);// /HH.m1.v[j][0];  // exp(-(HH.r.v[j][0]*HH.r.v[j][0]/w0/w0)-t[i]*t[i])/HH.m1.v[j][0];  //f.v[j][0]/HH.m1.v[j][0];
		F.v[i][1]=0.;//f.v[j][1]/HH.m1.v[j][0];
	}
	
	
	
	
	///////////////////////////////////////////////////////////////////////////
	//Loop
	///////////////////////////////////////////////////////////////////////////
	
	for (int gg=0;gg<Nz;gg++)
	{
		
		printf("gg =%d",gg);
		
		fftw_execute(pf);
		
		for(int i=0;i<Nt;i++)
		{				
			double aux0r=F.v[i][0];
			double aux0i=F.v[i][1];
			
			double psi=w[i]*w[i]*dz;
			
			F.v[i][0]=( aux0r*cos(psi)-aux0i*sin(psi) )/Nt;
			F.v[i][1]=( aux0i*cos(psi)+aux0r*sin(psi) )/Nt;
		}
		
		fftw_execute(pb);		
		
		//Print the output of the function
		int skiper=1;
		if( (gg%1)==0)
		{
			
			for(int i=0;i<Nt;i++)
				fprintf(sout0,"%e\n",F.v[i][0]*F.v[i][0]+F.v[i][1]*F.v[i][1]);
			//fprintf(sout0,"%e\n",Fretrieved.v[j*Nt+i][0]*Fretrieved.v[j*Nt+i][0]+Fretrieved.v[j*Nt+i][1]*Fretrieved.v[j*Nt+i][1]);
			//fprintf(sout0,"%e\n",F2.v[j*Nt+i][0]*F2.v[j*Nt+i][0]+F2.v[j*Nt+i][1]*F2.v[j*Nt+i][1]);
			//fprintf(sout0,"%d %e %e\n",j*Nt+i,F2.v[j*Nt+i][0],F2.v[j*Nt+i][1]);			
			
		}
		
				
	}//End loop
	
	
	
	
}