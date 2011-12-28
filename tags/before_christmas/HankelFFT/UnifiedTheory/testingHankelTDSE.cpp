/*
 *  testingHankelMatCreation.cpp
 *  
 *
 *  Created by Camilo Ruiz Méndez on 01/01/10.
 *  Copyright 2010 USAL. All rights reserved.
 *
 */


#include <iostream>
//#include <gsl/gsl_cblas.h>
#include <Accelerate/Accelerate.h> //contains data types used
#include <math.h>
#include "HankelMatrix.h"
//#include <gsl/gsl_sf_bessel.h>
//#include <complex>
//#define complex complex<double>

//using std;
int main()
{
	
	
	//complex I=complex(0.,1.);
	int Nr=200;//50;
	int Nz=400;//50;
	
	double R=30.;
	printf("Nr=%d\n",Nr);
		
	HankelMatrix HH(Nr,R);
	
	double dz=0.4;
	double nisq=pi/dz;
	double dw=dospi/Nz/dz;
	double wfactor=dz/sqrt(dospi);
	
	double *z;
	z = (double*) malloc(sizeof(double) * Nz);

	double *w;
	w = (double*) malloc(sizeof(double) * Nz);

	for(int i=0;i<Nz;i++)
        z[i]=(-Nz/2.+i)*dz;
	
	for(int i=0;i<Nz/2;i++)
        w[i]=i*dw;
	
	for(int i=Nz/2;i<Nz;i++)
        w[i]=-nisq+(i-Nz/2)*dw;
	
	
	FILE *out0,*out1,*out2,*out3,*out4,*out5;
	FILE *sout0;
	
	sout0=fopen("spout0.txt","w");
	out0=fopen("pout0.txt","w");
	out1=fopen("pout1.txt","w");
	out2=fopen("pout2.txt","w");
	out3=fopen("pout3.txt","w");
	out4=fopen("pout4.txt","w");
	out5=fopen("pout5.txt","w");
	
	
	//Transformation

	arrai F(Nr,Nz,0.);
	
	arrai f2(Nr,Nz,0.);
	arrai F2(Nr,Nz,0.);
	
	arrai Fretrieved(Nr,Nz,0.);
	arrai fretrieved(Nr,Nz,0.);
	
	arrai V(Nr,Nz,0.);
	
	fftw_plan pf[Nr];
	for(int j=0;j<Nr;j++)
		pf[j] = fftw_plan_dft_1d(Nz, F2.v+j*Nz, F2.v+j*Nz, FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_plan pb[Nr];
	for(int j=0;j<Nr;j++)
		pb[j] = fftw_plan_dft_1d(Nz, F2.v+j*Nz, F2.v+j*Nz, FFTW_BACKWARD, FFTW_ESTIMATE);
	
		
	double t_max = 1.;		//%	Maximum propagation distance
	int Nt = 600;//int(z_max/dz); 
	double dt =t_max/(Nt-1);
	
	
	///////////////////////////////////////////////////////////////////////////
	//Read the initial wave function
	///////////////////////////////////////////////////////////////////////////
	
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
		{			
			F.v[j*Nz+i][0]=exp(-z[i]*z[i]/9-HH.r.v[j][0]*HH.r.v[j][0]/9)/HH.m1.v[j][0];  // exp(-(HH.r.v[j][0]*HH.r.v[j][0]/w0/w0)-t[i]*t[i])/HH.m1.v[j][0];  //f.v[j][0]/HH.m1.v[j][0];
			F.v[j*Nz+i][1]=0.;//exp(-z[i]*z[i]-HH.r.v[j][0]*HH.r.v[j][0])*f.v[j][1]/HH.m1.v[j][0];
			V.v[j*Nz+i][0]=-4./sqrt(1.+z[i]*z[i]+HH.r.v[j][0]*HH.r.v[j][0]); 
		}
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			fprintf(out0,"%e %e\n",F.v[j*Nz+i][0],F.v[j*Nz+i][1]);	
	
	
	
	///////////////////////////////////////////////////////////////////////////
	//Read the initial wave function
	///////////////////////////////////////////////////////////////////////////
	
	
	__CLPK_doublecomplex alpha;
	__CLPK_doublecomplex beta;
	__CLPK_doublecomplex gamma;
	
	alpha.r=1;
	alpha.i=0;
	
	gamma.r=1;
	gamma.i=0;
	

	///////////////////////////////////////////////////////////////////////////
	//Loop
	///////////////////////////////////////////////////////////////////////////

	for (int gg=0;gg<Nz;gg++)
	{
		
		printf("gg =%d",gg);
		//Make F2 and Fretrieved variables to zero		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				F2.v[j*Nz+i][0]=0.;
				F2.v[j*Nz+i][1]=0.;
				
				Fretrieved.v[j*Nz+i][0]=0.;
				Fretrieved.v[j*Nz+i][1]=0.;
			}
		
	
		///////////////////////////////////////////////////////////////////////////
		//Forward
		///////////////////////////////////////////////////////////////////////////
		
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nz, Nr,
					 &alpha, HH.C.v, Nr,F.v,Nz,&gamma, F2.v, Nz);
		
		for(int j=0;j<Nr;j++)
			fftw_execute(pf[j]);
				
		//Multiply by the phase to make the propagation.
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{				
				double aux0r=F2.v[j*Nz+i][0];
				double aux0i=F2.v[j*Nz+i][1];
				double psi=(- 2*pi*HH.v.v[j][0]*2*pi*HH.v.v[j][0]-w[i]*w[i])*gg*dz;
				
				F2.v[j*Nz+i][0]=aux0r*exp(-psi)/Nz;//( aux0r*cos(psi)-aux0i*sin(psi) )/Nz;
				F2.v[j*Nz+i][1]=0.;//( aux0i*cos(psi)+aux0r*sin(psi) )/Nz;
			}
		
		
		for(int j=0;j<Nr;j++)
			fftw_execute(pb[j]);
		
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nz, Nr,
					 &alpha, HH.C.v, Nr,F2.v,Nz,&gamma, Fretrieved.v, Nz);
				
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{		
				f2.v[j*Nz+i][0]=F2.v[j*Nz+i][0]*HH.m2.v[j][0];
				f2.v[j*Nz+i][1]=F2.v[j*Nz+i][1]*HH.m2.v[j][0];
				
				fretrieved.v[j*Nz+i][0]=Fretrieved.v[j*Nz+i][0]*HH.m1.v[j][0];
				fretrieved.v[j*Nz+i][1]=Fretrieved.v[j*Nz+i][1]*HH.m1.v[j][0];
			}
		
		///////////////////////////////////////////////////////////////////////////
		//Forward
		///////////////////////////////////////////////////////////////////////////

		//Multiply by the phase to make the propagation.
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{				
				double aux0r=fretrieved.v[j*Nz+i][0];
				double aux0i=fretrieved.v[j*Nz+i][1];
				double psi=(V.v[j*Nz+i][0])*gg*dz;
				
				fretrieved.v[j*Nz+i][0]=aux0r*exp(-psi)/Nz;//( aux0r*cos(psi)-aux0i*sin(psi) )/Nz;
				fretrieved.v[j*Nz+i][1]=0.;//( aux0i*cos(psi)+aux0r*sin(psi) )/Nz;
			}
		
		
		
		///////////////////////////////////////////////////////////////////////////
		//Forward
		///////////////////////////////////////////////////////////////////////////
		
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nz, Nr,
					 &alpha, HH.C.v, Nr,F.v,Nz,&gamma, F2.v, Nz);
		
		for(int j=0;j<Nr;j++)
			fftw_execute(pf[j]);
		
		//Multiply by the phase to make the propagation.
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{				
				double aux0r=F2.v[j*Nz+i][0];
				double aux0i=F2.v[j*Nz+i][1];
				double psi=(- 2*pi*HH.v.v[j][0]*2*pi*HH.v.v[j][0]-w[i]*w[i])*gg*dz;
				
				F2.v[j*Nz+i][0]=aux0r*exp(psi)/Nz;//( aux0r*cos(psi)-aux0i*sin(psi) )/Nz;
				F2.v[j*Nz+i][1]=0.;//( aux0i*cos(psi)+aux0r*sin(psi) )/Nz;
			}
		
		
		for(int j=0;j<Nr;j++)
			fftw_execute(pb[j]);
		
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nz, Nr,
					 &alpha, HH.C.v, Nr,F2.v,Nz,&gamma, Fretrieved.v, Nz);
		
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{		
				f2.v[j*Nz+i][0]=F2.v[j*Nz+i][0]*HH.m2.v[j][0];
				f2.v[j*Nz+i][1]=F2.v[j*Nz+i][1]*HH.m2.v[j][0];
				
				fretrieved.v[j*Nz+i][0]=Fretrieved.v[j*Nz+i][0]*HH.m1.v[j][0];
				fretrieved.v[j*Nz+i][1]=Fretrieved.v[j*Nz+i][1]*HH.m1.v[j][0];
			}
		

		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				F.v[j*Nz+i][0]=0.;
				F.v[j*Nz+i][1]=0.;
			}
		
		
		//Print the output of the function
		int skiper=1;
		if( (gg%6)==0)
		{
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					fprintf(sout0,"%e\n",fretrieved.v[j*Nz+i][0]*fretrieved.v[j*Nz+i][0]+fretrieved.v[j*Nz+i][1]*fretrieved.v[j*Nz+i][1]);
		}
		
	}//End loop
	
}