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
//#include <complex.h>
//#define complex complex<double>

//using std;
int main()
{
	
	
	//complex I=complex(0.,1.);
	int Nr=1024;//50;
	int Nt=12;//50;
	
	int sampler=11;
	
	double R=.05;
	
	//Parametros para blas
	int lda=Nr;
	int ldb=1;
	int ldc=1;
	
	HankelMatrix HH(Nr,R);
	
	double dr = R/(Nr-1);	//%	Radial spacing	
	printf("Nr=%d\n",Nr);
	
	
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
	
	
	//Transformation
	
	//arrai phi;
	arrai f;
	arrai F;
	
	arrai f2;
	arrai F2;
	
	arrai Fretrieved;
	arrai fretrieved;
	
	//phi.init(N,0.);
	f.init(Nr,1,0.);
	F.init(Nr,Nt,0.);
	
	f2.init(Nr,Nt,0.);
	F2.init(Nr,Nt,0.);

	fretrieved.init(Nr,Nt,0.);
	Fretrieved.init(Nr,Nt,0.);
		
	double K=6.4387e+04;
	double z_max = .25;		//%	Maximum propagation distance
	int Nz = 200;//int(z_max/dz); 
	double dz = z_max/(Nz-1);
	
	arrai phi;
	phi.init(Nr,1,0.);
		
	
	///////////////////////////////////////////////////////////////////////////
	//Read the initial wave function
	///////////////////////////////////////////////////////////////////////////
	
	for(int j=0;j<Nr;j++)
	{
		//double r=i*dr;
		double read1,read2;
		fscanf(in1,"%lf %lf",&read1,&read2);
		
		f.v[j][0]= read1;  //exp(-(HH.r.v[i]*HH.r.v[i]/w0/w0))*cos(HH.v.v[*HH.v.v[i]*HH.r.v[i]*HH.r.v[i]/2./focus);
		f.v[j][1]= read2; // exp(-(HH.r.v[i]*HH.r.v[i]/w0/w0))*sin(HH.v.v[*HH.v.v[i]*HH.r.v[i]*HH.r.v[i]/2./focus);
	}
	
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nt;i++)
		{			
			F.v[j*Nt+i][0]=f.v[j][0]/HH.m1.v[j][0];
			F.v[j*Nt+i][1]=f.v[j][1]/HH.m1.v[j][0];
		}
	
	
	for(int j=0;j<Nr;j++)
		fprintf(out0,"%e %e %e %e\n",HH.r.v[j][0],HH.r.v[j][0],f.v[j][0],f.v[j][1]);
	
	
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
			for(int i=0;i<Nt;i++)
			{
				F2.v[j*Nt+i][0]=0.;
				F2.v[j*Nt+i][1]=0.;
				
				Fretrieved.v[j*Nt+i][0]=0.;
				Fretrieved.v[j*Nt+i][1]=0.;
			}
		
	
		///////////////////////////////////////////////////////////////////////////
		//Forward
		///////////////////////////////////////////////////////////////////////////
		
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,
					 &alpha, HH.C.v, Nr,F.v,Nt,&gamma, F2.v, Nt);
		

		//Multiply by the phase to make the propagation.
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nt;i++)
			{				
				double aux0r=F2.v[j*Nt+i][0];
				double aux0i=F2.v[j*Nt+i][1];
				
				phi.v[j][0]=(sqrt(K*K - 2*pi*HH.v.v[j][0]*2*pi*HH.v.v[j][0]) - K)*gg*dz; //The phase
				F2.v[j*Nt+i][0]=( aux0r*cos(phi.v[j][0])-aux0i*sin(phi.v[j][0]) );
				F2.v[j*Nt+i][1]=( aux0i*cos(phi.v[j][0])+aux0r*sin(phi.v[j][0]) );
			}
		
		
		///////////////////////////////////////////////////////////////////////////
		//Backward
		///////////////////////////////////////////////////////////////////////////
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,
					 &alpha, HH.C.v, Nr,F2.v,Nt,&gamma, Fretrieved.v, Nt);
		
		
		
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nt;i++)
			{		
				f2.v[j*Nt+i][0]=F2.v[j*Nt+i][0]*HH.m2.v[j][0];
				f2.v[j*Nt+i][1]=F2.v[j*Nt+i][1]*HH.m2.v[j][0];
				
				fretrieved.v[j*Nt+i][0]=Fretrieved.v[j*Nt+i][0]*HH.m1.v[j][0];
				fretrieved.v[j*Nt+i][1]=Fretrieved.v[j*Nt+i][1]*HH.m1.v[j][0];
			}
		

		
		//Print the output of the function
		int skiper=1;
		if(gg==0)
		{
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nt;i++)
					fprintf(sout0,"%d %e %e\n",j*Nt+i,F2.v[j*Nt+i][0],F2.v[j*Nt+i][1]);			
			
			for(int j=0;j<Nr;j++)
				for(int i=sampler;i<sampler+1;i++)
				{
					fprintf(out1,"%e %e %e\n",HH.v.v[j][0],f2.v[j*Nt+i][0],f2.v[j*Nt+i][1]);
					fprintf(out2,"%e %e %e\n",HH.v.v[j][0],F2.v[j*Nt+i][0],F2.v[j*Nt+i][1]);
				}
		}
		
		//Movie
		for(int j=0;j<Nr;j++)
			for(int i=sampler;i<sampler+1;i++)
				fprintf(out5,"%e\n",
						(Fretrieved.v[j*skiper*Nt+i][0]*Fretrieved.v[j*skiper*Nt+i][0]+Fretrieved.v[j*skiper*Nt+i][1]*Fretrieved.v[j*skiper*Nt+i][1]) );
		
		
		
	}//End loop
	
	
	//Print the output of the function
	
	for(int j=0;j<Nr;j++)
		for(int i=sampler;i<sampler+1;i++)
		{
			fprintf(out3,"%e %e\n", F2.v[j*Nt+i][0],F2.v[j*Nt+i][1]);
			fprintf(out4,"%e %e\n",	Fretrieved.v[j*Nt+i][0],Fretrieved.v[j*Nt+i][1] );
		}

}