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
	arrai f0,f1;
	arrai F0,F1;
	arrai f20,f21;
	arrai F20,F21;
	arrai Fretrieved0,Fretrieved1;
	arrai fretrieved0,fretrieved1; 
	arrai F0mat,F1mat,F20mat,F21mat;
	
	arrai FCmat,F2Cmat,Fretrieved;
	arrai f2C;
	
	//phi.init(N,0.);
	f0.init(Nr,1,0.);
	f1.init(Nr,1,0.);
	F0.init(Nr,1,0.);
	F1.init(Nr,1,0.);
	f20.init(Nr,1,0.);
	f21.init(Nr,1,0.);
	
	f2Cmat.init(Nr,Nt,0.);
	
	F20.init(Nr,1,0.);
	F21.init(Nr,1,0.);
	
	F0mat.init(Nr,Nt,0.);
	F1mat.init(Nr,Nt,0.);
	
	FCmat.init(Nr,Nt,0.); //This is the complex version
	
	F20mat.init(Nr,Nt,0.);
	F21mat.init(Nr,Nt,0.);
	
	F2Cmat.init(Nr,Nt,0.);//This is the complez version
	
	fretrieved0.init(Nr,1,0.);
	fretrieved1.init(Nr,1,0.);
	Fretrieved0.init(Nr,1,0.);
	Fretrieved1.init(Nr,1,0.);
	
	Fretrieved.init(Nr,Nt,0.);//This is the complez version
	
	double K=6.4387e+04;
	
	double z_max = .25;		//%	Maximum propagation distance
	int Nz = 200;//int(z_max/dz); 
	double dz = z_max/(Nz-1);
	
	arrai phi;
	phi.init(Nr,1,0.);
	
	//vector<complex> pp;
	//pp.init(N,0.);
	
	for(int i=0;i<Nr;i++)
	{
		//double r=i*dr;
		double read1,read2;
		fscanf(in1,"%lf %lf",&read1,&read2);
		
		f0.v[i][0]= read1;  //exp(-(HH.r.v[i]*HH.r.v[i]/w0/w0))*cos(HH.v.v[*HH.v.v[i]*HH.r.v[i]*HH.r.v[i]/2./focus);
		f1.v[i][0]= read2; // exp(-(HH.r.v[i]*HH.r.v[i]/w0/w0))*sin(HH.v.v[*HH.v.v[i]*HH.r.v[i]*HH.r.v[i]/2./focus);
	}
	
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nt;i++)
		{			
			FCmat.v[j*Nt+i][0]=f0.v[j][0]/HH.m1.v[j][0];
			FCmat.v[j*Nt+i][1]=f1.v[j][0]/HH.m1.v[j][0];
		}
	
	
	for(int i=0;i<Nr;i++)
	{
		fprintf(out0,"%e %e %e %e\n",HH.r.v[i][0],HH.r.v[i][0],f0.v[i][0],f1.v[i][0]);
	}
	
	
	//Loop
	
	
	//Loop

	for (int gg=0;gg<Nz;gg++)
	{
		
		printf("gg =%d",gg);
		//Make F2 and Fretrieved variables to zero		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nt;i++)
			{
				F2Cmat.v[j*Nt+i][0]=0.;
				F2Cmat.v[j*Nt+i][1]=0.;
				
				FCretrieved.v[i][0]=0.;
				FCretrieved.v[i][1]=0.;
			}
		
		
		/*
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, 1, Nr,
					 1.0, HH.C.v, lda, F0.v, ldb, 0.0, F20.v, ldc);
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, 1, Nr,
					 1.0, HH.C.v, lda, F1.v, ldb, 0.0, F21.v, ldc);
		
		
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,
					 1.0, HH.C.v, lda, F0mat.v, Nt, 0.0, F20mat.v, Nt);
		
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,
					 1.0, HH.C.v, lda, F1mat.v, Nt, 0.0, F21mat.v, Nt);
		*/
		__CLPK_doublecomplex alpha;
		__CLPK_doublecomplex beta;
		__CLPK_doublecomplex gamma;
		
		alpha.r=1;
		alpha.i=0;
		
		gamma.r=1;
		gamma.i=0;
		
		cblas_zgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,
					 &alpha, HH.C.v, Nr,FCmat.v,Nt,&gamma, F2Cmat.v, Nt);
		
		
		
		 //Apply the matrix  Forward Hankel
		// for(int j=0;j<Nr;j++)
		// for(int i=0;i<Nr;i++)
		// {
		// F20.v[j]+=HH.C.v[j*Nt+i]*F0.v[i];
		// F21.v[j]+=HH.C.v[j*Nt+i]*F1.v[i];
		// }
		 
		
		
		
		
		
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int i=0;i<Nr;i++)
		{
			double temp0=HH.m2.v[i][0];
			f2Cmat.v[i]=F2Cmat.v[i]*temp0;
			
//			f20.v[i][0]=F20.v[i][0]*HH.m2.v[i][0];
//			f21.v[i][0]=F21.v[i][0]*HH.m2.v[i][0];
		}
		
		
		//Multiply by the phase to make the propagation.
		for(int i=0;i<Nr;i++)
		{
			double aux0r=F20.v[i][0];
			double aux0i=F21.v[i][0];
			
			phi.v[i][0]=(sqrt(K*K - 2*pi*HH.v.v[i][0]*2*pi*HH.v.v[i][0]) - K)*gg*dz; //The phase
			//	phi.v[i]=0.;//(sqrt(K*K - 2*pi*HH.v.v[i]*2*pi*HH.v.v[i]) - K)*gg*dz; //The phase
			
			F20.v[i][0]=( aux0r*cos(phi.v[i][0])-aux0i*sin(phi.v[i][0]) );
			F21.v[i][0]=( aux0i*cos(phi.v[i][0])+aux0r*sin(phi.v[i][0]) );
			
		}
		
		
		
		
		//Multiply by the phase to make the propagation.
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nt;i++)
			{
//				double aux0r=F20mat.v[j*Nt+i][0];
//				double aux0i=F21mat.v[j*Nt+i][0];
				
				double aux0r=F2Cmat.v[j*Nt+i][0];
				double aux0i=F2Cmat.v[j*Nt+i][1];
				
				phi.v[j][0]=(sqrt(K*K - 2*pi*HH.v.v[j][0]*2*pi*HH.v.v[j][0]) - K)*gg*dz; //The phase
				//	phi.v[i]=0.;//(sqrt(K*K - 2*pi*HH.v.v[i]*2*pi*HH.v.v[i]) - K)*gg*dz; //The phase
				
//				F20mat.v[j*Nt+i][0]=( aux0r*cos(phi.v[j][0])-aux0i*sin(phi.v[j][0]) );
//				F21mat.v[j*Nt+i][0]=( aux0i*cos(phi.v[j][0])+aux0r*sin(phi.v[j][0]) );
				
				F2Cmat.v[j*Nt+i][0]=( aux0r*cos(phi.v[j][0])-aux0i*sin(phi.v[j][0]) );
				F2Cmat.v[j*Nt+i][1]=( aux0i*cos(phi.v[j][0])+aux0r*sin(phi.v[j][0]) );
				
			}
		
		
		
		/*
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, N, 1, Nr,
					 1.0, HH.C.v, lda, F20.v, ldb, 0.0,Fretrieved0.v, ldc);
		cblas_dgemm (CblasRowMajor, 
					 CblasNroTrans, CblasNoTrans, Nr, 1, Nr,
					 1.0, HH.C.v, lda, F21.v, ldb, 0.0, Fretrieved1.v, ldc);
		*/
		
		/*
		 //Apply the matrix Backward Hankel Transform over F2
		 for(int j=0;j<Nr;j++)
		 for(int i=0;i<Nr;i++)
		 {
		 Fretrieved0.v[j]+=HH.C.v[j*Nr+i]*F20.v[i];
		 Fretrieved1.v[j]+=HH.C.v[j*Nr+i]*F21.v[i];
		 }
		 */
		
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform
		for(int i=0;i<Nr;i++)
		{
			f20.v[i][0]=F20.v[i][0]*HH.m2.v[i][0];
			f21.v[i][0]=F21.v[i][0]*HH.m2.v[i][0];
			
			fretrieved0.v[i][0]=Fretrieved0.v[i][0]*HH.m1.v[i][0];
			fretrieved1.v[i][0]=Fretrieved1.v[i][0]*HH.m1.v[i][0];
		}
		
		
		
		//Print the output of the function
		int skiper=1;
		
		if(gg==20)
		{
			
			/*
			 for(int j=0;j<Nt;j++)
			 for(int i=0;i<Nr;i++)
			 fprintf(sout0,"%d %e %e\n",j*Nr+i,F0mat.v[j*Nr+i],F1mat.v[j*Nr+i]);
			 */
			
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nt;i++)
					fprintf(sout0,"%d %e %e\n",j*Nt+i,F2Cmat.v[j*Nt+i][0],F2Cmat.v[j*Nt+i][1]);
			
			//	for(int i=0;i<Nt*Nr;i++)
			//		fprintf(sout0,"%d %e %e\n",i,F20mat.v[i],F21mat.v[i]);
		}
		
		
		for(int i=0;i<Nr;i++)
		{
			if(gg==20)
			{
				fprintf(out1,"%e %e %e\n",HH.v.v[i][0],f20.v[i][0],f21.v[i][0]);
			//	fprintf(out2,"%e %e %e\n",HH.v.v[i][0],F20.v[i][0],F21.v[i][0]);
			}		
			
			//Movie
			fprintf(out5,"%e\n", // Fretrieved1.v[i*skiper]);//
					(Fretrieved0.v[i*skiper][0]*Fretrieved0.v[i*skiper][0]+Fretrieved1.v[i*skiper][0]*Fretrieved1.v[i*skiper][0]) );
			
		}
		
	}//End loop
	
	/*
	//Print the output of the function
	int skiper=1;
	for(int i=0;i<Nr;i++)
	{
		fprintf(out3,"%e %e\n",
				F20.v[i],F21.v[i]);
		fprintf(out4,"%e %e\n",
				Fretrieved0.v[i*skiper],Fretrieved1.v[i*skiper] );
	}
	*/
	
	
	
	
}