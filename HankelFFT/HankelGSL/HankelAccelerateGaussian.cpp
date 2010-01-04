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
	int N=2000;//50;
	double R=40.;//.05;
	
	//Parametros para blas
	int lda=N;
	int ldb=1;
	int ldc=1;

	HankelMatrix HH(N,R);
	
	//double r_max = .05;		//%	Maximum radius (5cm)
	double dr = R/(N-1);	//%	Radial spacing
	//r = nr*dr;
	
	printf("N=%d\n",N);
	
	
	FILE *in0,*in1;
	in0=fopen("besselCeros.txt","r");
	in1=fopen("InputFunction.txt","r");
	
	FILE *out0,*out1,*out2,*out3,*out4,*out5;
	out0=fopen("pout0.txt","w");
	out1=fopen("pout1.txt","w");
	out2=fopen("pout2.txt","w");
	out3=fopen("pout3.txt","w");
	out4=fopen("pout4.txt","w");
	out5=fopen("pout5.txt","w");
	
	
	/***********************************/
	//Transformation
	
	//arrai phi;
	arrai f0,f1;
	arrai F0,F1;
	arrai f20,f21;
	arrai F20,F21;
	arrai Fretrieved0,Fretrieved1;
	arrai fretrieved0,fretrieved1; 
	//phi.init(N,0.);
	f0.init(N,0.);
	f1.init(N,0.);
	F0.init(N,0.);
	F1.init(N,0.);
	f20.init(N,0.);
	f21.init(N,0.);
	F20.init(N,0.);
	F21.init(N,0.);
	fretrieved0.init(N,0.);
	fretrieved1.init(N,0.);
	Fretrieved0.init(N,0.);
	Fretrieved1.init(N,0.);
	
	
	double K=6.4387e+04;
	
	
	double z_max =2000.;		//%	Maximum propagation distance
	double dz = 0.005;//z_max/(Nz-1);
	int Nz = int(z_max/dz); 
	
	
	arrai phi;
	phi.init(N,0.);
	
	//vector<complex> pp;
	//pp.init(N,0.);
	
	double r0=0.;
	double FWHM=5.e-3;
	double Kr=5000.;
	double K=6.4387e+04;
	
	for(int i=0;i<N;i++)
	{
		
		//double r=i*dr;
		double read1,read2;
		fscanf(in1,"%lf %lf",&read1,&read2);
		
		f0.v[i]= exp(-(HH.r.v[i]*HH.r.v[i]/w0/w0))*cos(HH.v.v[*HH.v.v[i]*HH.r.v[i]*HH.r.v[i]/2./focus);
		f1.v[i]= exp(-(HH.r.v[i]*HH.r.v[i]/w0/w0))*sin(HH.v.v[*HH.v.v[i]*HH.r.v[i]*HH.r.v[i]/2./focus);
		
		F0.v[i]=f0.v[i]/HH.m1.v[i];
		F1.v[i]=f1.v[i]/HH.m1.v[i];
		
	}
	
	for(int i=0;i<N;i++)
	{
		fprintf(out0,"%e %e %e %e\n",HH.r.v[i],HH.r.v[i],f0.v[i],f1.v[i]);
		
		
	}
	
	/****************************/
	//Loop
	
	for (int gg=0;gg<Nz;gg++)
	{
		
		printf("gg =%d",gg);
		//Make F2 and Fretrieved variables to zero
		for(int i=0;i<N;i++)
		{
			F20.v[i]=0.;
			F21.v[i]=0.;
			
			Fretrieved0.v[i]=0.;
			Fretrieved1.v[i]=0.;
		}
		
		
		
		
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, N, 1, N,
					 1.0, HH.C.v, lda, F0.v, ldb, 0.0, F20.v, ldc);
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, N, 1, N,
					 1.0, HH.C.v, lda, F1.v, ldb, 0.0, F21.v, ldc);
		
		
		/*
		//Apply the matrix  Forward Hankel
		for(int j=0;j<N;j++)
			for(int i=0;i<N;i++)
			{
				F20.v[j]+=HH.C.v[j*N+i]*F0.v[i];
				F21.v[j]+=HH.C.v[j*N+i]*F1.v[i];
			}
		
		 */
		
		
		
		
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int i=0;i<N;i++)
		{
			
			f20.v[i]=F20.v[i]*HH.m2.v[i];
			f21.v[i]=F21.v[i]*HH.m2.v[i];
		}
		
		
		//Multiply by the phase to make the propagation.
		for(int i=0;i<N;i++)
		{
			double aux0r=F20.v[i];
			double aux0i=F21.v[i];
			
			phi.v[i]=(sqrt(K*K - 2*pi*HH.v.v[i]*2*pi*HH.v.v[i]) - K)*gg*dz; //The phase
			//	phi.v[i]=0.;//(sqrt(K*K - 2*pi*HH.v.v[i]*2*pi*HH.v.v[i]) - K)*gg*dz; //The phase
			
			F20.v[i]=( aux0r*cos(phi.v[i])-aux0i*sin(phi.v[i]) );
			F21.v[i]=( aux0i*cos(phi.v[i])+aux0r*sin(phi.v[i]) );
			
		}
		
		
		
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, N, 1, N,
					 1.0, HH.C.v, lda, F20.v, ldb, 0.0,Fretrieved0.v, ldc);
		cblas_dgemm (CblasRowMajor, 
					 CblasNoTrans, CblasNoTrans, N, 1, N,
					 1.0, HH.C.v, lda, F21.v, ldb, 0.0, Fretrieved1.v, ldc);

		
		/*
		//Apply the matrix Backward Hankel Transform over F2
		for(int j=0;j<N;j++)
			for(int i=0;i<N;i++)
			{
				Fretrieved0.v[j]+=HH.C.v[j*N+i]*F20.v[i];
				Fretrieved1.v[j]+=HH.C.v[j*N+i]*F21.v[i];
			}
		*/
		 
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform
		for(int i=0;i<N;i++)
		{
			f20.v[i]=F20.v[i]*HH.m2.v[i];
			f21.v[i]=F21.v[i]*HH.m2.v[i];
			
			fretrieved0.v[i]=Fretrieved0.v[i]*HH.m1.v[i];
			fretrieved1.v[i]=Fretrieved1.v[i]*HH.m1.v[i];
		}
		
		
		
		//Print the output of the function
		int skiper=1;
		for(int i=0;i<N;i++)
		{
			if(gg==0)
			{
				fprintf(out1,"%e %e %e\n",HH.v.v[i],f20.v[i],f21.v[i]);
				fprintf(out2,"%e %e %e\n",HH.v.v[i],F20.v[i],F21.v[i]);
			}		
	
			//Movie
			if((gg%100)==0)
				fprintf(out5,"%e\n", // Fretrieved1.v[i*skiper]);//
					(Fretrieved0.v[i*skiper]*Fretrieved0.v[i*skiper]+Fretrieved1.v[i*skiper]*Fretrieved1.v[i*skiper]) );
			
		}
		
	}//End loop
		
	
	//Print the output of the function
	int skiper=1;
	for(int i=0;i<N;i++)
	{
		fprintf(out3,"%e %e\n",
				F20.v[i],F21.v[i]);
		fprintf(out4,"%e %e\n",
				Fretrieved0.v[i*skiper],Fretrieved1.v[i*skiper] );
	}
		
		
		
	
	
	
	
}
