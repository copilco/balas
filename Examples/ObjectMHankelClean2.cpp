#include <iostream>
#include <math.h>
#include "grid.h"
#include "wavef3.h"
#include "HankelMatrix.h"
#include <gsl/gsl_sf_bessel.h>
//#include <complex.h>
//#define complex complex<double>

int main()
{
	

	complex I=complex(0.,1.);
	int N=1024;//50;
	double R=.05;
	int ord=0;

	
	HankelMatrix HH(N,R);
	
	double r_max = .05;		//%	Maximum radius (5cm)
	double dr = r_max/(N-1);	//%	Radial spacing
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
	
	grid g1;                                    //Define a grid object
	g1.set_grid(N,1,1,1.,1.,1.);
	
	field f,F,f2,F2,Fretrieved,fretrieved; 
	f.put_on_grid(g1);
	F.put_on_grid(g1);
	f2.put_on_grid(g1);
	F2.put_on_grid(g1);
	fretrieved.put_on_grid(g1);
	Fretrieved.put_on_grid(g1);
	
	
	
	//Define the function
	double r0=0.;
	double FWHM=5.e-3;
	double Kr=5000.;
	double K=6.4387e+04;
	
	int Nz = 200;			//%	Number of z positions
	double z_max = .25;		//%	Maximum propagation distance
	double dz = z_max/(Nz-1);
	
	vector<double> phi;
	phi.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		//f.w[i][0]=exp(-2*log(2)*((r[i]-r0)/FWHM)*((r[i]-r0)/FWHM)  )*cos(Kr*r[i]);
		//f.w[i][1]=exp(-2*log(2)*((r[i]-r0)/FWHM)*((r[i]-r0)/FWHM)  )*sin(Kr*r[i]);
		
		double f1,f2;
		fscanf(in1,"%lf %lf",&f1,&f2);
		
		f.w[i][0]=f1;//exp(-2*log(2)*((canonic_r[i]-r0)/FWHM)*((canonic_r[i]-r0)/FWHM)  )*cos(Kr*canonic_r[i]);
		f.w[i][1]=f2;//exp(-2*log(2)*((canonic_r[i]-r0)/FWHM)*((canonic_r[i]-r0)/FWHM)  )*sin(Kr*canonic_r[i]);
		
		F.w[i][0]=f.w[i][0]/HH.m1[i];
		F.w[i][1]=f.w[i][1]/HH.m1[i];

	}
	
	
	/****************************/
	//Loop
	
	for (int gg=0;gg<Nz;gg++)
	{
		
		printf("gg =%d",gg);
		//Make F2 and Fretrieved variables to zero
		for(int i=0;i<N;i++)
		{
			F2.w[i][0]=0.;
			F2.w[i][1]=0.;
			
			Fretrieved.w[i][0]=0.;
			Fretrieved.w[i][1]=0.;
		}
		
		//Apply the matrix  Forward Hankel
		for(int j=0;j<N;j++)
			for(int i=0;i<N;i++)
			{
				//F2.w[j][0]+=C[j*N+i]*F.w[i][0];
				//F2.w[j][1]+=C[j*N+i]*F.w[i][1];

				F2.w[j][0]+=HH.C[j*N+i]*F.w[i][0];
				F2.w[j][1]+=HH.C[j*N+i]*F.w[i][1];
			}
		
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int i=0;i<N;i++)
		{
			//f2.w[i][0]=F2.w[i][0]*m2[i];
			//f2.w[i][1]=F2.w[i][1]*m2[i];
			
			f2.w[i][0]=F2.w[i][0]*HH.m2[i];
			f2.w[i][1]=F2.w[i][1]*HH.m2[i];
		}
		
		
		//Multiply by the phase to make the propagation.
		for(int i=0;i<N;i++)
		{
			double aux0r=F2.w[i][0];
			double aux0i=F2.w[i][1];
			
			phi[i]=(sqrt(K*K - 2*pi*HH.v[i]*2*pi*HH.v[i]) - K)*gg*dz; //The phase
			
			F2.w[i][0]=( aux0r*cos(phi[i])-aux0i*sin(phi[i]) );
			F2.w[i][1]=( aux0i*cos(phi[i])+aux0r*sin(phi[i]) );
			
		}
		
		//Apply the matrix Backward Hankel Transform over F2
		for(int j=0;j<N;j++)
			for(int i=0;i<N;i++)
			{
				Fretrieved.w[j][0]+=HH.C[j*N+i]*F2.w[i][0];
				Fretrieved.w[j][1]+=HH.C[j*N+i]*F2.w[i][1];
			}
		
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform
		for(int i=0;i<N;i++)
		{
			f2.w[i][0]=F2.w[i][0]*HH.m2[i];
			f2.w[i][1]=F2.w[i][1]*HH.m2[i];
			
			fretrieved.w[i][0]=Fretrieved.w[i][0]*HH.m1[i];
			fretrieved.w[i][1]=Fretrieved.w[i][1]*HH.m1[i];
		}
		
		//Print the output of the function
		for(int i=0;i<N;i++)
		{
			fprintf(out3,"%e\n",
					(f2.w[i][0]*f2.w[i][0]+f2.w[i][1]*f2.w[i][1]) );
			fprintf(out4,"%e\n",
					(fretrieved.w[i][0]*fretrieved.w[i][0]+fretrieved.w[i][1]*fretrieved.w[i][1]) );
			fprintf(out5,"%e\n", // Fretrieved.w[i][1]);//
					(Fretrieved.w[i][0]*Fretrieved.w[i][0]+Fretrieved.w[i][1]*Fretrieved.w[i][1]) );
			
		}
		
	}
	
	
	
}
