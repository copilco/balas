#include <iostream>
#include <math.h>
#include "grid.h"
#include "wavef3.h"
#include "HankelMatrix.h"
#include "material.h"
#include "propagator.h"
#include <gsl/gsl_sf_bessel.h>
//#include <complex.h>
//#define complex complex<double>

int main()
{
	

	complex I=complex(0.,1.);

	int N=100;//24;//50;
	double R=.05;
	int ord=0;
	double dr=R/(N-1);
	
	vector<double> canonic_r;
	canonic_r.resize(N,0.);
	
	for(int i=0;i<N;i++)
		canonic_r[i] = i*dr; 
	
	HankelMatrix HH(N,R);
	
	int N1=100;//24;//50;
	double dt=.05;

	
	//double r_max = .05;		//%	Maximum radius (5cm)
	//double dr = r_max/(N-1);	//%	Radial spacing
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
	
	grid g1;    //Define a grid object
	g1.set_grid(N1,N,1,dt,1.,1.);  //ORDER IS VERY IMPORTANT DIM1 is for TIME dim2 is for SPACE
	cout << "n1 " << g1.n1;
	cout << "n2 " << g1.n2;
	cout << "n3 " << g1.n3;
	
	field f,F,f2,F2,Fretrieved,fretrieved; 
	f.put_on_grid(g1);
	F.put_on_grid(g1);
	f2.put_on_grid(g1);
	F2.put_on_grid(g1);
	fretrieved.put_on_grid(g1);
	Fretrieved.put_on_grid(g1);
	
	material h;
	h.put_on_grid(g1);
	h.set_v_one_over_rho();
	
	
	
	//Define the function
	double r0=0.;
	double FWHM=5.e-3;
	double Kr=5000.;
	double K=6.4387e+04;
	
	int Nz = 40;			//%	Number of z positions
	double z_max = .25;		//%	Maximum propagation distance
	double dz = z_max/(Nz-1);

	vector<double> phi;
	phi.resize(N,0.);
	
	for(int k=0;k<f.n3;k++)
		for(int j=0;j<f.n2;j++)
			for(int i=0;i<f.n1;i++)
			{
				//double f1,f2;
				//fscanf(in1,"%lf %lf",&f1,&f2);
		
				f.w[f.index(k,j,i)][0]=exp(-2*log(2)*((canonic_r[j]-r0)/FWHM)*((canonic_r[j]-r0)/FWHM)  )*cos(Kr*canonic_r[j])*exp(-2*log(2)*(f.x1[i]*f.x1[i])/.05 );
				f.w[f.index(k,j,i)][1]=exp(-2*log(2)*((canonic_r[j]-r0)/FWHM)*((canonic_r[j]-r0)/FWHM)  )*sin(Kr*canonic_r[j])*exp(-2*log(2)*(f.x1[i]*f.x1[i])/.05 );
		
				F.w[f.index(k,j,i)][0]=f.w[f.index(k,j,i)][0]/HH.m1[j];
				F.w[f.index(k,j,i)][1]=f.w[f.index(k,j,i)][1]/HH.m1[j];
			}
	
	
	/****************************/
	//Loop
	
	
	
	
	
	for (int gg=0;gg<Nz;gg++)
	{
		printf("gg =%d",gg);
		
		/**********************//**********************//**********************//**********************//**********************/
		//Dispersion
		/**********************//**********************//**********************//**********************//**********************/
		//prop_dispersion_1(f, h, dz);
		
		double factor=1./f.n1;
		/**********************/
		f.fft_over1_F();
		/**********************/
		
		for(int k=0;k<f.n3;k++)
			for(int j=0;j<f.n2;j++)
				for(int i=0;i<f.n1;i++)
				{
					
					double dispersion=f.q1[i]*f.q1[i]*h.a1;
					
					//The factor sould be: //(dx*dy*dz*dx*dy*dz/dospi/dospi/dospi);
					
					complex dispersion_operator=exp(-I*dz*dispersion );
					
					double aux0r=f.w[f.index(k,j,i)][0];
					double aux0i=f.w[f.index(k,j,i)][1];
					
					f.w[f.index(k,j,i)][0]=( aux0r*real(dispersion_operator)+aux0i*imag(dispersion_operator) )*factor;
					f.w[f.index(k,j,i)][1]=( aux0i*real(dispersion_operator)-aux0r*imag(dispersion_operator) )*factor;
					
					//w.w[w.index(k,j,i)][0]*=factor;
					//w.w[w.index(k,j,i)][1]*=factor;
				}
		
		/**********************/
		f.fft_over1_B();
		/**********************/
		
		for(int k=0;k<f.n3;k++)
			for(int j=0;j<f.n2;j++)
				for(int i=0;i<f.n1;i++)
				{
					
					F.w[f.index(k,j,i)][0]=f.w[f.index(k,j,i)][0]/HH.m1[j];
					F.w[f.index(k,j,i)][1]=f.w[f.index(k,j,i)][1]/HH.m1[j];
				}
		
		/**********************//**********************//**********************//**********************//**********************/
		//Dispersion
		/**********************//**********************//**********************//**********************//**********************/

		//Make F2 and Fretrieved variables to zero
		for(int i=0;i<f.n1*f.n2*f.n3;i++)
		{
			F2.w[i][0]=0.;
			F2.w[i][1]=0.;
			
			Fretrieved.w[i][0]=0.;
			Fretrieved.w[i][1]=0.;
		}
		
		//Apply the matrix  Forward Hankel
		
		for(int k=0;k<f.n3;k++)
			for(int jj=0;jj<HH.N;jj++)
				for(int j=0;j<f.n2;j++)
					for(int i=0;i<f.n1;i++)
					{
						F2.w[f.index(k,jj,i)][0]+=HH.C[jj *HH.N+j]*F.w[f.index(k,j,i)][0];
						F2.w[f.index(k,jj,i)][1]+=HH.C[jj *HH.N+j]*F.w[f.index(k,j,i)][1];
					}
		
		//Prepare the real transformed function. In this case just to display
		//This gives f2 previous to the transformation
		for(int k=0;k<f.n3;k++)
			for(int j=0;j<f.n2;j++)
				for(int i=0;i<f.n1;i++)
				{
					f2.w[f.index(k,j,i)][0]=F2.w[f.index(k,j,i)][0]*HH.m2[j];
					f2.w[f.index(k,j,i)][1]=F2.w[f.index(k,j,i)][1]*HH.m2[j];
				}
		
		
		//Multiply by the phase to make the propagation.
		for(int k=0;k<f.n3;k++)
			for(int j=0;j<f.n2;j++)
				for(int i=0;i<f.n1;i++)
				{
					double aux0r=F2.w[f.index(k,j,i)][0];
					double aux0i=F2.w[f.index(k,j,i)][1];
					
					phi[j]=(sqrt(K*K - 2*pi*HH.v[j]*2*pi*HH.v[j]) - K)*gg*dz; //The phase
					
					F2.w[f.index(k,j,i)][0]=( aux0r*cos(phi[j])-aux0i*sin(phi[j]) );
					F2.w[f.index(k,j,i)][1]=( aux0i*cos(phi[j])+aux0r*sin(phi[j]) );
					
				}
		
		
		//Apply the matrix Backward Hankel Transform over F2
		//for(int j=0;j<N;j++)
		//	for(int i=0;i<N;i++)
		//	{
		for(int k=0;k<f.n3;k++)
			for(int jj=0;jj<HH.N;jj++)
				for(int j=0;j<f.n2;j++)
					for(int i=0;i<f.n1;i++)
					{
								
						Fretrieved.w[f.index(k,jj,i)][0]+=HH.C[jj*HH.N+j]*F2.w[f.index(k,j,i)][0];
						Fretrieved.w[f.index(k,jj,i)][1]+=HH.C[jj*HH.N+j]*F2.w[f.index(k,j,i)][1];
					}
		
		//Prepare the scaled functions f2 and fretrieved after the propagator and the back Hankel transform

		for(int k=0;k<f.n3;k++)
			for(int j=0;j<f.n2;j++)
				for(int i=0;i<f.n1;i++)
				{
					f2.w[f2.index(k,j,i)][0]=F2.w[f2.index(k,j,i)][0]*HH.m2[j];
					f2.w[f2.index(k,j,i)][1]=F2.w[f2.index(k,j,i)][1]*HH.m2[j];
					
					fretrieved.w[f2.index(k,j,i)][0]=Fretrieved.w[f2.index(k,j,i)][0]*HH.m1[j];
					fretrieved.w[f2.index(k,j,i)][1]=Fretrieved.w[f2.index(k,j,i)][1]*HH.m1[j];
				}
		
		/*
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
		 */
		
		snapshot(out3,f2,1,1,1);
		snapshot(out4,fretrieved,1,1,1);
		snapshot(out5,Fretrieved,1,1,1);
		
		
	}
	
	
	
}
