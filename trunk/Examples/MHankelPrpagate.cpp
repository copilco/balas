#include <iostream>
#include <math.h>
#include "grid.h"
#include "wavef3.h"
#include <gsl/gsl_sf_bessel.h>
//#include <complex.h>
//#define complex complex<double>

int main()
{
	complex I=complex(0.,1.);
	int N=1024;//50;
	double R=.05;
	int ord=0;
	
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
	
	int zer=5*3001;
	int nn=5;//3001;
	int mm=3001;
	
	vector<double> zeros;
	zeros.resize(zer,0.);
	
	
	for(int j=0;j<mm;j++)
	{
		for(int i=0;i<nn;i++)
		{
			double f;
			fscanf(in0,"%lf",&f);
			zeros[j*nn+i]=f;
		}
		fscanf(in0,"\n");
	}
	
	
	
	vector<double> c;
	c.resize(N+1,0.);
	
	
	for(int j=0;j<N+1;j++)
	{
		c[j]=zeros[j*nn+(ord)];
		//fprintf(out0,"%e ",c[j]);
	}
	
	double	V = c[N]/(2*pi*R);    //% Maximum frequency
	printf("V %e",V);
	
	vector<double> r;
	r.resize(N,0.);	
	
	vector<double> canonic_r;
	canonic_r.resize(N,0.);	
	
	vector<double> v;
	v.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		r[i] = c[i]*R/c[N];   //% Radius vector
		canonic_r[i] = i*dr;   //% Radius vector
		v[i] = c[i]/(2*pi*R);
	}
	
	vector<double> jnm;
	jnm.resize(N*N,0.);
	vector<double> C;
	C.resize(N*N,0.);
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
			jnm[j*N+i]=c[i];
	
	
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N]);
			C[j*N+i]=(2./c[N])*gsl_sf_bessel_J0 (jnm[j*N+i]*jnm[i*N+j]/c[N])/abs(gsl_sf_bessel_J1(jnm[i*N+j]))/abs(gsl_sf_bessel_J1(jnm[j*N+i]));
		}
	
	vector<double> m1;
	m1.resize(N,0.);
	vector<double> m2;
	m2.resize(N,0.);
	
	for(int i=0;i<N;i++)
	{
		m1[i]=abs(gsl_sf_bessel_J1(c[i]))/R;
		m2[i]=m1[i]*R/V;
	}
	
	/********************/
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
	
	/*
	 vector<double> f;
	 f.resize(N,0.);
	 vector<double> F;
	 F.resize(N,0.);
	 vector<double> f2;
	 f2.resize(N,0.);
	 vector<double> F2;
	 F2.resize(N,0.);
	 vector<double> Fretrieved;
	 Fretrieved.resize(N,0.);
	 vector<double> fretrieved;
	 fretrieved.resize(N,0.);
	 */
	
	
	//Define the function
	double r0=0.;
	double FWHM=5.e-3;
	double Kr=5000.;
	double K=6.4387e+04;
	
	int Nz = 200;			//%	Number of z positions
	double z_max = .25;		//%	Maximum propagation distance
	double dz = z_max/(Nz-1);
	//z = (0:Nz-1)'*dz;
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
		
		F.w[i][0]=f.w[i][0]/m1[i];
		F.w[i][1]=f.w[i][1]/m1[i];
		
		F2.w[i][0]=0.;
		F2.w[i][1]=0.;
		
		Fretrieved.w[i][0]=0.;
		Fretrieved.w[i][0]=0.;
		
	}
	
	//Apply the matrix  Forward Hankel
	for(int j=0;j<N;j++)
		for(int i=0;i<N;i++)
		{
			F2.w[j][0]+=C[j*N+i]*F.w[i][0];
			F2.w[j][1]+=C[j*N+i]*F.w[i][1];
		}
	
	for(int i=0;i<N;i++)
	{
		f2.w[i][0]=F2.w[i][0]*m2[i];
		f2.w[i][1]=F2.w[i][1]*m2[i];
	}
	
	for(int i=0;i<N;i++)
	{
		 //f2.w[i][]=F2[i]*m2[i];
		 //fprintf(out3,"%e\n",phi[i]);
		 fprintf(out0,"%e %e %e %e\n",r[i],canonic_r[i],f.w[i][0],f.w[i][1]);
		 fprintf(out1,"%e %e %e\n",v[i],f2.w[i][0],f2.w[i][1]);
		 fprintf(out2,"%e %e %e\n",v[i],F2.w[i][0],F2.w[i][1]);
	}
	
	/****************************/
	//Loop
	
	for (int gg=0;gg<Nz;gg++)
	{
		
		printf("gg =%d",gg);
		
		
		for(int i=0;i<N;i++)
		{
			//f.w[i][0]=exp(-2*log(2)*((r[i]-r0)/FWHM)*((r[i]-r0)/FWHM)  )*cos(Kr*r[i]);
			//f.w[i][1]=exp(-2*log(2)*((r[i]-r0)/FWHM)*((r[i]-r0)/FWHM)  )*sin(Kr*r[i]);
			
			// canonic_r[i]
			/*
			double f1,f2;
			fscanf(in1,"%lf %lf",&f1,&f2);
			
			
			f.w[i][0]=f1;//exp(-2*log(2)*((canonic_r[i]-r0)/FWHM)*((canonic_r[i]-r0)/FWHM)  )*cos(Kr*canonic_r[i]);
			f.w[i][1]=f2;//exp(-2*log(2)*((canonic_r[i]-r0)/FWHM)*((canonic_r[i]-r0)/FWHM)  )*sin(Kr*canonic_r[i]);
			
			F.w[i][0]=f.w[i][0]/m1[i];
			F.w[i][1]=f.w[i][1]/m1[i];
			*/
			
			F2.w[i][0]=0.;
			F2.w[i][1]=0.;
			
			Fretrieved.w[i][0]=0.;
			Fretrieved.w[i][1]=0.;
		}
		
		//Apply the matrix  Forward Hankel
		for(int j=0;j<N;j++)
			for(int i=0;i<N;i++)
			{
				F2.w[j][0]+=C[j*N+i]*F.w[i][0];
				F2.w[j][1]+=C[j*N+i]*F.w[i][1];
			}
		
		for(int i=0;i<N;i++)
		{
			f2.w[i][0]=F2.w[i][0]*m2[i];
			f2.w[i][1]=F2.w[i][1]*m2[i];
		}
		
		//complex a;
		//complex b;
		for(int i=0;i<N;i++)
		{
			//complex kinoperator=exp(-I*dt*kinetic );
			
			double aux0r=F2.w[i][0];
			double aux0i=F2.w[i][1];
			
			phi[i]=(sqrt(K*K - 2*pi*v[i]*2*pi*v[i]) - K)*gg*dz;//  (sqrt(K*K - 2*pi*v[i]) - K)*gg*dz;
			
			F2.w[i][0]=( aux0r*cos(phi[i])-aux0i*sin(phi[i]) );
			F2.w[i][1]=( aux0i*cos(phi[i])+aux0r*sin(phi[i]) );
			
		}
		
		//Apply the matrix  Forward Hankel
		for(int j=0;j<N;j++)
			for(int i=0;i<N;i++)
			{
				Fretrieved.w[j][0]+=C[j*N+i]*F2.w[i][0];
				Fretrieved.w[j][1]+=C[j*N+i]*F2.w[i][1];
			}
		
		for(int i=0;i<N;i++)
		{
			f2.w[i][0]=F2.w[i][0]*m2[i];
			f2.w[i][1]=F2.w[i][1]*m2[i];
			
			fretrieved.w[i][0]=Fretrieved.w[i][0]*m1[i];
			fretrieved.w[i][1]=Fretrieved.w[i][1]*m1[i];
		}
		
		
		//Fretrieved = C*F2;           %% Inverse hankel transform
		
		//fretrieved = Fretrieved.*m1
		
		
		for(int i=0;i<N;i++)
		{
			//fprintf(out3,"%e %e %e\n",v[i],F2.w[i][0],F2.w[i][1]);
			//fprintf(out4,"%e %e %e\n",r[i],Fretrieved.w[i][0],Fretrieved.w[i][1]);
			
			/*
			
			//f2.w[i][]=F2[i]*m2[i];
			//fprintf(out3,"%e\n",phi[i]);
			fprintf(out0,"%e %e %e %e\n",r[i],canonic_r[i],f.w[i][0],f.w[i][1]);
			fprintf(out1,"%e %e %e\n",v[i],f2.w[i][0],f2.w[i][1]);
			fprintf(out2,"%e %e %e\n",v[i],F2.w[i][0],F2.w[i][1]);
			
			fprintf(out3,"%e %e %e\n",v[i],F2.w[i][0],F2.w[i][1]);
					//sqrt( F2.w[i][0]*F2.w[i][0]+F2.w[i][1]*F2.w[i][1]) , sqrt(f2.w[i][0]*f2.w[i][0]+f2.w[i][1]*f2.w[i][1] ));
			
			*/
			
			
			
//					sqrt( f.w[i][0]*f.w[i][0]+f.w[i][1]*f.w[i][1]) , 
//					sqrt(fretrieved.w[i][0]*fretrieved.w[i][0]+fretrieved.w[i][1]*fretrieved.w[i][1]) );
			fprintf(out4,"%e\n",
					(fretrieved.w[i][0]*fretrieved.w[i][0]+fretrieved.w[i][1]*fretrieved.w[i][1]) );
			fprintf(out5,"%e\n", // Fretrieved.w[i][1]);//
					(Fretrieved.w[i][0]*Fretrieved.w[i][0]+Fretrieved.w[i][1]*Fretrieved.w[i][1]) );
			
		}
		
	}
	
	
	
}
