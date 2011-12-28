/*
 *  testingCranck1D.cpp
 *  
 *
 *  Created by alexis on 02/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream.h>
#include <complex.h>
#include <alloc.h>
#include <stdio.h>
#include "HankelMatrix.h"

#define complex complex<double>

using namespace std;

#define ndimz 3000
#define dr 0.1
#define dt 0.01
#define nex 0
#define ssw 1	//Swicht to do uniform ssw=0, and non uniform ssw=1 grid


double rho0  = 100.;
double R     = 200.;

int Ntime    = 1000;
int snap     = 50;

complex *gam1=(complex *)malloc((ndimz+nex)*sizeof(complex));
void nrerror (char error_text[]);
void tridag(complex a[], complex b[], complex c[], complex r[],complex u[], int n);


int main()
{
	
	FILE *out0;
	FILE *out1;	
	FILE *out2;		
	FILE *out3;		
	
	out0 = fopen("out0.txt","w");
	out1 = fopen("out1.txt","w");	
	out2 = fopen("out2.txt","w");	
	out3 = fopen("out3.txt","w");		
	
	complex *phi,*rv,*phil;
	complex *a, *b, *c;
	double *v, *rho, *drr;

	
	phi = (complex *)malloc((ndimz+nex)*sizeof(complex));
	if(!phi) cout<<"\nerrror de asignacion en complex";
	
	rv  = (complex *)malloc((ndimz+nex)*sizeof(complex));
	if(!rv) cout<<"\nerrror de asignacion en complex";
	
	phil= (complex *)malloc((ndimz+nex)*sizeof(complex));
	if(!phil) cout<<"\nerrror de asignacion en complex";
	
	a= (complex *)malloc((ndimz+nex)*sizeof(complex));
	if(!a) cout<<"\nerrror de asignacion en complex";	

	b  = (complex *)malloc((ndimz+nex)*sizeof(complex));
	if(!b) cout<<"\nerrror de asignacion en complex";	
	
	c= (complex *)malloc((ndimz+nex)*sizeof(complex));
	if(!c) cout<<"\nerrror de asignacion en complex";	
	
	v   = (double *)malloc(((ndimz+nex))*sizeof(double));
	if(!v) cout<<"\nerrror de asignacion en complex";
	
	rho   = (double *)malloc(((ndimz+nex))*sizeof(double));
	if(!rho) cout<<"\nerrror de asignacion en complex";	
	
	drr   = (double *)malloc(((ndimz+nex))*sizeof(double));
	if(!drr) cout<<"\nerrror de asignacion en complex";		
	
	
	
	//Initializing at zero
	for(int i=0; i<(ndimz+nex); i++)
    {
		phi[i]	= complex(0.,0.); 
		a[i]    = phi[i];
		c[i]    = a[i];		
		v[i]	=  0.;
		rho[i]  =  0.;			
    }//End initializing
	
	
	
	//Refinition grid
	HankelMatrix H(ndimz,R);	
	for (int i=0;i<ndimz;i++)
		rho[i]=H.r[i];
	
	
	//NonUniform spacing
	for (int i=0; i<ndimz-1; i++)
		drr[i]=rho[i+1]-rho[i];
	drr[ndimz-1]=R - rho[ndimz-1];	
	
	
	cout << "\nR="<<R;
	cout << "   rmax="<<rho[ndimz-1];	
	cout << endl;
	
	
	//Parameters
	double rho1  = 100.;	
	double with0 = 15.;
	double h0    = 0.;	
	double v0    = 0.;
	double sig0    = 0.5;	
		
	//Gaussian initial condition
	for(int i=0;i<ndimz+nex;i++)
	{
		phi[i]     =  exp( - (rho[i] - rho0)*(rho[i] - rho0)/sig0/sig0 )*
						complex( cos( v0*(rho[i] - rho0) ),
								 sin( v0*(rho[i] - rho0) ) );//*/
		
		
		if (rho[i] >= rho1-with0/2. && rho[i] <= rho1+with0/2. )
			v[i]       = h0;  //1./sqrt(1+(rho[i] - rho1)*(rho[i] - rho1));//
		
		fprintf(out0," %e %e %12.10e \n", rho[i], v[i], drr[i]);
		
	}//End Gaussian
	
	
	//Making norm
	double norm0=0.0;
	for(int i=0;i<ndimz;i++)	
		norm0+= drr[i]*rho[i]*abs(phi[i])*abs(phi[i]);//dr*rho[i]*abs(phi[i])*abs(phi[i]);
	
	cout << "\nNorm= "<<norm0;
	cout << "   drr= "<<drr[0];
	
	
	//Making normalization
	for(int i=0;i<(ndimz+nex);i++)	
		phi[i]=phi[i]/sqrt(norm0);
		
	
	//Making norm
	norm0=0.0;
	for(int i=0;i<ndimz;i++)	
		norm0+= drr[i]*rho[i]*abs(phi[i])*abs(phi[i]);//dr*rho[i]*abs(phi[i])*abs(phi[i]);
	
	
    printf("\nNorm =%e\n",norm0);

	
	
	//Evaluating crank-nicholson formula	
	double hh0 = 0.;
	double dh0 = 0.;	
	double hh  = 0.;
	
	
	//Evaluation of diagonal up and diagonal down
	for (int i=0; i< ndimz-1; i++) {
		
		hh  = drr[i+1] + drr[i];
		dh0 = drr[i+1] - drr[i];		
		hh0 = drr[i+1]*drr[i+1] + drr[i]*drr[i];
		
		a[i]     =  ( complex( 0., - 1./hh0  ) +
				      complex( 0., - 1.*dh0/hh/hh0  ) +
				      complex( 0., + 1./2./hh/rho[i]  ) 
					 )*dt/2.;
		
		c[i]     =  ( complex( 0., - 1./hh0  ) +
				      complex( 0., + 1.*dh0/hh/hh0  ) +
					  complex( 0., - 1./2./hh/rho[i]  )
					 )*dt/2. ;
	}//*/

	int nend=ndimz-1;
	hh  = drr[nend] + drr[nend];
	dh0 = drr[nend] - drr[nend];		
	hh0 = drr[nend]*drr[nend] + drr[nend]*drr[nend];
	
	a[nend]     =  ( complex( 0., - 1./hh0  ) +
					 complex( 0., - 1.*dh0/hh/hh0  ) +
					 complex( 0., + 1./2./hh/rho[nend]  ) 
				   )*dt/2.;
	
	c[nend]     =  ( complex( 0., - 1./hh0  ) +
					  complex( 0., + 1.*dh0/hh/hh0  ) +
					  complex( 0., - 1./2./hh/rho[nend]  )
				   )*dt/2. ;
	//End evaluation of the diagonal up and down
	
	
	
	//Starting Left main diagonal
	for(int i=0; i<ndimz-1; i++){
		hh0 = drr[i+1]*drr[i+1] + drr[i]*drr[i];		
		b[i] =  complex( 1.0, 0.)+
				complex( 0.0, 2./hh0  + v[i] )*dt/2.;
	}
	hh0 = drr[nend]*drr[nend] + drr[nend]*drr[nend];		
	b[nend] =  complex( 1.0, 0.)+
				complex( 0.0, 2./hh0  + v[nend] )*dt/2.;
	//End of the main diagonal	
	//End left and Evaluation of crank's formula
	
	
	
	complex zero=complex(0.,0.);		
	//Time propagation loop
	for (int ktime=0; ktime<Ntime; ktime++) {
		
		
		for (int i=0; i<ndimz+nex; i++)
			phil[i]=complex(0.,0.);
						
		
		//Starting Right
		rv[0] =   conj( a[0] )*phi[1]
		        + conj( b[0] )*phi[0]  
				+ conj( c[0] )*phi[1];	
				
		for (int i=1; i<ndimz-1; i++) 	
			rv[i] = conj( a[i] )*phi[i-1]  +
					conj( b[i] )*phi[i]    + 
					conj( c[i] )*phi[i+1];
					
		rv[ndimz-1] = conj( a[ndimz-1]  )*phi[ndimz-2]  +
				      conj( b[ndimz-1]  )*phi[ndimz-1]  +
					  conj( c[ndimz-1]  )*phi[ndimz-2]*zero;
		
		//Finishing right
		
		
		//Solving Triagonal Matrix
		tridag(a,b,c,rv,phil,ndimz);
		

		for(int i=0;i<ndimz;i++)
			phi[i]=phil[i];
		
		
		norm0=0.;
		for(int i=0;i<ndimz;i++)
			norm0+= drr[i]*rho[i]*abs(phi[i])*abs(phi[i]);//Calculing norm
		
		//printf("Norm =%e    ErrorNorm =%e\n",norm0,abs(norm0-1.));	//Printing norm and error	
		fprintf(out3,"%10.17e\n",log10(fabs(1.-norm0)));  
		
		if((ktime%(Ntime/(snap-1)))==0){
			
			norm0=0.;
			for(int i=0;i<ndimz;i++)
				norm0+= drr[i]*rho[i]*abs(phi[i])*abs(phi[i]);//Calculing norm

			printf("Norm =%e    ErrorNorm =%e\n",norm0,abs(norm0-1.));	//Printing norm and error	
			//fprintf(out3,"%10.17e\n",log10(fabs(1.-norm0)));  
			
			
			for (int i=0; i<ndimz+nex; i++)
				fprintf(out1,"%10.17e \n", rho[i]*abs(phi[i]) ); //Save wave function multiply by rho axis
			
			
			fprintf(out2,"%10.17e \n", ktime*dt);
		}
		
		
	}//End time loop
	
	fprintf(out2,"%d \n %10.17e\n",Ntime,dt);
	
}




void tridag( complex a[], complex b[],complex c[], complex r[], complex u[], int n)
{
	int j;

	complex bet=b[0];//, *gam;
	
	
	if(b[0]==0.0) nrerror ("error 1 en tridag");
	u[0]=r[0]/bet;
	
	for(j=1;j<n;j++)
	{
		gam1[j] = c[j-1]/bet;
		bet     = b[j] - a[j]*gam1[j];
		
		if (bet==0.0) nrerror("error 2 en tridag");
		u[j]  = (r[j]-a[j]*u[j-1])/bet;
	}
	
	for(j=(n-2);j>=0;j--)
		u[j]-=gam1[j+1]*u[j+1];
	//free(gam);
}



void nrerror (char error_text[])
{
	fprintf(stderr, "numerical recipes runtime error..\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr,"... now exiting system..\n");
}

