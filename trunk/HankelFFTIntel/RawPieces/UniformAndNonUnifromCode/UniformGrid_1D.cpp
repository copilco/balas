/*
 *  testingCranck1D.cpp
 *
 *  Created by alexis on 02/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream.h>
#include <complex.h>
#include <alloc.h>
#include <stdio.h>

#define complex complex<double>

using namespace std;

#define ndimz 4000
#define dr 0.1
#define dt 0.01
#define nex 0

double rho0  = 100.;
int Ntime = 1000;
int snap  = 50;

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
	
	complex *phi,*rv,*phil,*a,*b,*c;
	double *v, *rho;

	
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
	

	
	//Doing zeros
	for(int i=0; i<(ndimz+nex); i++)
    {
		phi[i]	= complex(0.,0.); 
		a[i]    = phi[i];
		c[i]    = a[i];		
		v[i]	=  0.;
		rho[i]  =  0.;		
		
    }//zeros
	

	
	//Uniform Grid
	rho[0]=dr/2.;
	for (int i=1; i<ndimz; i++)
		rho[i] =rho[i-1]+dr ;
	//End uniform grid
	
		
	//Parameters
	double rho1  = 200.;		
	double with0 = 15.;
	double h0    = 0.;	
	double v0    = 0.;
	
	
	
	//Gaussian initial condition
	for(int i=0;i<ndimz+nex;i++)
	{
		phi[i]     =  exp( - (rho[i] - rho0)*(rho[i] - rho0) )*
						complex( cos( v0*(rho[i] - rho0) ),
								 sin( v0*(rho[i] - rho0) ) );	
		
		if (rho[i] >= rho1-with0/2. && rho[i] <= rho1+with0/2. )
			v[i]       = h0;  //1./sqrt(1+(rho[i] - rho1)*(rho[i] - rho1));//
		
		fprintf(out0," %e %e %12.10e \n", rho[i], v[i], dr);
		
	}//End Gaussian
	
	
	
	//Making norm
	double norm0=0.0;
	for(int i=0;i<ndimz;i++)	
		norm0+= dr*rho[i]*abs(phi[i])*abs(phi[i]);//dr*rho[i]*abs(phi[i])*abs(phi[i]);
	
	cout << "\nNorm: "<<norm0;
	cout << "\ndrr: "<<dr;
	
	
	//Making normalization
	for(int i=0;i<(ndimz+nex);i++)	
		phi[i]=phi[i]/sqrt(norm0);
	//End normalization
	
	
	
	//Making norm
	norm0=0.0;
	for(int i=0;i<ndimz;i++)	
		norm0+= dr*rho[i]*abs(phi[i])*abs(phi[i]);//dr*rho[i]*abs(phi[i])*abs(phi[i]);
	
	
	
    printf("\nNorm =%e\n",norm0);
	cout << endl;
	
	
	//Evaluating crank-nicholson formula	
	//Evaluation of up diagonal, main diagonal and down diagonal
	for (int i=0; i< ndimz; i++) {
		
		a[i]     =    complex( 0., - dt/4./dr/dr  ) +
				      complex( 0., + dt/8./dr/rho[i]  );
		
		b[i]	 =	  complex( 1.0, 0.)+
					  complex( 0.0, dt/2./dr/dr  + v[i]*dt/2. );
		
		c[i]     =    complex( 0., - dt/4./dr/dr  ) +
					  complex( 0., - dt/8./dr/rho[i]  );
	}//*/
	//End evaluation of the diagonal up and down
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
					  conj( c[ndimz-1]  )*phi[ndimz-1]*zero;		
		//Finishing right
		
		
		//Solving Triagonal Matrix
		tridag(a,b,c,rv,phil,ndimz);
		
		for(int i=0;i<ndimz;i++)
			phi[i]=phil[i];
		
		
		
		norm0=0.;
		for(int i=0;i<ndimz;i++)
			norm0+= dr*rho[i]*abs(phi[i])*abs(phi[i]);//Doing the norm
		
		//printf("Norm =%e    ErrorNorm =%e\n",norm0,abs(norm0-1.)); //Printing the norm and error			
		fprintf(out3,"%10.17e\n",log10(abs(1.-norm0)));
		
		
		if((ktime%(Ntime/snap))==0){
			
			norm0=0.;
			for(int i=0;i<ndimz;i++)
				norm0+= dr*rho[i]*abs(phi[i])*abs(phi[i]);//Doing the norm

			printf("Norm =%e    ErrorNorm =%e\n",norm0,abs(norm0-1.)); //Printing the norm and error			
			//fprintf(out3,"%10.17e\n",norm0);
			
			
			for (int i=0; i<ndimz+nex; i++)
				fprintf(out1,"%10.17e \n", abs(phi[i])*rho[i] );
			
			
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

