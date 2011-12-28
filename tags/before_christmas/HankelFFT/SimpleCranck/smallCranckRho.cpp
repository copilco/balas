/*
 *  smallCranck.cpp
 *  
 *
 *  Created by camilo Ruiz Méndez on 01/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream.h>
#include <complex.h>
#include <alloc.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define complex complex<double>
#include <fstream.h>
using namespace std;

#define ndimx 2000//2000 //4000


#define lightc 137.0
#define dospi 6.28318530717959
#define pi 3.1415926536

//void tridag(complex a,complex b[],complex c,complex r[],complex u[], int n);
void tridager(complex a00[], complex b[], complex c00[], complex r[], complex u[],int n);
void nrerror (char error_text[]);
int index(int ia,int ib);

int main(int arg)
{
	FILE *out0;
	
	out0 =fopen("movie.txt","w");
	
	/**********************************************************************/
    /****** Arrays declaration *************/
	
	complex *phic,*philx,*vlx,*sx,*ax,*cx;
	double *vr,*v;
	
	v=(double *)malloc(((ndimx+2) )*sizeof(double));
	if(!v) printf("\nerrror de asignacion en double");
	
	phic=(complex *)malloc(((ndimx+2) )*sizeof(complex));
	if(!phic) printf("\nerror de colocacion en la asignacion complex");
	
	philx=(complex *)malloc((ndimx+2)*sizeof(complex));
	if(!philx) printf("\nerrror de asignacion en double");
	
	vlx=(complex *)malloc((ndimx+2)*sizeof(complex));
	if(!vlx) printf("\nerrror de asignacion en double");
	
	sx=(complex *)malloc((ndimx+2)*sizeof(complex));
	if(!sx) printf("\nerror de colocacion en la asignacion complex");
	
	ax=(complex *)malloc((ndimx+2)*sizeof(complex));
	if(!sx) printf("\nerror de colocacion en la asignacion complex");
	
	cx=(complex *)malloc((ndimx+2)*sizeof(complex));
	if(!sx) printf("\nerror de colocacion en la asignacion complex");
	
	double dx=0.025;//.05/1024;//0.013;
	double dxx=dx*dx;
	double dt=0.005; //time step
	/****************************************************************************/
	
	for(int i=0;i<(ndimx+2);i++)
    {
		phic[i]=complex(0.,0.);
		v[i]=0.;
    }
	
	for(int i=1;i<=ndimx;i++)
    {
		double x=(i-0.5)*dx;// (-(ndimx+1)/2.+i)*dx;
		//v[i]=-2./sqrt(1.+x*x);
		phic[i]=exp(-2.*log(2)*x*x/36.);//*complex(3.*cos(x),3.*sin(x)) ;
		//exp(-2*log(2)*((x-x0)/FWHM).^2);
    }
	
	/*****************************************************************************/
	double norm1=0.0;
	for(int i=1;i<(ndimx+2);i++)
    {
		double x=(i-0.5)*dx;
		norm1+=x*dx*real((phic[i]*conj(phic[i])));
    }
	printf("\n norm readed of the ground state n=%10.9e",norm1);
	
	for(int i=0;i<(ndimx+2);i++)
		phic[i]/=sqrt(norm1);
	
	norm1=0.0;
	for(int i=1;i<(ndimx+2);i++)
    {
		double x=(i-0.5)*dx;
		norm1+=x*dx*real((phic[i]*conj(phic[i])));
    }
	printf("\n norm readed of the ground state n=%10.9e",norm1);
	
	/*****************************************************************************/
	//Exponential of the hamiltonian in x
	
	for (int h=1;h<5000;h++)
    {
				
		for(int i=1;i<=ndimx;i++)
		{
			double x=(i-0.5)*dx;
			
			ax[i]=complex(0,-dt/2./dxx+dt/4./dx/x);
			cx[i]=complex(0,-dt/2./dxx-dt/4./dx/x);
			
			vlx[i]=1.+complex(0.,dt/dxx);//+v[index(i,i)]*dt/4.;
			sx[i]=-ax[i]*phic[i-1]+conj(vlx[i])*phic[i]-cx[i]*phic[i+1];
		}
		
		
		ax[1]=complex(0.,-dt/dxx/1.);
		cx[1]=complex(0.,-dt/dxx/1.);
		vlx[1]= 1.+complex(0.,dt/dxx);//+v[index(1,i)]*dt/4.;
		sx[1]=conj(vlx[1])*phic[1]-cx[1]*phic[1+1];
		
		
		tridager(ax,vlx,cx,sx,philx,ndimx);
		
		
		//      zgtsv_(&nn,&nn,lower, middle, upp, b,&nn, &info);
		
		for(int i=1;i<=ndimx;i++)
			phic[i]=philx[i];
		
		
		norm1=0.0;
		for(int i=1;i<(ndimx+2);i++)
		{
			double x=(i-0.5)*dx;
			norm1+=x*dx*real((phic[i]*conj(phic[i])));
		}
		printf("\n norm readed of the ground state n=%10.9e",norm1);
		
		int skiper=10;
		if((h%100)==0)
			for(int i=0;i<(ndimx)/skiper;i++)
				fprintf(out0,"%e\n",dx*real(phic[i*skiper]*conj(phic[i*skiper])));
    }
	/***************************************************************************/
}

void tridager(complex a00[], complex b[], complex c00[], complex r[], complex u[],int n)
/*
 Solves for a vector u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n],
 b[1..n], c[1..n], and r[1..n] are input vectors and are not modiﬁed. */
{
	int j;
	complex bet,*gam;
	gam=new complex[(n+2)];
	if (b[1] == 0.0) nrerror("Error 1 in tridager **");
	u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++)
    {
		gam[j]=c00[j-1]/bet;
		bet=b[j]-a00[j]*gam[j];
		if (bet == 0.0) nrerror("Error 2 in tridager **"); //Algorithm fails; see below.
		u[j]=(r[j]-a00[j]*u[j-1])/bet;
    }
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
	free(gam);
}

void  nrerror (char error_text[])
{
    fprintf(stderr, "numerical recipes runtime error..\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr,"... now exiting system..\n");
}

