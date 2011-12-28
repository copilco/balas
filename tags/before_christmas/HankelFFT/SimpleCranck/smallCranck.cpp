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

#define ndimxs 2000//2000 //4000


#define lightc 137.0
#define dospi 6.28318530717959
#define pi 3.1415926536

void tridag(complex a,complex b[],complex c,complex r[],complex u[], int n);
void nrerror (char error_text[]);
int index(int ia,int ib);

int main(int arg)
{
  FILE *out0;
  
  out0 =fopen("movie.txt","w");
  
  /**********************************************************************/
    /****** Arrays declaration *************/
  
  complex *phic1,*philxs,*vlxs,*sxs;
  double *vr,*v;
  
  v=(double *)malloc(((ndimxs+2) )*sizeof(double));
  if(!v) printf("\nerrror de asignacion en double");
  
  phic1=(complex *)malloc(((ndimxs+2) )*sizeof(complex));
  if(!phic1) printf("\nerror de colocacion en la asignacion complex");
  
  philxs=(complex *)malloc((ndimxs+2)*sizeof(complex));
  if(!philxs) printf("\nerrror de asignacion en double");
  
  vlxs=(complex *)malloc((ndimxs+2)*sizeof(complex));
  if(!vlxs) printf("\nerrror de asignacion en double");
  
  sxs=(complex *)malloc((ndimxs+2)*sizeof(complex));
  if(!sxs) printf("\nerror de colocacion en la asignacion complex");
  
  double dxs=0.013;
  double dt=0.005; //time step
  /****************************************************************************/
  
  for(int i=0;i<(ndimxs+2);i++)
    {
      phic1[i]=complex(0.,0.);
      v[i]=0.;
    }
  
  for(int i=1;i<=ndimxs;i++)
    {
      double xs=(-(ndimxs+1)/2.+i)*dxs;
      v[i]=-2./sqrt(1.+xs*xs);
      phic1[i]=exp(-xs*xs);//*complex(3.*cos(xs),3.*sin(xs)) ;
    }
	
  /*****************************************************************************/
  double norm1=0.0;
  for(int j=0;j<(ndimxs+2);j++)
    norm1+=dxs*real((phic1[j]*conj(phic1[j])));
  printf("\n norm readed of the ground state n=%10.9e",norm1);
  
  for(int j=0;j<(ndimxs+2);j++)
    phic1[j]/=sqrt(norm1);
  
  norm1=0.0;
  for(int j=0;j<(ndimxs+2);j++)
    norm1+=dxs*real((phic1[j]*conj(phic1[j])));
  printf("\n norm readed of the ground state n=%10.9e",norm1);
  
  /*****************************************************************************/
  //Exponential of the hamiltonian in xs
  double campoxs=0.;
  
  complex axs=complex(-dt*campoxs/4./lightc/dxs,-dt/4./dxs/dxs);
  complex cxs=complex( dt*campoxs/4./lightc/dxs,-dt/4./dxs/dxs);
  
  for (int h=1;h<5000;h++)
    {
      
      for(int i=1;i<=ndimxs;i++)
	{
	  vlxs[i]=complex(1.,dt/2./dxs/dxs)+complex(0.,v[i])*dt/8.;
	  sxs[i]=-(axs)*phic1[i-1]+conj(vlxs[i])*phic1[i]-(cxs)*phic1[i+1];
	}
      
      tridag(axs,vlxs,cxs,sxs,philxs,ndimxs);
      int nn=ndimxs+2;
      
      //      zgtsv_(&nn,&nn,lower, middle, upp, b,&nn, &info);
      
      for(int i=1;i<=ndimxs;i++)
	phic1[i]=philxs[i];
      
      
	norm1=0.0;
	for(int j=0;j<(ndimxs+2);j++)
	  norm1+=dxs*real((phic1[j]*conj(phic1[j])));
	printf("\n norm readed of the ground state n=%10.9e",norm1);
	
	int skiper=10;
	if((h%100)==0)
	for(int j=0;j<(ndimxs)/skiper;j++)
	  fprintf(out0,"%e\n",dxs*real(phic1[j*skiper]*conj(phic1[j*skiper])));
    }
  /***************************************************************************/
}


