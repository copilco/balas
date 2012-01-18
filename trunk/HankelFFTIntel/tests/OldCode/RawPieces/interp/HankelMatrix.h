/*
 *  HankelMatrix.h
 *  
 *
 *  Created by Alejandro de la Calle on 26/09/11.
 *  
 *
 */

#ifndef HANKELMATRIX_H
#define HANKELMATRIX_H

#include <iostream>
#include <math.h>
#include <new>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "arrai.h"
#include "constant.h"


using namespace std;


class HankelMatrix
{
	
public:
	
	int Nr;
	double R;
	int ord;
	
	int nn;
	int mm;
	int zer;
	
	double	V;
	
	double *zeros;
	double *c;
	double *r;
	double *v;
	complex *C;
	double *m1;
	double *m2;
    
    double *dr;
    double *dv;
	
	FILE *in0;
		
	
	//////////////////////
	// Function members //
	//////////////////////
	
	// Constructor
	
	HankelMatrix(int _Nr,double _R)
	{
		Nr=_Nr;
		R=_R;
		ord=0.; //For this application we will only use the zeros order solutions
		
		nn=5;
		mm=3001;
		zer=nn*mm;
		
		//Read the zeros of the Bessel function
		in0=fopen("besselCeros.txt","r+");
		zeros= (double*)mkl_malloc(1*zer*sizeof(double),16);
		double f;
		
		//Read the zeros in the matrix called zeros, useful for higher order transformations
		for(int j=0;j<mm;j++)
		{
			for(int i=0;i<nn;i++)
			{
				fscanf(in0,"%lf",&f);
				zeros[j*nn+i]=f;
			}
			
			fscanf(in0,"\n");
			
		}
		
		//Fill the zeros onto the c vector
		c=(double*)mkl_malloc(1*(Nr+1)*sizeof(double),16);
		
		for(int j=0;j<Nr+1;j++)
		{
			c[j]=zeros[j*nn+(ord)];
		}
		
		
		//Maximum frequency
		double	V = c[Nr]/(2*pi*R);   
		
		//The axes for the rho and the krho coordinate
		r=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		v=(double*)mkl_malloc(1*Nr*sizeof(double),16);

		
		for(int i=0;i<Nr;i++)
		{
			r[i] = c[i]*R/c[Nr];   // Radius vector
			v[i] = c[i]/(2*pi*R);
			
		}
		
        dr=(double*)mkl_malloc(Nr*sizeof(double),16);
        dv=(double*)mkl_malloc(Nr*sizeof(double),16);
        
        for(int i=0;i<Nr-1;i++)
        {
            dr[i]=r[i+1]-r[i];
            dv[i]=v[i+1]-v[i];
        }
        
        dr[Nr-1]=R-r[Nr-1];
        dv[Nr-1]=c[Nr]/(2*pi*R);
        
		
		//The bessel functions and the matrix itself
		C=(complex*)mkl_malloc(Nr*Nr*sizeof(complex),16);
	
		for(int j=0;j<Nr;j++)
		{
			for(int i=0;i<Nr;i++)
			{
				real(C[j*Nr+i])=(2./c[Nr])*j0(c[j]*c[i]/c[Nr])/fabs(j1(c[i]))/fabs(j1(c[j]));
			}
		}
		
		m1=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		m2=(double*)mkl_malloc(1*Nr*sizeof(double),16);

		
		for(int i=0;i<Nr;i++)
		{
			m1[i]=fabs(j1(c[i]))/R;
			m2[i]=fabs(j1(c[i]))/V;
		}
		
	}
	
	// Destructor
	
	~HankelMatrix()
	{
		mkl_free(zeros);
		mkl_free(c);
		mkl_free(r);
		mkl_free(v);
		mkl_free(C);
		mkl_free(m1);
		mkl_free(m2);
		mkl_free(dr);
		mkl_free(dv);
		
	}
	
	
	
};

#endif //  HANKELMATRIX_H