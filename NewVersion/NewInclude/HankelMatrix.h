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
#include <vector>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "arrai.h"
#include "constant.h"

#include <iostream>
#include <fstream>


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
	double *kr;
	complex *C;
	double *Creal;
	double *m1;
	double *m2;
    
    double *dr;
    double *dv;
	double *dkr;
	
	FILE *in0;
		
	
	//////////////////////
	// Function members //
	//////////////////////
	
	// Constructor
	
	HankelMatrix(int _Nr,double _R)
	{
				
		initialize( _Nr,  _R);
		
		/*
		Nr=_Nr;
		R=_R;
		ord=0.; //For this application we will only use the zeros order solutions
		
		nn=5;
		mm=3001;
		zer=nn*mm;
		zeros= (double*)mkl_malloc(1*zer*sizeof(double),16);
		
		/****************************************************************/
		// Binary file
		/*
		fstream myfile ("besselZeros.bin", ios::in | ios::binary);
		int size = nn*mm;// myfile.tellg();		
		myfile.read ((char*)&zeros[0], sizeof (double)*(size) );		
		myfile.close();
		
		/****************************************************************/
		/*
		
		//Fill the zeros onto the c vector
		c=(double*)mkl_malloc(1*(Nr+1)*sizeof(double),16);
		
		for(int j=0;j<Nr+1;j++)
		{
			c[j]=zeros[j*nn+(ord)];
		};
		
		/****************************************************************
		for (int i=0;i<31;i++)
		{
			printf("\n");
			for (int j=0;j<nn;j++)
				printf("%4.25e ", zeros[j*nn+i]);
			
		}
		/****************************************************************/
		/*
		//Maximum frequency
		double	V = c[Nr]/(2.*pi*R);   
		
		//The axes for the rho, the v coordinate (spatial frequency), and the kr coordinate (wavevector)
		r=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		v=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		kr=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		
		
		for(int i=0;i<Nr;i++)
		{
			r[i] = c[i]*R/c[Nr];   // Radius vector
			v[i] = c[i]/(2.*pi*R);
			kr[i] = dospi*v[i];
		};
		
        dr	= (double*)mkl_malloc(Nr*sizeof(double),16);
        dv	= (double*)mkl_malloc(Nr*sizeof(double),16);
        dkr = (double*)mkl_malloc(Nr*sizeof(double),16);
		
		//The definition of the diferential of r is dr[0]=r[i+1]-r[i]
		//This mean that we need to provide a last number by hand, this is done in the next two lines.
		
        for(int i=0;i<Nr-1;i++)
        {
            dr[i]=r[i+1]-r[i];
            dv[i]=v[i+1]-v[i];
        };
        
        dr[Nr-1]=R-r[Nr-1];       //This way the array has Nr elements
        //dv[Nr-1]=c[Nr]/(2*pi*R);  
		dv[Nr-1]=dv[Nr-2];
		
		for(int i=0;i<Nr;i++)
			dkr[i]=dospi*dv[i];
		
		//The bessel functions and the matrix itself
		C=(complex*)mkl_malloc(Nr*Nr*sizeof(complex),16);
		for(int j=0;j<Nr*Nr;j++)
			C[j]=complex(0.,0.);
		

		//Auxiliar array for exproting to binary file
		Creal=(double*)mkl_malloc(Nr*Nr*sizeof(double),16);

	
		for(int j=0;j<Nr;j++)
		{
			for(int i=0;i<Nr;i++)
			{
				real(C[j*Nr+i])=(2./c[Nr])*j0(c[j]*c[i]/c[Nr])/fabs(j1(c[i]))/fabs(j1(c[j]));
				Creal[j*Nr+i]=real(C[j*Nr+i]);
			};
		};
		
	
		
		m1=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		m2=(double*)mkl_malloc(1*Nr*sizeof(double),16);

		
		for(int i=0;i<Nr;i++)
		{
			m1[i]=fabs(j1(c[i]))/R;
			m2[i]=fabs(j1(c[i]))/V;
		};//*/
		
		
	};//End Constructor 
	
	
	// Destructor	
	~HankelMatrix()
	{
		free_memory();
		
	};
	
	
	void initialize( int _Nr, double _R)	
	{
	
		Nr=_Nr;
		R=_R;
		ord=0.; //For this application we will only use the zeros order solutions
		
		nn=5;
		mm=3001;
		zer=nn*mm;
		zeros= (double*)mkl_malloc(1*zer*sizeof(double),16);
		
		/****************************************************************/
		// Binary file
		
		fstream myfile ("besselZeros.bin", ios::in | ios::binary);
		int size = nn*mm;// myfile.tellg();		
		myfile.read ((char*)&zeros[0], sizeof (double)*(size) );		
		myfile.close();
		
		/****************************************************************/
		
		
		//Fill the zeros onto the c vector
		c=(double*)mkl_malloc(1*(Nr+1)*sizeof(double),16);
		
		for(int j=0;j<Nr+1;j++)
			c[j]=zeros[j*nn+(ord)];
		
		/****************************************************************
		 for (int i=0;i<31;i++)
		 {
		 printf("\n");
		 for (int j=0;j<nn;j++)
		 printf("%4.25e ", zeros[j*nn+i]);
		 
		 }
		 /****************************************************************/
		
		//Maximum frequency
		double	V = c[Nr]/(2.*pi*R);   
		
		//The axes for the rho, the v coordinate (spatial frequency), and the kr coordinate (wavevector)
		r=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		v=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		kr=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		
		
		for(int i=0;i<Nr;i++)
		{
			r[i] = c[i]*R/c[Nr];   // Radius vector
			v[i] = c[i]/(2.*pi*R);
			kr[i] = dospi*v[i];
		};
		
        dr	= (double*)mkl_malloc(Nr*sizeof(double),16);
        dv	= (double*)mkl_malloc(Nr*sizeof(double),16);
        dkr = (double*)mkl_malloc(Nr*sizeof(double),16);
		
		//The definition of the diferential of r is dr[0]=r[i+1]-r[i]
		//This mean that we need to provide a last number by hand, this is done in the next two lines.
		
        for(int i=0;i<Nr-1;i++)
        {
            dr[i]=r[i+1]-r[i];
            dv[i]=v[i+1]-v[i];
        };
        
        dr[Nr-1]=R-r[Nr-1];       //This way the array has Nr elements
        //dv[Nr-1]=c[Nr]/(2*pi*R);  
		dv[Nr-1]=dv[Nr-2];
		
		for(int i=0;i<Nr;i++)
			dkr[i]=dospi*dv[i];
		
		//The bessel functions and the matrix itself
		C=(complex*)mkl_malloc(Nr*Nr*sizeof(complex),16);
		for(int j=0;j<Nr*Nr;j++)
			C[j]=complex(0.,0.);
			
			
			//Auxiliar array for exproting to binary file
			Creal=(double*)mkl_malloc(Nr*Nr*sizeof(double),16);
			
			
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nr;i++)
				{
					real(C[j*Nr+i])=(2./c[Nr])*j0(c[j]*c[i]/c[Nr])/fabs(j1(c[i]))/fabs(j1(c[j]));
					Creal[j*Nr+i]=real(C[j*Nr+i]);
				};
		
		
		
		m1=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		m2=(double*)mkl_malloc(1*Nr*sizeof(double),16);
		
		
		for(int i=0;i<Nr;i++)
		{
			m1[i]=fabs(j1(c[i]))/R;
			m2[i]=fabs(j1(c[i]))/V;
		};
		
	};//End of initialize
	
	
	void getAxisBin()
	{
		
		//////////////////////////////////////////////////////////
		//Exporting to a bin file the values of the Hankel axis
		//////////////////////////////////////////////////////////
		
		
		fstream raxis ("HankelRhoaxis.bin", ios::out | ios::binary);
		raxis.write ((char*)&r[0], sizeof (double)*(Nr) );
		raxis.close();
		
		fstream vaxis ("HankelVaxis.bin", ios::out | ios::binary);
		vaxis.write ((char*)&v[0], sizeof (double)*(Nr) );
		vaxis.close();
		
		
		
	};
	
	void getHankelMatrixBin()
	{
		//////////////////////////////////////////////////////////
		//Exporting to a bin file the values of the Hankel matrix.
		//////////////////////////////////////////////////////////
		
		fstream beselMat ("BesselMatrix.bin", ios::out | ios::binary);
		beselMat.write ((char*)&Creal[0], sizeof (double)*(Nr*Nr) );
		beselMat.close();
		
		
		//////////////////////////////////////////////////////////
		//Exporting to a bin file the values of the Hankel axis
		//////////////////////////////////////////////////////////
		
		
		fstream m1file ("Hankelm1.bin", ios::out | ios::binary);
		m1file.write ((char*)&m1[0], sizeof (double)*(Nr) );
		m1file.close();
		
		fstream m2file ("Hankelm2.bin", ios::out | ios::binary);
		m2file.write ((char*)&m2[0], sizeof (double)*(Nr) );
		m2file.close();
		
		
	};
	
	
	void resize_HankelMatrix(int _Nr, double _R)
	{
		free_memory();
		initialize(_Nr, _R);
	};
	
	
	void free_memory()
	{
		mkl_free(zeros);
		mkl_free(c);
		mkl_free(r);
		mkl_free(v);
		mkl_free(kr);
		mkl_free(C);
		mkl_free(Creal);
		mkl_free(m1);
		mkl_free(m2);
		mkl_free(dr);
		mkl_free(dv);
		mkl_free(dkr);			
				
	};
	
	
	
};

#endif //  HANKELMATRIX_H