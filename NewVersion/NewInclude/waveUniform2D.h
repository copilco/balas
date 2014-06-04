//
//  waveUniform2D.h
//  
//
//  Created by Alexis Chacón on 29/04/2014.
//  
//
//

#ifndef WAVEUNIFORM2D_H
#define WAVEUNIFORM2D_H


#include <complex>
#include <math.h>
#include <stdio.h>
#include "silo.h"
#include <new>
#include "omp.h"
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "tools.h"


class waveUniform2D
{
	
public:
	
	int Nr,Nz;
	double *r,dr;
	double *z,dz;
	double R;
	double e_charge;
	
	complex *phi;
	double *pot;
	int Nodes;
	complex az, cz, *ar, *cr;
	complex *rv_r,*rv_z;
	complex *phi_r;
	complex *br, *gamr;
	complex *bz, *gamz;
	
	//double **coords;
	double *xcoord;//		= ( double* ) malloc( NX*NY*NZ*sizeof(double) );
	double *ycoord;//		= ( double* ) malloc( NX*NY*NZ*sizeof(double) );
	double *zcoord;//		= ( double* ) malloc( NX*NY*NZ*sizeof(double) );		
	double *edensity;//      = ( double* ) malloc( NX*NY*NZ*sizeof(double) );	
	
	//Constructor 
	waveUniform2D()
	{
		r		= NULL;
		z		= NULL;
		
		xcoord	= NULL;
		ycoord	= NULL;
		zcoord	= NULL;
		edensity= NULL;
		
		phi		= NULL;
		pot		= NULL;		
		
		rv_r	= NULL;
		rv_z	= NULL;
		
		phi_r	= NULL;
		br		= NULL; 
		gamr	= NULL;
		bz		= NULL;
		gamz	= NULL;
		
	};
	
	/***********************************************************/	
	//	Fill the arrays and copy the grids in rho and k_rho
	/***********************************************************/	
	inline int index(int j,int i)
	{
		int indexer=j*Nz+i;
		return indexer;
	};
	
	int index3D(int k, int j, int i, int NY, int NZ)
	{
		return k*NY*NZ+j*NZ+i;
	};	
	
	void initialize(HankelMatrix &HH,int _Nz,double _dz)
	{
		Nr	= HH.Nr;
		R	= HH.R;
		dr	= R/( (double)(Nr));
		
		Nz	=_Nz;
		dz	=_dz;
		
		e_charge=-1.;
		Nodes   = 1;
#pragma omp parallel
		Nodes	= omp_get_num_threads();
		
		
		r		= (double*)mkl_malloc(  Nr*sizeof(double),16);
		z		= (double*)mkl_malloc(  Nz*sizeof(double),16);
		
		//*coords	= (double*)mkl_malloc( 3*sizeof(double),16 );
		
		
		phi		= (complex*)mkl_malloc( Nz*Nr*sizeof(complex),16);
		pot		= (double*)mkl_malloc(  Nz*Nr*sizeof(double),16);
		 		
		
		ar		= (complex*)mkl_malloc(Nr*sizeof(complex),16);
		cr		= (complex*)mkl_malloc(Nr*sizeof(complex),16);
		
		
		rv_z	= (complex*)mkl_malloc(Nodes*Nz*sizeof(complex),16);
		bz		= (complex*)mkl_malloc(Nodes*Nz*sizeof(complex),16);
		gamz	= (complex*)mkl_malloc(Nodes*Nz*sizeof(complex),16);
		
		
		
		rv_r	= (complex*)mkl_malloc(Nodes*Nr*sizeof(complex),16);
		phi_r	= (complex*)mkl_malloc(Nodes*Nr*sizeof(complex),16);
		br		= (complex*)mkl_malloc(Nodes*Nr*sizeof(complex),16);
		gamr	= (complex*)mkl_malloc(Nodes*Nr*sizeof(complex),16);		

	
		memset(r,		0,	sizeof(double)*Nr);
		memset(z,		0,	sizeof(double)*Nz);
				
		memset(phi,		0,	sizeof(complex)*Nr*Nz);
		memset(pot,		0,	sizeof(double)*Nr*Nz );
		
		memset(ar,		0,	sizeof(complex)*Nr);
		memset(cr,		0,	sizeof(complex)*Nr);
		
		memset(rv_r,	0, sizeof(complex)*Nr*Nodes);
		memset(phi_r,	0, sizeof(complex)*Nr*Nodes);
		memset(br,		0, sizeof(complex)*Nr*Nodes);	
		memset(gamr,	0, sizeof(complex)*Nr*Nodes);			
		
		memset(rv_z,	0, sizeof(complex)*Nz*Nodes);
		memset(bz,		0, sizeof(complex)*Nz*Nodes);	
		memset(gamz,	0, sizeof(complex)*Nz*Nodes);	
		
		
		
		r[0]=dr/2.;
		
		for (int j=1;j<Nr;j++)
			r[j] =r[j-1]+dr;
		

		for (int i=0;i<Nz;i++)
		{			
			//z[i] =(-double(Nz/2.) + i)*dz;
			z[i] =(-double(Nz-1)/2. + i)*dz;
		};
	
		
	};
	
	
	~waveUniform2D()
	{
		free_memory();
	};//*/
	
	/***********************************************************/	
	//	      	   Norms in both spaces
	/***********************************************************/	
	
	
	double norm()
	{
		double norm=0.0;
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				norm+=dz*dr*r[j]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return norm;
		
	};
	
	
	
	
	//***********************************************************//	
	//		   	Normalize function
	//***********************************************************//
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr*Nz;i++)
			phi[i]=phi[i]/sqrt(norm1);
		
	};
	
	
	//***********************************************************
	//	       	Preparing the arrays fof the Crank
	//************************************************************
	void PrepareCrankArrays(complex dt )
	{
		int chunk = Nodes;
									
		
//#pragma omp parallel for shared(chunk,dt) schedule(dynamic,chunk)
		for( int i=0; i<Nr; i++ )
		{						
			ar[i]     =    complex( 0., - 1./dr/dr/4. + 1./dr/r[i]/8. )*dt;				
			cr[i]     =    complex( 0., - 1./dr/dr/4. - 1./dr/r[i]/8. )*dt;
		};//End left part
		
		

	};
	
	
	
	
	
	
	//***********************************************************
	//	       	Preparing the arrays fof the Crank
	//************************************************************
	void PrepareCrankArraysOnRho(complex dt )
	{
		for( int i=0; i<Nr; i++ )
		{						
			ar[i]     =    complex( 0., - 1./dr/dr/4. + 1./dr/r[i]/8. )*dt;				
			cr[i]     =    complex( 0., - 1./dr/dr/4. - 1./dr/r[i]/8. )*dt;
		};//End left part		
				
	};	
	
	
	
	
	
	void PrepareCrankArraysOnZ(complex dt )
	{
		az = complex(0.,-1./dz/dz/4.)*dt;
		cz = complex(0.,-1./dz/dz/4.)*dt;
		
	};		
	
	
	
	
	//***************************************************
	//    	    Z propagator imaginary time
	//***************************************************	
	void Zprop( complex dt )	
	{	

#pragma omp parallel
		{
			int number_threads	= omp_get_num_threads();
			int id_thread		= omp_get_thread_num();
			int Nr_thread		= ceil(Nr / number_threads);
			int Nr_i			= id_thread * Nr_thread;
			int Nr_f			= Nr_i + Nr_thread;
			if (id_thread+1 == number_threads) Nr_f = Nr;
			
			for (int j=Nr_i; j<Nr_f; j++ )
			{      			
				
				rv_z[id_thread*Nz]	= ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt )*phi[index(j,0)] + 
										- cz*phi[index(j,1)]; //Big vector which contains all the nodes paralel problem
				
				bz[id_thread*Nz]	=  1. + complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt;			
				
				
				for (int i=1; i<Nz-1; i++)
				{
					bz[i+id_thread*Nz]		=  1. + complex( 0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt;
					rv_z[i+id_thread*Nz]	= - az*phi[index(j,i-1)]+ 
											  + ( 1. - complex(0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt )*phi[index(j,i)] 
					                          - cz*phi[index(j,i+1)];
				};
				
				
				
				bz[Nz-1+id_thread*Nz]		=  1. + complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt;								
				rv_z[Nz-1+id_thread*Nz]		=	- az*phi[index(j,Nz-2)]  
												+ ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix for Z
				trid_simple( az, &bz[id_thread*Nz], cz, 
							&rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz );
				

			};//End loop Nr and Nz
		
		};//End of the pragma
		
	};//Z propagator imaginary time	
	
	//*/
	
	
	
	
	
	
	
	//*******************************************************
	//	       RHO_Operator imaginary time
	//*******************************************************
	
	
	void Rprop( complex dt )
	{		

		
#pragma omp parallel
		{
			int number_threads	= omp_get_num_threads();
			int id_thread		= omp_get_thread_num();
			int Nz_thread		= ceil(Nz / number_threads);
			int Nz_i			= id_thread * Nz_thread;
			int Nz_f			= Nz_i + Nz_thread;
			if (id_thread+1 == number_threads) Nz_f = Nz;
			
			//  RHO_Operator dt	
			for ( int i=Nz_i; i<Nz_f; i++ )
			{
				
				//tid = omp_get_thread_num();
				
				//printf("Thread number %d is doing loop number %d.\n",tid,i);
				
				rv_r[id_thread*Nr]		=	-ar[0]*phi[index(1,i)] +
											(1.  -  complex(0., 1./dr/dr/2. + pot[index(0,i)]/4. )*dt)*phi[index(0,i)]
											-cr[0]*phi[index(1,i)];			
				
				br[id_thread*Nr]		= 1. + complex(0., 1./2./dr/dr + pot[index(0,i)]/4. )*dt;			
				
				
				for (int j=1; j<Nr-1; j++)
				{ 
					br[id_thread*Nr+j]		= 1. + complex(0., 1./2./dr/dr + pot[index(j,i)]/4. )*dt;	
					
					rv_r[id_thread*Nr+j]	=  -ar[j]*phi[index(j-1,i)] + 
												(1. -  complex( 0., 1./dr/dr/2. + pot[index(j,i)]/4. )*dt)*phi[index(j,i)]
											   -cr[j]*phi[index(j+1,i)];
				};
				
				br[Nr-1+id_thread*Nr]   =  1. + complex(0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt;			
				
				rv_r[Nr-1+id_thread*Nr] =  -ar[Nr-1]*phi[index(Nr-2,i)]  +
											(1.  -  complex( 0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt)*phi[index(Nr-1,i)];
				
				
				//Solving Triagonal Matrix
				tridag(&ar[0], &br[id_thread*Nr], &cr[0], &rv_r[id_thread*Nr], &phi_r[id_thread*Nr], &gamr[id_thread*Nr], Nr);
				//tridag(ar, &br[id_thread*Nr], cr, &rv_r[id_thread*Nr], &phi_r[id_thread*Nr], &gamr[id_thread*Nr], Nr);
				
				for (int j=0; j<Nr; j++)
					phi[index(j,i)] = phi_r[id_thread*Nr+j];
				
				
			};//End loop on Nz and Nr
									
		};//End of the pragma
		
	};//End Rho_propagator imaginary time
	
	//*/
	
	
	
	
	
	//**********************************************************
	//	     Z propagator in velocity gauge or p·A gauge
	//**********************************************************
	void Zprop_PAG( complex dt, double av_z )
	{
		
		double avsquare = av_z*av_z;		
		
		az = complex( + e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;
		cz = complex( - e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;		
		
		
#pragma omp parallel
		{
			int number_threads	= omp_get_num_threads();
			int id_thread		= omp_get_thread_num();
			int Nr_thread		= ceil(Nr / number_threads);
			int Nr_i			= id_thread * Nr_thread;
			int Nr_f			= Nr_i + Nr_thread;
			if (id_thread+1 == number_threads) Nr_f = Nr;
			
			for (int j=Nr_i; j<Nr_f; j++ )
			{
				
				bz[id_thread*Nz]     =   1. + complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare)/4. )*dt;					
				rv_z[id_thread*Nz]	 =  (1. - complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare )/4. )*dt)*phi[index(j,0)] 
										- cz*phi[index(j,1)];	
				
				
				for (int i=1; i<Nz-1; i++) 
				{
					
					bz[i+id_thread*Nz]		=   1. + complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt;				
					
					rv_z[i+id_thread*Nz]	=	-az*phi[index(j,i-1)]  + 
												(1. - complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt )*phi[index(j,i)] 
												-cz*phi[index(j,i+1)];
					
				};
				
				bz[Nz-1+id_thread*Nz]    =   1. + complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4. )*dt;				
				rv_z[Nz-1+id_thread*Nz]	=	-az*phi[index(j,Nz-2)]  +
											(1. - complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4.)*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, &bz[id_thread*Nz], cz, &rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz);				
			};
			
		}; //End pragma 
	
	};//Z propagator imaginary time	
	
	

	

	//****************************************************************//
	//          Z propagator in length gauge or r·E gauge             
	//****************************************************************//	
	void Zprop_REG( complex dt, double efield_z )
	{
		
		
#pragma omp parallel
		{
			int number_threads	= omp_get_num_threads();
			int id_thread		= omp_get_thread_num();
			int Nr_thread		= ceil(Nr / number_threads);
			int Nr_i			= id_thread * Nr_thread;
			int Nr_f			= Nr_i + Nr_thread;
			if (id_thread+1 == number_threads) Nr_f = Nr;
			

			for (int j=Nr_i; j<Nr_f; j++ )
			{
				
				bz[id_thread*Nz]     =    1. + complex( 0., 1./2./dz/dz + (pot[index(j,0)] - 2.*e_charge*z[0]*efield_z )/4. )*dt;							
				rv_z[id_thread*Nz]   =   (1. - complex( 0., 1./2./dz/dz + (pot[index(j,0)] - 2.*e_charge*z[0]*efield_z )/4. )*dt )*phi[index(j,0)] 
										- cz*phi[index(j,1)];
				
				
				for (int i=1; i<Nz-1; i++)
				{
					bz[i+id_thread*Nz]	 = 1. + complex( 0., 1./2./dz/dz + (pot[index(j,i)] - 2.*e_charge*z[i]*efield_z )/4. )*dt;				
					
					rv_z[i+id_thread*Nz] =	- az*phi[index(j,i-1)]  + 
											+ (1. - complex( 0., 1./2./dz/dz + (pot[index(j,i)] - 2.*e_charge*z[i]*efield_z )/4. )*dt )*phi[index(j,i)] 
											- cz*phi[index(j,i+1)];
				};
				
				
				bz[Nz-1+id_thread*Nz]    =   1. + complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - 2.*e_charge*z[Nz-1]*efield_z )/4. )*dt;							
				rv_z[Nz-1+id_thread*Nz]	 =  - az*phi[index(j,Nz-2)]  
											+ (1. - complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - 2.*e_charge*z[Nz-1]*efield_z )/4. )*dt )*phi[index(j,Nz-1)];
				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, &bz[id_thread*Nz], cz, &rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz);
				
				
			};//end loop Nr and Nz
			
		};//End of the pragma 
		
		
		
	};//Z propagator imaginary time		
	
	
	
	
	
	void parallel_evolution( complex dtev)
	{
		complex dt=dtev;
		
#pragma omp parallel
		{		
			int j, i;
				
			int id_thread= omp_get_thread_num();
			
			dt=dtev/2.;			
			//Evolution on Z				
#pragma omp for private(j)
			//int id_thread= omp_get_thread_num();
			for (j=0; j<Nr; j++ )
			{      			
			
				
				rv_z[id_thread*Nz]	= ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt )*phi[index(j,0)] + 
									       - cz*phi[index(j,1)]; //Big vector which contains all the nodes paralel problem
				
				bz[id_thread*Nz]	=  1. + complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt;			
				
				for ( i=1; i<Nz-1; i++)
				{
					bz[i+id_thread*Nz]		=  1. + complex( 0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt;
					rv_z[i+id_thread*Nz]	= - az*phi[index(j,i-1)]+ 
											+ ( 1. - complex(0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt )*phi[index(j,i)] 
											- cz*phi[index(j,i+1)];
				};
				
				
				
				bz[Nz-1+id_thread*Nz]		=  1. + complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt;								
				rv_z[Nz-1+id_thread*Nz]		=	- az*phi[index(j,Nz-2)]  
											+ ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix for Z
				trid_simple( az, &bz[id_thread*Nz], cz, 
							 &rv_z[id_thread*Nz],  &phi[index(j,0)],  &gamz[id_thread*Nz], Nz );
				
				
			};//End loop Nr and Nz
			//}
			
			dt=dtev;
			//Evolution on RHO
#pragma omp for private(i)			
			//  RHO_Operator dt	
			for ( i=0; i<Nz; i++ )
			{
				
				//tid = omp_get_thread_num();
				
				//printf("Thread number %d is doing loop number %d.\n",tid,i);
				
				rv_r[id_thread*Nr]		=	-ar[0]*phi[index(1,i)] +
											+(1.  -  complex(0., 1./dr/dr/2. + pot[index(0,i)]/4. )*dt)*phi[index(0,i)]
											-cr[0]*phi[index(1,i)];			
				
				br[id_thread*Nr]		= 1. + complex(0., 1./2./dr/dr + pot[index(0,i)]/4. )*dt;			
				
				
				for ( j=1; j<Nr-1; j++)
				{ 
					br[j+id_thread*Nr]		= 1. + complex(0., 1./2./dr/dr + pot[index(j,i)]/4. )*dt;	
					
					rv_r[j+id_thread*Nr]		=   -ar[j]*phi[index(j-1,i)] + 
													(1. -  complex( 0., 1./dr/dr/2. + pot[index(j,i)]/4. )*dt)*phi[index(j,i)]
												-cr[j]*phi[index(j+1,i)];
				};
				
				br[Nr-1+id_thread*Nr]   =  1. + complex(0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt;			
				
				rv_r[Nr-1+id_thread*Nr] =  -ar[Nr-1]*phi[index(Nr-2,i)]  +
											(1.  -  complex( 0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt)*phi[index(Nr-1,i)];
				
				
				//Solving Triagonal Matrix
				tridag(&ar[0], &br[id_thread*Nr], &cr[0], &rv_r[id_thread*Nr], &phi_r[id_thread*Nr], &gamr[id_thread*Nr], Nr);
				
				
				for (j=0; j<Nr; j++)
					phi[index(j,i)] = phi_r[j+id_thread*Nr];
				
				
			};//End loop on Nz and Nr			
			
			
			
			dt=dtev/2.;
			//Evolution on Z
#pragma omp for private(j)		
			for (j=0; j<Nr; j++ )
			{      			
				
				rv_z[id_thread*Nz]	= ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt )*phi[index(j,0)] + 
				- cz*phi[index(j,1)]; //Big vector which contains all the nodes paralel problem
				
				bz[id_thread*Nz]	=  1. + complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt;			
				
				for ( i=1; i<Nz-1; i++)
				{
					bz[i+id_thread*Nz]		=  1. + complex( 0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt;
					rv_z[i+id_thread*Nz]	=	- az*phi[index(j,i-1)]+ 
												+ ( 1. - complex(0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt )*phi[index(j,i)] 
												- cz*phi[index(j,i+1)];
				};
				
				
				
				bz[Nz-1+id_thread*Nz]		=  1. + complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt;								
				rv_z[Nz-1+id_thread*Nz]		=	- az*phi[index(j,Nz-2)]  
												+ ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix for Z
				trid_simple( az, &bz[id_thread*Nz], cz, 
							&rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz );
				
				
			};//End loop Nr and Nz	//*/
			
			
		};//End of pragma
	};

	
	
	
	
	
	void parallel_evolution_laser_PA( complex dtev, double av_z)
	{
		complex dt=dtev/2.;
		
		double avsquare = av_z*av_z;		
		
		az = complex( + e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;
		cz = complex( - e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;			
		
#pragma omp parallel
		{		
			int j, i;
			
					
			int id_thread= omp_get_thread_num();
			
			
			//dt=dtev/2.;			
			//Evolution on Z				
#pragma omp for private(j)
				for (j=0; j<Nr; j++ )
				{
					
					bz[id_thread*Nz]     =   1. + complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare)/4. )*dt;	//Big vector which contains all the nodes paralel problem
					
					rv_z[id_thread*Nz]	 =  (1. - complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare )/4. )*dt)*phi[index(j,0)] 
											- cz*phi[index(j,1)]; //Big vector which contains all the nodes paralel problem	
					
					
					for (i=1; i<Nz-1; i++) 
					{
						
						bz[i+id_thread*Nz]		=   1. + complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt;				
						
						rv_z[i+id_thread*Nz]	=	-az*phi[index(j,i-1)]  + 
													(1. - complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt )*phi[index(j,i)] 
													-cz*phi[index(j,i+1)];
						
					};
					
					bz[Nz-1+id_thread*Nz]    =   1. + complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4. )*dt;			
					
					rv_z[Nz-1+id_thread*Nz]	=	-az*phi[index(j,Nz-2)]  
												+(1. - complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4.)*dt )*phi[index(j,Nz-1)];				
					
					//Solving Triagonal Matrix	for Z
					trid_simple(az, &bz[id_thread*Nz], cz, &rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz);
					
				};
			
			
			
			dt=dtev;
			//Evolution on RHO
#pragma omp for private(i)			
			//  RHO_Operator dt	
			for (i=0; i<Nz; i++ )
			{
				
				//tid = omp_get_thread_num();
				
				//printf("Thread number %d is doing loop number %d.\n",tid,i);
				
				rv_r[id_thread*Nr]			=	-ar[0]*phi[index(1,i)] +
												+(1.  -  complex(0., 1./dr/dr/2. + pot[index(0,i)]/4. )*dt)*phi[index(0,i)]
												-cr[0]*phi[index(1,i)];			
				
				br[id_thread*Nr]			=  1. + complex(0., 1./2./dr/dr + pot[index(0,i)]/4. )*dt;			
				
				
				for (j=1; j<Nr-1; j++)
				{ 
					br[j+id_thread*Nr]		= 1. + complex(0., 1./2./dr/dr + pot[index(j,i)]/4. )*dt;	
					
					rv_r[j+id_thread*Nr]	=   -ar[j]*phi[index(j-1,i)]  
												+(1. -  complex( 0., 1./dr/dr/2. + pot[index(j,i)]/4. )*dt)*phi[index(j,i)]
												-cr[j]*phi[index(j+1,i)];
				};
				
				br[Nr-1+id_thread*Nr]   =  1. + complex(0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt;			
				
				rv_r[Nr-1+id_thread*Nr] =  -ar[Nr-1]*phi[index(Nr-2,i)] 
											+(1.  -  complex( 0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt)*phi[index(Nr-1,i)];
				
				
				//Solving Triagonal Matrix
				tridag(&ar[0], &br[id_thread*Nr], &cr[0], &rv_r[id_thread*Nr], &phi_r[id_thread*Nr], &gamr[id_thread*Nr], Nr);
				
				
				for (j=0; j<Nr; j++)
					phi[index(j,i)] = phi_r[j+id_thread*Nr];
				
				
			};//End loop on Nz and Nr			
			
			
			
			dt=dtev/2.;
			//Evolution on Z
#pragma omp for private(j)		
			for ( j=0; j<Nr; j++ )
			{
				
				bz[id_thread*Nz]     =   1. + complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare)/4. )*dt;	
				
				rv_z[id_thread*Nz]	 =  (1. - complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare )/4. )*dt)*phi[index(j,0)] 
										- cz*phi[index(j,1)];	
				
				
				for ( i=1; i<Nz-1; i++) 
				{
					
					bz[i+id_thread*Nz]		=   1. + complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt;				
					
					rv_z[i+id_thread*Nz]	=	-az*phi[index(j,i-1)]   
												+ (1. - complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt )*phi[index(j,i)] 
												-cz*phi[index(j,i+1)];
					
				};
				
				bz[Nz-1+id_thread*Nz]     =   1. + complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4. )*dt;			
				
				rv_z[Nz-1+id_thread*Nz]	  =	-az*phi[index(j,Nz-2)] 
											+(1. - complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4.)*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, &bz[id_thread*Nz], cz, &rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz);
				
			};
			
			
		};//End of pragma
	};	
	
	
	
	
	
	
	void parallel_evolution_laser_RE( complex dtev, double efield_z )
	{
		complex dt=dtev/2.;
		
		//az = complex(0.,-1./dz/dz/4.)*dt;
		//cz = complex(0.,-1./dz/dz/4.)*dt;
			
		
#pragma omp parallel
		{		
			int j, i;
						
			int id_thread= omp_get_thread_num();
			
			
			//Evolution on Z				
#pragma omp for private(j)
			for (j=0; j<Nr; j++ )
			{
				
				bz[id_thread*Nz]     =    1. + complex( 0., 1./2./dz/dz + (pot[index(j,0)] - 2.*e_charge*z[0]*efield_z )/4. )*dt;				
				rv_z[id_thread*Nz]   =   (1. - complex( 0., 1./2./dz/dz + (pot[index(j,0)] - 2.*e_charge*z[0]*efield_z )/4. )*dt )*phi[index(j,0)] 
										- cz*phi[index(j,1)];
				
				
				for (i=1; i<Nz-1; i++)
				{
					bz[i+id_thread*Nz]	 = 1. + complex( 0., 1./2./dz/dz + (pot[index(j,i)] - 2.*e_charge*z[i]*efield_z )/4. )*dt;				
					
					rv_z[i+id_thread*Nz] =	- az*phi[index(j,i-1)]   
											+ (1. - complex( 0., 1./2./dz/dz + (pot[index(j,i)] - 2.*e_charge*z[i]*efield_z )/4. )*dt )*phi[index(j,i)] 
											- cz*phi[index(j,i+1)];
				};
				
				
				
				bz[Nz-1+id_thread*Nz]    =   1. + complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - 2.*e_charge*z[Nz-1]*efield_z )/4. )*dt;				
				rv_z[Nz-1+id_thread*Nz]	 =  - az*phi[index(j,Nz-2)]  
											+ (1. - complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - 2.*e_charge*z[Nz-1]*efield_z )/4. )*dt )*phi[index(j,Nz-1)];
				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, &bz[id_thread*Nz], cz, &rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz);
				
				
			};//end loop Nr and Nz
			
			
			
			dt=dtev;
			//Evolution on RHO
#pragma omp for private(i)	
			for (i=0; i<Nz; i++ )
			{
				rv_r[id_thread*Nr]			=	-ar[0]*phi[index(1,i)] 
												+ (1.  -  complex(0., 1./dr/dr/2. + (pot[index(0,i)] )/4. )*dt)*phi[index(0,i)]
												-cr[0]*phi[index(1,i)];	
				
				br[id_thread*Nr]			=  1. + complex(0., 1./2./dr/dr + (pot[index(0,i)] )/4. )*dt;			
				
				
				for (j=1; j<Nr-1; j++)
				{ 
					br[j+id_thread*Nr]		=    1. + complex(0., 1./2./dr/dr + (pot[index(j,i)] )/4. )*dt;	
					
					
					
					rv_r[j+id_thread*Nr]	=   -ar[j]*phi[index(j-1,i)] 
												+(1. -  complex( 0., 1./dr/dr/2. + (pot[index(j,i)] )/4. )*dt)*phi[index(j,i)]
												-cr[j]*phi[index(j+1,i)];
				};
				
				br[Nr-1+id_thread*Nr]		=  1. + complex(0., 1./2./dr/dr + (pot[index(Nr-1,i)] )/4. )*dt;			
				
				rv_r[Nr-1+id_thread*Nr]		=  -ar[Nr-1]*phi[index(Nr-2,i)]  
												+(1.  -  complex( 0., 1./2./dr/dr + (pot[index(Nr-1,i)] )/4. )*dt)*phi[index(Nr-1,i)];
				
				
				//Solving Triagonal Matrix
				tridag(&ar[0], &br[id_thread*Nr], &cr[0], &rv_r[id_thread*Nr], &phi_r[id_thread*Nr], &gamr[id_thread*Nr], Nr);
				
				
				for (j=0; j<Nr; j++)
					phi[index(j,i)] = phi_r[j+id_thread*Nr];
			};//End loop on Nz and Nr
			
			
			
			
			dt=dtev/2.;
			//Evolution on Z
#pragma omp for private(j)		
			for (j=0; j<Nr; j++ )
			{
				
				bz[id_thread*Nz]     =    1. + complex( 0., 1./2./dz/dz + ( pot[index(j,0)] - 2.*e_charge*z[0]*efield_z )/4. )*dt;			
				
				rv_z[id_thread*Nz]   =   (1. - complex( 0., 1./2./dz/dz + ( pot[index(j,0)] - 2.*e_charge*z[0]*efield_z )/4.  )*dt )*phi[index(j,0)] 
				- cz*phi[index(j,1)];
				
				
				for (i=1; i<Nz-1; i++)
				{
					bz[i+id_thread*Nz]	 =  1. + complex( 0., 1./2./dz/dz +  (pot[index(j,i)] - 2.*e_charge*z[i]*efield_z )/4. )*dt;				
					
					rv_z[i+id_thread*Nz] =	- az*phi[index(j,i-1)]  + 
											+ (1. - complex( 0., 1./2./dz/dz + (pot[index(j,i)] - 2.*e_charge*z[i]*efield_z )/4. )*dt )*phi[index(j,i)] 
											- cz*phi[index(j,i+1)];
				};
				
				
				bz[Nz-1+id_thread*Nz]    =   1. + complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - 2.*e_charge*z[Nz-1]*efield_z )/4. )*dt;			
				
				rv_z[Nz-1+id_thread*Nz]	 =  - az*phi[index(j,Nz-2)]  
											+ (1. - complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - 2.*e_charge*z[Nz-1]*efield_z )/4. )*dt )*phi[index(j,Nz-1)];
				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, &bz[id_thread*Nz], cz, &rv_z[id_thread*Nz], &phi[index(j,0)], &gamz[id_thread*Nz], Nz);
				
				
			};//end loop Nr and Nz
			
			
		};//End of pragma
	};		
	

	
	
		
	//************************************************//
	//		       Hydrogen  potential                //
	//************************************************//
	void set_potential_hlike2D(double _charge_nuclei, double _soft_core )
	{		
		double soft_core=_soft_core;		
		
		for(int j=0;j<Nr;j++)
		  for(int i=0;i<Nz;i++)
		       pot[index(j,i)]= e_charge*_charge_nuclei/sqrt(soft_core + r[j]*r[j]+z[i]*z[i]);
		
	}//End Hydrogen potential

	
	
	//************************************************//
	//		       Poschl Teller Potential             //
	//************************************************//	
	void set_poschl_teller2D(double V0, double alpha )
	{		
		double rprime   =0. ;
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				rprime			= sqrt(r[j]*r[j] + z[i]*z[i]); 
				pot[index(j,i)] = -V0*(V0+1.)/2.*alpha*alpha*pow(cosh(alpha*rprime) , -2. );
			}
	}//End Poschl Teller Potential	
	
	
	//************************************************//
	//		       Tong Potential well                //
	//************************************************//
	void set_tong_potential_2D(double _charge_nuclei, double _soft_core, double *alpha )
	{		
		double soft_core=_soft_core;		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				pot[index(j,i)]= e_charge*( _charge_nuclei 
										    + alpha[0]*exp( - alpha[1]*sqrt( r[j]*r[j] + z[i]*z[i] ) )
										    + alpha[2]*sqrt( r[j]*r[j] + z[i]*z[i] )*exp( - alpha[3]*sqrt( r[j]*r[j] + z[i]*z[i] ) )
										    + alpha[4]*exp( - alpha[5]*sqrt( r[j]*r[j] + z[i]*z[i] ) )
										   )/sqrt( soft_core + r[j]*r[j] + z[i]*z[i] );
		
	}//End potential	
	
	

	//************************************************//
	//		     Argon  Muller Potential well         //
	//************************************************//
	void set_Armuller_potential_2D(double _charge_nuclei, double _soft_core, double *alpha )
	{		
		double soft_core=_soft_core;		
		double F  = 2.5;
		double G  = 2.01785;
		double Rx = 3.0;
		
		double rr=0.;
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				rr=sqrt( r[j]*r[j] + z[i]*z[i] );
				pot[index(j,i)]= e_charge*( _charge_nuclei 
										   + alpha[0]*exp( - alpha[1]*rr )
										   + (17.-alpha[0])*exp( - alpha[2]*rr )
										   )/sqrt( soft_core + r[j]*r[j] + z[i]*z[i] );
				
				/*if(rr<=Rx)
				{
					pot[index(j,i)]+= 0*e_charge*( pow((Rx-rr)/G,5)-pow((Rx-rr)/G,4) )*F;
				}*/
			};
	}//End potential		
	
	
	
	//************************************************//
	//		       Hydrogen molecule potential H2+     //
	//************************************************//
	void set_potential_hplus2D( double zmol1, double zmol2, double soft_core1, double soft_core2, double R, double shift )
	{		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				pot[index(j,i)]= e_charge*zmol1/sqrt( soft_core1 + r[j]*r[j] + (z[i]+shift-R)*(z[i]+shift-R) ) +
								 e_charge*zmol2/sqrt( soft_core2 + r[j]*r[j] + (z[i]+shift)*(z[i]+shift) );
		
	}//End Hydrogen potential	
	

	
	//*******************************************************************//
	//    CO or CO2 Molecular Hydrogen Model  PRL 111, 073902 (2013) //
	//************************************************//
	void set_potential_CO( double *Zmol_i, double *Zmol_o, double *R, double *soft_core, double *sigma )
	{		
		double sqr0	 =0;
		double sqr1	 =0;
		double temp0 =0;
		double temp1 =0;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++){
				sqr0	= r[j]*r[j] + (z[i]-R[0])*(z[i]-R[0]);
				sqr1	= r[j]*r[j] + (z[i]-R[1])*(z[i]-R[1]);
				
				temp0	=(Zmol_i[0] - Zmol_o[0])*exp(-sqr0/sigma[0]) + Zmol_o[0];
				temp1	=(Zmol_i[1] - Zmol_o[1])*exp(-sqr1/sigma[1]) + Zmol_o[1];
				
				pot[index(j,i)]= e_charge*( temp0/sqrt(soft_core[0]+sqr0)+
										    temp1/sqrt(soft_core[1]+sqr1));
			}
	}//End Hydrogen potential		
	
	
	
	
	
	//"Static Dipole or expectation position value"	
	double static_dipole2D(){
		
		double dip=0;
		
		for (int j=0; j<Nr; j++) 
			for(int i=0;i<Nz;i++)
				dip+= abs( phi[index(j,i)] )*abs( phi[index(j,i)] )*z[i]*r[j]*dz*dr;
		
		return dip;
		
	}
	
	
	
	//****************************************************************************
	//          Expected Kinetic energy value finite difference method
	//****************************************************************************
	double kinetic_energy_finite()
	{
		
		double me=1.;
		double a1=1./2./me;
		double dos=2.;
		
		complex kinetic=complex(0.,0.);
		
				
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				complex psi = phi[index(j,i)];
				
				complex psi_ip = complex(0.,0.);
				complex psi_im = complex(0.,0.);
				
				complex psi_jp = complex(0.,0.);
				complex psi_jm = complex(0.,0.);
				
				
				/******************************************************************/
				if (i-1>=0)
					psi_ip=phi[index(j,i-1)];
				
				if (i+1<Nz)
					psi_im=phi[index(j,i+1)];
				
				/******************************************************************/
				if (j-1>=0)
					psi_jp=phi[index(j-1,i)];
				
				if (j+1<Nr)
					psi_jm=phi[index(j+1,i)];
				
				/******************************************************************/
				
				
				
				complex ksquare=				
				-a1*(  psi_ip - dos*psi + psi_im )/dz/dz 
				-a1*( (psi_jp - dos*psi + psi_jm )/dr/dr + 
					 ( psi_jp - psi_jm )/r[j]/dr/dos );
				
				
				kinetic+= dz*dr*r[j]*(conj(psi)*ksquare);
				
			}
		return real(kinetic);	
	}//End kinetic energy expected value	
	
	
		
		
	//**********************************************************
	//	    Expected Potential energy value
	//**********************************************************
	double potential_energy()
	{		
		double potE=0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				potE+= dz*dr*r[j]*pot[index(j,i)]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return potE;
		
	} 
	// End potential energy expected value	
	
	
	
	//**********************************************************
	//	    z-position expectation value
	//**********************************************************
	double z_position_value(double z0, double zf, double rf)
	{		
		double z_exp  = 0.;
		double z_norm = 0.;
		double edensity=0.;
		
		int i_z0      = floor(abs(z0-z[0])/dz);
		int i_zf      = floor(abs(zf-z[0])/dz);
		int j_rf      = floor(abs(rf-r[0])/dr);		
		
		if(i_zf==0)
			i_zf=1;

		if(i_zf>=Nz)
			i_zf=Nz-1;		
		
		if(j_rf==0)
			j_rf=1;
		
		if(j_rf>=Nr)
			j_rf=Nr-1;
		

		
		for(int j=0;j<j_rf;j++)
			for(int i=i_z0;i<i_zf;i++)
			{
				edensity = pow(abs(phi[index(j,i)]), 2.);
				z_exp+= dz*dr*r[j]*z[i]*edensity;
				z_norm+=dz*dr*r[j]*edensity; 
				
			}
		
		return z_exp/z_norm;
		
	} 
	// End potential energy expected value		
	
	
	
	
	//*****************************************
	//     Gaussian initial function
	//*****************************************
	void rho_gaussian( double r0, double z0, double rsigma, double zsigma )
		{
		for(int j=0; j<Nr; j++)
			for(int i=0; i<Nz; i++)
				phi[index(j,i)] = exp( -( (r[j]-r0)*(r[j]-r0)/rsigma/rsigma + (z[i]-z0)*(z[i]-z0)/zsigma/zsigma ) );
		}//End Gaussian
	 
	 
	 
	 
	 
	//****************************************                                                                            
	//    Gaussian with initial velocity                                                                                          
	//****************************************                                                                                             
	void velocity_gaussian( double r0, double z0, double vr0, double vz0, double rsigma, double zsigma )
	{	    
		for(int j=0; j<Nr; j++)
			for(int i=0; i<Nz; i++)
				phi[index(j,i)] = exp( -((r[j]-r0)*(r[j]-r0)/rsigma/rsigma + (z[i]-z0)*(z[i]-z0)/zsigma/zsigma ) )*exp(I*(vr0*(r[j]-r0) +vz0*(z[i]-z0) ) );
			
	}//End Gaussian Velocity 
        
	
	//***************************************// 
	//               Save axes				 //
	//***************************************//	
	
	void saveAxes(fstream &file,int skiper1=1,int skiper2=1)
	{
		for(int j=0;j<Nr/skiper1;j++)
			file << r[j*skiper1] << endl;
		
		
		for(int i=0;i<Nz/skiper2;i++)
			file << z[i*skiper2] << endl;
		
			file.flush();
	}//End save axes
	
	
	//***************************************
	//           Mask function 2D
	//***************************************
	void mask_function2D(waveUniform2D &wp,const double &x0,const double &x1,const double &sigma,double rho0=0.0,double z0=0.0)
	{
		double x;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				x = sqrt( (r[j]-rho0) * (r[j]-rho0) + (z[i]-z0) * (z[i]-z0) );
				
				if(x<=x0)
					wp.phi[index(j,i)]=0.;
				
				if(x>x0 && x<=x1)
					wp.phi[index(j,i)]=phi[index(j,i)]*(1.-exp(-(x-x0)*(x-x0)/sigma/sigma));
				
				if(x>x1)
					wp.phi[index(j,i)]=phi[index(j,i)];
				
			};
		
	};//End mask function 2D
	
	
	
	//***************************************
	//        Inverse Mask function 2D
	//***************************************
	void inverse_mask_function2D(waveUniform2D &wp,const double &x0,const double &x1,const double &sigma,double rho0=0.0,double z0=0.0)
	{
		double x;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				x = sqrt( (r[j]-rho0) * (r[j]-rho0) + (z[i]-z0) * (z[i]-z0) );
				
				if(x<=x0)
					wp.phi[index(j,i)]=phi[index(j,i)];
				
				if(x>x0 && x<=x1)
					wp.phi[index(j,i)]=phi[index(j,i)]*(1.-exp(-(x-x1)*(x-x1)/sigma/sigma));
				
				if(x>x1)
					wp.phi[index(j,i)]=complex(.0,.0);
				
			};
		
	};//End mask function 2D		
	
	
	
	
	//******************************************// 
	//                absorber
	//******************************************// 	
	void absorber(const double &frac_r_left,const double &frac_z_left,const double &frac_r_right,const double &frac_z_right,const double &exponent)
	{
		
		// Try 1./6. in the exponent.
		
		double mask_start_z_right=z[int(Nz*(1.-frac_z_right))];
		double mask_start_r_right=r[int(Nr*(1.-frac_r_right))];
		
		double mask_start_z_left=z[int(Nz*frac_z_left)+1];
		double mask_start_r_left=r[int(Nr*frac_r_left)+1];
		
		
		double argument_z_right;
		double argument_r_right;
		
		
		double argument_z_left;
		double argument_r_left;	
		
		
		double mask_z_right;
		double mask_r_right;
		
		double mask_z_left;
		double mask_r_left;
		
		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				argument_z_right=(pi/2.)*(z[i]-mask_start_z_right)/(z[Nz-1]-mask_start_z_right+1.e-20);
				argument_r_right=(pi/2.)*(r[j]-mask_start_r_right)/(r[Nr-1]-mask_start_r_right+1.e-20);
				
				
				argument_z_left=(pi/2.)*(z[i]-mask_start_z_left)/(z[0]-mask_start_z_left+1.e-20);
				argument_r_left=(pi/2.)*(r[j]-mask_start_r_left)/(r[0]-mask_start_r_left+1.e-20);
				
				
				mask_z_right=pow(fabs(cos(argument_z_right)),exponent);
				mask_r_right=pow(fabs(cos(argument_r_right)),exponent);
				
				
				mask_z_left=pow(fabs(cos(argument_z_left)),exponent);
				mask_r_left=pow(fabs(cos(argument_r_left)),exponent);
				
				
				
				if (i< int(Nz*frac_z_left))
				{		  
					real(phi[index(j,i)])*=mask_z_left;	
					imag(phi[index(j,i)])*=mask_z_left;	
				}
				
				if (i> int(Nz*(1.-frac_z_right)))
				{		  
					real(phi[index(j,i)])*=mask_z_right;	
					imag(phi[index(j,i)])*=mask_z_right;	
				}
				
				if (j< int(Nr*frac_r_left))
				{		  
					real(phi[index(j,i)])*=mask_r_left;	
					imag(phi[index(j,i)])*=mask_r_left;	
				}
				
				if (j> int(Nr*(1.-frac_r_right)))
				{		  
					real(phi[index(j,i)])*=mask_r_right;	
					imag(phi[index(j,i)])*=mask_r_right;	
				};
				
			};//End the loop on ij
	};//End absorver
	
	
	
	//***********************************************
	//                  Binary writer
	//***********************************************
	void binwrite(FILE *afile )
	{		
		double wreal;
		double wimag;
		
		for(int j=0; j<Nr; j++)
			for(int i=0; i<Nz; i++)
			{
				wreal = double(real(phi[index(j,i)]));
				wimag = double(imag(phi[index(j,i)]));
				
				fwrite(&wreal , 1 , sizeof(wreal) , afile );
				fwrite(&wimag , 1 , sizeof(wimag) , afile );
				
				//fwrite(&wreal , sizeof(wreal) , 1 , afile );
				//fwrite(&wimag , sizeof(wimag) , 1 , afile );
			
			};
		
		fflush(afile);
	};//End Writer
	
	
	
	
	
	//***********************************************
	//				  Binary reader
	//***********************************************
	void binread(FILE *afile)
	{	
		double wreal;
		double wimag;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fread(&wreal , 1 , sizeof(wreal) , afile );
				fread(&wimag , 1 , sizeof(wimag) , afile );
				
				//fread(&wreal , sizeof(wreal) , 1 , afile );
				//fread(&wimag , sizeof(wimag) , 1 , afile );
				
				real(phi[index(j,i)])	= wreal;
				imag(phi[index(j,i)])	= wimag;		
			};	
		fflush(afile);
	};//End reader	
	

	//***********************************************
	//				  Snapshot
	//***********************************************	
	void snapshot(fstream &file,int skiper1=1,int skiper2=1)
	{
		
		double norm=0.0;
		
		for(int j=0;j<Nr/skiper2;j++)
			for(int i=0;i<Nz/skiper1;i++)
			{
				norm=dz*dr*r[j*skiper2]*pow( abs(phi[index(j*skiper2,i*skiper1)] ), 2.);
				file << norm << endl;
			};
		
		file.flush();
	};
	

	
	void silo_axes(double mintheta, double maxtheta, int Ntheta, int skiper1, int skiper2)
	{
		
		if(mintheta < 0.)
		{
			cout << "\n***********************\n Take care \n";			
			cout << "\nmintheta must be > 0 \n";
			exit(1);
		};
		
		if(maxtheta > dospi)
		{
			cout << "\n***********************\n Take care \n";
			cout << "\n\nmaxtheta must be < 2*pi \n\n";
			exit(1);
		};		
		
		int NX			= Nz/skiper1;
		int NY			= Nr/skiper2;
		int NZ			= Ntheta;
		
		double theta	= 0.;		
		double dtheta	= (maxtheta-mintheta)/((double)(Ntheta-1));
		
		
		xcoord			= ( double* ) mkl_malloc( NX*NY*NZ*sizeof(double),16 );
		ycoord			= ( double* ) mkl_malloc( NX*NY*NZ*sizeof(double),16 );
		zcoord			= ( double* ) mkl_malloc( NX*NY*NZ*sizeof(double),16 );		
		edensity		= ( double* ) mkl_malloc( NX*NY*NZ*sizeof(double),16 );	
		
				
		
		for(int k=0;k<NZ;k++)
		{
			theta = mintheta + k*dtheta;
			for (int j=0; j<NY; j++) 
			{
				for (int i=0; i<NX; i++) 
				{

					xcoord[index3D(k,j,i,NY,NX)] = r[j*skiper2]*cos(theta);
					ycoord[index3D(k,j,i,NY,NX)] = r[j*skiper2]*sin(theta);
					zcoord[index3D(k,j,i,NY,NX)] = z[i*skiper1];
					
				};
			};
		};
		/*
		
		coords[0] = xcoord;
		coords[1] = ycoord;
		coords[2] = ycoord;		*/
		

	};
	
	
	//***********************************************
	//			  Snapshot On Silo File 
	//***********************************************	
	void snapshot_silo_file(double mintheta, double maxtheta, int Ntheta, int ksnapshot, double t, int skiper1=1,int skiper2=1)
	{

		DBfile *dbfile = NULL; // Open the silo file
		DBoptlist *optList;		
		
		double *coords[3];
		
		coords[0] = xcoord;
		coords[1] = ycoord;
		coords[2] = zcoord;					//*/
		
		
		int NX			= Nz/skiper1;
		int NY			= Nr/skiper2;
		int NZ			= Ntheta;
		
		char filename[80];		
		
		double eDensity = 0.0;
		double dtheta   = (maxtheta-mintheta)/( (double)(Ntheta-1) );
		
		for(int k=0; k<NZ; k++)
			for(int j=0;j<NY;j++)
				for(int i=0;i<NX;i++)
				{
					eDensity  = log10(dtheta*dz*dr*r[j*skiper2]*pow( abs( phi[index(j*skiper2,i*skiper1)] ), 2. )+1.e-14);
					edensity[index3D(k,j,i,NY,NX)]  = eDensity;
				};
		
		
		
		//Creating a file that contain the mesh and variable//
		sprintf(filename, "edensity%.4d.silo", ksnapshot );
		printf("\n....\nCreating Electron Density File \"%s\".\n\n", filename  );	
		dbfile = DBCreate(filename, 
						  DB_CLOBBER, 
						  DB_LOCAL, 
						  "Laser-Atom Interaction ", 
						  DB_HDF5);
		
		
		int ndims			= 3;
		int dims[3];
		char *coordnames[3] ={"x","y","z"};
		
		
		dims[0]				= NX ;
		dims[1]				= NY ;
		dims[2]				= NZ ;
		
		
		optList				= DBMakeOptlist(10);
		DBAddOption(optList, DBOPT_DTIME, &t);
		
		
		
		DBPutQuadmesh(dbfile, "quadmesh", coordnames, 
					  coords, dims, ndims,
					  DB_DOUBLE, DB_NONCOLLINEAR, optList );	
		
		
		DBPutQuadvar1(dbfile, "EDensity", "quadmesh", 
					  edensity, dims, ndims, NULL,
					  0, DB_DOUBLE, DB_NODECENT, optList );

		DBPutQuadvar1(dbfile, "EDensity2", "quadmesh", 
					  edensity, dims, ndims, NULL,
					  0, DB_DOUBLE, DB_NODECENT, optList );		
		
		DBFreeOptlist(optList);	
		DBClose(dbfile);
	};//*/
	
	
	
	//***********************************************
	//			  Snapshot On Silo File 
	//***********************************************	
	void save_potential_silo_file(double mintheta, double maxtheta, int Ntheta, int ksnapshot, double t, int skiper1=1,int skiper2=1)
	{
		
		
		DBfile *dbfile = NULL; // Open the silo file
		DBoptlist *optList;		
		
		double *coords[3];
		
		coords[0] = xcoord;
		coords[1] = ycoord;
		coords[2] = zcoord;	//*/
		
		
		int NX			= Nz/skiper1;
		int NY			= Nr/skiper2;
		int NZ			= Ntheta;
		
		char filename[80];		
		
		double eDensity = 0.0;
		double dtheta   = (maxtheta-mintheta)/( (double)(Ntheta-1) );
		
		for(int k=0; k<NZ; k++)
			for(int j=0;j<NY;j++)
				for(int i=0;i<NX;i++)
				{
					eDensity  = pot[index(j*skiper2,i*skiper1)];
					edensity[index3D(k,j,i,NY,NX)]  = eDensity;
				};
		
		
		
		//Creating a file that contain the mesh and variable//
		sprintf(filename, "Potential%.4d.silo", ksnapshot );
		printf("\n....\nCreating Potential file \"%s\".\n\n", filename  );	
		dbfile = DBCreate(filename, 
						  DB_CLOBBER, 
						  DB_LOCAL, 
						  "Potential System ", 
						  DB_HDF5);
		
		
		int ndims			= 3;
		int dims[3];
		char *coordnames[3] ={"x","y","z"};
		
		
		dims[0]				= NX ;
		dims[1]				= NY ;
		dims[2]				= NZ ;
		
		
		optList				= DBMakeOptlist(10);
		DBAddOption(optList, DBOPT_DTIME, &t);
		
		
		
		DBPutQuadmesh(dbfile, "quadmesh", coordnames, 
					  coords, dims, ndims,
					  DB_DOUBLE, DB_NONCOLLINEAR, optList );	
		
		
		DBPutQuadvar1(dbfile, "EDensity", "quadmesh", 
					  edensity, dims, ndims, NULL,
					  0, DB_DOUBLE, DB_NODECENT, optList );
		
		DBFreeOptlist(optList);	
		DBClose(dbfile);
	};	//*/
	
	
	void snapshot_square_wave(fstream &file,int skiper1=1,int skiper2=1)
	{
		double norm=0.0;
		
		for(int j=0;j<Nr/skiper2;j++)
			for(int i=0;i<Nz/skiper1;i++)
			{
				norm=pow(abs(phi[index(j*skiper2,i*skiper1)]),2.);
				file << norm << endl;
			};
		file.flush();		
	};
	
	
	
	void fplaceWF( complex *wsmall, int Nr_Small, int Nz_Small)
	{	
		int iplacer = floor((Nz - Nz_Small)/2);
		memset(phi, 0, sizeof(complex)*Nr*Nz );
		
		
		for(int j=0; j<Nr_Small; j++)
			for(int i=0; i<Nz_Small; i++)
				phi[index(j,iplacer+i)] = wsmall[j*Nz_Small+i];
		
	};//End replacer of wavefunction
	
	
	
	void resize_waveUniform2D(HankelMatrix &HH_resize, int _Nz_rezise, double _dz )
	{
		
		int NrSmall = Nr;
		int NzSmall = Nz;
		
		complex *temp_phi;
		temp_phi = (complex*)mkl_malloc( NrSmall*NzSmall*sizeof(complex), 16);
		
		
		for (int i=0; i<NrSmall*NzSmall;i++)
			temp_phi[i]=phi[i];
		
		
		//Free memory
		free_memory();
		
		
		//Initializing new arrays
		initialize(HH_resize, _Nz_rezise, _dz);
		
		
		//Replacing data on larger wavefunction 
		fplaceWF(temp_phi, NrSmall, NzSmall);
		
		
		mkl_free(temp_phi);
		
	};//End resize of the function //
	
	
	
	void free_memory()
	{
		if(r!=NULL)
		{
			mkl_free(r);
			mkl_free(z);
		
			mkl_free(phi);	
			mkl_free(pot);
		};
		
		if (ar!=NULL)	mkl_free(ar);
		
		if (br!=NULL)	mkl_free(br);
		if (cr!=NULL)	mkl_free(cr);
		
		if (rv_r!=NULL)	mkl_free(rv_r);	
		if (phi_r!=NULL)mkl_free(phi_r);
		if (gamr!=NULL)	mkl_free(gamr);			
		
		
		if (rv_z!=NULL)	mkl_free(rv_z);		
		if (bz!=NULL)	mkl_free(bz);
		if (bz!=NULL)	mkl_free(gamz);	
		
		
		if(xcoord!=NULL)
		{
			mkl_free(xcoord);
			mkl_free(ycoord);
			mkl_free(zcoord);
			mkl_free(edensity);
		};//*/
	};
	
	
};
#endif //
