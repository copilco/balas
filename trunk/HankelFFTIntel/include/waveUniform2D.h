/*
 *  wave.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 29/11/11.
 *  Copyright 2011 USAL. All rights reserved.
 *
 */

#ifndef WAVEUNIFORM_H
#define WAVEUNIFORM_H

#include "tools.h"
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

class waveUniform2D
{
	
public:

	int Nr,Nz;
	double *r,dr;
	double *z,dz;
	double R;
	
	complex *phi;
	double *v;
	
	complex *rv_r,*rv_z;
	complex *phi_r,*phi_z;
	complex *ar, *br, *cr,*gamr;
	complex az,*bz,cz,*gamz;
	
	
	/***********************************************************/	
	//Fill the arrays and copy the grids in rho and k_rho
	/***********************************************************/	
	inline int index(int j,int i)
	{
		int indexer=j*Nz+i;
		return indexer;
	}
	
	void initialize(HankelMatrix &HH,int _Nz,double _dz)
	{
		Nr=HH.Nr;
		R=HH.R;
		dr=R/Nr;
		
		Nz=_Nz;
		dz=_dz;
		
		/*
		phi=new complex[Nr*Nz];		
		v=new double[Nr*Nz];		

		r=new double[Nr];
		z=new double[Nz];
		
		rv_r=    new complex[Nr];
		phi_r=   new complex[Nr];
		ar=      new complex[Nr];
		br=      new complex[Nr];
		cr=      new complex[Nr];
		gamr=    new complex[Nr];
		
		rv_z=    new complex[Nz];
		phi_z=   new complex[Nz];
		bz=      new complex[Nz];
		gamz=    new complex[Nz];	
		*/
		
		phi=(complex*)mkl_malloc(Nz*Nr*sizeof(complex),16);
		v=(double*)mkl_malloc(Nz*Nr*sizeof(double),16);

		r=(double*)mkl_malloc(Nr*sizeof(double),16);
		z=(double*)mkl_malloc(Nr*sizeof(double),16);

		
		rv_r=(complex*)mkl_malloc(Nr*sizeof(complex),16);
		phi_r=(complex*)mkl_malloc(Nr*sizeof(complex),16);
		ar=(complex*)mkl_malloc(Nr*sizeof(complex),16);
		br=(complex*)mkl_malloc(Nr*sizeof(complex),16);
		cr=(complex*)mkl_malloc(Nr*sizeof(complex),16);
		gamr=(complex*)mkl_malloc(Nr*sizeof(complex),16);
				
		rv_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
		phi_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
		bz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
		gamz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	
		
		/***********************************************************/	
		//Fill the arrays with zeros.
		/***********************************************************/
	
		for (int i=0; i<Nr*Nz; i++)
		{
			phi[i]=complex(0.,0.);
			v[i]=0.;
		}
		
		
		for(int i=0;i<Nr;i++)
		{
			ar[i]   =complex(0.,0.);
			cr[i]   =complex(0.,0.);			
			br[i]   =complex(0.,0.);	
			rv_r[i] =complex(0.,0.);
			phi_r[i]=complex(0.,0.);
			gamr[i]=complex(0.,0.);
		}
		
		
		for(int i=0;i<Nz;i++)
		{
			rv_z[i] =complex(0.,0.);
			phi_z[i]=complex(0.,0.);
			bz[i]=complex(0.,0.);
			gamz[i]=complex(0.,0.);
		}
		
		/***********************************************************/	
		//Fill the grids onto the wavefunction
		/***********************************************************/
		
		r[0]=dr/2.;
		for (int i=1; i<Nr; i++)
			r[i] =r[i-1]+dr ;
		
		for (int i=0; i<Nz; i++)
			z[i] =(-(Nz-1)/2. + i)*dz;
		
	}
	
	/***********************************************************/	
	//Norms in both spaces
	/***********************************************************/	
	
	
	double norm()
	{
		double norm=0.0;
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nr;i++)
				norm+=dz*dr*r[j]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return norm;
	}
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr*Nz;i++)
			phi[i]=phi[i]/sqrt(norm1);
	}
	
	
	
	/***********************************************************/	
	//Preparing the arrays fof the Crank
	/***********************************************************/
	
	void PrepareCrankArrays(complex dt)
	{
		for( int i=0; i<Nr; i++ )
		{						
			ar[i]     =   (complex( 0., - 1./4./dr/dr  ) +
						   complex( 0., + 1./8./dr/r[i]  ))*dt;
			
			//a2		  =   complex(0., 1./2./dr/dr  + v[index(i,j,nz)]/4.)*dt ;
			//br[i]	  =   a1 + a2;
			
			cr[i]     =   (complex( 0., - 1./4./dr/dr  ) +
						   complex( 0., - 1./8./dr/r[i]  ))*dt;
		}//End left part
		
		
		az = complex( 0., -1./dz/dz/4. )*dt;
		cz = complex( 0., -1./dz/dz/4. )*dt;
		
		/*
		 for (int i=0; i< Nr; i++) 
		 {
		 a[i]     =    complex( 0., - dt/4./dr/dr  ) +
		 complex( 0., + dt/8./dr/r[i]  );
		 
		 b[i]	 =	  complex( 1.0, 0.)+
		 complex( 0.0, dt/2./dr/dr     );
		 //	complex( 0.0, dt/2./dr/dr  + v[i]*dt/2. );
		 
		 c[i]     =    complex( 0., - dt/4./dr/dr  ) +
		 complex( 0., - dt/8./dr/r[i]  );
		 }
		 */
	}
	
	
	void Zprop1( complex dt )
	{
	
		
		complex a2 =complex(0.,0.);		
		
		//az = complex( 0., -1./dz/dz/4. )*dt;
		//cz = complex( 0., -1./dz/dz/4. )*dt;
		
		//*******************************
		//===============================//
		//Z_operator 1/2*dt
		for (int j=0; j<Nr; j++ )
		{
			
			//Left part Z
			for(int i=0; i<Nz; i++)	
			{
				a2= complex(0.,1./2./dz/dz  + v[index(j,i)]/4.) ;
				bz[i]	=	complex(1.,0.) + a2*dt;
			}
			
			
			//Right part Z
			rv_z[0]		=   conj(bz[0])*phi[index(j,0)]   + 
			conj( cz  )*phi[index(j,1)];	
			
			
			for (int i=1; i<Nz-1; i++) 
				rv_z[i] =	conj( az    )*phi[index(j,i-1)] + 
				            conj( bz[i] )*phi[index(j,i)]   + 
				            conj( cz    )*phi[index(j,i+1)];
			
			
			rv_z[Nz-1]	=	conj( az  )*phi[index(j,Nz-2)]  +
		       	            conj( bz[Nz-1]  )*phi[index(j,Nz-1)];	
			//Finishing right part Z
			
			//Zeros on Z
			//zeros(phi_z);		
			//for(int i=0;i<Nz;i++)
			//	phi_z[i]=complex(0.,0.);
			
			
			//Solving Triagonal Matrix	for Z
			trid_simple(az,bz,cz,rv_z,phi_z,gamz,Nz);
			//tridagS(az,bz,cz,rv_z,phi_z, Nz,gamz);
			
			//Save function 
			for (int i=0; i<Nz; i++)
				phi[j*Nz+i] =phi_z[i];		//psi0[index(i,j,nz)] = psi0[index(i,j,nz)];//
		}	
	}
	
		
	
	
	void Zprop( complex dt )
	{
		
		/*
		 for (int i=0; i<Nr*Nz; i++)
		 {
		 v[i]=0.;
		 }*/
		
		
		//complex a1 =complex(1.,0.);
		complex a2 =complex(0.,0.);		
		
		//bz.resize(    nz, 0. );
		//rv_z.resize(  nz, 0. );
		//phi_z.resize( nz, 0. );
		
		//az = complex( 0., -1./dz/dz/4. )*dt;
		//cz = complex( 0., -1./dz/dz/4. )*dt;
		
		//*******************************
		//===============================//
		//Z_operator 1/2*dt
		for (int j=0; j<Nr; j++ )
		{
			
			//Left part Z
			for(int i=0; i<Nz; i++)	
			{
				a2= complex(0.,1./2./dz/dz  + v[j*Nz+i]/4.) ;
				bz[i]	=	complex(1.,0.) + a2*dt;
			}
			
			
			//Right part Z
			rv_z[0]		=   conj(bz[0])*phi[j*Nz+0]   + 
			conj( cz  )*phi[j*Nz+1];	
			
			
			for (int i=1; i<Nz-1; i++) 
				rv_z[i] =	conj( az    )*phi[ j*Nz+(i-1)] + 
				conj( bz[i] )*phi[j*Nz+i]   + 
				conj( cz    )*phi[j*Nz+i+1];
			
			
			rv_z[Nz-1]	=	conj( az  )*phi[j*Nz+Nz-2]  +
			conj( bz[(Nz-1)]  )*phi[j*Nz+(Nz-1)];	
			//Finishing right part Z
			
			//Zeros on Z
			//zeros(phi_z);		
			//for(int i=0;i<Nz;i++)
			//	phi_z[i]=complex(0.,0.);
			
			
			//Solving Triagonal Matrix	for Z
			//trid_simple(az,bz,cz,rv_z,phi_z,gamz,Nz);
			tridagS(az,bz,cz,rv_z,phi_z, Nz,gamz);
			
			//Save function 
			for (int i=0; i<Nz; i++)
				phi[j*Nz+i] =phi_z[i];		//psi0[index(i,j,nz)] = psi0[index(i,j,nz)];//
		}	
	}
	

	//**********************************
	//==================================//
	//RHO_Operator dt	Complete Version 
	
	void Rprop( complex dt )
	{	
		
		complex a1 = complex(1.,0.);
		complex a2 = complex(0.,0.);	
		
		
		//**********************************
		//**********************************
		//==================================//
		//RHO_Operator dt	
		for ( int j=0; j<Nz; j++ )
		{
			
			//Left part Rho
			for( int i=0; i<Nr; i++ )
			{						
				ar[i]     =   (complex( 0., - 1./4./dr/dr  ) +
							   complex( 0., + 1./8./dr/r[i]  ))*dt;
				
				a2		  =   complex(0., 1./2./dr/dr  + v[index(i,j)]/4.)*dt ;
				br[i]	  =   a1 + a2;
				
				cr[i]     =   (complex( 0., - 1./4./dr/dr  ) +
							   complex( 0., - 1./8./dr/r[i]  ))*dt;
			}//End left part
			
			
			
			//Right part Rho 
			//complex phi0 = phi[index(1,j)];//complex(0.,0.);//
			
			rv_r[0]		 =	conj( ar[0]  )*phi[index(1,j)] +
			                conj( br[0]  )*phi[index(0,j)] + 
			                conj( cr[0]  )*phi[index(1,j)];			
			
			for (int i=1; i<Nr-1; i++) 
				rv_r[i]	=	conj( ar[i]  )*phi[index(i-1,j)] + 
				            conj( br[i]  )*phi[index(i,j)]   + 
				            conj( cr[i]  )*phi[index(i+1,j)];
			
			rv_r[Nr-1]	 =	conj( ar[Nr-1]  )*phi[index(Nr-2,j)]  +
			                conj( br[Nr-1]  )*phi[index(Nr-1,j)];
			//Finishing right
			
			
			
			//Zeros
			//zeros(phi_r);
			
			
			//Solving Triagonal Matrix			
			tridagAlexis(ar, br, cr, rv_r, phi_r,gamr, Nr);
			
			
			for (int i=0; i<Nr; i++)
				phi[index(i,j)] =	phi_r[i];//phi[index(i,j,nz)];
		}
	}

	//******************************************// End Rho_propagator
	//==========================================//
		
		
	/*
	void KineticPropCrankUniform(double dt)
	{
		complex zero=complex(0.,0.);	
		//Starting Right
		rv[0] =   conj( a[0] )*phi[1]
		+ conj( b[0] )*phi[0]  
		+ conj( c[0] )*phi[1];	
		
		for (int i=1; i<Nr-1; i++) 	
			rv[i] = conj( a[i] )*phi[i-1]  +
			conj( b[i] )*phi[i]    + 
			conj( c[i] )*phi[i+1];
		
		rv[Nr-1] = conj( a[Nr-1]  )*phi[Nr-2]  +
		conj( b[Nr-1]  )*phi[Nr-1]  +
		conj( c[Nr-1]  )*phi[Nr-1]*zero;
		
		
		//Solving Triagonal Matrix
		tridag(a,b,c,rv,phil,Nr,gam1);
		
		for(int i=0;i<Nr;i++)
			phi[i]=phil[i];
	}
	 */
		
		
};
#endif // TOOLS_H
