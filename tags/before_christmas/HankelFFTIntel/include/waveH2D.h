/*
 *  wave.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 29/11/11.
 *  Copyright 2011 USAL. All rights reserved.
 *
 */

#ifndef WAVE_H
#define WAVE_H

#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "tools.h"

class waveH2D
{
	
public:

	int Nr,Nz;
	int space_indicator;  //Defines in which space you are.
	double *r,*v,*dr,*dv;
	double *z,dz,*q,dq,qmax;

	complex *F,*phiHank;
	complex *F2,*f2;
	double *pot;
	
	complex *rv,*phil;
	complex *rv_z;
	complex *phi_z;
	complex *a, *b, *c,*gam1;
	complex az,*bz,cz,*gamz;
	
	MKL_LONG status;
	
	DFTI_DESCRIPTOR_HANDLE hand;
	double fscale;// = dz/sqrt(dospi);// 1.;//dz/sqrt(dospi);
	double bscale;// = dw/sqrt(dospi);//1./Nz;//dw/sqrt(dospi);
	
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
		
		Nz=_Nz;
		dz=_dz;
		
		
		//space_indicator=0;//We are in position rho space
		
		
		phiHank=new complex[Nr*Nz];		
		F=new complex[Nr*Nz];
		
		F2=new complex[Nr*Nz];
		f2=new complex[Nr*Nz];
		
		pot=new double[Nr*Nz];
				
		r=new double[Nr];
		v=new double[Nr];
		dr=new double[Nr];
		dv=new double[Nr];
		
		z=new double[Nz];
		q=new double[Nz];
		
		
		rv=    new complex[Nr];
		phil=  new complex[Nr];
		a=     new complex[Nr];
		b=     new complex[Nr];
		c=     new complex[Nr];
		gam1=     new complex[Nr];
		
		rv_z=    new complex[Nz];
		phi_z=   new complex[Nz];
		bz=      new complex[Nz];
		gamz=    new complex[Nz];	
		
		
		/***********************************************************/	
		//Fill the arrays with zeros.
		/***********************************************************/
		
		for(int i=0;i<Nr*Nz;i++)
		{
			phiHank[i]=complex(0.,0.);
			F[i]=complex(0.,0.);
			
			F2[i]=complex(0.,0.);
			f2[i]=complex(0.,0.);
			
			pot[i]=0.;
			
		}
		
		/***********************************************************/	
		//Fill the grids onto the wavefunction
		/***********************************************************/
		
		for(int i=0;i<Nr;i++)
		{
			r[i]=HH.r[i];
			v[i]=HH.v[i];
			dr[i]=HH.dr[i];
			dv[i]=HH.dv[i];
		}
		
		for (int i=0; i<Nz; i++)
			z[i] =(-(Nz-1)/2. + i)*dz;
		
		double dq=dospi/Nz/dz;
		double qmax=pi/dz;
		
		for(int k=0;k<Nz/2;k++)
			q[k]=(k)*dq;
		
		for(int k=Nz/2;k<Nz;k++)
			q[k]=-qmax+(k-Nz/2)*dq;
		
		
		//DFTI_DESCRIPTOR_HANDLE
		hand=0;
		fscale = dz/sqrt(dospi);// 1.;//dz/sqrt(dospi);
		bscale = dq/sqrt(dospi);//1./Nz;//dw/sqrt(dospi);
		
		status = DftiCreateDescriptor(&hand,DFTI_DOUBLE, DFTI_COMPLEX,1,(MKL_LONG)Nz);
		status = DftiSetValue(hand, DFTI_FORWARD_SCALE, fscale);
		status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, bscale);
		status = DftiCommitDescriptor(hand);


		/***********************************************************/	
		//Preparing the arrays fof the Crank
		/***********************************************************/
		
		
		//Evaluation of diagonal up and diagonal down
		for (int i=0; i< Nr; i++)
		{
			a[i]=complex(0.,0.);
			c[i]=complex(0.,0.);			
			b[i]=complex(0.,0.);	
			rv[i]=complex(0.,0.);
			phil[i]=complex(0.,0.);
			gam1[i]=complex(0.,0.);
		}
		
		for(int i=0;i<Nz;i++)
		{
			rv_z[i] =complex(0.,0.);
			phi_z[i]=complex(0.,0.);
			bz[i]=complex(0.,0.);
			gamz[i]=complex(0.,0.);
		}
		
				
		
	}
	
	/***********************************************************/	
	//Norms in both spaces
	/***********************************************************/	
	
	
	double norm()
	{
		double norm=0.0;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				norm+=dz*dr[j]*r[j]*real(conj(phiHank[index(j,i)])*phiHank[index(j,i)]);
		
		return norm;
	}
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr*Nz;i++)
		{
			phiHank[i]=phiHank[i]/sqrt(norm1);
		}
	}
	
	
	double vnorm()
	{
		double vnorm=0.0;
		for(int i=0;i<Nr;i++)
			vnorm+=dv[i]*v[i]*real(conj(f2[i])*f2[i]);
		
		return vnorm;
	}
	
	/***********************************************************/	
	// Scaling the axis
	/***********************************************************/	
	
	
	void phi2F(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				F[index(j,i)]=phiHank[index(j,i)]/(HH.m1[j]);
		
		
	}
	
	void F2phi(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				phiHank[index(j,i)]=F[index(j,i)]*HH.m1[j];
	}
	
	void F22f2(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				f2[index(j,i)]=F2[index(j,i)]*HH.m2[j];
		
	}

	
	void f22F2(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				F2[index(j,i)]=f2[index(j,i)]/HH.m2[j];
		
	}
	
	/***********************************************************/	
	//Hankel Transforms.
	/***********************************************************/	
	
	
	void HankelTransform(HankelMatrix &HH)
	{
		complex alpha=complex(1.,0.);
		complex gamma=complex(0.,0.);
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans, Nr, Nz, Nr,&alpha, HH.C, Nr,F,Nz,&gamma, F2, Nz);
	}
	
	void HankelTransformBack(HankelMatrix &HH)
	{
		complex alpha=complex(1.,0.);
		complex gamma=complex(0.,0.);
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans,Nr,Nz,Nr,&alpha,HH.C,Nr,F2,Nz,&gamma,F, Nz);
	}
	
	
	
	void FFTFor()
	{
		for(int j=0;j<Nr;j++)
			status=DftiComputeForward(hand,&phiHank[j*Nz]);
	}

	void FFTBack()
	{
		for(int j=0;j<Nr;j++)
			status=DftiComputeBackward(hand,&phiHank[j*Nz]);
	}
	
	void PropKinHFFT(double dt)
	{
		double fase=0.;
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=dospi*dospi*q[i]*q[i]*dt/2.;
				phiHank[index(j,i)]*=exp(I*fase);
			}
		
	}
	
	
	
	/***********************************************************/	
	//Preparing the arrays fof the Crank
	/***********************************************************/
	
	void PrepareCrankArrays(double dt)
	{
	
		//Evaluating crank-nicholson formula	
		double hh0 = 0.;
		double dh0 = 0.;	
		double hh  = 0.;
		
		
		//Evaluation of diagonal up and diagonal down
		for (int i=0; i< Nr-1; i++) 
		{
			
			hh  = dr[i+1] + dr[i];
			dh0 = dr[i+1] - dr[i];		
			hh0 = dr[i+1]*dr[i+1] + dr[i]*dr[i];
			
			a[i]     =  ( complex( 0., - 1./hh0  ) +
						 complex( 0., - 1.*dh0/hh/hh0  ) +
						 complex( 0., + 1./2./hh/r[i]  ) 
						 )*dt/2.;
			
			c[i]     =  ( complex( 0., - 1./hh0  ) +
						 complex( 0., + 1.*dh0/hh/hh0  ) +
						 complex( 0., - 1./2./hh/r[i]  )
						 )*dt/2. ;
			
			b[i] =  complex( 1.0, 0.)+ complex( 0.0, 2./hh0 )*dt/2.;
			//b[i] =  complex( 1.0, 0.)+ complex( 0.0, 2./hh0  + v[i] )*dt/2.;
		}
		
		int nend=Nr-1;
		hh  = dr[nend] + dr[nend];
		dh0 = dr[nend] - dr[nend];		
		hh0 = dr[nend]*dr[nend] + dr[nend]*dr[nend];
		
		a[nend]     =  ( complex( 0., - 1./hh0  ) +
						complex( 0., - 1.*dh0/hh/hh0  ) +
						complex( 0., + 1./2./hh/r[nend]  ) 
						)*dt/2.;
		
		c[nend]     =  ( complex( 0., - 1./hh0  ) +
						complex( 0., + 1.*dh0/hh/hh0  ) +
						complex( 0., - 1./2./hh/r[nend]  )
						)*dt/2. ;
		
		b[nend] =  complex( 1.0, 0.)+complex( 0.0, 2./hh0   )*dt/2.;
		//b[nend] =  complex( 1.0, 0.)+complex( 0.0, 2./hh0  + v[nend] )*dt/2.;
		
		
		az = complex( 0., -1./dz/dz/4. )*dt;
		cz = complex( 0., -1./dz/dz/4. )*dt;
		
	}
	
	void KineticPropCrankNonUniform(double dt)
	{
				
		for (int i=0;i<Nz;i++)
		{
			
			//Starting Right
			rv[0] =   conj(a[0])*phiHank[index(1,i)]+conj(b[0])*phiHank[index(0,i)]+conj( c[0] )*phiHank[index(1,i)];	

			for (int j=1; j<Nr-1; j++) 	
				rv[j] = conj(a[j])*phiHank[index(j-1,i)]+conj( b[j] )*phiHank[index(j,i)]+conj(c[j])*phiHank[index(j+1,i)];

			//Finishing right
			rv[Nr-1] = conj(a[Nr-1])*phiHank[index(Nr-2,i)]+conj(b[Nr-1])*phiHank[index(Nr-1,i)];//+conj(c[Nr-1]);//*phiHank[index(Nr-2,i)]*zero;
			
			
			//Solving Triagonal Matrix
			tridag(a,b,c,rv,phil,Nr,gam1);
			
			for(int j=0;j<Nr;j++)
				phiHank[index(j,i)]=phil[j];
	
		}		  
	}
	
	
	void Zprop1( double dt )
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
				a2= complex(0.,1./2./dz/dz  + pot[index(j,i)]/4.) ;
				bz[i]	=	complex(1.,0.) + a2*dt;
			}
			
			
			//Right part Z
			rv_z[0]		=   conj(bz[0])*phiHank[index(j,0)]   + 
			conj( cz  )*phiHank[index(j,1)];	
			
			
			for (int i=1; i<Nz-1; i++) 
				rv_z[i] =	conj( az    )*phiHank[index(j,i-1)] + 
				conj( bz[i] )*phiHank[index(j,i)]   + 
				conj( cz    )*phiHank[index(j,i+1)];
			
			
			rv_z[Nz-1]	=	conj( az  )*phiHank[index(j,Nz-2)]  +
			conj( bz[Nz-1]  )*phiHank[index(j,Nz-1)];	
			//Finishing right part Z
			
			//Zeros on Z
			//zeros(phiHank_z);		
			//for(int i=0;i<Nz;i++)
			//	phiHank_z[i]=complex(0.,0.);
			
			
			//Solving Triagonal Matrix	for Z
			trid_simple(az,bz,cz,rv_z,phi_z,gamz,Nz);
			//tridagS(az,bz,cz,rv_z,phi_z, Nz,gamz);
			
			//Save function 
			for (int i=0; i<Nz; i++)
				phiHank[j*Nz+i] =phi_z[i];		//psi0[index(i,j,nz)] = psi0[index(i,j,nz)];//
		}	
	}
	
		
};

#endif // TOOLS_H
