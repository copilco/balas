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

#include "tools.h"

class wave
{
	
public:

	int Nr;
	int space_indicator;  //Defines in which space you are.
	double *r,*v,*dr,*dv;

	complex *F,*phiHank;
	complex *F2,*f2;
	
	complex *rv,*phil;
	complex *a, *b, *c,*gam1;
	
	
	/***********************************************************/	
	//Fill the arrays and copy the grids in rho and k_rho
	/***********************************************************/	
	
	void initialize(HankelMatrix &HH)
	{
		Nr=HH.Nr;
		space_indicator=0;//We are in position rho space
		
		
		phiHank=new complex[Nr];		
		F=new complex[Nr];
		
		F2=new complex[Nr];
		f2=new complex[Nr];
				
		r=new double[Nr];
		v=new double[Nr];
		dr=new double[Nr];
		dv=new double[Nr];
		
		
		rv=    new complex[Nr];
		phil=  new complex[Nr];
		a=     new complex[Nr];
		b=     new complex[Nr];
		c=     new complex[Nr];
		gam1=     new complex[Nr];
		
		
		/***********************************************************/	
		//Fill the arrays with zeros.
		/***********************************************************/
		
		for(int i=0;i<Nr;i++)
		{
			phiHank[i]=complex(0.,0.);
			F[i]=complex(0.,0.);
			
			F2[i]=complex(0.,0.);
			f2[i]=complex(0.,0.);
			
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
				
		
		
	}
	
	/***********************************************************/	
	//Norms in both spaces
	/***********************************************************/	
	
	
	double norm()
	{
		double norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=dr[i]*r[i]*real(conj(phiHank[i])*phiHank[i]);
		}
		return norm;
	}
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr;i++)
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
		for(int i=0;i<Nr;i++)
			F[i]=phiHank[i]/(HH.m1[i]);
		
		
	}
	
	void F2phi(HankelMatrix &HH)
	{
		for(int i=0;i<Nr;i++)
			phiHank[i]=F[i]*HH.m1[i];
	}
	
	void F22f2(HankelMatrix &HH)
	{
		for(int i=0;i<Nr;i++)
			f2[i]=F2[i]*HH.m2[i];
		
	}

	
	void f22F2(HankelMatrix &HH)
	{
		for(int i=0;i<Nr;i++)
			F2[i]=f2[i]/HH.m2[i];
		
	}
	
	/***********************************************************/	
	//Hankel Transforms.
	/***********************************************************/	
	
	
	void HankelTransform(HankelMatrix &HH,int Nt)
	{
		complex alpha=complex(1.,0.);
		complex gamma=complex(0.,0.);
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,&alpha, HH.C, Nr,F,Nt,&gamma, F2, Nt);
	}
	
	void HankelTransformBack(HankelMatrix &HH, int Nt)
	{
		complex alpha=complex(1.,0.);
		complex gamma=complex(0.,0.);
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans,Nr,Nt,Nr,&alpha,HH.C,Nr,F2,Nt,&gamma,F, Nt);
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
	}
	
	void KineticPropCrankNonUniform(double dt)
	{
		//Starting Right
		rv[0] =   conj(a[0])*phiHank[1]+conj(b[0])*phiHank[0]+conj( c[0] )*phiHank[1];	
		
		for (int i=1; i<Nr-1; i++) 	
			rv[i] = conj(a[i])*phiHank[i-1]+conj( b[i] )*phiHank[i]+conj(c[i])*phiHank[i+1];

		//Finishing right
		complex zero=complex(0.,0.);
		rv[Nr-1] = conj(a[Nr-1])*phiHank[Nr-2]+conj(b[Nr-1])*phiHank[Nr-1]+conj(c[Nr-1])*phiHank[Nr-2]*zero;
		
		
		
		//Solving Triagonal Matrix
		tridag(a,b,c,rv,phil,Nr,gam1);
		
		for(int i=0;i<Nr;i++)
			phiHank[i]=phil[i];
	}
		
		
};

#endif // TOOLS_H
