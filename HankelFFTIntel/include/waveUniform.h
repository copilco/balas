/*
 *  wave.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 29/11/11.
 *  Copyright 2011 USAL. All rights reserved.
 *
 */

#include "tools.h"

class waveUniform
{
	
public:

	int Nr;
	double *r,dr;
	double R;
	
	complex *phi;
	
	complex *rv,*phil;
	complex *a, *b, *c,*gam1;
	
	
	/***********************************************************/	
	//Fill the arrays and copy the grids in rho and k_rho
	/***********************************************************/	
	
	void initialize(HankelMatrix &HH)
	{
		Nr=HH.Nr;
		R=HH.R;
		dr=R/(Nr-1);
		
		phi=new complex[Nr];		

		r=new double[Nr];
		
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
			phi[i]=complex(0.,0.);
			
			a[i]=complex(0.,0.);
			c[i]=complex(0.,0.);			
			b[i]=complex(0.,0.);	
			rv[i]=complex(0.,0.);
			phil[i]=complex(0.,0.);
			gam1[i]=complex(0.,0.);
		}
		
		/***********************************************************/	
		//Fill the grids onto the wavefunction
		/***********************************************************/
		
		r[0]=dr/2.;
		for (int i=1; i<Nr; i++)
			r[i] =r[i-1]+dr ;
		
	}
	
	/***********************************************************/	
	//Norms in both spaces
	/***********************************************************/	
	
	
	double norm()
	{
		double norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=dr*r[i]*real(conj(phi[i])*phi[i]);
		}
		return norm;
	}
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr;i++)
		{
			phi[i]=phi[i]/sqrt(norm1);
		}
	}
	
		
	
	/***********************************************************/	
	//Preparing the arrays fof the Crank
	/***********************************************************/
	
	void PrepareCrankArrays(double dt)
	{
	
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
	}
	
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
		
		
};
