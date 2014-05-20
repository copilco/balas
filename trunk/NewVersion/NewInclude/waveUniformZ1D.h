/*
 *  wave.h
 *  
 *
 *  Created by Alexis Chacón on 26/03/2012.
 *  Copyright 2011 USAL. All rights reserved.
 *
 */

#ifndef WAVEUNIFORMZ1D_H
#define WAVEUNIFORMZ1D_H

#include <complex>
#include <math.h>
#include <new>
#include "omp.h"
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "constant.h"
#include "tools.h"

class waveUniformZ1D
{
	
public:
	
	int n;
	int Nz;
	double *z; 
	double dz;
	double e_charge;
	double *pot;
	
	
	
	double kmomenta;
	double fn_m;
	double fn;
	double fn_p;		
	double kinetic_value;
	double me;
	double a1;

	
	
	complex *phi;
	complex *rv_z;
	complex *phi_z;
	
	complex az;
	complex *bz;
	complex cz;	
	complex *gamz;

	complex psin_1;
	complex Dpsin_1;	
	
	
	
	
	/***********************************************************/	
	//Fill the arrays and copy the grids in rho and k_rho
	/***********************************************************/	
	
	void initialize(int _nz, double _dz)
	{
		Nz    = _nz;
		dz    = _dz;
		
		z     = (double*)mkl_malloc(Nz*sizeof(double),16);		
		phi   = (complex*)mkl_malloc(Nz*sizeof(complex),16);
		pot   = (double*)mkl_malloc(Nz*sizeof(double),16);		
		
		rv_z  = (complex*)mkl_malloc(Nz*sizeof(complex),16);
		phi_z = (complex*)mkl_malloc(Nz*sizeof(complex),16);
		bz    = (complex*)mkl_malloc(Nz*sizeof(complex),16);
		gamz  = (complex*)mkl_malloc(Nz*sizeof(complex),16);		
		
		e_charge	= -1.;
		me			=  1.;
		a1			=  1./2./me;
		
		/***********************************************************/	
		//Fill the arrays with zeros.
		/***********************************************************/
	
		for(int i=0;i<Nz;i++)
		{
			phi[i]    = complex(0.,0.);
			phi_z[i]  = complex(0.,0.);
			pot[i]    = 0.;
			
			bz[i]     = complex(0.,0.);	
			rv_z[i]   = complex(0.,0.);
			gamz[i]   = complex(0.,0.);				
		
		}		
		
		
		
		
		/***********************************************************/	
		          //Fill the grids onto the wavefunction
		/***********************************************************/
		for (int i=0;i<Nz;i++)
			z[i] =(-(Nz-1)/2. + i)*dz;
		
	}
	
	
	
	
	
	/***********************************************************/	
	                 //Norms in both spaces
	/***********************************************************/	
	double norm()
	{
		double norm0=0.0;
		for(int i=0;i<Nz;i++)
			norm0+= dz*abs(phi[i])*abs(phi[i]);
		
		return norm0;
	}

	
	
	
	
	
	/***********************************************************/	
	             //  normalize in both spaces
	/***********************************************************/	
	void normalize()
	{
		double norm0=norm();
		for(int i=0;i<Nz;i++)
			phi[i]=phi[i]/sqrt(norm0);
		
	}
	
		
	
	
	
	//***************************************************
	//    	    Z propagator imaginary time
	//***************************************************	
	void Zprop( complex dt )	
	{	
		
		
		az = complex(0., -1./dz/dz/4.)*dt;
		cz = complex(0., -1./dz/dz/4.)*dt;
		
		
		
		
		bz[0]	=   1. + complex( 0., 1./2./dz/dz + pot[0]/2. )*dt;	
		rv_z[0]	= ( 1. - complex( 0., 1./2./dz/dz + pot[0]/2. )*dt )*phi[0] - cz*phi[1];
		
		
		
		
		for (int i=1; i<Nz-1; i++)
			{
				bz[i]	=  1. + complex( 0., 1./2./dz/dz + pot[i]/2. )*dt;
				rv_z[i] = -az*phi[i-1] + ( 1. - complex(0., 1./2./dz/dz + pot[i]/2. )*dt )*phi[i] - cz*phi[i+1];
			}
		
		
		bz[Nz-1]	=  1. + complex( 0., 1./2./dz/dz + pot[Nz-1]/2. )*dt;				
		rv_z[Nz-1]  = - az*phi[Nz-2] + ( 1. - complex( 0., 1./2./dz/dz + pot[Nz-1]/2. )*dt )*phi[Nz-1];	
				
		
		
		
		//Solving Triagonal Matrix for Z
		trid_simple( az, bz, cz, rv_z, phi_z, gamz, Nz );
				
				
		
		//Save function 
		for (int i=0; i<Nz; i++)
			phi[i] = phi_z[i];
						
			
	}
		
	
	
	
	
	
	
	//**********************************************************
	//	     Z propagator in velocity gauge or p·A gauge
	//**********************************************************
	void Zprop_PAG( complex dt, double av_z )
	{
		
		double avsquare = av_z*av_z;
		az = complex( + e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;
		cz = complex( - e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;
		
			
				
		bz[0]    =   1. + complex(0., 1./2./dz/dz + (pot[0] + e_charge*e_charge*avsquare )/2. )*dt;				
		rv_z[0]	 =  (1. - complex(0., 1./2./dz/dz + (pot[0] + e_charge*e_charge*avsquare )/2. )*dt)*phi[0] - cz*phi[1];	
				
				
		for (int i=1; i<Nz-1; i++) 
		  {
					
			bz[i]   =   1. + complex(0., 1./2./dz/dz + (pot[i] + e_charge*e_charge*avsquare)/2. )*dt;
			rv_z[i] =  -az*phi[i-1]  + 
			           (1. - complex(0., 1./2./dz/dz + (pot[i] + e_charge*e_charge*avsquare)/2. )*dt )*phi[i] - cz*phi[i+1];
					
		  }
				
		bz[Nz-1]    = 1. + complex(0., 1./2./dz/dz + (pot[Nz-1] + e_charge*e_charge*avsquare)/2. )*dt;
		rv_z[Nz-1]	= -az*phi[Nz-2] + (1. - complex(0., 1./2./dz/dz + (pot[Nz-1] + e_charge*e_charge*avsquare)/2.)*dt )*phi[Nz-1];
				
		
		
		//Solving Triagonal Matrix	for Z
		trid_simple(az, bz, cz, rv_z, phi_z, gamz, Nz);
				
				
				
		//Save function 
		for (int i=0; i<Nz; i++)
			phi[i] = phi_z[i];
					
		
	}//Z propagator imaginary time	
	
	
	
	
	
	//****************************************************************//
	//          Z propagator in length gauge or r·E gauge             
	//****************************************************************//	
	void Zprop_REG( complex dt, double efield_z )
	{
		
		az = complex(0.,-1./dz/dz/4.)*dt;
		cz = complex(0.,-1./dz/dz/4.)*dt;		
		
		
		
		
		bz[0]     =    1. + complex( 0., 1./2./dz/dz + (pot[0] - e_charge*z[0]*efield_z )/2. )*dt;							
		rv_z[0]   =   (1. - complex( 0., 1./2./dz/dz + (pot[0] - e_charge*z[0]*efield_z )/2. )*dt )*phi[0] - cz*phi[1];
				
				
		for (int i=1; i<Nz-1; i++)
			{
				bz[i]   =  1. + complex( 0., 1./2./dz/dz + (pot[i] - e_charge*z[i]*efield_z )/2. )*dt;				
				rv_z[i] =  - az*phi[i-1]  + 
					(      1. - complex( 0., 1./2./dz/dz + (pot[i] - e_charge*z[i]*efield_z )/2. )*dt )*phi[i] - cz*phi[i+1];
			}
				
				
		bz[Nz-1]    =   1. + complex( 0., 1./2./dz/dz + (pot[Nz-1] - e_charge*z[Nz-1]*efield_z )/2. )*dt;				
		rv_z[Nz-1]	=  -az*phi[Nz-2]  +
				       (1. - complex( 0., 1./2./dz/dz + (pot[Nz-1] - e_charge*z[Nz-1]*efield_z )/2. )*dt )*phi[Nz-1];
				
		
		
		
		//Solving Triagonal Matrix	for Z
		trid_simple(az, bz, cz, rv_z, phi_z, gamz, Nz);
				
				
				
		//Save function 
		for (int i=0; i<Nz; i++)
			phi[i] = phi_z[i];

		
	}//Z propagator imaginary time			

	
	
	
	
	
	//************************************************//
	//		       Hydrogen  potential                //
	//************************************************//
	void set_potential_hlike1D(double _charge_nuclei, double _soft_core )
	{		
		double soft_core=_soft_core;		
		
		for(int i=0;i<Nz;i++)
			pot[i]= e_charge*_charge_nuclei/sqrt(soft_core + z[i]*z[i]);
		
	}//End Hydrogen potential	
	
	
	
	
	
	
	
	//****************************************************************************
	//          Expected Kinetic energy value finite difference method
	//****************************************************************************
	double kinetic_energy_finite()
	{
		

		
		complex kinetic=complex(0.,0.);
		
		
		for(int i=0;i<Nz;i++)
			{
				
				complex psi		= phi[i];				
				complex psi_ip	= complex(0.,0.);
				complex psi_im	= complex(0.,0.);
				
				
				/******************************************************************/
				if (i-1>=0)
					psi_ip=phi[i-1];
				
				if (i+1<Nz)
					psi_im=phi[i+1];
								
				/******************************************************************/
				
				
				complex ksquare = -a1*(  psi_ip - 2.*psi + psi_im )/dz/dz;

				kinetic+= dz*conj(psi)*ksquare;
				
			}
		return real(kinetic);	
	}//End kinetic energy expected value	
	
	
	
	
	
	
	
	//**********************************************************
	//	    Expected Potential energy value
	//**********************************************************
	double potential_energy()
	{		
		double potE=0.;
		
		for(int i=0;i<Nz;i++)
			potE+= dz*pot[i]*real(conj(phi[i])*phi[i]);
		
		return potE;
		
	} 
	// End potential energy expected value	
	
	
	
	
	
	
	
	//*****************************************
	//     Gaussian initial function
	//*****************************************
	void gaussian1D( double z0, double zsigma )
	{
		for(int i=0; i<Nz; i++)
			phi[i] = exp( - (z[i]-z0)*(z[i]-z0)/zsigma/zsigma  );
	}//End Gaussian	
	
	
	
	
	
	
	
	
	
	
	//***************************************
	//           Mask function 1D
	//***************************************
	void mask_function1D(waveUniformZ1D &wp,const double &x0,const double &x1,const double &sigma, double z0=0.0)
	{
		double x;
		

		for(int i=0;i<Nz;i++)
			{
				x = sqrt( (z[i]-z0) * (z[i]-z0) );
				
				if(x<=x0)
					wp.phi[i]=0.;
				
				if(x>x0 && x<=x1)
					wp.phi[i]= phi[i]*(1.-exp(-(x-x0)*(x-x0)/sigma/sigma));
				
				if(x>x1)
					wp.phi[i]=phi[i];
				
			}
		
	}//End mask function 1D
	
	
	
	
	
	
	
	
	
	//******************************************// 
	//                absorber
	//******************************************// 	
	void absorber1D(const double &frac_z_left, const double &frac_z_right,const double &exponent)
	{
		
		// Try 1./6. in the exponent.
		
		double mask_start_z_right=z[int(Nz*(1.-frac_z_right))];		
		double mask_start_z_left=z[int(Nz*frac_z_left)+1];
		
		
		double argument_z_right;
		double argument_z_left;	
		
		double mask_z_right;		
		double mask_z_left;
		
		
		
		for(int i=0;i<Nz;i++)
			{
				argument_z_right=(pi/2.)*(z[i]-mask_start_z_right)/(z[Nz-1]-mask_start_z_right+1.e-20);
				
				argument_z_left = (pi/2.)*(z[i]-mask_start_z_left)/(z[0]-mask_start_z_left+1.e-20);				
				
				mask_z_right=pow(fabs(cos(argument_z_right)),exponent);	
			
				mask_z_left=pow(fabs(cos(argument_z_left)),exponent);
				
				
				
				if (i< int(Nz*frac_z_left))
				{		  
					real(phi[i])*=mask_z_left;	
					imag(phi[i])*=mask_z_left;	
				}
				
				if (i> int(Nz*(1.-frac_z_right)))
				{		  
					real(phi[i])*=mask_z_right;	
					imag(phi[i])*=mask_z_right;	
				}
				
				
			}//End the loop on ij
	}//End absorver	
	
	
	
	
	
	
	
	
	
	//***********************************************
	//                  Binary writer
	//***********************************************
	void binwrite(FILE *afile )
	{		
		double wreal;
		double wimag;
		
		for(int i=0; i<Nz; i++)
		{
			wreal = double(real(phi[i]));
			wimag = double(imag(phi[i]));
				
			fwrite(&wreal , 1 , sizeof(wreal) , afile );
			fwrite(&wimag , 1 , sizeof(wimag) , afile );
			
		}
	}//End Writer
	
	
	
	
	
	
	
	
	//***********************************************
	//				  Binary reader
	//***********************************************
	void binread(FILE *afile)
	{	
		double wreal;
		double wimag;
		
		for(int i=0;i<Nz;i++)
		{
			fread(&wreal , 1 , sizeof(wreal) , afile );
			fread(&wimag , 1 , sizeof(wimag) , afile );
				
			//fread(&wreal , sizeof(wreal) , 1 , afile );
			//fread(&wimag , sizeof(wimag) , 1 , afile );
				
			real(phi[i])	= wreal;
			imag(phi[i])	= wimag;		
		}
	}//End reader	
	
	
	
	void saveAxes(FILE *file,int skiper1=1)
	{
		for(int i=0;i<Nz/skiper1;i++)
			fprintf(file,"%e\n",z[i*skiper1]);
	
	}
	
	
	
	//***********************************************
	//				  Snapshot
	//***********************************************	
	void snapshot(FILE *file,int skiper1=1)
	{		
		for(int i=0;i<Nz/skiper1;i++)
				fprintf(file,"%12.16e \n",dz*abs(phi[i*skiper1])*abs(phi[i*skiper1]));
	}		
	

	
	void wavesnapshot(FILE *file,int skiper1=1)
	{		
		for(int i=0;i<Nz/skiper1;i++)
			fprintf(file,"%12.16e %12.16e\n",real(phi[i*skiper1]),imag(phi[i*skiper1]));
	}			
	
	
	

	
	/*************************************************
            Two conditions for Numerov method
	 **************************************************/
	void Numerov1D() 
	{
		//Kinetic energy for momentum kmomenta
		kinetic_value = kmomenta*kmomenta/2.;
		fn_m = 0.;
		fn	 = 0.;
		fn_p = 0.;
		
		complex phim;
		complex phin;
		complex phip;
		
		
		if (kmomenta>=0){ //Positive Momentum condition
			for (int i=1; i<Nz-1; i++) {
				//Evaluating variable to Numerov's formula				
				fn_m	    = 2.*me*(kinetic_value-pot[i-1]);
				fn			= 2.*me*(kinetic_value-pot[i]);
				fn_p	    = 2.*me*(kinetic_value-pot[i+1]);
				
				
				phim = phi[i-1];
				phin = phi[i];	
				
				
				//Doing Numerov's Formula from left to Right
				phip = (  2.*( 1.  -  5.*dz*dz*fn/12.   )*phin 
								   - ( 1.    +   dz*dz*fn_m/12.  )*phim  )
				                    /( 1.    +   dz*dz*fn_p/12.  );
				
				
				phi[i+1] = phip;
				
			}//End loop positive momentum
		}else { //Negative Momentum condition
			 for (int i=Nz-2; i>0; i--) {
				//Evaluating variable to Numerov's formula
				fn_m	    = 2.*me*(kinetic_value-pot[i-1]);
				fn			= 2.*me*(kinetic_value-pot[i]);
				fn_p	    = 2.*me*(kinetic_value-pot[i+1]);
				
				
				phip		= phi[i+1];
				phin		= phi[i];	
				
				//Doing Numerov's Formula from right to left
				phim		= (  2.*( 1.  -  5.*dz*dz*fn/12. )*phin 
								  -  ( 1.  +   dz*dz*fn_p/12. )*phip  )
				                     /( 1.  +  dz*dz*fn_m/12.  );
				
				
				phi[i-1]    = phim;
				
			}//End loop for	
		}// end Else condition for negative momentum*/ 
		
		
		
		/*
		for (int i=1; i<Nz-1; i++) {
			//Evaluating variable to Numerov's formula				
			fn_m	    = 2.*me*(kinetic_value-pot[i-1]);
			fn			= 2.*me*(kinetic_value-pot[i]);
			fn_p	    = 2.*me*(kinetic_value-pot[i+1]);
			
			
			phim = phi[i-1];
			phin = phi[i];	
			
			
			//Doing Numerov's Formula from left to Right
			phip = (  2.*( 1.  -  5.*dz*dz*fn/12.   )*phin 
					- ( 1.    +   dz*dz*fn_m/12.  )*phim  )
			/( 1.    +   dz*dz*fn_p/12.  );
			
			
			phi[i+1] = phip;
			
		}//End loop positive momentum	*/
		
		
	} //End Numerov TwoWay
	
	
	
	
	
	
	
	
	
	/*******************************************
	  //Derivative finite method third steps
	*******************************************/
	inline complex diff_th( complex yn1, complex yn_1, double h){
		return (yn1-yn_1)/2./h;
	}
	
	
	
	
	
	
	
	
	
	/********************************************************
	//Complex product of the wavefunction by a scalar
	 ****************************************************/
	inline void complexproduct(complex a){
		//This function does the product of a wavefunction vector by an scalar
		//complex factor=1./sqrt(2.*pi)*a;//sqrt((w.n1-1)*w.dx1);//sqrt(w.dx1);
		
		for (int i=0; i<Nz; i++) 
			phi[i]=phi[i]*a;
	}
	
	
	
	
	
	
	
	
	/*******************************************
	//Norm factor from scattering wave
	 ***************************************/
	inline complex norm_factor( complex psi, complex Dpsi, double k, double xn)
	{
		return 2.*k*exp( I*k*xn )/( k*psi - I*Dpsi)/sqrt(2.*pi);
		//2.*k*complex( cos(k*xn) , sin(k*xn) )/( k*psi - I*Dpsi);
	}
	
	
	
	
	
	
	
	
	/******************************************************************
	//New continuum wave function building each of the momentum k
	 ************************************************************/
	void Continuum_WF( )
	{
		
		
		/*if (kmomenta>=0)
			n = Nz-2; //point to derivative the wave function solution
		else 
			n = 1;
		
		
		//Set intial condition
		real(phi[0]) = cos( kmomenta*z[0] );
		imag(phi[0]) = sin( kmomenta*z[0] );
		
		real(phi[1]) = cos( kmomenta*z[1] );
		imag(phi[1]) = sin( kmomenta*z[1] );
		
		
		//cout << "\nPhase To Positive Condition";
		//cout << "\nkmomenta*z[0]=" <<kmomenta*z[0]<<endl;
		//cout << "\nkmomenta*z[1]=" <<kmomenta*z[1];	*/
		
		if (kmomenta>=0){//POSITIVE MOMENTUM

			
			n = Nz-2; //point to derivative the wave function solution
			
			//Set intial condition
			real(phi[0]) = cos( kmomenta*z[0] );
			imag(phi[0]) = sin( kmomenta*z[0] );
			
			real(phi[1]) = cos( kmomenta*z[1] );
			imag(phi[1]) = sin( kmomenta*z[1] );
			
			
			//cout << "\nPhase To Positive Condition";
			//cout << "\nkmomenta*z[0]=" <<kmomenta*z[0]<<endl;
			//cout << "\nkmomenta*z[1]=" <<kmomenta*z[1];			
			
			
		}
		else{//NEGATIVE MOMENTUM
			n = 1;
			
			//Set intial condition			
			real(phi[Nz-1]) = cos( kmomenta*z[Nz-1] );
			imag(phi[Nz-1]) = sin( kmomenta*z[Nz-1] );			
			
			real(phi[Nz-2]) = cos( kmomenta*z[Nz-2] );
			imag(phi[Nz-2]) = sin( kmomenta*z[Nz-2] );			

			//cout << "\nPhase To Negative Condition";
			//cout << "\nkmomenta*z[Nz-1]=" <<kmomenta*z[Nz-1]<<endl;
			//cout << "\nkmomenta*z[1]=" <<kmomenta*z[1];						
			
						
		} //*/
		
		
		
		//Applying Numerov function
		Numerov1D();
		
		//Building variable to normalization factor
		psin_1			=  phi[n];
		Dpsin_1			=  diff_th(phi[n+1],phi[n-1],dz);		
		
		
		
		//calculing normalization factor
		complex factor	= norm_factor( psin_1, Dpsin_1, kmomenta, z[n]);	
		complexproduct(factor);
		
		
	}//End continuum wavefunction
	
	
	
	
	
	
	
	//*****************************************************************************
	//      Reemplace wavefunction wsmall in other bigger wavefunction wlarge
	//*****************************************************************************
	void placeWF( waveUniformZ1D &wlarge, waveUniformZ1D wsmall)
	{	
		int iplacer=floor((wlarge.Nz - wsmall.Nz)/2);
		
		for(int i=0; i<wlarge.Nz; i++)
			wlarge.phi[i]=complex(0.,0.);
		

		for(int i=0; i< wsmall.Nz; i++)
			wlarge.phi[iplacer+i] = wsmall.phi[i];
		
	}//End reemplace wavefunction	
	

	
	
	
	//**********************************************************
	//                 Projection function
	//**********************************************************
	complex projection(waveUniformZ1D w0, waveUniformZ1D w )
	{
		
		complex proj=complex(0.,0.);
		
		for (int i=0; i<w.Nz; i++)
			proj+= conj(w0.phi[i])*w.phi[i]*w.dz;
		
		return proj;
	}//End project function	
	
	
	
	
	
	
	
	//*********************************************************************
	//        Project out removes wavefunction w0 of wavefunction w
	//*********************************************************************
	void project_out(waveUniformZ1D w0, waveUniformZ1D &w )
	{
		
		complex proj= w.projection(w0, w);
		
		for (int i=0; i<w.Nz; i++)
			w.phi[i]-= proj*w0.phi[i];
		
	}//End of project out	
	
	
	
	
	
	
};
#endif // TOOLS_H
