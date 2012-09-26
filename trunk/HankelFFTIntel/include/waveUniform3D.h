//
//  waveUniform3D.h
//  
//
//  Created by Camilo Ruiz Méndez on 29/11/11.
//  
//
//

#ifndef WAVEUNIFORM3D_H
#define WAVEUNIFORM3D_H


#include <complex>
#include <math.h>
#include <new>
#include "omp.h"
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "tools.h"


class waveUniform3D
{
	
public:
  
  int Nr,Nzr,Nzs;
  double *r,dr;
  double *zr,dzr;
  double *zs,dzs;
	
  double R;
  double e_charge;
  
  complex *phi;
  double *pot;
  
  complex azs,czs, azr,czr, *ar, *cr;
	
	
  //complex *rv_r,*rv_z;
  //complex *phi_r,*phi_z;
  //complex *ar, *br, *cr,*gamr;
  //complex az,*bz,cz,*gamz;
  
  
  /***********************************************************/	
  //	Fill the arrays and copy the grids in rho and k_rho
  /***********************************************************/	
  inline int index(int k,int j,int i)
    {
      int indexer=k*Nzr*Nzs+j*Nzs+i;
      return indexer;
    }
  
  void initialize(HankelMatrix &HH,int _Nzr,double _dzr, int _Nzs, double _dzs)
  {
    Nr=HH.Nr;
    R=HH.R;
    dr=R/Nr;
	  dr=0.3;
    
    Nzs=_Nzs;
    dzs=_dzs;
    
    Nzr=_Nzr;
    dzr=_dzr;

    e_charge=-1.;
    
    
    phi=(complex*)mkl_malloc(Nzs*Nzr*Nr*sizeof(complex),16);
    pot=(double*)mkl_malloc(Nzs*Nzr*Nr*sizeof(double),16);
    
    r=(double*)mkl_malloc(Nr*sizeof(double),16);
    zr=(double*)mkl_malloc(Nzr*sizeof(double),16);
    zs=(double*)mkl_malloc(Nzs*sizeof(double),16);
		
    ar = (complex*)mkl_malloc(Nr*sizeof(complex),16);
    cr = (complex*)mkl_malloc(Nr*sizeof(complex),16);
    
    /*
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
    */
    
    /***********************************************************/	
    //	      	Fill the arrays with zeros.
    /***********************************************************/
    
    for (int i=0; i<Nr*Nzr*Nzs; i++)
      {
	phi[i]=complex(0.,0.);
	pot[i]=0.;
      }
    
    /*
      for(int i=0;i<Nr;i++)
      {
      ar[i]    = complex(0.,0.);
      cr[i]    = complex(0.,0.);			
      br[i]    = complex(0.,0.);	
      rv_r[i]  = complex(0.,0.);
      phi_r[i] = complex(0.,0.);
      gamr[i]  = complex(0.,0.);
      }
      
      
      for(int i=0;i<Nz;i++)
      {
      rv_z[i] =	complex(0.,0.);
      phi_z[i]=	complex(0.,0.);
      bz[i]	=	complex(0.,0.);
      gamz[i]	=	complex(0.,0.);
      }
    */
    //***********************************************************	
    //    	  Fill the grids onto the wavefunction
    //***********************************************************
    
    //r[0]=dr/2.;
    
    for (int k=0;k<Nr;k++)
      r[k]=(k+0.5)*dr;
    
    for (int i=0;i<Nzs;i++)
      zs[i]=(-(Nzs-1)/2.+i)*dzs;
    
    for (int j=0;j<Nzr;j++)
      zr[j]=(-(Nzs-1)/2.+j)*dzr;
    
  }
  
  /***********************************************************/	
  //	      	   Norms in both spaces
  /***********************************************************/	
  
	
  double norm()
  {
    double norm=0.0;
    for(int k=0;k<Nr;k++)
      for(int j=0;j<Nzr;j++)
	for(int i=0;i<Nzs;i++)
	  norm+=2.*2.*pi*dzs*dzr*dr*r[k]*real(conj(phi[index(k,j,i)])*phi[index(k,j,i)]);
    
    return norm;
  }
	
  
  //***********************************************************//	
  //		   	Normalize function
  //***********************************************************//
  
  void normalize()
  {
    double norm1=norm();
    for(int i=0;i<Nr*Nzs*Nzr;i++)
      phi[i]=phi[i]/sqrt(norm1);
  }
  
  
  //***********************************************************
  //	       	Preparing the arrays fof the Crank
  //************************************************************
  void PrepareCrankArrays(complex dt )
  {
    int chunk = 8;
    
#pragma omp parallel for shared(chunk,dt) schedule(dynamic,chunk)
    for( int i=0; i<Nr; i++ )
      {						
	ar[i]     =    complex( 0., - 1./dr/dr/4. + 1./dr/r[i]/8. )*dt;			
	
	cr[i]     =    complex( 0., - 1./dr/dr/4. - 1./dr/r[i]/8. )*dt;
      }//End left part
    
    
  }
  
	
  
  
  
  
	//***************************************************
	//    	    Z propagator imaginary time
	//***************************************************	
	/*

	void Zprop( complex dt )	
	{	
		
		az = complex(0.,-1./dz/dz/4.)*dt;
		cz = complex(0.,-1./dz/dz/4.)*dt;
		
		int tid;
		int chunk = 8;
		
		
#pragma omp parallel shared(chunk,dt) private(tid)
		{
	
			complex *rv_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *phi_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *bz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *gamz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			
			
			//Start loop to Nr and Nz 	
#pragma omp for schedule(dynamic,chunk)
			for (int j=0; j<Nr; j++ )
			{      			
				//tid = omp_get_thread_num();
				
				//printf("Thread number %d is doing loop number %d.\n",tid,j);
				
				rv_z[0]	= ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt )*phi[index(j,0)] + 
				- cz*phi[index(j,1)];
				
				bz[0]	=  1. + complex( 0., 1./2./dz/dz + pot[index(j,0)]/4. )*dt;			
				
				for (int i=1; i<Nz-1; i++)
				{
					bz[i]	=  1. + complex( 0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt;
					rv_z[i] = -az*phi[index(j,i-1)]+( 1. - complex(0., 1./2./dz/dz + pot[index(j,i)]/4. )*dt )*phi[index(j,i)] - cz*phi[index(j,i+1)];
				}
				
				bz[Nz-1]	=  1. + complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt;				
				
				rv_z[Nz-1]  = - az*phi[index(j,Nz-2)]  
				+ ( 1. - complex( 0., 1./2./dz/dz + pot[index(j,Nz-1)]/4. )*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix for Z
				trid_simple( az, bz, cz, rv_z, phi_z, gamz, Nz );
				
				
				//Save function 
				for (int i=0; i<Nz; i++)
					phi[index(j,i)] = phi_z[i];
			}//End loop Nr and Nz
			
			
			mkl_free(rv_z);
			mkl_free(bz);
			mkl_free(phi_z);
			mkl_free(gamz);
			
		}
		
	}//Z propagator imaginary time	
	
	
	
	
	
	//*******************************************************
	//	       RHO_Operator imaginary time
	//*******************************************************
	void Rprop( complex dt )
	{
		
		int tid;
		int chunk = 8;
		
		
#pragma omp parallel shared(chunk,dt) private(tid)
		{
			
			complex *rv_r=(complex*)mkl_malloc(Nr*sizeof(complex),16);
			complex *phi_r=(complex*)mkl_malloc(Nr*sizeof(complex),16);
			complex *br=(complex*)mkl_malloc(Nr*sizeof(complex),16);
			complex *gamr=(complex*)mkl_malloc(Nr*sizeof(complex),16);
			
			
#pragma omp for schedule(dynamic,chunk)
			//  RHO_Operator dt	
			for ( int i=0; i<Nz; i++ )
			{
				
				//tid = omp_get_thread_num();
				
				//printf("Thread number %d is doing loop number %d.\n",tid,i);
				
				rv_r[0] = -ar[0]*phi[index(1,i)] +
				(1.  -  complex(0., 1./dr/dr/2. + pot[index(0,i)]/4. )*dt)*phi[index(0,i)]
				-cr[0]*phi[index(1,i)];			
				
				br[0]   = 1. + complex(0., 1./2./dr/dr + pot[index(0,i)]/4. )*dt;			
				
				
				for (int j=1; j<Nr-1; j++)
				{ 
					br[j] = 1. + complex(0., 1./2./dr/dr + pot[index(j,i)]/4. )*dt;	
					
					rv_r[j]	= -ar[j]*phi[index(j-1,i)] + 
					(1. -  complex( 0., 1./dr/dr/2. + pot[index(j,i)]/4. )*dt)*phi[index(j,i)]
					-cr[j]*phi[index(j+1,i)];
				}
				
				br[Nr-1]   =  1. + complex(0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt;			
				
				rv_r[Nr-1] =  -ar[Nr-1]*phi[index(Nr-2,i)]  +
				(1.  -  complex( 0., 1./2./dr/dr + pot[index(Nr-1,i)]/4. )*dt)*phi[index(Nr-1,i)];
				
				
				//Solving Triagonal Matrix
				tridag(ar, br, cr, rv_r, phi_r, gamr, Nr);
				
				
				for (int j=0; j<Nr; j++)
					phi[index(j,i)] = phi_r[j];
			}//End loop on Nz and Nr
			
			
			mkl_free(br);
			mkl_free(rv_r);
			mkl_free(phi_r);
			mkl_free(gamr);
			
		}
		
	}//End Rho_propagator imaginary time
	
	
  	
	//**********************************************************
	//	     Z propagator in velocity gauge or p·A gauge
	//**********************************************************
	void Zprop_PAG( complex dt, double av_z )
	{
		
		int tid;
		int chunk = 8;		
		double avsquare = av_z*av_z;
		
		
		az = complex( + e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;
		cz = complex( - e_charge*av_z/dz/4., -1./dz/dz/4. )*dt;
		
		
#pragma omp parallel shared(chunk,dt) private(tid)
		{
			
			complex *rv_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *phi_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *bz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *gamz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			
			
			//Start loop to Nr and Nz
#pragma omp for schedule(dynamic,chunk)
			for (int j=0; j<Nr; j++ )
			{
				
				bz[0]    =   1. + complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare)/4. )*dt;	
				
				rv_z[0]	 =  (1. - complex(0., 1./2./dz/dz + (pot[index(j,0)] + e_charge*e_charge*avsquare )/4. )*dt)*phi[index(j,0)] 
				- cz*phi[index(j,1)];	
				
				
				for (int i=1; i<Nz-1; i++) {
					
					bz[i]  =   1. + complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt;				
					
					rv_z[i] =	-az*phi[index(j,i-1)]  + 
					(1. - complex(0., 1./2./dz/dz + (pot[index(j,i)] + e_charge*e_charge*avsquare)/4. )*dt )*phi[index(j,i)] 
					-cz*phi[index(j,i+1)];
					
				}
				
				bz[Nz-1]    =   1. + complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4. )*dt;			
				
				rv_z[Nz-1]	=	-az*phi[index(j,Nz-2)]  +
				(1. - complex(0., 1./2./dz/dz + (pot[index(j,Nz-1)] + e_charge*e_charge*avsquare)/4.)*dt )*phi[index(j,Nz-1)];				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, bz, cz, rv_z, phi_z, gamz, Nz);
				
				
				
				//Save function 
				for (int i=0; i<Nz; i++)
					phi[index(j,i)] = phi_z[i];
			}
			
			mkl_free(rv_z);
			mkl_free(bz);
			mkl_free(phi_z);
			mkl_free(gamz);
			
		}
		
	}//Z propagator imaginary time	
	
	
	*/
	
	
	// ****************************************************************//
	//          Z propagator in length gauge or r·E gauge             
	// ****************************************************************//	
	  /*
	void Zprop_REG( complex dt, double efield_z )
	{
		
		int tid;
		int chunk = 8;	
		
		az = complex(0.,-1./dz/dz/4.)*dt;
		cz = complex(0.,-1./dz/dz/4.)*dt;		
		
#pragma omp parallel shared(chunk,dt) private(tid)
		{
			
			complex *rv_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *phi_z=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *bz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			complex *gamz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
			
			
			
			//Start loop to Nr and Nz 
#pragma omp for schedule(dynamic,chunk)
			for (int j=0; j<Nr; j++ )
			{
				
				bz[0]     =    1. + complex( 0., 1./2./dz/dz + (pot[index(j,0)] - e_charge*z[0]*efield_z )/4. )*dt;			
				
				rv_z[0]   =   (1. - complex( 0., 1./2./dz/dz + (pot[index(j,0)] - e_charge*z[0]*efield_z )/4. )*dt )*phi[index(j,0)] 
				- cz*phi[index(j,1)];
				
				
				for (int i=1; i<Nz-1; i++)
				{
					bz[i] = 1. + complex( 0., 1./2./dz/dz + (pot[index(j,i)] - e_charge*z[i]*efield_z )/4. )*dt;				
					
					rv_z[i] =	- az*phi[index(j,i-1)]  + 
					(1. - complex( 0., 1./2./dz/dz + (pot[index(j,i)] - e_charge*z[i]*efield_z )/4. )*dt )*phi[index(j,i)] 
					- cz*phi[index(j,i+1)];
				}
				
				
				bz[Nz-1]    =   1. + complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - e_charge*z[Nz-1]*efield_z )/4. )*dt;			
				
				rv_z[Nz-1]	=  -az*phi[index(j,Nz-2)]  +
				(1. - complex( 0., 1./2./dz/dz + (pot[index(j,Nz-1)] - e_charge*z[Nz-1]*efield_z )/4. )*dt )*phi[index(j,Nz-1)];
				
				
				//Solving Triagonal Matrix	for Z
				trid_simple(az, bz, cz, rv_z, phi_z, gamz, Nz);
				
				
				
				//Save function 
				for (int i=0; i<Nz; i++)
					phi[index(j,i)] = phi_z[i];
			}//end loop Nr and Nz
			
			mkl_free(rv_z);
			mkl_free(bz);
			mkl_free(phi_z);
			mkl_free(gamz);
			
		}
		
	}//Z propagator imaginary time		
	
	
	  */
		
  //************************************************//
  //		       Hydrogen  potential                //
  //************************************************//
  void set_potential_Helium(double _asquare, double _bsquare )
  {		
	  //double asquare=_asquare;		
	  //double bsquare=_bsquare;
	  
	  double asquare=0.135;//0.138;//0.135;//0.01365;//0.135;
	  double bsquare=0.;
	  
	  for(int k=0;k<Nr;k++)
		  for(int j=0;j<Nzr;j++)
			  for(int i=0;i<Nzs;i++)
			  {
				  pot[index(k,j,i)]=
				  -2./sqrt(asquare+r[k]*r[k]/4.+(zs[i]+zr[j]/2.)*(zs[i]+zr[j]/2.))
				  -2./sqrt(asquare+r[k]*r[k]/4.+(zs[i]-zr[j]/2.)*(zs[i]-zr[j]/2.))
				  +1./sqrt(bsquare+r[k]*r[k]+zr[j]*zr[j]);
			  }
	  //= e_charge*_charge_nuclei/sqrt(soft_core + r[j]*r[j]+z[i]*z[i]);
	  
  }//End Hydrogen potential
	
	
  //****************************************************************************
  //          Expected Kinetic energy value finite difference method
  //****************************************************************************
  double kinetic_energy_finite()
  {

    double ham=0.;

    double kinetic_x1=4.;
    double kinetic_x2=1.;
    double kinetic_x3=1.;
    double momentum=1.;
    double split_ion=3.;
    double split_ee=2.;

    complex pzs,pzr,pr;
    
	  for(int k=0;k<Nr;k++)
      {
		  for(int j=0;j<Nzr;j++)
		  {
			  for(int i=0;i<Nzs;i++)
			  {
				  
				  
				  if(i==0)
				  {
					  pzs= ( phi[index(k,j,i+1)]-
							2.*phi[index(k,j,i)])/dzs/dzs/kinetic_x1;
				  }
				  else if(i==(Nzs-1))
				  {
					  pzs=(	phi[index(k,j,i-1)]-
						   2.*phi[index(k,j,i)])/dzs/dzs/kinetic_x1;
				  }
				  else
				  {
					  pzs=(	phi[index(k,j,i-1)]+
						   phi[index(k,j,i+1)]-
						   2.*phi[index(k,j,i)])/dzs/dzs/kinetic_x1;
				  }
				  
				  
				  if(j==0)
				  {
					  pzr=(phi[index(k,j+1,i)]-
						   2.*phi[index(k,j,i)])/dzr/dzr/kinetic_x2;
				  }
				  else if(j==(Nzr-1))
				  {
					  pzr=(phi[index(k,j-1,i )]-1.*phi[index(k,j,i)])/dzr/dzr/kinetic_x2;
				  }
				  else
				  {
					  pzr=(phi[index(k,j-1,i)]+
						   phi[index(k,j+1,i)]-
						   2.*phi[index(k,j,i)])/dzr/dzr/kinetic_x2;
				  }
				  
				  
				  if(k==0)
				  {
					  pr=( 	  phi[index(k+1,j,i)]*(1./dr/dr/kinetic_x3+momentum/r[k]/dr/2.)-
						  2.*phi[index(k,j,i)]/dr/dr/kinetic_x3);
				  }
				  else if(k==(Nr-1))
				  {
					  pr=( phi[index(k-1,j,i)]*(1./dr/dr/kinetic_x3-momentum/r[k]/dr/2.)-
						  2.*phi[index(k,j,i)]/dr/dr/kinetic_x3);
				  }
				  else
				  {
					  pr=( phi[index(k-1,j,i)]*(1./dr/dr/kinetic_x3-momentum/r[k]/dr/2.)+
						  phi[index(k+1,j,i)]*(1./dr/dr/kinetic_x3+momentum/r[k]/dr/2.)-
						  2.*phi[index(k,j,i)]/dr/dr/kinetic_x3);
				  }
				  
				  
				  ham+=real(conj(phi[index(k,j,i)])*(-pzs-pr-pzr+pot[index(k,j,i)]*phi[index(k,j,i)])*2.*2.*pi*r[k]*dzs*dr*dzr);
				  
				  
				  
				  //ham+=abs(conj(phi[index(k,j,i)])*(pot[index(k,j,i)]*phi[index(k,j,i)]))*2.*2.*pi*r[k]*dzs*dr*dzr;

				  
				  //ham+=real(conj(phi[index(k,j,i)])*(pot[index(k,j,i)]*phi[index(k,j,i)])*2.*2.*pi*r[k]*dzs*dr*dzr);
				  
				  //ham+=real(conj(phi[index(k,j,i)])*(-pzs-pr-pzr))*2.*2.*pi*r[k]*dzs*dr*dzr;
				  
				  //+(v[index(k,j,i)]+vr[k*(ndimzr+2)+j] )
				  
				  
				  
			  }
		  }
      }
	  
  return ham;

    /*
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
    /*
      if (i-1>=0)
      psi_ip=phi[index(j,i-1)];
      
      if (i+1<Nz)
      psi_im=phi[index(j,i+1)];
      
      /******************************************************************/
    /*
      if (j-1>=0)
      psi_jp=phi[index(j-1,i)];
      
      if (j+1<Nr)
      psi_jm=phi[index(j+1,i)];
      
      /******************************************************************/
    
    /*		
      
    complex ksquare=				
    -a1*(  psi_ip - dos*psi + psi_im )/dz/dz 
    -a1*( (psi_jp - dos*psi + psi_jm )/dr/dr + 
    ( psi_jp - psi_jm )/r[j]/dr/dos );
    
    
    kinetic+= dz*dr*r[j]*(conj(psi)*ksquare);
    
    }
    return real(kinetic);
    
    */	
  }//End kinetic energy expected value	
  
	
  
  
  //**********************************************************
  //	    Expected Potential energy value
  //**********************************************************
/*	
double potential_energy()
	{		
		double potE=0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				potE+= dz*dr*r[j]*pot[index(j,i)]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return potE;
		
	} 
	// End potential energy expected value	
	
	
*/
	
	
	//*****************************************
	//     Gaussian initial function
	//*****************************************
/*
	void rho_gaussian( double r0, double z0, double rsigma, double zsigma )
		{
		for(int j=0; j<Nr; j++)
			for(int i=0; i<Nz; i++)
				phi[index(j,i)] = exp( -( (r[j]-r0)*(r[j]-r0)/rsigma/rsigma + (z[i]-z0)*(z[i]-z0)/zsigma/zsigma ) );
		}//End Gaussian
	 
*/ 
	 
	 
	 
	//****************************************                                                                            
	//    Gaussian with initial velocity                                                                                          
	//****************************************                                                                                             
/*
	void velocity_gaussian( double r0, double z0, double vr0, double vz0, double rsigma, double zsigma )
	{	    
		for(int j=0; j<Nr; j++)
			for(int i=0; i<Nz; i++)
				phi[index(j,i)] = exp( -((r[j]-r0)*(r[j]-r0)/rsigma/rsigma + (z[i]-z0)*(z[i]-z0)/zsigma/zsigma ) )*exp(I*(vr0*(r[j]-r0) +vz0*(z[i]-z0) ) );
			
	}//End Gaussian Velocity 
*/      
	
	//***************************************
	//               Save axes
	//***************************************	
	/*
	void saveAxes(fstream &file)
	{
		for(int j=0;j<Nr;j++)
			file << r[j] << endl;
		
		
		for(int i=0;i<Nz;i++)
			file << z[i] << endl;
		
		
	}//End save axes
	*/
	
	//***************************************
	//           Mask function 2D
	//***************************************

	
	//******************************************// 
	//                absorber
	//******************************************//

/*
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
					real(phi[index(j,i)])*=mask_r_left;	
					imag(phi[index(j,i)])*=mask_r_left;	
				}
				
			}//End the loop on ij
	}//End absorver
	
*/	
	
	//***********************************************
	//                  Binary writer
	//***********************************************

/*
	void binwrite(FILE *afile )
	{		
		double wreal;
		double wimag;
		
		for(int j=0; j<Nr; j++)
			for(int i=0; i<Nz; i++)
			{
				wreal = real(phi[index(j,i)]);
				wimag = imag(phi[index(j,i)]);
				
				//fwrite(&wreal , 1 , sizeof(wreal) , afile );
				//fwrite(&wimag , 1 , sizeof(wimag) , afile );
				
				fwrite(&wreal , sizeof(wreal) , 1 , afile );
				fwrite(&wimag , sizeof(wimag) , 1 , afile );
			
			}
	}//End Writer
	
	
	
*/	
	
	//***********************************************
	//				  Binary reader
	//***********************************************
/*
	void binread(FILE *afile)
	{	
		double wreal;
		double wimag;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				//fread(&wreal , 1 , sizeof(wreal) , afile );
				//fread(&wimag , 1 , sizeof(wimag) , afile );
				
				fread(&wreal , sizeof(wreal) , 1 , afile );
				fread(&wimag , sizeof(wimag) , 1 , afile );
				
				real(phi[index(j,i)])	= wreal;
				imag(phi[index(j,i)])	= wimag;		
			}	
	}//End reader	
	
*/
	//***********************************************
	//				  Snapshot
	//***********************************************	

/*
	void snapshot(fstream &file,int skiper1=1,int skiper2=1)
	{
		
		double norm=0.0;
		
		for(int j=0;j<Nr/skiper2;j++)
			for(int i=0;i<Nz/skiper1;i++)
			{
				norm=dz*dr*r[j]*real(conj(phi[index(j*skiper2,i*skiper1)])*phi[index(j*skiper2,i*skiper1)]);
				file << norm << endl;
			}
		
	}
*/	
	
	
};
#endif //
