//
//  wave.h
//  
//
//  Created by Camilo Ruiz Méndez on 29/11/11.
//  
//
//
#ifndef WAVEH2D_H
#define WAVEH2D_H

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
	double *r, *dr;
	double *v, *dv, *kr, *dkr;
	double *z,dz,*q,dq,qmax;
	
	complex *F,*phi;
	complex *G,*phik;
	
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
	};
	
	void initialize(HankelMatrix &HH,int _Nz,double _dz)
	{
		Nr=HH.Nr;
		
		Nz=_Nz;
		dz=_dz;
		
		
		//space_indicator=0; // We are in position rho space
		
		
		phi		= new complex[Nr*Nz];		
		F		= new complex[Nr*Nz];
		
		G		= new complex[Nr*Nz];
		phik	= new complex[Nr*Nz];
		
		pot		= new double[Nr*Nz];
				
		r		= new double[Nr];
		dr		= new double[Nr];
		
		v		= new double[Nr];
		dv		= new double[Nr];
		kr		= new double[Nr];
		dkr		= new double[Nr];
		
		
		z		= new double[Nz];
		q		= new double[Nz];
		
		
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
		
		
		/*******************************/	
		//  Fill the arrays with zeros //
		/*******************************/
		
		for(int i=0;i<Nr*Nz;i++)
		{
			phi[i]=complex(0.,0.);
			F[i]=complex(0.,0.);
			
			G[i]=complex(0.,0.);
			phik[i]=complex(0.,0.);
			
			pot[i]=0.;
			
		};
		
		/*****************************************/	
		// Fill the grids onto the wavefunction  //
		/*****************************************/
		
		for(int i=0;i<Nr;i++)
		{
			r[i]=HH.r[i];
			dr[i]=HH.dr[i];
			
			v[i]=HH.v[i];
			dv[i]=HH.dv[i];

			kr[i]=HH.kr[i];
			dkr[i]=HH.dkr[i];
			
		};
		
		for (int i=0; i<Nz; i++){
			//z[i] =(-double(Nz/2.) + i)*dz;
			z[i] =(-double(Nz-1)/2. + i)*dz;
		};
	
		dq=dospi/( (double)(Nz) )/dz;
		
		qmax=pi/dz;
		
		for(int k=0;k<Nz/2;k++)
			q[k]=(k)*dq;
		
		for(int k=Nz/2;k<Nz;k++)
			q[k]=-qmax+(k-double(Nz/2.))*dq;
		
		
		//DFTI_DESCRIPTOR_HANDLE
		hand=0;
		fscale = dz/sqrt(dospi);// 1.;//dz/sqrt(dospi);
		bscale = dq/sqrt(dospi);//1./Nz;//dw/sqrt(dospi);
		
		status = DftiCreateDescriptor(&hand,DFTI_DOUBLE, DFTI_COMPLEX,1,(MKL_LONG)Nz);
		status = DftiSetValue(hand, DFTI_FORWARD_SCALE, fscale);
		status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, bscale);
		status = DftiCommitDescriptor(hand);


		/****************************/	
		//  Defining Crank arrays   //
		/****************************/
		
		
		//Evaluation of diagonal up and diagonal down
		for (int i=0; i< Nr; i++)
		{
			a[i]=complex(0.,0.);
			c[i]=complex(0.,0.);			
			b[i]=complex(0.,0.);	
			rv[i]=complex(0.,0.);
			phil[i]=complex(0.,0.);
			gam1[i]=complex(0.,0.);
		};
		
		for(int i=0;i<Nz;i++)
		{
			rv_z[i] =complex(0.,0.);
			phi_z[i]=complex(0.,0.);
			bz[i]=complex(0.,0.);
			gamz[i]=complex(0.,0.);
		};
		
				
		
	};
	
	
	
	~waveH2D()
	{
		free_memory();
				
	};
	
	
	
	void free_memory()
	{
		
		delete[] phi;		
		delete[] F;
		
		delete[] G;
		delete[] phik;
		
		delete[] pot;
		
		delete[] r;
		delete[] dr;
		
		delete[] v;
		delete[] dv;
		delete[] kr;
		delete[] dkr;
		
		
		delete[] z;
		delete[] q;
		
		
		delete[] rv ;
		delete[] phil;
		
		delete[] a ;
		delete[] b ;
		delete[] c;
		
		delete[] gam1;
		
		delete[] rv_z;
		
		delete[] phi_z;
		delete[] bz ;
		delete[] gamz;		
	};
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/****************************************/	
	//   Set potencial of the hamiltonian   //
	/****************************************/

	
	//************************************************//
	//		       Hydrogen  potential                //
	//************************************************//	
	void set_potential_hlike2D(double charge_nuclei=1., double soft_core=0.)
	{
		for(int j=0; j<Nr; j++)
			for (int i=0; i<Nz; i++)
				pot[index(j,i)]=charge_nuclei*charge_electron_au/sqrt(soft_core+(r[j])*(r[j])+(z[i])*(z[i]));
	
	};
	
	
	
	//************************************************//
	//		       Poschl Teller Potential             //
	//************************************************//	
	void set_poschl_teller2D(double V0, double alpha )
	{		
		double rprime   =0. ;
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				rprime   = sqrt(r[j]*r[j] + z[i]*z[i]); 
				pot[index(j,i)]= -V0*(V0+1.)*alpha*alpha*pow(cosh(alpha*rprime) , -2. );
			}
	};//End Poschl Teller Potential		
	
	//************************************************//
	//		       Hydrogen molecule potential                //
	//************************************************//
	void set_potential_hplus2D( double zmol1, double zmol2, double soft_core1, double soft_core2, double R, double shift )
	{		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				pot[index(j,i)]= charge_electron_au*zmol1/sqrt( soft_core1 + r[j]*r[j] + (z[i]+shift-R)*(z[i]+shift-R) ) +
				charge_electron_au*zmol2/sqrt( soft_core2 + r[j]*r[j] + (z[i]+shift)*(z[i]+shift) );
		
	};//End Hydrogen potential	
	
	
	
	/***************************/	
	//   Norms in both spaces  //
	/***************************/	
	
	
	double norm()
	{
		double norm=0.0;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				norm+=dz*dr[j]*r[j]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return norm;
	};
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr*Nz;i++)
		{
			phi[i]=phi[i]/sqrt(norm1);
		};
	};
	
	
	double qnorm(HankelMatrix &HH)
	{
		double qnorm=0.0;
		
		FFTFor();
		phi2F(HH);
		HankelTransform(HH);
		G2phik(HH);
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				qnorm+=dq*dv[j]*v[j]*real(conj(phik[index(j,i)])*phik[index(j,i)]);
		
		
		phik2G(HH);
		HankelTransformBack(HH);
		F2phi(HH);
		FFTBack();
		
		return qnorm;
	
	};
	
	/*********************/	
	// Scaling the axis  //
	/*********************/	
	
	
	void phi2F(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				F[index(j,i)]=phi[index(j,i)]/(HH.m1[j]);
		
		
	};
	
	void F2phi(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				phi[index(j,i)]=F[index(j,i)]*HH.m1[j];
	};
	
	void G2phik(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				phik[index(j,i)]=G[index(j,i)]*HH.m2[j];
		
	};

	
	void phik2G(HankelMatrix &HH)
	{
		for(int i=0;i<Nz;i++)
			for(int j=0;j<Nr;j++)
				G[index(j,i)]=phik[index(j,i)]/HH.m2[j];
		
	};
	
	/*************************/	
	//  Hankel Transforms    //
	/*************************/	
	
	
	void HankelTransform(HankelMatrix &HH)
	{
		complex alpha=complex(1.,0.);
		complex gamma=complex(0.,0.);
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans, Nr, Nz, Nr,&alpha, HH.C, Nr,F,Nz,&gamma, G, Nz);
	}
	
	void HankelTransformBack(HankelMatrix &HH)
	{
		complex alpha=complex(1.,0.);
		complex gamma=complex(0.,0.);
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans,Nr,Nz,Nr,&alpha,HH.C,Nr,G,Nz,&gamma,F, Nz);
	}
	
	
	
	void FFTFor()
	{
		for(int j=0;j<Nr;j++)
			status=DftiComputeForward(hand,&phi[j*Nz]);
	}

	void FFTBack()
	{
		for(int j=0;j<Nr;j++)
			status=DftiComputeBackward(hand,&phi[j*Nz]);
	}
	
	
	
	/***************/	
	// Propagators //
	/***************/
	
	void prop_kinetic(HankelMatrix &HH, complex dt,double vec_pot1=0.0, double vec_pot2=0.0)
	{
		
		complex fase=complex(0.,0.);
		
		
		FFTFor();        
		phi2F(HH);
		HankelTransform(HH);
		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=1./2.*(dospi*v[j]+charge_electron_au*vec_pot2/lightC_au)*(dospi*v[j]+charge_electron_au*vec_pot2/lightC_au)*dt+1./2.*(q[i]+charge_electron_au*vec_pot1/lightC_au)*(q[i]+charge_electron_au*vec_pot1/lightC_au)*dt;
				G[index(j,i)]*=exp(-I*fase);
			}
		
		
		HankelTransformBack(HH);
		F2phi(HH);
		FFTBack();       
		
		
	}
	
	void prop_potencial(complex &dt)
	{
		
		complex fase=complex(0.,0.);
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=1./2.*pot[index(j,i)]*dt;
				phi[index(j,i)]*=exp(-I*fase);
			}
		
		
	}
	
	
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
				
				pot[index(j,i)]= -( temp0/sqrt(soft_core[0]+sqr0)+
										   temp1/sqrt(soft_core[1]+sqr1));
			}
	}// CO or CO2 Molecular Hydrogen Model  PRL 111, 073902 (2013) 
	
	
	
		
	
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
					phi[index(j,i)].real()*=mask_z_left;	
					phi[index(j,i)].imag()*=mask_z_left;	
				}
				
				if (i> int(Nz*(1.-frac_z_right)))
				{		  
					phi[index(j,i)].real()*=mask_z_right;	
					phi[index(j,i)].imag()*=mask_z_right;	
				}
				
				if (j< int(Nr*frac_r_left))
				{		  
					phi[index(j,i)].real()*=mask_r_left;	
					phi[index(j,i)].imag()*=mask_r_left;	
				}
				
				if (j> int(Nr*(1.-frac_r_right)))
				{		  
					phi[index(j,i)].real()*=mask_r_right;	
					phi[index(j,i)].imag()*=mask_r_right;	
				}
				
			}//End the loop on ij
		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/***********************************************/	
	//   Operators in Schrödinger representation   //
	/***********************************************/
	

	double expectedZ()
	{
		double zExpected=0.0;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				zExpected+=dz*dr[j]*r[j]*z[i]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return zExpected;
	}
	
	double expectedRHO()
	{
		double RhoExpected=0.0;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				RhoExpected+=dz*dr[j]*r[j]*r[j]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		return RhoExpected;
	}
	
	
	double kinetic_energy(HankelMatrix &HH)
	{
		
		FFTFor();
		phi2F(HH);
		HankelTransform(HH);
		G2phik(HH);
		
		
		
		double energy = 0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				energy+=dq*dv[j]*v[j]*((q[i]*q[i])/2.+dospi*dospi*(v[j]*v[j])/2.)*real(conj(phik[index(j,i)])*phik[index(j,i)]);
				//energy+=dq*dkr[j]*kr[j]*((q[i]*q[i])/2.+(kr[j]*kr[j])/2.)/dospi/dospi*real(conj(phik[index(j,i)])*phik[index(j,i)]);
			}
		
		phik2G(HH);
		HankelTransformBack(HH);
		F2phi(HH);
		FFTBack();
		
		
		return energy;
		
	}
	
	
	double pot_energy()
	{
		double potE=0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				potE+=dz*dr[j]*r[j]*pot[index(j,i)]*real(conj(phi[index(j,i)])*phi[index(j,i)]);
		
		
		return potE;
		
	}
	
	
	double projection(waveH2D &wave1)
	{
		double proj = 0.0;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				proj+=dz*dr[j]*r[j]*abs(conj(wave1.phi[index(j,i)])*phi[index(j,i)]);
		
		
		return proj;
	}
	
	
	/****************************/	
	//   Utility functions      //
	/****************************/
	
	void mask_function2D(waveH2D &mom_w,const double &x0,const double &x1,const double &sigma,double rho0=0.0,double z0=0.0)
	{
		
		double x;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				x = sqrt( (r[j]-rho0) * (r[j]-rho0) + (z[i]-z0) * (z[i]-z0) );
				
				if(x<=x0)
					mom_w.phi[index(j,i)]*=0.;
				
				if(x>x0 && x<=x1)
					mom_w.phi[index(j,i)]=phi[index(j,i)]*(1.-exp(-(x-x0)*(x-x0)/16.));
				
				if(x>x1)
					mom_w.phi[index(j,i)]=phi[index(j,i)];
				
			}
		
	}
	
	
	void saveAxes(fstream &file, int skiper1, int skiper2)
	{
		for(int j=0;j<Nr/skiper1;j++)
			file << r[j*skiper1] << endl;
		
		for(int j=0;j<Nz/skiper2;j++)
			file << z[j*skiper2] << endl;
		file.flush();
	}


	
	void saveQAxes(fstream &file, int skiper1, int skiper2)
	{
		for(int j=0;j<Nr/skiper1;j++)
			file << dospi*v[j*skiper1] << endl;
		
		for(int j=0;j<Nz/skiper2;j++)
			file << q[j*skiper2] << endl;
		file.flush();		
	}

	
	void snapshot(fstream &file,int skiper1=1,int skiper2=1)
	{
		
		double norm=0.0;
		
		for(int j=0;j<Nr/skiper2;j++)
			for(int i=0;i<Nz/skiper1;i++)
			{
				norm=dz*dr[j*skiper2]*r[j*skiper2]*real(conj(phi[index(j*skiper2,i*skiper1)])*phi[index(j*skiper2,i*skiper1)]);
				file << norm << endl;
			}
		
		file.flush();		
		
	}
	
	void qsnapshot(fstream &file,HankelMatrix &HH,int skiper1=1,int skiper2=1)
	{
		double qnorm=0.0;
		
		FFTFor();
		phi2F(HH);
		HankelTransform(HH);
		G2phik(HH);
		
		for(int j=0;j<Nr/skiper2;j++)
			for(int i=0;i<Nz/skiper1;i++)
			{
				qnorm=dq*dv[j*skiper2]*v[j*skiper2]*real(conj(phik[index(j*skiper2,i*skiper1)])*phik[index(j*skiper2,i*skiper1)]);
				file << qnorm << endl;
			}
		
		phik2G(HH);
		HankelTransformBack(HH);
		F2phi(HH);
		FFTBack();
		
		file.flush();		

	}

	
	void qsnapshot_wave(fstream &file,HankelMatrix &HH,int skiper1=1,int skiper2=1)
	{
		double qnorm=0.0;
		
		FFTFor();
		phi2F(HH);
		HankelTransform(HH);
		G2phik(HH);
		
		for(int j=0;j<Nr/skiper2;j++)
			for(int i=0;i<Nz/skiper1;i++)
			{
				//qnorm=dq*dv[j]*v[j]*real(conj(phik[index(j*skiper2,i*skiper1)])*phik[index(j*skiper2,i*skiper1)]);
				file << real(phik[index(j*skiper2,i*skiper1)]) << " " << imag(phik[index(j*skiper2,i*skiper1)]) << endl;
			}
		
		phik2G(HH);
		HankelTransformBack(HH);
		F2phi(HH);
		FFTBack();
		
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
	
	
	
	void resize_waveH2D(HankelMatrix &HH, int _Nz, double _dz)
	{
		int Nr_Small = Nr;
		int Nz_Small = Nz;
		
		complex *temp_phi = (complex*)malloc(Nr_Small*Nz_Small*sizeof(complex));
		for (int i=0; i<Nr_Small*Nz_Small; i++)
			temp_phi[i]=phi[i];
		

		free_memory();
		
		
		initialize( HH, _Nz, _dz);	
		
		
		fplaceWF( temp_phi, Nr_Small, Nz_Small);
		
				
		free(temp_phi);
	};
	
	
};

#endif // TOOLS_H
