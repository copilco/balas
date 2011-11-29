/*
 *  Hankel1D.cpp
 *  
 *
 *  Created by Alejandro de la Calle on 28/09/11.
 *  
 *
 */

#include <iostream>
#include <math.h>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "constant.h"
#include "arrai.h"
#include "HankelMatrix.h"


int main()
{
	// Parameters!!
		
	int Nr=1000;
	double dr=0.03;
	//double r;
	double r0=10.;
	int Nt=1;
	
	
	double dt=0.1;
	
	arrai phi(Nr,1);
	arrai F(Nr,1);
	arrai F2(Nr,1);
	arrai f2(Nr,1);
	arrai Fp2(Nr,1);
	arrai fp2(Nr,1);
	arrai phi2(Nr,1);
	arrai Fretrieved(Nr,1);
	arrai phiRetrieved(Nr,1);
	double phase[Nr];
	
	double aux0r;
	double aux0i;
	double K=6.4387e+04;
	
	
	HankelMatrix HH(Nr,200.);
	
	
	
	
	// BLAS parameters
	
	complex alpha=complex(1.,0.);
	complex gamma=complex(0.,0.);
	
	
	FILE *outPhi, *outPhiRetrieved;
	outPhi=fopen("outPhi.txt","w+");
	outPhiRetrieved=fopen("outPhiRetrieved.txt","w+");
	
	// Fill the function
	
	for(int i=0;i<Nr;i++)
	{
		//r=(-(Nr+1)/2.+i)*dr;
		//rho=dr*i/2.;
		phi.v[i]=exp(-(HH.r[i]-r0)*(HH.r[i]-r0)/0.5/0.5);
	}
	
	// Check the initial normalization
	
	double norm=0.0;
	for(int i=0;i<Nr;i++)
	{
		norm+=HH.dr[i]*HH.r[i]*real(conj(phi.v[i])*phi.v[i]);
	}
	
	for(int i=0;i<Nr;i++)
	{
		phi.v[i]=phi.v[i]/sqrt(norm);
	}
	
	norm=0.0;
	for(int i=0;i<Nr;i++)
	{
		norm+=HH.dr[i]*HH.r[i]*real(conj(phi.v[i])*phi.v[i]);
	}
	
	printf("Norm coordinate space =%e\n",norm);
	double norm_orig=norm;
	
	// Scale F
	
	
	
	int Ng=100;
	//double ddr=0.05;
	int Nsnap=20;
    
	// Add the phase
	for (int gg=0;gg<Ng;gg++)
	{
		
		
		for(int i=0;i<Nr;i++)
		{			
			F.v[i]=phi.v[i]/(HH.m1[i]);
			//fprintf(outPhi,"%12.5e %12.5e\n",HH.r[i],abs(phi.v[i]));
			
		}
	
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans, Nr, Nt, Nr,&alpha, HH.C, Nr,F.v,Nt,&gamma, F2.v, Nt);
		
		
		for(int i=0;i<Nr;i++)
		{		
			f2.v[i]=F2.v[i]*HH.m2[i];
		}
		
		
		norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=HH.dv[i]*HH.v[i]*real(conj(f2.v[i])*f2.v[i]);
		}
		
		printf("Norm frequency space =%e\n",norm);
		
		for(int i=0;i<Nr;i++)
		{
			Fp2.v[i]=complex(0.,0.);
			Fretrieved.v[i]=complex(0.,0.);
			phiRetrieved.v[i]=complex(0.,0.);
		}
		
		
		for(int i=0;i<Nr;i++)
		{
			
			phase[i]=HH.v[i]*HH.v[i]*dt/2.; //The phase
			Fp2.v[i]=F2.v[i]*exp(I*phase[i]);
		}
		
		for(int i=0;i<Nr;i++)
		{		
			fp2.v[i]=Fp2.v[i]*HH.m2[i];
			fprintf(outPhi,"%e %e\n",f2.v[i].real(),f2.v[i].imag());
		}
		
		norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=HH.dv[i]*HH.v[i]*real(conj(fp2.v[i])*fp2.v[i]);
		}
		
		//printf("Norm frequency space (after phase) =%e\n",norm);
		
		// Backwards
		
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans,Nr,Nt,Nr,&alpha,HH.C,Nr,Fp2.v,Nt,&gamma,Fretrieved.v, Nt);
		
		
		
		for(int i=0;i<Nr;i++)
		{		
			phi.v[i]=Fretrieved.v[i]*HH.m1[i];
			//fprintf(outPhiRetrieved,"%e %e\n",HH.r[i],abs(phi.v[i]));
		}
		
        
        if (gg%(Ng/(Nsnap-1))==0 )
            for(int i=0;i<Nr;i++)
            {
                fprintf(outPhiRetrieved,"%12.5e\n",abs(phi.v[i]));
            }
		
		
		norm=0.;
		for(int i=0;i<Nr;i++)
		{
			norm+=HH.dr[i]*HH.r[i]*real(conj(phi.v[i])*phi.v[i]);
		}
		printf("Norm retrieved (coordinate space) =%e\n",norm);
		
		printf("Error norm =%e\n",norm-norm_orig);
		
		
		
	}
	 
	 
	
	fclose(outPhi);
	fclose(outPhiRetrieved);



}