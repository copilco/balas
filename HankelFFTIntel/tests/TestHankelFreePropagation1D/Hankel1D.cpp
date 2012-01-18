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
	int Nt=1;
	
	
	
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
	
	
	HankelMatrix HH(Nr,200.);
	
	
	
	
	// BLAS parameters
	
	complex alpha=complex(1.,0.);
	complex gamma=complex(0.,0.);
	
	
	FILE *outPhi, *outPhiRetrieved,*normrecord,*timefile,*rhofile;
	outPhi=fopen("outPhi.txt","w+");
	outPhiRetrieved=fopen("outPhiRetrieved.txt","w+");
	normrecord=fopen("normrecord.txt","w+");
	timefile=fopen("timefile.txt","w+");
	rhofile=fopen("rhofile.txt","w+");
	
	// Fill the function
	
	double r0=100.;
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
	
	
	
	int Ng=1000;
	double dt=0.005;
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
		
		printf("Norm frequency space =%e ",norm);
		
		for(int i=0;i<Nr;i++)
		{
			Fp2.v[i]=complex(0.,0.);
			Fretrieved.v[i]=complex(0.,0.);
			phiRetrieved.v[i]=complex(0.,0.);
		}
		
		
		for(int i=0;i<Nr;i++)
		{
			
			phase[i]=dospi*dospi*HH.v[i]*HH.v[i]*dt/2.; //The phase
			//double fase=HH.v[i]*HH.v[i]*dt/2.;
			//Fp2.v[i]=F2.v[i]*complex(cos(fase),-sin(fase) ); //exp(-I*fase);//F2.v[i];//*exp(-I*fase);//phase[i]);
			Fp2.v[i]=F2.v[i]*exp(-I*phase[i]);//F2.v[i];//*exp(-I*fase);//phase[i]);
			//printf("%e %e %e\n",abs(Fp2.v[i]), abs(F2.v[i]),log10(abs(F2.v[i])-abs(Fp2.v[i]) )  );
		}
		 
		
		for(int i=0;i<Nr;i++)
		{		
			fp2.v[i]=Fp2.v[i]*HH.m2[i];
			//fprintf(outPhi,"%e %e\n",f2.v[i].real(),f2.v[i].imag());
		}
		
		norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=HH.dv[i]*HH.v[i]*real(conj(fp2.v[i])*fp2.v[i]);
		}
		
		printf("Norm frequency space (after phase) =%e\n",1.-norm);
		
		// Backwards
		
		cblas_zgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans,Nr,Nt,Nr,&alpha,HH.C,Nr,Fp2.v,Nt,&gamma,Fretrieved.v, Nt);
		
		
		
		for(int i=0;i<Nr;i++)
		{		
			phi.v[i]=Fretrieved.v[i]*HH.m1[i];
			//fprintf(outPhiRetrieved,"%e %e\n",HH.r[i],abs(phi.v[i]));
		}
		
        
		if((gg%10)==0)
		{
			
            for(int i=0;i<Nr;i++)
            {
                fprintf(outPhiRetrieved,"%12.5e\n",HH.r[i]*abs(phi.v[i]));
            }
			//fprintf(timefile,"%e\n",gg*dt);
		}
		
		norm=0.;
		for(int i=0;i<Nr;i++)
		{
			norm+=HH.dr[i]*HH.r[i]*real(conj(phi.v[i])*phi.v[i]);
		}
		printf("Norm retrieved (coordinate space) =%e ",norm);
		printf("Error norm =%e\n",norm-norm_orig);
		fprintf(normrecord,"%e\n",log10(abs(norm-norm_orig)));
		
		
	}
	
	
	for(int i=0;i<Nr;i++)
		fprintf(rhofile,"%12.5e\n",HH.r[i] );
	
	
	fclose(outPhi);
	fclose(outPhiRetrieved);



}