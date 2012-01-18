/*
 *  mainTestWaveH.cpp
 *  
 *
 *  Created by camilo on 29/11/11.
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
#include "wave.h"
#include "waveUniform.h"
#include "interp.h"


int main()
{
	
	FILE *out1;
	out1=fopen("out1.txt","w");
	
	FILE *out2;
	out2=fopen("out2.txt","w");
	
	// Parameters!!
	
	int Nr=1000;
	int Nt=1;
	
	HankelMatrix HH(Nr,200.);
	
	wave wHank;
	wHank.initialize(HH);
	
	printf("%d\n",wHank.Nr);
	
	double r0=100.;
	double sigma=10.;
	for(int i=0;i<HH.Nr;i++)
	{
		wHank.phiHank[i]=exp(-(wHank.r[i]-r0)*(wHank.r[i]-r0)/sigma/sigma);
	}
	
	printf("%e\n",wHank.norm());
	wHank.normalize();
	printf("%e\n",wHank.norm());
	
	
	waveUniform w;
	w.initialize(HH);
	
	
	double dt=0.5;
	double *fase=new double[Nr];	


	
	//w.PrepareCrankArrays(dt);
	for (int ktime=0; ktime<500; ktime++)
	{
		//w.KineticPropCrankUniform(dt);

		wHank.phi2F(HH);
		wHank.HankelTransform(HH, Nt);
		
		for (int i=0;i<Nr;i++)
		{
			fase[i]=dospi*dospi*HH.v[i]*HH.v[i]*dt/2.;
			wHank.F2[i]=wHank.F2[i]*exp(I*fase[i]);
		}
		
		wHank.HankelTransformBack(HH, Nt);

		wHank.F2phi(HH);
		
		interpH2U(wHank, w);
		
		
		for (int i=0; i<HH.Nr; i++)
			fprintf(out1,"%10.17e %10.17e %10.17e %10.17e \n", wHank.r[i],abs(wHank.phiHank[i]), w.r[i], abs(w.phi[i])); //Save wave function multiply by rho axis
		
		
		for (int i=0; i<HH.Nr; i++)
			fprintf(out2,"%10.17e %10.17e %10.17e %10.17e %10.17e %10.17e \n", wHank.r[i],real(wHank.phiHank[i]),imag(wHank.phiHank[i]), w.r[i], real(w.phi[i]), imag(w.phi[i]) ); //Save wave function multiply by rho axis
		
		printf("ENorm of the Hankel %e\n",1.-wHank.norm());
		printf("ENorm of the Uniform wave %e\n",1-w.norm());
		
		
		//printf("%e\n",1.-w.norm());
	}
	
	
	/*
	 
	 interpU2H(w, wHank);
	 
	 
	 for (int i=0; i<HH.Nr; i++)
	 fprintf(out1,"%10.17e %10.17e %10.17e %10.17e \n", wHank.r[i],abs(wHank.phiHank[i]), w.r[i], abs(w.phi[i])); //Save wave function multiply by rho axis
	 
	
	for (int i=0; i<HH.Nr; i++)
		fprintf(out2,"%10.17e %10.17e %10.17e %10.17e %10.17e %10.17e \n", wHank.r[i],real(wHank.phiHank[i]),imag(wHank.phiHank[i]), w.r[i], real(w.phi[i]), imag(w.phi[i]) ); //Save wave function multiply by rho axis

	*/
	printf("ENorm of the Uniform wave %e\n",1-w.norm());
	printf("ENorm of the Hankel %e\n",1.-wHank.norm());
	
	
	
}
