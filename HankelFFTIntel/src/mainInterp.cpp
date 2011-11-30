/*
 *  mainTestWaveH.cpp
 *  
 *
 *  Created by camilo on 29/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
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
	
	waveUniform w;
	w.initialize(HH);
	
	printf("%d\n",w.Nr);
	
	double r0=100.;
	double sigma=10.;
	for(int i=0;i<HH.Nr;i++)
	{
		w.phi[i]=exp(-(w.r[i]-r0)*(w.r[i]-r0)/sigma/sigma);
	}
	
	printf("%e\n",w.norm());
	w.normalize();
	printf("%e\n",w.norm());
	
	
	wave wHank;
	wHank.initialize(HH);
	
	
	double dt=0.005;
	w.PrepareCrankArrays(dt);
	for (int ktime=0; ktime<100; ktime++)
	{
		w.KineticPropCrankUniform(dt);
		
		interpU2H(w, wHank);
		
		
		for (int i=0; i<HH.Nr; i++)
			fprintf(out1,"%10.17e %10.17e %10.17e %10.17e \n", wHank.r[i],abs(wHank.phiHank[i]), w.r[i], abs(w.phi[i])); //Save wave function multiply by rho axis
		
		
		for (int i=0; i<HH.Nr; i++)
			fprintf(out2,"%10.17e %10.17e %10.17e %10.17e %10.17e %10.17e \n", wHank.r[i],real(wHank.phiHank[i]),imag(wHank.phiHank[i]), w.r[i], real(w.phi[i]), imag(w.phi[i]) ); //Save wave function multiply by rho axis
		
		printf("ENorm of the Uniform wave %e\n",1-w.norm());
		printf("ENorm of the Hankel %e\n",1.-wHank.norm());
		
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
