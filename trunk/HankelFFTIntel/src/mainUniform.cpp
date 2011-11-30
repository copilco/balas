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
#include "waveUniform.h"


int main()
{
	
	FILE *out1;
	out1=fopen("out1.txt","w");
	
	// Parameters!!
	
	int Nr=1000;
	int Nt=1;
	
	
	HankelMatrix HH(Nr,200.);
	
	waveUniform w;
	w.initialize(HH);
	
	printf("%d\n",w.Nr);
	
	double r0=100.;
	for(int i=0;i<HH.Nr;i++)
	{
		w.phi[i]=exp(-(w.r[i]-r0)*(w.r[i]-r0)/0.5/0.5);
	}
	
	printf("%e\n",w.norm());
	w.normalize();
	printf("%e\n",w.norm());
	
	
	double dt=0.005;
	w.PrepareCrankArrays(dt);
	for (int ktime=0; ktime<10000; ktime++)
	{
		w.KineticPropCrankUniform(dt);
		
		if((ktime%100)==0)
		{
			for (int i=0; i<HH.Nr; i++)
				fprintf(out1,"%10.17e \n", w.r[i]*abs(w.phi[i]) ); //Save wave function multiply by rho axis
		}
		printf("%e\n",1.-w.norm());
	}
	
	
	
}
