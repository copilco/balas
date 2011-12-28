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


int main()
{
	
	FILE *out1;
	out1=fopen("out1.txt","w");
	
	// Parameters!!
	
	int Nr=1000;
	int Nt=1;
	
	
	HankelMatrix HH(Nr,200.);
	
	wave wHank;
	wHank.initialize(HH);
	
	printf("%d\n",wHank.Nr);
	
	double r0=100.;
	for(int i=0;i<HH.Nr;i++)
	{
		wHank.phiHank[i]=exp(-(HH.r[i]-r0)*(HH.r[i]-r0)/0.5/0.5);
	}
	
	printf("%e\n",wHank.norm());
	wHank.normalize();
	printf("%e\n",wHank.norm());
	
	
	double dt=0.005;
	wHank.PrepareCrankArrays(dt);
	for (int ktime=0; ktime<1000; ktime++)
	{
		wHank.KineticPropCrankNonUniform(dt);
		
		if((ktime%100)==0)
		{
			for (int i=0; i<HH.Nr; i++)
				fprintf(out1,"%10.17e \n", wHank.r[i]*abs(wHank.phiHank[i]) ); //Save wave function multiply by rho axis
		}
		printf("%e\n",1.-wHank.norm());
	}
	
	
	
}
