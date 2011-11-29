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
	
	wHank.phi2F(HH);
	wHank.HankelTransform(HH,1);
	wHank.F22f2(HH);
	
	printf("%e\n",wHank.vnorm());
	
	wHank.HankelTransformBack(HH,1);
	
	
	wHank.F2phi(HH);	
	printf("%e\n",wHank.norm());
	
}
