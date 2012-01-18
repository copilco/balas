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
#include "waveUniform2D.h"
#include "tools.h"


int main()
{
	
	FILE *out2,*out0;

	out2 = fopen("out2.txt","w");
	out0 = fopen("out0.txt","w");
	// Parameters!!
	
	int Nr=100;
	int Nz=200;
	
	HankelMatrix HH(Nr,20.);
	
	waveUniform2D w;
	w.initialize(HH,Nz,0.1);
	printf("\n Nr=%d\n",w.Nr);
	printf("   Nz=%d\n",w.Nz);
	
	for(int j=0; j<w.Nr; j++)
		fprintf(out0,"%e \n",w.r[j]);
	
	for (int i=0; i<w.Nz; i++)
		fprintf(out0,"%e \n",w.z[i]);		
	
	double r0=10.;
	double sigmar=.7;
	for(int j=0;j<HH.Nr;j++)
		for(int i=0;i<Nz;i++)
		{
			w.phi[w.index(j,i)]=exp(-(w.r[j]-r0)*(w.r[j]-r0)/sigmar/sigmar-(w.z[i]*w.z[i]));
			//w.v[w.index(j,i)]=-exp(-(w.r[j]-r0)*(w.r[j]-r0)/sigmar/sigmar-(w.z[i]*w.z[i]));
		}

	printf("%e\n",w.norm());
	w.normalize();
	printf("%e\n",w.norm());
	
	/*
	for(int j=0;j<w.Nr;j++)
		for(int i=0;i<w.Nz;i++)
		{
			fprintf(out1,"%10.17e \n", abs(w.phi[w.index(j,i)]) ); //Save wave function multiply by rho axis
		}
*/	
	

	
	
	
	complex dt=complex(0.01,0.);
	//double dt=0.01;
	w.PrepareCrankArrays(dt);
	
	complex a2;
	
	for (int ktime=0; ktime<1200; ktime++)
	{
		w.Zprop1(dt);
		
		if((ktime%100)==0)
		{
			
			for(int j=0;j<HH.Nr*Nz;j++)
				fprintf(out2,"%e \n", abs(w.phi[j])); //Save wave function multiply by rho axis
			
		}
		 
		//fprintf(out1,"%e \n", real(w.phi[1]));
		printf("%d %e\n",ktime,1.-w.norm());
	}
	 

	for(int l=0;l<Nr*Nz;l++)
		fprintf(out2,"%e\n",abs(w.phi[l]) );
	
	
	
	
	
}
