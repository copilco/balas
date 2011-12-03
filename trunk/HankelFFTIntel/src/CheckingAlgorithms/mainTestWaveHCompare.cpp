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
//#define MKL_Complex16 std::complex<double>
//#include "mkl.h"
//#include "mkl_dfti.h"
#include "arrai.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "wave.h"
#include "interp.h"

int main()
{
	FILE *out1,*out1U,*ejes;
	out1=fopen("out1.txt","w+");
	out1U=fopen("out1U.txt","w+");
	ejes=fopen("ejes.txt","w+");
	
	// Parameters!!
	
	int Nr=1000;
	int Nt=1;
	
	
	HankelMatrix HH(Nr,200.);
	
	/**********************/
	//  Prepare Hankel
	/**********************/
	wave wHank;
	wHank.initialize(HH);
	
	/**********************/
	//  Prepare the  Uniform
	/**********************/

	waveUniform w;
	w.initialize(HH);

	
	cout << "mainTestWaveH. Running example...\n" ;
	printf("Number of points: %d\n",wHank.Nr);
	
	
	
	/*******************************/
	//  Start with Hankel definition
	/*******************************
	
	double r0=100.;
	double sigma=10.;
	for(int i=0;i<HH.Nr;i++)
	{
		wHank.phiHank[i]=exp(-(HH.r[i]-r0)*(HH.r[i]-r0)/sigma/sigma);
	}
	
	printf("Initial normalization: %e \n",wHank.norm());
	wHank.normalize();
	printf("Normalization: %e\n",wHank.norm());
	
	interpH2U(wHank, w);
	
	printf("%e\n",1.-w.norm());
	w.normalize();
	printf("%e\n",1.-w.norm());
	
	/******************************/
	//Start with Uniform definition
	/******************************/
	
	double r0=100.;
	double sigma=10.;
	for(int i=0;i<HH.Nr;i++)
	{
		w.phi[i]=exp(-(w.r[i]-r0)*(w.r[i]-r0)/sigma/sigma)*sin(w.r[i]*0.5);
		wHank.phiHank[i]=exp(-(wHank.r[i]-r0)*(wHank.r[i]-r0)/sigma/sigma)*sin(wHank.r[i]*0.5);
	}
	
	
	
	/******************************
	//Interp Alejandro
	interpU2H(w,wHank);
	/******************************/

	/******************************
	//Read interp Hankel Matlab
	 
	double *temp=new double[w.Nr];
	
	fstream raxis ("phiHankelmatlab.bin", ios::in | ios::binary);
	raxis.read ((char*)&temp[0], sizeof (double)*(w.Nr) );
	raxis.close();
	
	printf("ss %e,\n",temp[500]);
	
	for(int i=0;i<HH.Nr;i++)
	{
		wHank.phiHank[i]=complex(temp[i],0.);
	}
	/******************************/
	
	
	printf("Initial normalization: %e \n",w.norm());
	w.normalize();
	printf("Normalization: %e\n",w.norm());
	
	printf("%e\n",1.-wHank.norm());
	wHank.normalize();
	printf("%e\n",1.-wHank.norm());
	
	/******************************/
	//Start the loop
	/******************************/
	
	
	
	double dt=0.1;
	double *fase=new double[Nr];	

	/**********************/
	//  Prepare the  Uniform
	/**********************/
	
	w.PrepareCrankArrays(dt);

	/**********************/
	//  Start the time loop
	/**********************/

	
	for (int i=0; i<HH.Nr; i++)
	{
		fprintf(ejes,"%10.17e %10.17e\n", wHank.r[i],w.r[i]); 
	
	}
	
	
	for (int ktime=0; ktime<1000; ktime++)
	{
		
		wHank.phi2F(HH);
		wHank.HankelTransform(HH, Nt);
		wHank.F22f2(HH);
		

		/**********************/
		//  Prepare the  Uniform
		w.KineticPropCrankUniform(dt);
		/**********************/
		
		printf("Nomalization in space (before phase): %e\n",1.-wHank.norm());
		printf("UUUUNomalization in space (before phase): %e\n",1.-w.norm());
		
		for (int i=0;i<Nr;i++)
		{
			fase[i]=dospi*dospi*HH.v[i]*HH.v[i]*dt/2.; ;//HH.v[i]*HH.v[i]/2.*ktime*dt;
			wHank.F2[i]=wHank.F2[i]*exp(-I*fase[i]);
		}
		
		printf("Nomalization in frequency space (after phase): %e\n",1.-wHank.vnorm());
		
		wHank.HankelTransformBack(HH, Nt);
		wHank.F2phi(HH);	
		printf("Final normalization: %e\n",wHank.norm());

		
		
		if((ktime%10)==0)
		{
			for (int i=0; i<HH.Nr; i++)
			{
				//fprintf(out1,"%10.17e\n", wHank.r[i]*abs(wHank.phiHank[i])); 
				//fprintf(out1U,"%10.17e\n", w.r[i]*abs(w.phi[i])); 
				
				fprintf(out1,"%10.17e\n", abs(wHank.phiHank[i])); 
				fprintf(out1U,"%10.17e\n", abs(w.phi[i])); 

			}
		}
		
		//printf("ENorm of the Hankel %e\n",1.-wHank.norm());
	
	}
	
	
	
	wHank.phi2F(HH);
	wHank.HankelTransform(HH,1);
	wHank.F22f2(HH);
	
	printf("Nomalization in frequency space: %e\n",wHank.vnorm());
	
	wHank.HankelTransformBack(HH,1);
	
	
	wHank.F2phi(HH);	
	printf("Final normalization: %e\n",wHank.norm());
	
	
	
	
	fclose(out1);
	
}
