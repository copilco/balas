/*
 *  testingHankelMatCreation.cpp
 *  
 *
 *  Created by Camilo Ruiz Méndez on 01/01/10.
 *  Copyright 2010 USAL. All rights reserved.
 *
 */


#include <iostream>
//#include <gsl/gsl_cblas.h>
#include <Accelerate/Accelerate.h> //contains data types used
#include <math.h>
#include "HankelMatrix.h"
//#include <gsl/gsl_sf_bessel.h>
//#include <complex.h>
//#define complex complex<double>

//using std;
int main()
{
	
	
	//complex I=complex(0.,1.);
	int Nr=1024;//50;
	int Nt=12;//50;
	
	double R=.05;
	
	//Parametros para blas
	int lda=Nr;
	int ldb=1;
	int ldc=1;
	
	HankelMatrix HH(Nr,R);
	
	double dr = R/(Nr-1);	//%	Radial spacing	
	printf("Nr=%d\n",Nr);
	
	
	FILE *in0,*in1;
	in0=fopen("besselCeros.txt","r");
	in1=fopen("InputFunction.txt","r");
	
	FILE *out0,*out1,*out2,*out3,*out4,*out5;
	FILE *sout0;
	
	sout0=fopen("spout0.txt","w");
	out0=fopen("pout0.txt","w");
	out1=fopen("pout1.txt","w");
	out2=fopen("pout2.txt","w");
	out3=fopen("pout3.txt","w");
	out4=fopen("pout4.txt","w");
	out5=fopen("pout5.txt","w");
	
	
	
	//Transformation
	
	//arrai phi;
	arrai f0,f1;
	arrai F0,F1;
	arrai f20,f21;
	arrai F20,F21;
	arrai Fretrieved0,Fretrieved1;
	arrai fretrieved0,fretrieved1; 
	arrai F0mat,F1mat,F20mat,F21mat;
	
	//phi.init(N,0.);
	f0.init(Nr,1,0.);
	f1.init(Nr,1,0.);
	F0.init(Nr,1,0.);
	F1.init(Nr,1,0.);
	f20.init(Nr,1,0.);
	f21.init(Nr,1,0.);
	F20.init(Nr,1,0.);
	F21.init(Nr,1,0.);
	
	F0mat.init(Nr,Nt,0.);
	F1mat.init(Nr,Nt,0.);
	
	F20mat.init(Nr,Nt,0.);
	F21mat.init(Nr,Nt,0.);
	
	fretrieved0.init(Nr,1,0.);
	fretrieved1.init(Nr,1,0.);
	Fretrieved0.init(Nr,1,0.);
	Fretrieved1.init(Nr,1,0.);
	
	
	
	
}