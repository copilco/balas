#include <iostream.h>
#include <stdio.h>
#include <cstdlib>
#include "clapack.h"
#include "math.h"
#include <complex.h>
#define complex complex<double>

using namespace std;
int main (int argc, const char * argv[]) {
  
  FILE *out0;
  out0 =fopen("movie.txt","w");
  
  
  __CLPK_integer n=2000;	
  __CLPK_integer n1=1;
  
  __CLPK_doublecomplex *lower;
  lower = (__CLPK_doublecomplex*) malloc( (n-1)*sizeof(__CLPK_doublecomplex) );
  
  __CLPK_doublecomplex *middle;
  middle = (__CLPK_doublecomplex*) malloc( n*sizeof(__CLPK_doublecomplex) );
  
  __CLPK_doublecomplex *upp;
  upp = (__CLPK_doublecomplex*) malloc( (n-1)*sizeof(__CLPK_doublecomplex) );
  
  __CLPK_doublecomplex *b;
  b = (__CLPK_doublecomplex*) malloc( n*sizeof(__CLPK_doublecomplex) );
  
  __CLPK_integer info=0;
  
  __CLPK_doublecomplex *wave;
  wave = (__CLPK_doublecomplex*) malloc( n*sizeof(__CLPK_doublecomplex) );
  
  double dx=0.013;
  double dt=0.005;
  
  for (int i=0;i<n;i++)
    {
      double x=(-n/2+i)*dx;
      wave[i].r=exp(-x*x);//*cos(x);
      wave[i].i=0.;//exp(-x*x)*sin(x);
    }
	
	double norm=0.;
	for (int i=0;i<n;i++)
	  norm+=dx*(wave[i].r*wave[i].r+wave[i].i*wave[i].i);

	for (int i=0;i<n;i++)
	{
		wave[i].r/=sqrt(norm);
		wave[i].i/=sqrt(norm);
	}
	
	 norm=0.;
	for (int i=0;i<n;i++)
	   norm+=dx*(wave[i].r*wave[i].r+wave[i].i*wave[i].i);
	printf("NNorm=%e\n",norm);
	
	
	
	
	
	complex wav_m1;	
	complex wav;
	complex wav_p1;
	complex rightside;
	
	for(int h=1;h<5000;h++)
	{
		
		for(int i=0;i<n-1;i++)
		{
			lower[i].r=0.;
			lower[i].i=-dt/2./dx/dx;
			
			upp[i].r=0.;
			upp[i].i=-dt/2./dx/dx;
		}
		
		for(int i=0;i<n;i++)
		{
		middle[i].r=1.;
		middle[i].i=dt/dx/dx;
		}
		

		
		complex lo  = complex( 0., -dt/2./dx/dx);
		complex mi  = complex( 1.,  dt/dx/dx);
		complex up  = complex( 0., -dt/2./dx/dx);
		
		/*
		wav    = complex(  wave[0].r,  wave[0].i);
		wav_p1 = complex(  wave[1].r,  wave[1].i);
		
		rightside=-(up)*wav_p1+conj(mi)*wav;
		b[0].r=real(rightside);
		b[0].i=imag(rightside);
		*/
		
		for(int i=0;i<n;i++)
		{
			wav_m1 = complex(  wave[i-1].r,  wave[i-1].i);	
			wav    = complex(  wave[i].r,  wave[i].i);
			wav_p1 = complex(  wave[i+1].r,  wave[i+1].i);
						
			if(i==0)
				rightside=-(up)*wav_p1+conj(mi)*wav;
			if(i==(n-1))
				rightside=conj(mi)*wav-(lo)*wav_m1;
			else
				rightside=-(up)*wav_p1+conj(mi)*wav-(lo)*wav_m1;
			
			b[i].r=real(rightside);
			b[i].i=imag(rightside);
		}
		
		/*
		wav    = complex(  wave[n-1].r,  wave[n-1].i);
		wav_m1 = complex(  wave[n-2].r,  wave[n-2].i);	
		
		rightside=conj(mi)*wav-(lo)*wav_m1;
		b[n-1].r=real(rightside);
		b[n-1].i=imag(rightside);
		 */
		
		int n2=1;
		zgtsv_(&n,&n1,lower, middle, upp,b,&n, &info);

		for (int i=0;i<n;i++)
		{
		  wave[i].r=b[i].r;
		  wave[i].i=b[i].i;
		}
		
		std::cout << "info = " << info << "\n";
		double norm=0.;
		for (int i=0;i<n;i++)
		  norm+=dx*(wave[i].r*wave[i].r+wave[i].i*wave[i].i);
		printf("norm=%e\n",norm);
		
		int skiper=10;
		if((h%100)==0)
		for (int i=0;i<n/skiper;i++)
		  fprintf(out0,"%e\n",dx*(wave[i*skiper].r*wave[i*skiper].r+wave[i*skiper].i*wave[i*skiper].i));
		
	}
	
		std::cout << "info = " << info << "\n";
}
