/*
 *  arrai.h
 *  
 *
 *  Created by Alejandro de la Calle on 26/09/11.
 *  
 *
 */

#ifndef ARRAI_H
#define ARRAI_H

#include <iostream>
#include <math.h>
#include <new>
#include <vector>
#include <complex>
//#define MKL_Complex16 std::complex<double>
#define complex complex<double>
//#include "mkl.h"
//#include "mkl_dfti.h"
#include "constant.h"

using namespace std;

class arrai
{
	friend class HankelMatrix;
public:
	
	int Nr, Nt;
	complex *v;
	
	
	
	//////////////////////
	// Function members //
	//////////////////////
	
	// Constructor
	
	arrai(int _Nr,int _Nt,double filler=0.0)
	{
		Nr=_Nr;
		Nt=_Nt;
		v= (complex*)malloc(Nr*Nt*sizeof(complex));
		
		
		for(int i=0;i<Nt*Nr;i++)
		{
			real(v[i])=filler;
			imag(v[i])=filler;
		}
		
	}
	
	
	
	// Destructor
	
	~arrai()
	{ 
		free(v);

	}
	
	// Overloaded operators
	/*
	static void* operator new(size_t sz)
	{
		void *p=(void*)mkl_malloc(sz,16);
		return p;
		
	}
	
	static void operator delete(void* ptr)
	{
		mkl_free(ptr);
	}
	*/
	
	
	void init(int _Nr,int _Nt,double filler=0.0)
	{
		Nr=_Nr;
		Nt=_Nt;
		v= (complex*)malloc(Nr*Nt*sizeof(complex));
		
		
		for(int i=0;i<Nt*Nr;i++)
		{
			real(v[i])=filler;
			imag(v[i])=filler;
		}
		
	}
	
	void set(double filler)
	{
		for(int i=0;i<Nt*Nr;i++)
		{
			real(v[i])=filler;
			imag(v[i])=filler;
		}
	
	}
	
	
};

#endif //  ARRAI_H
		