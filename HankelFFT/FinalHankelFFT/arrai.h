/*
 *  arrai.h
 *  
 *
 *  Created by camilo Ruiz Méndez on 26/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include "fftw3.h"

class arrai
{
public:
	
	int Nt,Nr;
	fftw_complex *v;
//	fftw_plan *p;
	
	
	arrai();
	arrai(int _Nr,int _Nt, double filler);
	void init(int _Nr, int _Nt, double filler);
};

arrai::arrai(int _Nr, int _Nt, double filler)
{
	Nt=_Nt;
	Nr=_Nr;
	v = (fftw_complex*) fftw_malloc(Nt*Nr * sizeof(fftw_complex));
	
	for(int i=0;i<Nt*Nr;i++)
	{
		v[i][0]=filler;
		v[i][1]=filler;
	}
}

void arrai::init(int _Nr, int _Nt, double filler)
{
	Nt=_Nt;
	Nr=_Nr;
	v = (fftw_complex*) fftw_malloc(Nt*Nr * sizeof(fftw_complex));
	
	for(int i=0;i<Nt*Nr;i++)
	{
		v[i][0]=filler;
		v[i][1]=filler;
	}
	
}

arrai::arrai()
{
}


/*
arrai::arrai(int _N, double filler)
{
	N=_N;
	v = (double *)calloc( N, sizeof( double ) );
	for(int i=0;i<N;i++)
		v[i]=filler;
}
*/