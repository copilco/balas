#ifndef CYLINDRICAL_H
#define CYLINDRICAL_H

#include <iostream>
#include <math.h>
#include <new>
#include "waveUniform2D.h"


using namespace std;


//*****************************************************************************
//      Reemplace wavefunction wsmall in other bigger wavefunction wlarge
//*****************************************************************************
void placeWF( waveUniform2D &wlarge, waveUniform2D &wsmall)
{	
	int iplacer=floor((wlarge.Nz - wsmall.Nz)/2);
	//int jplacer=floor((wlarge.Nr - wsmall.Nr)/2);
	
	//for(int i=0; i<wlarge.Nr*wlarge.Nz; i++)
	//	wlarge.phi[i]=complex(0.,0.);
	memset(wlarge.phi, 0, sizeof(complex)*wlarge.Nr*wlarge.Nz );
	
	for(int j=0; j<wsmall.Nr; j++)
		for(int i=0; i< wsmall.Nz; i++)
			wlarge.phi[wlarge.index(j,iplacer+i)] = wsmall.phi[wsmall.index(j,i)];
	
}//End reemplace wavefunction




//**********************************************************
//                 Projection function
//**********************************************************
/*
complex projection(waveUniform2D &w0, waveUniform2D &w )
{
	
	complex proj=complex(0.,0.);
	complex *p = NULL;	
	int number_threads=10;
	
	#pragma omp parallel
		number_threads	= omp_get_num_threads();
	
	
	p  = (complex *)mkl_malloc( number_threads*sizeof(complex), 16 );
	memset( p, 0, sizeof(complex)*number_threads);	
	
	
#pragma omp parallel
	{
		int id_thread		= omp_get_thread_num();
		int Nr_thread		= ceil(w.Nr / number_threads);
		int Nr_i			= id_thread * Nr_thread;
		int Nr_f			= Nr_i + Nr_thread;
		
		if (id_thread+1 == number_threads) Nr_f = w.Nr;
		
		for (int j=Nr_i; j<Nr_f; j++)
			for (int i=0; i<w.Nz; i++)		
				p[id_thread]+= w.r[j]*conj(w0.phi[w0.index(j,i)])*w.phi[w.index(j,i)]*w.dr*w.dz;
			
	};
	
	for (int i=0;i<number_threads;i++)
		proj+=p[i];
	
	mkl_free(p);
	
	return proj;
	
};//End project function*/




complex projection(waveUniform2D &w0, waveUniform2D &w )
{
	
	complex proj=complex(0.,0.);
			
		for (int j=0; j<w.Nr; j++)
			for (int i=0; i<w.Nz; i++)		
				proj+= w.r[j]*conj(w0.phi[w0.index(j,i)])*w.phi[w.index(j,i)]*w.dr*w.dz;
		
	
	return proj;
	
};//End project function*/





//*********************************************************************
//        Project out: removes wavefunction w0 of wavefunction w
//*********************************************************************
void project_out(waveUniform2D &w0, waveUniform2D &w )
{
	
	complex proj= projection(w0, w);	
	
	for (int i=0; i<w.Nr*w.Nz; i++)
		w.phi[i]-= proj*w0.phi[i];
	
}//End of project out




#endif
