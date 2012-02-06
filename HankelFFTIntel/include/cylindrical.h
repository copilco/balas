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
	
	
	for(int i=0; i< wlarge.Nr*wlarge.Nz; i++)
		wlarge.phi[i]=complex(0.,0.);
	
	for(int j=0; j<wsmall.Nr; j++)
		for(int i=0; i< wsmall.Nz; i++)
			wlarge.phi[wlarge.index(j,iplacer+i)] = wsmall.phi[wsmall.index(j,i)];
	
}//End reemplace wavefunction




//**********************************************************
//                 Projection function
//**********************************************************
complex projection(waveUniform2D w0, waveUniform2D w )
{
	
	complex proj=complex(0.,0.);
	
	for (int j=0; j<w.Nr; j++)
		for (int i=0; i<w.Nz; i++)
			proj+= w.r[j]*conj(w0.phi[w0.index(j,i)])*w.phi[w.index(j,i)]*w.dr*w.dz;
	
	return proj;
}//End project function




//*********************************************************************
//        Project out removes wavefunction w0 of wavefunction w
//*********************************************************************
void project_out(waveUniform2D w0, waveUniform2D &w )
{
	
	complex proj= projection(w0, w);	
	
	for (int i=0; i<w.Nr*w.Nz; i++)
		w.phi[i]-= proj*w0.phi[i];
	
}//End of project out




#endif
