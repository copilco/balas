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

class arrai
{
public:
	
	int N;
	double *v;
	arrai();
	arrai(int _N, double filler);
	void init(int _N, double filler);
};

arrai::arrai(int _N, double filler)
{
	N=_N;
	v = (double *)calloc( N, sizeof( double ) );
	for(int i=0;i<N;i++)
		v[i]=filler;
}

void arrai::init(int _N, double filler)
{
	N=_N;
	v = (double *)calloc( N, sizeof( double ) );
	for(int i=0;i<N;i++)
		v[i]=filler;
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