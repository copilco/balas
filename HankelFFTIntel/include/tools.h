/*
 *  tools.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 30/11/11.
 *  
 *
 */
#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <math.h>
#include <new>

using namespace std;


void nrerror (char error_text[])
{
	fprintf(stderr, "Runtime error..\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr,"... now exiting system..\n");
}



void tridag(complex *a,complex *b,complex *c,complex *r,complex *u, int n,complex *gam1)
{
	
	complex bet=b[0]; // Declare and define auxiliar array
	
	if(b[0]==0.0) nrerror ("error 1 en tridag");
	u[0]=r[0]/bet;
	
	for(int j=1;j<n;j++)
	{
		gam1[j] = c[j-1]/bet;
		bet     = b[j] - a[j]*gam1[j];
		
		if (bet==0.0) nrerror("error 2 en tridag");
		u[j]  = (r[j]-a[j]*u[j-1])/bet;
	}
	
	for(int j=(n-2);j>=0;j--)
		u[j]-=gam1[j+1]*u[j+1];

}






#endif // TOOLS_H