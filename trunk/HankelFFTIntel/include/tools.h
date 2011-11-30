/*
 *  tools.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 30/11/11.
 *  Copyright 2011 USAL. All rights reserved.
 *
 */

#include <stdio.h>

void nrerror (char error_text[])
{
	fprintf(stderr, "numerical recipes runtime error..\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr,"... now exiting system..\n");
}



void tridag( complex a[], complex b[],complex c[], complex r[], complex u[], int n,complex gam1[])
{
	int j;
	
	complex bet=b[0];//, *gam;
	
	
	if(b[0]==0.0) nrerror ("error 1 en tridag");
	u[0]=r[0]/bet;
	
	for(j=1;j<n;j++)
	{
		gam1[j] = c[j-1]/bet;
		bet     = b[j] - a[j]*gam1[j];
		
		if (bet==0.0) nrerror("error 2 en tridag");
		u[j]  = (r[j]-a[j]*u[j-1])/bet;
	}
	
	for(j=(n-2);j>=0;j--)
		u[j]-=gam1[j+1]*u[j+1];
	//free(gam);
}



