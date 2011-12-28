/*
 *  testingCranck2D.cpp
 *  
 *
 *  Created by alexis on 08/11/2011
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream.h>
#include <complex.h>
#include <vector>
#include <math.h>
#include "cylindrical.h"
#define complex complex<double>


using namespace std;

#define nr 520
#define nz 680

#define dr 0.05
#define dz 0.05

#define dt 0.01

int Ntime = 200;
int snap  = 30;


int main()
{	
	
	FILE *out0;
	FILE *out1;	
	FILE *out2;		
	FILE *out3;	
	FILE *out4;	
	
	
	out0 = fopen("out0.txt","w");
	out1 = fopen("out1.txt","w");	
	out2 = fopen("out2.txt","w");	
	out3 = fopen("out3.txt","w");	
	out4 = fopen("out4.txt","w");	
	
	vector<double>  rho (nr, 0.);       //Rho axis
	vector<double>  z   (nz, 0.);	    //z axis
	vector<double>  v   (nr*nz, 0.);    //energy potential
	
	
	vector<complex> psi0   ( nr*nz, 0. );  //initial wavefunction
	vector<complex> psi1   ( nr*nz, 0. );  //propagating array
	
	vector<complex> phi_r  ( nr, 0. );	   //auxiliar rho
	vector<complex> phi_z  ( nz, 0. );	   //auxiliar z
	
	vector<complex> rv_r   ( nr, 0. );     //auxiliar right rho
	vector<complex> rv_z   ( nz, 0. );	   //auxiliar right z
	
	
	//Diagonal for rho and z
	complex ar;
	vector<complex> br (nr, 0.);
	complex cr;
	
	complex az;
	vector<complex> bz (nz, 0.);
	complex cz;	
	
	
	
	
	//Definition grid
	rho[0] = dr/2.;	
	fprintf(out0,"%e \n",rho[0]);	
	for(int i=1; i<nr; i++){
		rho[i]  =  rho[i-1] + dr;
		fprintf(out0,"%e \n",rho[i]);
	}
	
	for (int k=0; k<nz; k++){
		z[k]    =  (-(nz-1)/2. + k)*dz;
		fprintf(out0,"%e \n",z[k]);		
	}
	
	
	
	//Parameters 
	double rho0  = 10.;
	double rho00 = 12.;
	double z0    = 0.;		
	double v0r   = 5.;
	double v0z   = 0.0;	
	double sigma = 0.1;
	
	
	
	//Initial test function	
	for(int i=0; i<nr; i++)
		for (int j=0; j<nz; j++){
			psi0[index(i,j,nz)]	= rho[i]*exp(  -(rho[i] - rho0)*(rho[i] - rho0)/sigma -(z[j] - z0)*(z[j] - z0)/sigma )*
								  complex( cos( v0r*(rho[i] - rho0) + v0z*(z[j] - z0) ),
										   sin( v0r*(rho[i] - rho0) + v0z*(z[j] - z0) ) ); 
			if (rho[i]>=rho00*0.)
				v[index(i,j,nz)]=-100.0/sqrt(1.+(rho[i]-rho00)*(rho[i]-rho00)+(z[j] - 2.*z0)*(z[j] - 2.*z0) )*1.;
		}
	//End test function
	
	cout<<"\nNorm: "<< cnorm( psi0, rho, nr, nz, dr ,dz ) <<endl;
	
	
	cnormalize(psi0, rho, nr, nz, dr, dz);
	printf("Norm: %e    ENorm: %10.4e\n",cnorm( psi0, rho, nr, nz, dr, dz ),abs(cnorm( psi0, rho, nr, nz, dr, dz )-1.));	
	cout << endl;
	
	
	
	//Saving initial function
	for (int k=0; k < nr; k++)
		for (int k2=0;k2<nz;k2++)
			fprintf(out1,"%e \n", abs(psi0[index(k,k2,nz)])*abs(psi0[index(k,k2,nz)])*rho[k] );
	fprintf(out2,"%e \n",0.);
	fprintf(out4,"%e \n",cnorm( psi0, rho, nr, nz, dr, dz ));
	//End initial function
	
	
	//********************************************************
	//================================================
	//Crank-Nicholson code	
	
	//Time propagation loop
	for (int ktime=0; ktime<Ntime; ktime++) {
		
		Zprop( psi0, v, complex( dt/2.,  0.), nz, nr, dz );
		Rprop( psi0, v, rho, complex( dt, 0.), nz, nr, dr );
		Zprop( psi0, v, complex( dt/2.,  0.), nz, nr, dz );
		
 		if(ktime%(Ntime/(snap-1)) == 0){			
			for (int k=0; k < nr; k++)
				for (int k2=0;k2<nz;k2++)
					fprintf(out1,"%e \n", abs(psi0[index(k,k2,nz)])*abs(psi0[index(k,k2,nz)])*rho[k] );

			fprintf(out2,"%e \n",ktime*dt);
			printf("Norm: %e    ENorm: %10.4e\n",cnorm( psi0, rho, nr, nz, dr, dz ),abs(cnorm( psi0, rho, nr, nz, dr, dz )-1.));
			fprintf(out4,"%e \n",cnorm( psi0, rho, nr, nz, dr, dz ));			
		}//
		//cout<<"\nNorm: "<< rnorm( psi0, nz, nr, dz, dr ) <<endl;		
	}//End time loop//*/


	fprintf(out3,"%d \n %d \n %d \n %d \n %e\n",nr, nz, Ntime, snap+1, dt);
} //
