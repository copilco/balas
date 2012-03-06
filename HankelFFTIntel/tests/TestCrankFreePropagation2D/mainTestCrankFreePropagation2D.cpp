//
//  mainCrankUniformProp.cpp
//  
//
//  Created by camilo on 21/12/11.
//
//
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include "arrai.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "waveUniform2D.h"
#include "tools.h"


int main()
{

	cout << "\n\n////////////////////////////////////////////////////" << endl;
    cout << "////////////////////////////////////////////////////" << endl;
    cout << "mainTestCrankFreePropagation2D. Running example..." << endl;
    cout << "///////////////////////////////////////////////////" << endl;
	cout << "///////////////////////////////////////////////////" << endl;
	
	
	fstream axis("axis.txt",ios::out);
	fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
	
	
	//////////////////
    //  Parameters  //
    //////////////////
	
	int Nr=520;
	int Nz=680;
	
	double dz=0.3;
	double dr=0.3;
	complex dt=complex(0.01,0.);
	
	
	int Ntime=100;
	int snap=10;
	
	//Gaussian parameters
    
	double Rmax  = Nr*dr;
    double rho0  = 0.;//Rmax/2.;
	double rho00 = 0.;
	double z0    = 0.;		
	double v0r   = 0.;//5.;
	double v0z   = 0.0;	
	double sigma = 5.;
	
	// Print out the information on the screen
	
    cout << "\nNr: " << Nr << endl;
    cout << "Nz: " << Nz << endl;
    cout << "dz: " << dz << endl;
    cout << "dt: " << dt << endl;
    cout << "Ntime: " << Ntime << endl;
    cout << "Snapshots: " << snap << endl;
    cout << "Rmax: " << Rmax << endl;
    cout << "rho0: " << rho0 << endl;
    cout << "rho00: " << rho00 << endl;
    cout << "z0: " << z0 << endl;
    cout << "v0r: " << v0r << endl;
    cout << "v0z: " << v0z << endl;
    cout << "Sigma: " << sigma << endl;

	
	////////////////////////
    //  Declare objects  //
    ///////////////////////
	
	
	HankelMatrix HH(Nr,Rmax);
	
	waveUniform2D w;
	w.initialize(HH,Nz,dz);

	
	//Initial test function	
	
	for(int j=0;j<HH.Nr;j++)
		for(int i=0;i<Nz;i++)
		{
			w.phi[w.index(j,i)]= exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma )*exp(I* (v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0)) );
		}
	
	
	// Check the norm
	
	cout << "\n\nOriginal norm: " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm: " << w.norm() << endl << endl;
	
	
	//////////////////////////////////////////
    //  Save initial wavefunction and axis  //
    //////////////////////////////////////////
    
	for(int j=0;j<Nr;j++)
		axis << w.r[j] << endl;
	
	for(int j=0;j<Nz;j++)
		axis << w.z[j] << endl;
	
    for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out0 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;
	
	
	////////////////////////////
    //  Prepare Crank Arrays  //
    ////////////////////////////

	
	w.PrepareCrankArrays(dt);
	
	
	////////////////////////////
    //  Start temporal loop   //
    ////////////////////////////
	
	// This variables is for measure the time elapsed in propagation
	
	time_t start, end;
	double diff;
	
	time(&start);
	
	for (int ktime=0; ktime<Ntime; ktime++)
	{
		
		cout << "Loop number: " << ktime << endl;
		
		
		///////////////
        //  Evolve   //
        ///////////////
		
		
		w.Zprop(dt/2.);
		w.Rprop(  dt );
		w.Zprop(dt/2.);
		
		
		//cout << "Error norm: " << 1.-w.norm() << endl << endl;
		
		
		if(ktime%(Ntime/(snap-1)) == 0)
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;	
		
		
		/*
		 if((ktime%(Ntime/snap))==0)
		 {				
		 for(int j=0;j<HH.Nr;j++)
		 for(int i=0;i<Nz;i++)
		 fprintf(out2,"%e \n",abs(w.phi[w.index(j,i)])*abs(w.phi[w.index(j,i)])*w.r[j]);     //Save wave function multiply by rho axis
		 
		 
		 double rexp=0.;
		 for(int j=0;j<HH.Nr;j++)
		 for(int i=0;i<Nz;i++)
		 rexp+= abs(w.phi[w.index(j,i)])*abs(w.phi[w.index(j,i)])*w.r[j]*w.r[j]*w.dz*w.dr;
		 
		 
		 
		 double zexp=0.;
		 for(int j=0;j<HH.Nr;j++)
		 for(int i=0;i<Nz;i++)
		 zexp+= abs(w.phi[w.index(j,i)])*abs(w.phi[w.index(j,i)])*w.z[i]*w.r[j]*w.dz*w.dr;
		 
		 fprintf(out2_,"%e %e %e\n",1.-w.norm(),rexp,zexp);
		 printf("%d %e\n",ktime,1.-w.norm()); 
		 }
		 */
		
		
		/////////////////////////
        //  Apply the absorber //
        /////////////////////////
		
		
		//w.absorber(0.1,0.,0.1,0.,1./6.);
		
		
	}
	
	time(&end);
	
	diff = difftime(end,start);
	
	cout << "Time took by propagation: " << diff << endl;
	
	axis.close();
	out0.close();
	out1.close();
	
}
