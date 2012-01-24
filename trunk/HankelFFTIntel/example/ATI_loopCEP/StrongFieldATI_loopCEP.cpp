//
//  StrongFieldATI_loopCEP.cpp
//  
//
//  Created by Alejandro de la Calle on 24/01/12.
//  
//
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include "HankelMatrix.h"
#include "constant.h"
#include "waveH2D.h"


int main()
{
	
	cout << "\n\n/////////////////////////////////////////////" << endl;
    cout << "/////////////////////////////////////////////" << endl;
    cout << "StrongFieldATI_loopCEP. Running test...  " << endl;
    cout << "/////////////////////////////////////////////" << endl;
	cout << "/////////////////////////////////////////////" << endl;
	
	fstream axis("axis.txt",ios::out);
	fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
	
	
	//////////////////
    //  Parameters  //
    //////////////////
    
    int Nr=520;
    int Nz=680;
    
    double dz=0.1;
    double absdt = 0.01;
	complex dt;
	complex t;
	
	int charge_nuclei=1.;
	double soft_core=2.;
	
	int Ntime=2500;
    int snap=20;

	//Gaussian parameters
    
	double Rmax  = 100.;//ceil(Nr*dr);
    double rho0  = Rmax/2.;
	double rho00 = 12.;
	double z0    = 0.;		
	double v0r   = 0.;//5.;
	double v0z   = 0.0;	
	double sigma = 30.;
	
	//Laser parameters
	
	double efield_z;
    double avect_z;
	
	double e0=sqrt(3.e14/3.5e16);
	double ww=0.057;
	double period=dospi/ww;
	double cycle_number=4.;
	double cep=0.;
	
	
	////////////////////////
    //  Declare objects  //
    ///////////////////////
    
    HankelMatrix HH(Nr,Rmax);
	
    waveH2D w;
	w.initialize(HH,Nz,dz);
    
	waveH2D w1;
	w1.initialize(HH,Nz,dz);
	
	
	w.rho0=rho0;
	w.z0=z0;
	
    // Print out the information on the screen
	
    cout << "\nNr: " << Nr << endl;
    cout << "Nz: " << Nz << endl;
    cout << "dz: " << dz << endl;
    cout << "dt: " << absdt << endl;
    cout << "Ntime: " << Ntime << endl;
    cout << "Snapshots: " << Ntime/snap << endl;
    cout << "Rmax: " << Rmax << endl;
    cout << "rho0: " << rho0 << endl;
    cout << "rho00: " << rho00 << endl;
    cout << "z0: " << z0 << endl;
    cout << "v0r: " << v0r << endl;
    cout << "v0z: " << v0z << endl;
    cout << "Sigma: " << sigma << endl;

	
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			w.phi[w.index(j,i)]= w.r[j]*exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma )*exp(I*(v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0)) );	
		}
    
	//Set potential
	
	w.set_potencial_hlike1D(charge_nuclei,soft_core);
    
    // Check the norm
    
    cout << "\n\nOriginal norm: " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm: " << w.norm() << endl << endl;
	
	
	//////////////////////////////////////////
    //  Save initial wavefunction and axis  //
    //////////////////////////////////////////
    
	
	 w.saveAxes(axis);
	 w.snapshot(out0,HH);
	
	
    ////////////////////////////
    //  Start temporal loop  //
    ///////////////////////////
	
	
	///////////////
	//  Evolve   //
	///////////////
	
	
	////////////////////////////////////////////
	// Take the snapshot of the wavefunction  //
	////////////////////////////////////////////
	/*
	if((ktime%snap) == 0)
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;	
	*/
	//w.absorber(0.1,0.1,0.1,0.1,1./6.);
	
	
	axis.close();
    out0.close();
    out1.close();
	
	
}


