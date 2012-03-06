//
//  mainTestUnitaryHankelProp2D.cpp
//  
//
//  Created by de la Calle Negro Alejandro on 18/12/11.
//  
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include "arrai.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "waveH2D.h"
#include "tools.h"

int main()
{
    cout << "\n\n////////////////////////////////////////////////////" << endl;
    cout << "////////////////////////////////////////////////////" << endl;
    cout << "mainTestHankelFreePropagation2D. Running example..." << endl;
    cout << "////////////////////////////////////////////////////" << endl;
	cout << "////////////////////////////////////////////////////" << endl;
	
    
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
    double dt=0.1;
    
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
	
	double fase;
    
    HankelMatrix HH(Nr,Rmax);

    waveH2D w;
	w.initialize(HH,Nz,dz);
    
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			w.phi[w.index(j,i)]= exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma )*exp(I* (v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0))  );
			
			//if (w.r[j]>=rho00*0.)
				//v[w.index(j,i)]=-100.0/sqrt(1.+(w.r[j]-rho00)*(w.r[j]-rho00)+(w.z[i] - 2.*z0)*(w.z[i] - 2.*z0) )*1.;
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
    //  Start temporal loop  //
    ///////////////////////////
    
    for (int ktime=0; ktime<Ntime; ktime++)
	{
       
        cout << "Loop number: " << ktime << endl;
        
        
        ///////////////
        //  Evolve   //
        ///////////////
        
        
		w.FFTFor();
        w.phi2F(HH);
        w.HankelTransform(HH);
		
		
		fase=0.;
		
        for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=dospi*dospi*HH.v[j]*HH.v[j]*dt/2.+w.q[i]*w.q[i]*dt/2.;
				w.G[w.index(j,i)]*=exp(-I*fase);
			}
        
        
        w.HankelTransformBack(HH);
        w.F2phi(HH);
		w.FFTBack();
        
        
        ////////////////////////////////////////////
        // Take the snapshot of the wavefunction  //
        ////////////////////////////////////////////
		
		cout << "Error norm: " << 1.-w.norm() << endl << endl;
        
		if(ktime%(Ntime/(snap-1)) == 0)
			
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;	
        
		
		
		/////////////////////////
        //  Apply the absorber //
        /////////////////////////
		
		
		//w.absorber(0.1,0.1,0.1,0.1,1./6.);
	
	
    }

    
	axis.close();
    out0.close();
    out1.close();
    
}

