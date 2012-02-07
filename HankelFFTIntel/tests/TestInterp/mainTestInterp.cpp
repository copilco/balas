//
//  mainTestInterp.cpp
//  
//
//  Created by de la Calle Negro Alejandro on 03/01/12.
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
#include "waveUniform2D.h"
#include "tools.h"
#include "interp.h"


int main()
{
    cout << "\n\n/////////////////////////////////////////////" << endl;
    cout << "/////////////////////////////////////////////" << endl;
    cout << "mainTestInterp. Running example..." << endl;
    cout << "/////////////////////////////////////////////" << endl;
	cout << "/////////////////////////////////////////////" << endl;
	
    
	fstream axis("axis.txt",ios::out);
    fstream out0("out0.txt",ios::out);
    //fstream out1("outH2U.txt",ios::out);
    //fstream out2("outU2H.txt",ios::out);
    fstream out3("outErrorH2U.txt",ios::out);
    fstream out4("outErrorU2H.txt",ios::out);
	//fstream out5("outDiffHankel.txt",ios::out);
	//fstream out6("outDiffCrank.txt",ios::out);
    
    
    //////////////////
    //  Parameters  //
    //////////////////
    
    int Nr=250;
    int Nz=500;
    
    double dz=0.3;
	double absdt = 0.1;
    complex dt = complex(0.1,0.);
    
    int Ntime=5000;
    int snap=20;
    
    //Gaussian parameters
    
	double Rmax  = 75.;
    double rho0  = 0.;
	double rho00 = 0.;
	double z0    = 0.;		
	double v0r   = 0.;
	double v0z   = 0.;	
	double sigma = 20.;
    
    // Print out the information on the screen
    
    cout << "\nNr: " << Nr << endl;
    cout << "Nz: " << Nz << endl;
    cout << "dr: " << Rmax/double(Nr) << endl;
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
    
    waveH2D wH;
	wH.initialize(HH,Nz,dz);
    
    waveH2D wHInterp;
	wHInterp.initialize(HH,Nz,dz);
    
    waveUniform2D wU;
	wU.initialize(HH,Nz,dz);
    
    waveUniform2D wUInterp;
	wUInterp.initialize(HH,Nz,dz);
    
    
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			wH.phi[wH.index(j,i)]= exp(  -(wH.r[j] - rho0)*(wH.r[j] - rho0)/sigma/sigma -(wH.z[i] - z0)*(wH.z[i] - z0)/sigma/sigma )*exp(I* v0r*(wH.r[j] - rho0) + v0z*(wH.z[i] - z0)  );
			wU.phi[wU.index(j,i)]= exp(  -(wU.r[j] - rho0)*(wU.r[j] - rho0)/sigma/sigma -(wU.z[i] - z0)*(wU.z[i] - z0)/sigma/sigma )*exp(I* v0r*(wU.r[j] - rho0) + v0z*(wU.z[i] - z0)  );

        }
    
    
    // Check the norm
    
    cout << "\n\nOriginal norm in Hankel: " << wH.norm() << endl;
    wH.normalize();
    cout << "Initial norm in Hankel: " << wH.norm() << endl;
    
	cout << "\n\nOriginal norm in Uniform: " << wU.norm() << endl;
    wU.normalize();
    cout << "Initial norm in Uniform: " << wU.norm() << endl;
    
    /////////////////////////////////
    //  Save initial wavefunction  //
    /////////////////////////////////
    
	for(int j=0;j<Nr;j++)
		axis << wU.r[j] << endl;
    
    for(int j=0;j<Nz;j++)
		axis << wU.z[j] << endl;
	
    for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out0 << abs(conj(wH.phi[wH.index(j,i)])*wH.phi[wH.index(j,i)])*wH.r[j] << "  " << abs(conj(wU.phi[wU.index(j,i)])*wU.phi[wU.index(j,i)])*wU.r[j] << endl;

	
    ////////////////////////////
    //  Prepare Crank arrays  //
    ////////////////////////////
	
    wU.PrepareCrankArrays(dt);
    
    double fase;
    
	    
    ////////////////////////////
    //  Start temporal loop  //
    ///////////////////////////
    
    for (int ktime=0; ktime<Ntime; ktime++)
	{
        
        cout << "Loop number: " << ktime << endl;
        
                
        //////////////////////////////////
        //  Evolve in Hankel Transform  //
        //////////////////////////////////
		
		 wH.prop_kinetic(HH,dt);
		 
        ////////////////////////////////
        //  Evolve in Crank-Nicolson  //
        ////////////////////////////////
        
        
        wU.Zprop1(dt/2.);
		wU.Rprop(dt);
		wU.Zprop1(dt/2.);
        
        
        /////////////////////////////////////////////////////
        //  Interpolate the wavefuntion (Hankel to Uniform) //
        //////////////////////////////////////////////////////        
        
		
		interpH2U(wH, wUInterp);
        interpU2H(wU, wHInterp);
        
		
        //cout << wHInterp.norm() << endl; 
        //cout << wUInterp.norm() << endl;
        
        
        ///////////////////////////////////
        //  Write the error in the norm  //
        ///////////////////////////////////
		/*
        for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
            {
				//out1 << ktime << "  " << abs(conj(wH.phi[wH.index(j,i)])*wH.phi[wH.index(j,i)])*wH.r[j] << "  " << abs(conj(wUInterp.phi[wUInterp.index(j,i)])*wUInterp.phi[wUInterp.index(j,i)])*wUInterp.r[j] << endl;
                //out2 << ktime << "  " << abs(conj(wU.phi[wU.index(j,i)])*wU.phi[wU.index(j,i)])*wU.r[j] << "  " << abs(conj(wHInterp.phi[wHInterp.index(j,i)])*wHInterp.phi[wHInterp.index(j,i)])*wHInterp.r[j] << endl;
				out5 << ( abs(conj(wH.phi[wH.index(j,i)])*wH.phi[wH.index(j,i)])*wH.r[j] ) - ( abs(conj(wHInterp.phi[wHInterp.index(j,i)])*wHInterp.phi[wHInterp.index(j,i)])*wHInterp.r[j] ) << endl;
				out6 << ( abs(conj(wU.phi[wU.index(j,i)])*wU.phi[wU.index(j,i)])*wU.r[j] ) - ( abs(conj(wUInterp.phi[wUInterp.index(j,i)])*wUInterp.phi[wUInterp.index(j,i)])*wUInterp.r[j] ) << endl;
			}
		*/
        
		cout << 1.-wH.norm() << "  " << 1.-wUInterp.norm() << endl;
        
        out3 << ktime << "  " << 1.-wH.norm() << "  " << 1.-wUInterp.norm() << endl;
        
        cout << 1.-wU.norm() << "  " << 1.-wHInterp.norm() << endl;

        out4 << ktime << "  " << 1.-wU.norm() << "  " << 1.-wHInterp.norm() << endl;
        
    }
    
    
	axis.close();
    out0.close();
    //out1.close();
    //out2.close();
    out3.close();
    out4.close();
	//out5.close();
	//out6.close();
    
}


