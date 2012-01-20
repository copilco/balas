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
	
    
    fstream out0("out0.txt",ios::out);
    //fstream out1("outH2U.txt",ios::out);
    //fstream out2("outU2H.txt",ios::out);
    fstream out3("outErrorH2U.txt",ios::out);
    fstream out4("outErrorU2H.txt",ios::out);
    
    
    //////////////////
    //  Parameters  //
    //////////////////
    
    int Nr=1000;//520;
    int Nz=1000;//680;
    
    double dr=0.05;
    double dz=0.1;
    double dt=0.005;
    
    int Ntime=200;
    int snap=30;
    
    //Gaussian parameters
    
	double Rmax  = ceil(Nr*dz);
    double rho0  = Rmax/2.;
	double rho00 = 12.;
	double z0    = 0.;		
	double v0r   = 0.;//5.;
	double v0z   = 0.0;	
	double sigma = 25.;
    
    // Print out the information on the screen
    
    cout << "\nNr: " << Nr << endl;
    cout << "Nz: " << Nz << endl;
    cout << "dr: " << dr << endl;
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
    
    waveH2D wHinterp;
	wHInterp.initialize(HH,Nz,dz);
    
    waveUniform2D wU;
	wU.initialize(HH,Nz,dz);
    
    waveUniform2D wUInterp;
	wUInterp.initialize(HH,Nz,dz);
    
    
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			//w.phi[w.index(j,i)]=exp(-(w.r[j]-rho0)*(w.r[j]-rho0)/sigma/sigma-(w.z[i]*w.z[i])/sigma/sigma);
			wH.phi[w.index(j,i)]= wH.r[j]*exp(  -(wH.r[j] - rho0)*(wH.r[j] - rho0)/sigma -(wH.z[i] - z0)*(wH.z[i] - z0)/sigma )*exp(I* v0r*(wH.r[j] - rho0) + v0z*(wH.z[i] - z0)  );
			//if (w.r[j]>=rho00*0.)
            //v[w.index(j,i)]=-100.0/sqrt(1.+(w.r[j]-rho00)*(w.r[j]-rho00)+(w.z[i] - 2.*z0)*(w.z[i] - 2.*z0) )*1.
        }
    
    
    // Check the norm
    
    cout << "\n\nOriginal norm: " << wH.norm() << endl;
    wH.normalize();
    cout << "Initial norm: " << wH.norm() << endl;
    
    
    /////////////////////////////////
    //  Save initial wavefunction  //
    /////////////////////////////////
    
    for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out0 << abs(conj(wH.phi[wH.index(j,i)])*wH.phi[wH.index(j,i)])*wH.r[j] << endl;
    
    
	// Copy the wavefunction to the others objects
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
		{
			wU.phi[wU.index(j,i)]=wH.phi[wH.index(j,i)];
		}
	
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
        
        
		wH.FFTFor();
        wH.phi2F(HH);
        wH.HankelTransform(HH);
		
		fase=0.;
		
        for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=dospi*dospi*HH.v[j]*HH.v[j]*dt/2.+wH.q[i]*wH.q[i]*dt/2.;
				wH.G[wH.index(j,i)]*=exp(I*fase);
			}
        
		
        wH.HankelTransformBack(HH);
        wH.F2phi(HH);
		wH.FFTBack();

        
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
        
		
        cout << wHInterp.norm() << endl; 
        cout << wUInterp.norm() << endl;
        
        
        ///////////////////////////////////
        //  Write the error in the norm  //
        ///////////////////////////////////
		
        for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
            {
				//out1 << ktime << "  " << abs(conj(wH.phi[wH.index(j,i)])*wH.phi[wH.index(j,i)])*wH.r[j] << "  " << abs(conj(wHInterp.phi[wHInterp.index(j,i)])*wHInterp.phi[wHInterp.index(j,i)])*wHInterp.r[j] << endl;
                //out2 << ktime << "  " << abs(conj(wU.phi[wU.index(j,i)])*wU.phi[wU.index(j,i)])*wU.r[j] << "  " << abs(conj(wUInterp.phi[wUInterp.index(j,i)])*wUInterp.phi[wUInterp.index(j,i)])*wUInterp.r[j] << endl;
            }
    
        
		cout << ktime << "  " << 1.-wH.norm() << "  " << 1.-wHInterp.norm() << endl;
        
        out3 << ktime << "  " << 1.-wH.norm() << "  " << 1.-wHInterp.norm() << endl;
        
        cout << ktime << "  " << 1.-wU.norm() << "  " << 1.-wUInterp.norm() << endl;

        out4 << ktime << "  " << 1.-wU.norm() << "  " << 1.-wUInterp.norm() << endl;
        
    }
    
    
    out0.close();
    //out1.close();
    //out2.close();
    out3.close();
    out4.close();
    
}


