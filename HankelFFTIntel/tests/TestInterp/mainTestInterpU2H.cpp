index(j,i)//
//  mainTestInterpU2H.cpp
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
#include "waveUniform2D.h"
#include "tools.h"
#include "interp.h"


int main()
{
    cout << "\n\n\n/////////////////////////////////////////" << endl;
    cout << "/////////////////////////////////////////" << endl;
    cout << "///mainTestInterpU2H. Running example..." << endl;
    cout << "/////////////////////////////////////////" << endl;
	cout << "/////////////////////////////////////////" << endl;
	
    
	fstream axes("axes.txt",ios::out);
    fstream out0("out0.txt",ios::out);
    fstream out1("outU2H.txt",ios::out);
	fstream out2("outErrorU2H.txt",ios::out);

    
    //////////////////
    //  Parameters  //
    //////////////////
    
    int Nr=520;
    int Nz=680;
    
    double dr=0.1;
    double dz=0.1;
    complex dt=complex(0.01,0.0);
    
    int Ntime=100;
    int snap=20;
    
    //Gaussian parameters
    
	double Rmax  = ceil(Nr*dz);
    double rho0  = Rmax/2.;
	double rho00 = 12.;
	double z0    = 0.;		
	double v0r   = 0.;//5.;
	double v0z   = 0.0;	
	double sigma = 10.;
    
    // Print the information on the screen
        
    cout << "\nNr: " << Nr << endl;
    cout << "Nz: " << Nz << endl;
    cout << "dr: " << dr << endl;
    cout << "dz: " << dz << endl;
    cout << "dt: " << dt.real() << " " << dt.imag() << endl;
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
    
    waveUniform2D w;
	w.initialize(HH,Nz,dz);
    
    
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			w.phi[w.index(j,i)]= w.r[j]*exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma )*
            exp(I* v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0)  );
			//if (wH.r[j]>=rho00*0.)
				//v[wH.index(j,i)]=-100.0/sqrt(1.+(wH.r[j]-rho00)*(wH.r[j]-rho00)+(wH.z[i] - 2.*z0)*(wH.z[i] - 2.*z0) )*1.;
		}
    
    
    // Check the norm
    
    cout << "\n\nOriginal norm: " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm: " << w.norm() << endl << endl;
    
    
    /////////////////////////////////
    //  Save initial wavefunction  //
    /////////////////////////////////
    
    for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out0 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;

	//Save axes
	
	for(int j=0;j<Nr;j++)
		axes << w.r[j] << "  " << wH.r[j] << endl;
	
	
	for(int i=0;i<Nz;i++)
		axes << w.z[i] << "  " << wH.z[i] << endl;

    ////////////////////////////
    //  Prepare Crank arrays  //
    ////////////////////////////
	
    w.PrepareCrankArrays(dt);

    ////////////////////////////
    //  Start temporal loop  //
    ///////////////////////////
    
    for (int ktime=0; ktime<1; ktime++)
	{
       
        cout << "Loop number: " << ktime << endl;
        
        /////////////////////////////////////////////////////
        //  Evolve with Split operator in Hankel and FFT   //
        /////////////////////////////////////////////////////
        
        
        w.Zprop1(dt/2.);
		w.Rprop(dt);
		w.Zprop1(dt/2.);

        
        //////////////////////////////////////////////////////
        //  Interpolate the wavefuntion (Hankel to Uniform) //
        //////////////////////////////////////////////////////        
			
		interpU2H(w, wH);
        
        ///////////////////////////////////
        //  Write the error in the norm  //
        ///////////////////////////////////

        for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				out1 << ktime << "  " << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << "  " << abs(conj(wH.phi[wH.index(j,i)])*wH.phi[wH.index(j,i)])*wH.r[j] << endl;
		
		
		out2 << ktime << "  " << 1.-wH.norm() << "  " << 1.-w.norm() << endl;
		
		
		cout << "Propagated (Crank) error norm: " << 1.-w.norm() << "\nInterpolated error norm: " << 1.-wH.norm() << endl;        
    }
	
    
	axes.close();
    out0.close();
    out1.close();
	out2.close();
    
}

