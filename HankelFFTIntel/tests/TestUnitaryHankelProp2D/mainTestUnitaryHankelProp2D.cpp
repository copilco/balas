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
    cout << "\n\n/////////////////////////////////////////////" << endl;
    cout << "/////////////////////////////////////////////" << endl;
    cout << "mainTestPropHankel2D. Running example..." << endl;
    cout << "/////////////////////////////////////////////" << endl;
	cout << "/////////////////////////////////////////////" << endl;
	
    
    fstream out0("out0.txt",ios::out);
    fstream out1("outHankelError.txt",ios::out);

    
    //////////////////
    //  Parameters  //
    //////////////////
    
    int Nr=520;
    int Nz=680;
    
    double dr=0.1;
    double dz=0.1;
    double dt=0.01;
    
    int Ntime=5000;
    int snap=30;
    
    //Gaussian parameters
    
	double Rmax  = ceil(Nr*dr);
    double rho0  = Rmax/2.;
	double rho00 = 12.;
	double z0    = 0.;		
	double v0r   = 0.0;//5.;
	double v0z   = 0.0;	
	double sigma = 10;
    
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

    waveH2D w;
	w.initialize(HH,Nz,dz);
    
    waveH2D wrho;
	wrho.initialize(HH,Nz,dz);
    
    waveH2D wz;
	wz.initialize(HH,Nz,dz);
    
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			//w.phi[w.index(j,i)]=exp(-(w.r[j]-rho0)*(w.r[j]-rho0)/sigma/sigma-(w.z[i]*w.z[i])/sigma/sigma);
			w.phi[w.index(j,i)]= w.r[j]*exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma )*
            exp(I* v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0)  );
			//if (w.r[j]>=rho00*0.)
				//v[w.index(j,i)]=-100.0/sqrt(1.+(w.r[j]-rho00)*(w.r[j]-rho00)+(w.z[i] - 2.*z0)*(w.z[i] - 2.*z0) )*1.;
		}
    
    
    // Check the norm
    
    cout << "\n\nOriginal norm: " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm: " << w.norm() << endl;
    
    
    /////////////////////////////////
    //  Save initial wavefunction  //
    /////////////////////////////////
    
    for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out0 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;

    
	// Copy the wavefunction to the others objects
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
		{
			wrho.phi[wrho.index(j,i)]=w.phi[w.index(j,i)];
			wz.phi[wz.index(j,i)]=w.phi[w.index(j,i)];
		}
	
	
    ///////////////////////////////////////////
    //  Define auxiliar variables for tests  //
    ///////////////////////////////////////////
    
    double fase;
    double faseHankel;
    double faseFFT;
    
	DFTI_DESCRIPTOR_HANDLE handle=0;
	MKL_LONG stat;
	
	
	stat = DftiCreateDescriptor(&handle,DFTI_DOUBLE, DFTI_COMPLEX,1,(MKL_LONG)Nz);
	stat = DftiSetValue(handle, DFTI_FORWARD_SCALE, wz.fscale);
	stat = DftiSetValue(handle, DFTI_BACKWARD_SCALE, wz.bscale);
	stat = DftiCommitDescriptor(handle);
    
    ////////////////////////////
    //  Start temporal loop  //
    ///////////////////////////
    
    for (int ktime=0; ktime<Ntime; ktime++)
	{
       
        cout << "Loop number: " << ktime << endl;
        
        /////////////////////////////////////////////////////
        //  Evolve only in rho axis with Hankel transform  //
        /////////////////////////////////////////////////////
        
        
		//wrho.FFTFor();
        wrho.phi2F(HH);
        wrho.HankelTransform(HH);
		
		
		faseHankel=0.;
        
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				faseHankel=dospi*dospi*HH.v[j]*HH.v[j]*dt/2.;
				wrho.G[wrho.index(j,i)]*=exp(I*faseHankel);
			}


        wrho.HankelTransformBack(HH);
        wrho.F2phi(HH);
		//wrho.FFTBack();
		
        
        //////////////////////////////////////
        //  Evolve only in z axis with FFT  //
        //////////////////////////////////////
        
        
		wz.FFTFor();
        //wz.phi2F(HH);
        //wz.HankelTransform(HH);
	
        
		faseFFT=0.;
        
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				faseFFT=wz.q[i]*wz.q[i]*dt/2.;
				wz.phi[wz.index(j,i)]*=exp(I*faseFFT);
			}
		

        //wz.HankelTransformBack(HH);
        //wz.F2phi(HH);
		wz.FFTBack();

		
        
        
        ///////////////////////////
        //  Evolve in both axis  //
        ///////////////////////////
        
        
		w.FFTFor();
        w.phi2F(HH);
        w.HankelTransform(HH);
		
		fase=0.;
		
        for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=dospi*dospi*HH.v[j]*HH.v[j]*dt/2.+w.q[i]*w.q[i]*dt/2.;
				w.G[w.index(j,i)]*=exp(I*fase);
			}

        
        w.HankelTransformBack(HH);
        w.F2phi(HH);
		w.FFTBack();

        
        
        ///////////////////////////////////
        //  Write the error in the norm  //
        ///////////////////////////////////
		
		cout << "Error norm in both axis: " << 1.-w.norm() << "\nError norm in HankelTransform: " << 1.-wrho.norm() << "\nError norm in FFT: " << 1.-wz.norm() << endl;
        
        out1 << ktime << "  " << 1.-w.norm() << "  " << 1.-wrho.norm() << "  " << 1.-wz.norm() << endl;
		
		
        
    }

    
    out0.close();
    out1.close();
    
}

