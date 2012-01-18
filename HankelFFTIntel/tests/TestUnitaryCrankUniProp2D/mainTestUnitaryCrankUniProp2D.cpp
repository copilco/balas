//
//  mainTestUnitaryCrankUniProp2D.cpp
//  
//
//  Created by de la Calle Negro Alejandro on 26/12/11.
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
	
	cout << "\n\n/////////////////////////////////////////////" << endl;
    cout << "/////////////////////////////////////////////" << endl;
    cout << "mainTestUnitaryCrankUniProp2D. Running example..." << endl;
    cout << "/////////////////////////////////////////////" << endl;
	cout << "/////////////////////////////////////////////" << endl;
	
	
	fstream axis("axis.txt",ios::out);
	fstream out0("out0.txt",ios::out);
    fstream out1("outCrankUniError.txt",ios::out);
	
	
	//////////////////
    //  Parameters  //
    //////////////////
	
	int Nr=1000;//520;
	int Nz=1000;//680;
	
	double dz=0.1;
	double dr=0.1;
	complex dt=complex(0.01,0.);
	
	
	int Ntime=5000;
	int snap=30;
	
	//Gaussian parameters
    
	double Rmax  = ceil(Nr*dr);
    double rho0  = Rmax/2.;
	double rho00 = 12.;
	double z0    = 0.;		
	double v0r   = 0.;//5.;
	double v0z   = 0.0;	
	double sigma = 25.;
	
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
    
    waveUniform2D wrho;
	wrho.initialize(HH,Nz,dz);
    
    waveUniform2D wz;
	wz.initialize(HH,Nz,dz);

	
	//Initial test function	
	
	for(int j=0;j<HH.Nr;j++)
		for(int i=0;i<Nz;i++)
		{
			//w.phi[w.index(j,i)]=exp(-(w.r[j]-r0)*(w.r[j]-r0)/sigmar/sigmar-(w.z[i]*w.z[i]))*complex(cos(v0r*(w.r[j]-r0)+v0z*w.z[i]), sin(v0r*(w.r[j]-r0)+v0z*w.z[i])  );
			w.phi[w.index(j,i)]= w.r[j]*exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma )*exp(I* v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0)  );
		}
			
	
	// Check the norm
	
	cout << "\n\nOriginal norm: " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm: " << w.norm() << endl;
	
	
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
    
    
    // Copy the wavefunction to the others objects
	
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
		{
			wrho.phi[wrho.index(j,i)]=w.phi[w.index(j,i)];
			wz.phi[wz.index(j,i)]=w.phi[w.index(j,i)];
		}
	
	
	////////////////////////////
    //  Prepare Crank Arrays  //
    ////////////////////////////

	
	w.PrepareCrankArrays(dt);
    
    wrho.PrepareCrankArrays(dt);
    
    wz.PrepareCrankArrays(dt);
	
	////////////////////////////
    //  Start temporal loop   //
    ////////////////////////////
	
	
	for (int ktime=0; ktime<Ntime; ktime++)
	{
		
		cout << "//////////////////////////////////////////////////////////////////////////////////////////" << endl;
        cout << "Loop number: " << ktime << endl;
		
		///////////////////////////
        //  Evolve in rho axis   //
        ///////////////////////////
		
		wrho.Rprop(dt);
    
        ////////////////////////
        //  Evolve in z axis  //
        ////////////////////////
		
		wz.Zprop1(dt);
		
        
		///////////////////////////
        //  Evolve in both axis  //
        ///////////////////////////
		
		
		w.Zprop1(dt/2.);
		w.Rprop(  dt );
		w.Zprop1(dt/2.);
		
    
        ///////////////////////////////////
        //  Write the error in the norm  //
        ///////////////////////////////////
		
		cout << "Both error  norm: " << 1.-w.norm() << "  Rho error norm: " << 1.-wrho.norm() << "  Z error norm: " << 1.-wz.norm() << endl;
        cout << "//////////////////////////////////////////////////////////////////////////////////////////" << endl;
        
        out1 << ktime << "  " << 1.-w.norm() << "  " << 1.-wrho.norm() << "  " << 1.-wz.norm() << endl;
        
        
        
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
	}
	
	
	
	axis.close();
	out0.close();
	out1.close();
	
}
