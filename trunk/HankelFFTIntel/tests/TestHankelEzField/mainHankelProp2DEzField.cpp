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
    cout << "mainHankelProp2DEzField. Running example..." << endl;
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
    
    double echarge=-1.;
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
	double sigma = 1.;
	
	//Laser parameters
	
	double efield_z;
    double avect_z;
	
	double e0=sqrt(3.e14/3.5e16);
	double ww=0.057;
	double period=dospi/ww;
	double cycle_number=4.;
	double cep=0.;
	

	
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
    
    
    
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
	
	complex fase;
    
    HankelMatrix HH(Nr,Rmax);

    waveH2D w;
	w.initialize(HH,Nz,dz);
    
    //Initial test function	
    
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			//w.phi[w.index(j,i)]=exp(-(w.r[j]-rho0)*(w.r[j]-rho0)/sigma/sigma-(w.z[i]*w.z[i])/sigma/sigma);
			w.phi[w.index(j,i)]= w.r[j]*exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma )*exp(I*(v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0)) );
			
			//Set potential
			w.pot[w.index(j,i)]=echarge/sqrt(soft_core+(w.r[j]-rho0)*(w.r[j]-rho0)+(w.z[i]-z0)*(w.z[i]-z0));
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
    
	dt=complex(absdt,0.);

    for (int ktime=0; ktime<Ntime; ktime++)
	{
       
        cout << "Loop number: " << ktime << endl;
        
		
		t = ktime*1.*dt;
		//cep=ktime*dospi/40;
		
		efield_z=e0*sin(ww*abs(t)/2./cycle_number)*sin(ww*abs(t)/2./cycle_number)*sin(ww*abs(t)+cep);
		
		avect_z+=-efield_z*abs(dt)*lightC_au;
		
		//cout << "Avector_z: " << avect_z << endl;
        
        ///////////////
        //  Evolve   //
        ///////////////
        
		
        w.prop_kinetic(HH,dt/2.,avect_z,0.0);
		w.prop_potencial(dt);
		w.prop_kinetic(HH,dt/2.,avect_z,0.);
		
		
        ////////////////////////////////////////////
        // Take the snapshot of the wavefunction  //
        ////////////////////////////////////////////
		
		cout << "Error norm: " << 1.-w.norm() << endl << endl;
        
		if((ktime%snap) == 0)
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;	
        
		
		
		/////////////////////////
        //  Apply the absorber //
        /////////////////////////
		
		
		w.absorber(0.1,0.1,0.1,0.1,1./6.);
		
		
        
    }

    
	axis.close();
    out0.close();
    out1.close();
	
	
}

