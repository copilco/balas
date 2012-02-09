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
    cout << "\n\n//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "mainTestImaginaryXUVIonization. Running example..." << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    
    
    fstream axis("axis.txt",ios::out);
    fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
	fstream out2("out2.txt",ios::out);
	fstream outEne("outEne.txt",ios::out);
	fstream outProj("outProj.txt",ios::out);
    
    //////////////////
    //  Parameters  //
    //////////////////
    
    int Nr=3000;
    int Nz=2000;
    
    double dz=0.3;
    //double dr=0.1;
    double absdt=0.05;
    complex t;
	
    double echarge=-1.;
    double soft_core=.789545;;//.472;//.789545;
    
    int Ntime=4000;
    int snap=20;
	
	
	//Laser parameters
	
	double efield_z;
    double avect_z;
	
	double e0=sqrt(3.e14/3.5e16);
	double ww=2.;
	double period=dospi/ww;
	double cycle_number=4.;
	double cep=0.;
	int Ntreal=floor(period*cycle_number/absdt);
	
    
    //Gaussian parameters
    
    double Rmax  = 200;//ceil(Nr*dr);
    double rho0  = 0.;//Rmax/2.;
    double rho00 = 0.;//12.;
    double z0    = 0.;		
    double v0r   = 0.;//5.;
    double v0z   = 0.0;	
    double sigma = 4.;
    
    // Print out the information on the screen
    
    cout << "\nNr: " << Nr << endl;
    cout << "Nz: " << Nz << endl;
    cout << "dz: " << dz << endl;
    cout << "dt: " << absdt << endl;
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
	
    complex fase;
    
    HankelMatrix HH(Nr,Rmax);
    
    waveH2D w;
    w.initialize(HH,Nz,dz);
    
    //Initial test function	
    
    
    for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
        {
			//w.phi[w.index(i,j)]=exp(-(w.r[j]-rho0)*(w.r[j]-rho0)/sigma/sigma-(w.z[i]*w.z[i])/sigma/sigma);
			w.phi[w.index(j,i)]= exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma );//*exp(I* (v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0))  );
			
			
			//Set potential
			w.pot[w.index(j,i)]=-2./sqrt(soft_core+(w.r[j]-rho0)*(w.r[j]-rho0)+(w.z[i]-z0)*(w.z[i]-z0));
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
    
	
	//w.SaveGroundState();
	
    
    //////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////
    
    
    complex dt=-I*absdt;
    //complex dt=complex(0.,-absdt);
	//complex dt=complex(absdt,0.);
	
	double ene1 = 0.;
	double ene2 = 0.;
	
    for (int ktime=0; ktime<Ntime; ktime++)
	{
		
        cout << "Loop number: " << ktime << endl;
        
		
		///////////////
		//  Evolve   //
		///////////////
		
				
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		
		w.prop_kinetic(HH,dt/2.);
		w.prop_potencial(dt);
		w.prop_kinetic(HH,dt/2.);
		
		
		w.normalize();
		
		//cout << "Norm in frequency space: " << w.qnorm(HH) << endl;
		
		ene1 = w.kinetic_energy(HH)+w.pot_energy();
		
		cout << "Energy (Expected value):  " << ene1 << "   Error: " << log10(abs(ene1-ene2)) << endl;
		cout << "Norm after Transform  " << w.norm() << endl;
		
		outEne << ene1 << "  " << log10(abs(ene1-ene2)) << endl;
		
		ene2=ene1;
		
		////////////////////////////////////////////
		// Take the snapshot of the wavefunction  //
		////////////////////////////////////////////
		/*
		if(ktime%(Ntime/(snap-1)) == 0)
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;	
		*/
	}
	
	
    //////////////////////////////////////////////////////
	//  Copy the ground state to another wavefunction   //
	//////////////////////////////////////////////////////
	
	waveH2D wClone;
    wClone.initialize(HH,Nz,dz);
	
	
	for(int j=0; j<Nr; j++)
		for (int i=0; i<Nz; i++)
		{
			wClone.phi[wClone.index(j,i)] = w.phi[w.index(j,i)];
			wClone.pot[wClone.index(j,i)] = w.pot[w.index(j,i)];
		}
	
	
	//////////////////////////////////////////////////////
	// Evolve again, this time the cloned wavefunction  //
	//////////////////////////////////////////////////////
	
	double proj;
	dt=complex(absdt,0.);
	
	
	for (int ktime=0; ktime<Ntreal; ktime++)
	{
		cout << "Loop number: " << ktime << endl;
	
		
		t = ktime*1.*dt;
		//cep=ktime*dospi/40;
		
		efield_z=e0*sin(ww*abs(t)/2./cycle_number)*sin(ww*abs(t)/2./cycle_number)*sin(ww*abs(t)+cep);
		
		avect_z+=-efield_z*abs(dt)*lightC_au;
		
		
		
		wClone.prop_kinetic(HH,dt/2.,avect_z);
		wClone.prop_potencial(dt);
		wClone.prop_kinetic(HH,dt/2.,avect_z);
		
	
		////////////////////////////////////////////
		// Take the snapshot of the wavefunction  //
		////////////////////////////////////////////
		
		if(ktime%(Ntreal/(snap-1)) == 0)
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out2 << abs(conj(wClone.phi[w.index(j,i)])*wClone.phi[w.index(j,i)])*wClone.r[j] << endl;	
		
		
		cout << "Norm in propagation: " << wClone.norm() << endl;
		
		// Projection of the two wavefunctions 
		
		proj = 0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				proj+=w.dz*w.dr[j]*w.r[j]*abs(conj(wClone.phi[w.index(j,i)])*w.phi[w.index(j,i)]);
		
		
		outProj << 1.-proj << endl;
		cout << "< wClone | w > = " << proj << "  Error in projection: " << 1.-proj << endl;
		
		w.absorber(0.0,0.1,0.1,0.1,1./6.);
	
	}
	
	
	for (int ktime=0; ktime<2*Ntreal; ktime++)
	{
		cout << "Loop number: " << ktime << endl;
		
		
		t = ktime*1.*dt;
		//cep=ktime*dospi/40;
		
		efield_z=e0*sin(ww*abs(t)/2./cycle_number)*sin(ww*abs(t)/2./cycle_number)*sin(ww*abs(t)+cep);
		
		avect_z=-efield_z*abs(dt)*lightC_au;
		
		
		
		wClone.prop_kinetic(HH,dt/2.);
		wClone.prop_potencial(dt);
		wClone.prop_kinetic(HH,dt/2.);
		
		
		////////////////////////////////////////////
		// Take the snapshot of the wavefunction  //
		////////////////////////////////////////////
		
		if(ktime%(Ntreal/(snap-1)) == 0)
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out2 << abs(conj(wClone.phi[w.index(j,i)])*wClone.phi[w.index(j,i)])*wClone.r[j] << endl;	
		
		
		cout << "Norm in propagation: " << wClone.norm() << endl;
		
		// Projection of the two wavefunctions 
		
		proj = 0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				proj+=w.dz*w.dr[j]*w.r[j]*abs(conj(wClone.phi[w.index(j,i)])*w.phi[w.index(j,i)]);
		
		
		outProj << 1.-proj << endl;
		cout << "< wClone | w > = " << proj << "  Error in projection: " << 1.-proj << endl;
		
		w.absorber(0.0,0.1,0.1,0.1,1./6.);
		
	}
	
	
	axis.close();
    out0.close();
    out1.close();
	out2.close();
    outEne.close();
	outProj.close();
	
    
}

