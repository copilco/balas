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
    cout << "mainTestImaginaryHankelProp2D. Running example..." << endl;
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
    
    int Nr=1000;
    int Nz=2000;
    
    double dz=.2;
    double dr=0.1;
    double absdt=0.1;
    
    double echarge=-1.;
    double soft_core=.472;//.789545;
    
    int Ntime=5000;
    int snap=20;
    
    //Gaussian parameters
    
    double Rmax  = 150;//ceil(Nr*dr);
    double rho0  = 0.;//Rmax/2.;
    double rho00 = 0.;//12.;
    double z0    = 0.;		
    double v0r   = 0.;//5.;
    double v0z   = 0.0;	
    double sigma = 20.;
    
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
			w.phi[w.index(j,i)]= exp(  -(w.r[j] - rho0)*(w.r[j] - rho0)/sigma/sigma -(w.z[i] - z0)*(w.z[i] - z0)/sigma/sigma );//*exp(I* (v0r*(w.r[j] - rho0) + v0z*(w.z[i] - z0))  );
			
			//Set potential
			w.pot[w.index(j,i)]=-1./sqrt(soft_core+(w.r[j]-rho0)*(w.r[j]-rho0)+(w.z[i]-z0)*(w.z[i]-z0));
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
		
		w.FFTFor();        
		w.phi2F(HH);
		w.HankelTransform(HH);
		
		
		fase=complex(0.,0.);
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=dospi*dospi*HH.v[j]*HH.v[j]*dt/4.+w.q[i]*w.q[i]*dt/4.;
				w.G[w.index(j,i)]*=exp(-I*fase);
			}
		
		
		w.HankelTransformBack(HH);
		w.F2phi(HH);
		w.FFTBack();        
		
		
		
		/////////////////////////////////////////////////////////////////////
		
		fase=complex(0.,0.);
		
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=w.pot[w.index(j,i)]*dt/2.;
				w.phi[w.index(j,i)]*=exp(-I*fase);
			}
		
		
		
		//////////////////////////////////////////////////////////////////////
		
		
		w.FFTFor();
		w.phi2F(HH);
		w.HankelTransform(HH);
		
		
		fase=complex(0.,0.);
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
			{
				fase=dospi*dospi*HH.v[j]*HH.v[j]*dt/4.+w.q[i]*w.q[i]*dt/4.;
				w.G[w.index(j,i)]*=exp(-I*fase);
			}
		
		
		w.HankelTransformBack(HH);
		w.F2phi(HH);
		w.FFTBack();
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		/*
		w.prop_kinetic(HH,dt/2.);
		w.prop_potencial(dt);
		w.prop_kinetic(HH,dt/2.);
		*/
		
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
	
	
	for (int ktime=0; ktime<Ntime; ktime++)
	{
		cout << "Loop number: " << ktime << endl;
		
		wClone.prop_kinetic(HH,dt/2.);
		wClone.prop_potencial(dt);
		wClone.prop_kinetic(HH,dt/2.);
		

	
		////////////////////////////////////////////
		// Take the snapshot of the wavefunction  //
		////////////////////////////////////////////
		/*
		if(ktime%(Ntime/(snap-1)) == 0)
			for(int j=0;j<Nr;j++)
				for(int i=0;i<Nz;i++)
					out2 << abs(conj(wClone.phi[w.index(j,i)])*wClone.phi[w.index(j,i)])*wClone.r[j] << endl;	
		*/
		
		
		// Projection of the two wavefunctions 
		
		proj = 0.;
		
		for(int j=0;j<Nr;j++)
			for(int i=0;i<Nz;i++)
				proj+=w.dz*w.dr[j]*w.r[j]*abs(conj(wClone.phi[w.index(j,i)])*w.phi[w.index(j,i)]);
		
		
		outProj << 1.-proj << endl;
		cout << "< wClone | w > = " << proj << "  Error in projection: " << 1.-proj << endl;
		
		//w.absorber(0.1,0.1,0.1,0.1,1./6.);
	
	}
	
	
	axis.close();
    out0.close();
    out1.close();
	out2.close();
    outEne.close();
	outProj.close();
	
    
}

