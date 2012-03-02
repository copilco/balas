//
//  mainTestUnitaryCrankUniProp2D.cpp
//  
//
//  Created by de la Calle Negro Alejandro on 26/12/11.
//
//
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include "arrai.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "waveH2D.h"
#include "waveUniform2D.h"
#include "tools.h"
#include "cylindrical.h"
#include "interp.h"

int main()
{   
	
	cout << "\n\n//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "mainTestXUVplusATTO. Running example..." << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
   	
	
	
	fstream axes("axes.txt",ios::out);
    fstream qaxes("qaxes.txt",ios::out);
	fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    fstream out3("out3.txt",ios::out);    
    fstream out4("out4.txt",ios::out);
    fstream out5("out5.txt",ios::out);	
    fstream out6("out6.txt",ios::out);	
	
	
    //////////////////
    //  Parameters  //
    //////////////////    
    
	int Nr=500;
    int Nz=1000;
    
    double dz = 0.3;
    double dr = 0.3;
    complex dt = complex(0.05,0.);
    
	
	// Hydrogen parameters 
    
	double charge_nuclei = 1.;
    double soft_core     = 0.0073685;//0.00087;
	
	
	//Laser parameters

	double efieldz;
	double avect_z;
    double wx      = 0.057;//1.1;
	double periodx = dospi/wx;
	int ncycle     = 2;
	double sigmax  = ncycle*periodx/(8.*log(2.));
    double efield0 = sqrt( 5.0e14/(3.5e16) );
	double cep   = 0.0  ;
	
	double t       = 0.;
	double tmin    = -1.2*periodx;
	double tmax    = ncycle*periodx;	
	
	int snap=500;
    int Ntime =  floor( (tmax-tmin)/abs(dt) ) +1;
	
	
	//Gaussian parameters
   
	double Rmax  = Nr*dr;
	
	
    //Print out the information on the screen 
	
	//cout.setf(ios::scientific);
	cout.precision(8);
	
	cout << "\nNtime= "<<Ntime <<endl;
	cout << "Nr= " << Nr << endl;
    cout << "Nz= " << Nz << endl;
	cout << "dz= " << dz << endl;
    cout << "dt= " << abs(dt) << endl;
	cout << "Rmax: " << Rmax << endl;
	
    
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
	
    HankelMatrix HH( Nr, Rmax );
	HankelMatrix HGround(250, 150);
	
    waveUniform2D w;
	waveUniform2D wGround;
	waveUniform2D w0;
	
	
	waveH2D w1;
	
	
    w.initialize(  HH, Nz, dz);
    w0.initialize( HH, Nz, dz);
	wGround.initialize( HGround, 400, dz);
	w1.initialize( HH, Nz, dz);
	
	
	
	FILE *bwave;
	bwave=fopen("./binwave.bin","rb") ;
	//w.binread(bwave);
	wGround.binread(bwave);

	placeWF( w, wGround);	
	
	
    // Set hydrogen potential
	
    w.set_potential_hlike2D( charge_nuclei, soft_core );
    w0.set_potential_hlike2D( charge_nuclei, soft_core ); 
	w1.set_potential_hlike2D( charge_nuclei, soft_core ); 
    
	// Save axes
	
    w.saveAxes(axes);
	w1.saveQAxes(qaxes);
    

	//Check the norm	
    cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm() << endl;
    
    

  // --------------------------> Time Evolution <---------------------------- //
    
	// Mask parameters
	
    double masc0 = 15.;
    double masc1 = 30.;
	double sigma_masc = 1.0 ;
	double Energy1 = 0.0;
	double Energy2 = 0.0;
	
	//////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////
	
	dt = complex(0.,-0.05);
	
	//Preparing Crank-Nicholson arrrays
    
	w.PrepareCrankArrays(dt);
	
	for (int ktime=0; ktime<1000; ktime++)
	{
		
		cout << "Loop number: " << ktime << endl;
		
		///////////////
		//  Evolve   //
		///////////////
		
		w.Zprop(  dt/2. );
		w.Rprop(  dt );
		w.Zprop(  dt/2. );
		w.normalize();	   
		
		//************************************
		//  Write the error in the norm  
		//************************************
		
		Energy1 = w.kinetic_energy_finite()+w.potential_energy();
		
		cout << "Energy (Expected value):  " << Energy1 << "   Error: " << log10(abs(Energy1-Energy2)) << endl;
		
		Energy2=Energy1;
		
		
	}//End propagation

	Energy1 = 0.0;
	dt = complex(0.05,0.);
	
	//Preparing Crank-Nicholson arrrays
    
	w.PrepareCrankArrays(dt);
	
	
	//Start Loop propagation on time
    
	for (int ktime=0; ktime<Ntime+1; ktime++)
	{
		
		t= tmin + ktime*abs(dt);
		
		// Gausian pulse
		efieldz = efield0*exp(-t*t/sigmax/sigmax)*sin( wx*t + cep );
		
		// Sin^2 pulse
		//efieldz = efield0*sin(wx*t/2./ncycle)*sin(wx*t/2./ncycle)*sin(wx*t+cep);
		
		
		avect_z+=-efieldz*abs(dt);
		
		
		//w.Zprop_PAG(dt/2., avect_z );
		//w.Rprop(dt);
		//w.Zprop_PAG(dt/2., avect_z );
		
		w.Zprop_REG(dt/2., efieldz );
		w.Rprop(dt);
		w.Zprop_REG(dt/2., efieldz );
		
		w.absorber(0.0, 0.1, 0.1, 0.1, 1./6.);	   
		
		
		//***************************************
	    //     Write the error in the norm      
	    //*************************************** 
		
	    if (ktime%snap==0)
		{
			
	     	Energy1 = w.kinetic_energy_finite()+w.potential_energy();
			
			
			cout << "Loop number = " << ktime << "  Energy = " << Energy1 << "  Norm = " << w.norm() << "  Error in the norm = " << 1.-w.norm() << endl;
		
			
			w.mask_function2D( w0, masc0, masc1, sigma_masc);	
			
	       	out1 << ktime << "  " << 1.-w.norm() << " " << Energy1 << endl;
			
			out2 << ktime << "  " << 1.-w.norm() <<" "<< w.norm()-w0.norm()<<" "<< 1.-abs(projection(w, wGround ))  << endl;
			
			
			//Snapshot full
			w.snapshot(out3,1,1);
			
			//Snapshot mask
			w0.snapshot(out4,1,1);
			
			interpU2H(w0,w1);
			w1.qsnapshot(out5,HH);
			
			
		}
		
		out6 << t << " " << efieldz <<endl;
		
	}//End loop propagation on time
	
	
	
	axes.close();
	qaxes.close();
	out0.close();
	out1.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	out6.close();
	
	
}
