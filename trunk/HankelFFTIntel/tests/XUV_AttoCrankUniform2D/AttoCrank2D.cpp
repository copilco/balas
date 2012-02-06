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
#include "cylindrical.h"

int main()
{   
   	
    fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    fstream out3("out3.txt",ios::out);    
    fstream out4("out4.txt",ios::out);
    fstream out5("out5.txt",ios::out);	
    fstream out6("out6.txt",ios::out);	
    fstream out7("out7.txt",ios::out);		
    
	
	
	
	
	
    //////////////////
    //  Parameters  //
    //////////////////    
    int Nr=250;//520;
    int Nz=400;//680;
    
    double dz = 0.3;
    double dr = 0.3;
    double dt = 0.05; //complex(0.01,0.);
    
	
	
	
	
	//Laser parameters
    double wx      = 1.1;
	double periodx = dospi/wx;
	int ncycle     = 5;
	double sigmax  = ncycle*periodx/(8.*log(2.));
    double efield0 = sqrt( 5.0e12/(3.5e16) );
	
	double t       = 0.;
	double tmin    = -2.5*periodx;
	double tmax    =  5.0*periodx;	
	
	
	
	
	//Snapers
	int snap=30;
    int Ntime = floor( (tmax-tmin)/dt ) +1;
    
	cout << "\nNtime= "<<Ntime <<endl;
	
	
	
    //Gaussian parameters
    double Rmax  = Nr*dr;
	
    
	
    
    //Print out the information on the screen    
    cout << "\nNr= " << Nr << endl;
    cout << "Nz= " << Nz << endl;
    cout << "dz= " << dz << endl;
    cout << "dt= " << dt << endl;
	
    
    
    
    
    
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
    HankelMatrix HH( Nr, Rmax );
    HankelMatrix HH0( Nr, Rmax );	
	
    waveUniform2D w;
	waveUniform2D w0;
	waveUniform2D w1;
	
    w.initialize(  HH, Nz, dz );
    w0.initialize( HH0, Nz, dz);
    w1.initialize(  HH, Nz, dz );	
    //w1.initialize( HH, 2*Nz, dz);    
	
    
	
	FILE *bwave;
	bwave=fopen("../ImaginaryPropCrank2D/binwave.bin","rb") ;
	w0.binread(bwave);
	placeWF( w, w0);	
	
	
	cout << "\nNorm= "<<w.norm() <<endl;
	
	
	//return 0;
	// Hydrogen parameters 
    double charge_nuclei = 1.;
    double soft_core     = 0.0073685;//0.00087;
    
    w.set_potential_hlike2D( charge_nuclei, soft_core );
    w0.set_potential_hlike2D( charge_nuclei, soft_core );    
    
    w.saveAxes( out0 );
    
    
	
	
	
	//Check the norm	
    cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm() << endl;
    
    
    
	
	

	
    //****************************************************    
    //Propagation on time
    
	//Preparing Crank-Nicholson
    w.PrepareCrankArrays(dt);
    
    

	
	double cep   = 0.0  ;
    double masc0 = 15.  ;
    double masc1 = 20.  ;
	double sigma_masc = 1.0 ;
	
	
	//Start Loop propagation on time
    for (int ktime=0; ktime<Ntime+1; ktime++)
	{
		t=tmin + ktime*dt;
		
		double efieldz = efield0*exp(-t*t/sigmax/sigmax)*sin( wx*t + cep );
		
		w.Zprop_REG( complex( dt/2., 0.), efieldz );
		w.Rprop(  complex( dt, 0.) );
		w.Zprop_REG( complex( dt/2., 0.), efieldz );
		
		/*		 w.Zprop_PAG(dt/2., efieldz );
		 w.Rprop(  complex(dt,0.) );
		 w.Zprop_PAG(dt/2., efieldz );		*/
		
		
		w.absorber(0.0, 0.1, 0.1, 0.1, 1./4.);	   
		
		
		//***************************************
	    //     Write the error in the norm      
	    //*************************************** 
		
	    if (ktime%snap==0){
			
			//placeWF(w1, w);	
	     	double Energy1=w.kinetic_energy_finite()+w.potential_energy();
			cout << " ktime= "<<ktime ;
	     	cout << "  ErrorNorm= "           << abs(1.-w.norm());
	       	cout << "  Energy= "                   << Energy1 <<endl;	
			
			
			w1.mask_function2D( w1, masc0, masc1, sigma_masc, 0., 0.);			
	       	out1 << ktime << "  " << abs(1.-w.norm()) << endl;
			
			//Snapshot full
			w.snapshot(out3,1,1);

			//Snapshot mask
			w1.snapshot(out6,1,1);
			
			//fprintf(out4, "%e %e %e \n ", ErrorEnergy1, abs(ErrorEnergy1-ErrorConvergence), abs(1.-w.norm()) );			
			out4 << abs(1.-w.norm()) << " " <<Energy1 << endl;
			
			
		}

		//fprintf(out5, "%e %e \n", t, efieldz);
		out5 << t << " " <<efieldz<<endl;
	}//End loop propagation on time
	//return 0;	
	
	
	out0.close();
	out1.close();
	out2.close();
	//out3.close();
	//out4.close();
}
