//
//  mainTestUnitaryCrankUniProp2D.cpp
//  
//
//  Created by de la Calle Negro Alejandro and Alexis Chacon on 26/12/11.
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
#include "laser.h"
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
	
	FILE *laserout;
	laserout=fopen("outa.txt","w");
	
	
	

	
	
	//Laser parameters
	double efieldz;
	double avect_z;
    double wfreq   = 0.25;     //1.1;
	int ncycle     = 2;
	double cep     = pi/2. ;
	
	
	//BUILDing LASER PULSE
	
	int npulses= 1;	                       // Number of pulses	
	
	string env_name ="rsin2";              /**** 
											IMPORTANT The Name of the envelop, may be: 
											rect, sin2, gauss or rsin2.
											*****/ 
	
	double tstart     = 0;					      // Start time of the first pulse
	double realdt     = 0.05;                     // Temporary increase
	double offset     = 5.;                       // Time before the pulse atomic unit
	double outset     = 50.;	                  // Time after the pulse atomic unit
	
	
	
	
	laser fpulse(npulses);                        // Creator of Laser Pulses
	
	
	
	
	
	//First pulse Laser
	fpulse.I0[0]      = 8.75e13;  		          // Intensity W/cm^2 
	fpulse.e[0]       = 0.00;  			          // Elliptical of the pulse
	fpulse.w0[0]      = wfreq;                       // Central frequency
   	fpulse.cycles0[0] = ncycle;                   // Cycles number
   	fpulse.cep0[0]    = cep;                      // Carrier Envelop Phase
	fpulse.phi_rel[0] = 0;                        // Relative phase between the polarization Ex and Ey
	
	
	
	
	double period0 = dospi/fpulse.w0[0];          // Period
	
	//fpulse.delay0[0]  = 0.0;                    // Delay 			
	
	
	
	
	fpulse.envelope[0]=env_name;                  // Envelop name
	cout << "\nEnvelope= "<<fpulse.envelope[0];
		
	
   	fpulse.laser_pulses(realdt, tstart,  offset,  outset);  // Making the linear polarization pulse Ey
	
	
	
	// Save the laser pulse
	for (int ktime=0; ktime<fpulse.Nmaxt; ktime++)
		fprintf(laserout,"%e %e %e %e %e\n",
				fpulse.g.t[ktime],
				real(fpulse.efield.f[ktime]), real(fpulse.avector.f[ktime]),
				imag(fpulse.efield.f[ktime]), imag(fpulse.avector.f[ktime]));
    //End save the laser pulse	
	
	

	
    //////////////////
    //  Parameters  //
    //////////////////    
    
	
	int smallNr = 250;
	int smallNz = 400;	
	
	int largeNr     = 1700;
    int largeNz     = 3000;
    
    double dz  = 0.3;
    double dr  = 0.3;
    complex imagdt = complex(0.05,0.);
    
	
	
	
	
	// Hydrogen parameters     
	double charge_nuclei = 1.;
    double soft_core     = 0.0073685;	
	
	
	
	
	
	//Rmax parameter and snaper
	double smallRmax  = smallNr*dr;	
	double largeRmax  = largeNr*dr;
	
	int snap          = floor(fpulse.Nmaxt/10);	
	int snaper1       = 4;
	int snaper2       = 4;
	
	
	
	
    //Print out the information on the screen 
	
	//cout.setf(ios::scientific);
	cout.precision(8);
	
	cout << "\nNtime= "     <<   fpulse.Nmaxt <<endl;
	cout << "smallNr= "     <<   smallNr      << endl;
	cout << "largeNr= "     <<   largeNr      << endl;	
    cout << "smallNz= "     <<   smallNz      << endl;
    cout << "largeNz= "     <<   largeNz      << endl;	
	cout << "dz= "          <<   dz           << endl;
    cout << "dt= "          <<   abs(imagdt)  << endl;
	cout << "smallRmax= "   <<   smallRmax    << endl;
	cout << "lageRmax= "    <<   largeRmax    << endl;	
    
	
	
	
	
	
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
	HankelMatrix HGround(smallNr, smallRmax);	
    HankelMatrix HH( largeNr, largeRmax );

	
	
    waveUniform2D w;
	waveUniform2D w0;
	waveUniform2D wGround;	
	
	
	//Hankel wave
	waveH2D w1;
	
	
	
	wGround.initialize( HGround, smallNz, dz);	
	
    w.initialize(  HH, largeNz, dz);
    w0.initialize( HH, largeNz, dz);
	w1.initialize( HH, largeNz, dz);


	
	FILE *bwave;
	bwave=fopen("binwave.bin","rb") ;
	
	wGround.binread(bwave);
	
		
    // Set hydrogen potential	
    wGround.set_potential_hlike2D( charge_nuclei, soft_core );	
    w.set_potential_hlike2D( charge_nuclei, soft_core );
    w0.set_potential_hlike2D( charge_nuclei, soft_core ); 
	w1.set_potential_hlike2D( charge_nuclei, soft_core ); 
    
	
	
	// Save axes	
    w.saveAxes(axes);
	w1.saveQAxes(qaxes);
    
	placeWF( w, wGround);

	
	//Check the norm	
    cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm() << endl;
    
    cout <<"\nEnergy= "<<w.kinetic_energy_finite()+w.potential_energy();


  // --------------------------> Time Evolution <---------------------------- //
    
	// Mask parameters
	
    double masc0       = 15.;
    double masc1       = 20.;
	double sigma_masc  = 1.0 ;
	double Energy1     = 0.0;
	double Energy2     = 0.0;
	
	
	
	
	
	
	//////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////
	
	imagdt = complex( 0., -0.05 );
	
	//Preparing Crank-Nicholson arrrays
    
	w.PrepareCrankArrays( imagdt );
	
	for (int ktime=0; ktime<500; ktime++)
	{
		
		
		
		///////////////
		//  Evolve   //
		///////////////
		
		w.Zprop(  imagdt/2. );
		w.Rprop(  imagdt );
		w.Zprop(  imagdt/2. );
		w.normalize();	   
		
		//************************************
		//  Write the error in the norm  
		//************************************
		if (ktime%50==0){
			cout << "Loop number: " << ktime << endl;
			Energy1 = w.kinetic_energy_finite()+w.potential_energy();
		
		     cout << "Energy (Expected value):  " << Energy1 << "   Error: " << log10(abs(Energy1-Energy2)) << endl;
		}
		Energy2=Energy1;
		
		
	}//End propagation

	//*/
	
	
	//Place wavefunction on bigger grid
	cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm() << endl;
    cout <<"\nEnergy= "<<w.kinetic_energy_finite()+w.potential_energy();	
		
	
	
	
	Energy1 = 0.0;
	imagdt = complex(0.05,0.);
	
	//Preparing Crank-Nicholson arrrays
    
	w.PrepareCrankArrays(imagdt);
	
	
	
	
	
	//Start Loop propagation on time
    
	for (int ktime=0; ktime<fpulse.Nmaxt; ktime++)
	{
		
		//w.Zprop_PAG(dt/2., avect_z );
		//w.Rprop(dt);
		//w.Zprop_PAG(dt/2., avect_z );
		
		
		
		efieldz=imag(fpulse.efield.f[ktime]);
		
		
		
		w.Zprop_REG(imagdt/2., efieldz );
		w.Rprop(imagdt);
		w.Zprop_REG(imagdt/2., efieldz );
		
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
			w.snapshot(out3, snaper1, snaper2);
			
			//Snapshot mask
			w0.snapshot(out4, snaper1, snaper2);
			
			interpU2H(w0,w1);
			w1.qsnapshot(out5,HH, snaper1, snaper2);
			
			out6 << ktime <<"  " << fpulse.g.t[ktime] << " " << efieldz <<endl;			
			
		}
		

		
	}//End loop propagation on time
	

	w.mask_function2D( w0, masc0, masc1, sigma_masc);	
	
	out1 << fpulse.Nmaxt << "  " << 1.-w.norm() << " " << Energy1 << endl;
	
	out2 << fpulse.Nmaxt << "  " << 1.-w.norm() <<" "<< w.norm()-w0.norm()<<" "<< 1.-abs(projection(w, wGround ))  << endl;
	
	
	//Snapshot full
	w.snapshot(out3,snaper1,snaper2);
	
	//Snapshot mask
	w0.snapshot(out4,snaper1,snaper2);
	
	interpU2H(w0,w1);
	w1.qsnapshot(out5,HH,snaper1,snaper2);	
	
	FILE *finalwave;
	finalwave=fopen("finalwave.bin","w+") ;
	
	w.binwrite(finalwave);
	
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
