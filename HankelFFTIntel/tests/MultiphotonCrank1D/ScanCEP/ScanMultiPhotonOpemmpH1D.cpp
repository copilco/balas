//
//  mainTestUnitaryCrankUniProp2D.cpp
//  
//
//  Created by de la Calle Negro Alejandro and Alexis Chacon on 26/12/11.
//   
//
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include "arrai.h"
#include "constant.h"
#include "waveUniformZ1D.h"
#include "tools.h"
#include "cylindrical.h"
#include "laser.h"
#include <sstream>

int iparam = 0;     // integer parameter
int main(int argc,char *argv[])
{
	
	mkl_set_num_threads(12);	
    iparam = atoi(argv[optind]);
    cout << "CEP parameter = " << iparam << endl;
	
    
	
	int Nfile = 8;
    FILE *out[Nfile];
	
    FILE *bwave0;
    FILE *bwave1;

	string ofile ="out";
	
	for(int i=0; i<Nfile; i++)
	{
		stringstream aout0;
		stringstream aout1;
		
		aout0 << ofile << i ;
		string aout_s0 = aout0.str() + "D" ;
		
		
		aout1 << aout_s0 << iparam ;
		
		string file_f = aout1.str()+".txt";
		out[i] = fopen(file_f.c_str(),"w");
	}
	
	
	
	
    bwave0 = fopen("/Users/alexis/Dropbox/BALAS/TestingUniformeZ1D/binwavel.bin","rb") ;
	
	
	
	
	
	FILE *laserout;
	laserout=fopen("outa.txt","w");
	
	
	
	
	
	
	
	
	
	
	/*********************************************************
	        Laser parameters and declare object
	**************************/
	
	
	//Laser parameters
	double efieldz;
	double avect_z;
    double wfreq   = 0.25;     //1.1;
	int ncycle     = 6;
	double cep     = pi/2.+iparam*dospi/20. ;
	
	
	//BUILDing LASER PULSE
	
	int npulses= 1;	                       // Number of pulses	
	
	string env_name ="rsin2";              /**** 
											IMPORTANT The Name of the envelop, may be: 
											rect, sin2, gauss or rsin2.
											*****/ 
	
	double tstart     = 0;					      // Start time of the first pulse
	double realdt     = 0.01;                     // Temporary increase
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
	
	

	
	
	
	
	
	
	
	
	
	
	/***************************************************
	   Definition parameter on grid and declare object
	 **************************/
    int Nz     = 150000;
	int lNz    = 250000;
    double dz  = 0.05;
		

    
	
	// Hydrogen parameters     
	double charge_nuclei = 1.;
    double soft_core     = 2.;	
	
	
	
	int snap          = floor(fpulse.Nmaxt/10);	
	int snaper1       = 1;
	
	
	
	
    //Print out the information on the screen 
	cout << "\ndz= " << dz;
	cout << "  Nz= " << Nz;
	
	
	
	
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
    waveUniformZ1D w;
    waveUniformZ1D wlarge;
    waveUniformZ1D wlarge1;	

	
	
    w.initialize(Nz,dz);
    wlarge.initialize(  lNz, dz);	
    wlarge1.initialize( lNz, dz);
	
	
	
	
	//Potential well
	w.set_potential_hlike1D( charge_nuclei, soft_core );	    
    wlarge.set_potential_hlike1D( charge_nuclei, soft_core );	    
    wlarge1.set_potential_hlike1D( charge_nuclei, soft_core );    
	
    
	
	
	
	
	// Initial test function	
    w.binread(bwave0);
	wlarge.placeWF(wlarge,w);    		// Place wavefunction
    
	
	
	
	
	//*******************************************
    //  Save initial wavefunction and axis  //
    //*******************************************    
    for(int i=0;i<wlarge.Nz;i++)
		fprintf(out[0],"%e \n", wlarge.z[i]);
	
	
	
	
	
	
	//Check the norm	
    cout << "\nOriginal norm= " << wlarge.norm();
    wlarge.normalize();
    cout << "    Initial norm= " << wlarge.norm();    
    cout <<"\nEnergy= "<<wlarge.kinetic_energy_finite()+wlarge.potential_energy();
	
	
	
	
	
	
	
	
	
	
  // --------------------------> Time Evolution <---------------------------- //
    
	// Mask parameters
	
    double masc0       = 15.;
    double masc1       = 20.;
	double sigma_masc  = 1.0 ;
	double Energy1     = 0.0;
	double Energy2     = 0.0;
	
	
	
	
	
	
	//************************************
	//  Imaginary time propagation
	//************************************
	complex dt= complex( 0., -realdt );
	
	for (int ktime=0; ktime<700; ktime++)
	{		
		wlarge.Zprop(  dt );
		wlarge.normalize();	   
		
		
		
		//************************************
		//  Write the error in the norm  
		//************************************
		if (ktime%50==0){
			cout << "\nLoop number= " << ktime;
			Energy1 = wlarge.kinetic_energy_finite()+wlarge.potential_energy();
		
		     cout << "Energy =  " << Energy1 << "   Error: " << log10(abs(Energy1-Energy2)) << endl;
		}
		Energy2 = wlarge.kinetic_energy_finite()+wlarge.potential_energy();		

		
	}//End propagation
	
	

	
	//Place wavefunction on bigger grid
	cout << "\nOriginal norm= " << wlarge.norm() << endl;
    wlarge.normalize();
    cout << "Initial norm= " << wlarge.norm() << endl;
    cout <<"\nEnergy= "<<wlarge.kinetic_energy_finite()+wlarge.potential_energy();	
		
	
	
	
	
	
	
	
	//*************************************************************
	//     Time real evolution of the wave function on the laser      
	//*************************************** 			
	
	//Time step
	dt = realdt;
	
	
	int kont = 1;
	int Nmom = 10000;
	
	double k    =  0.;
	double dk   =  0.001;
	double kmin = -dk*(Nmom-1)/2.;	
	
	
	Energy1 = 0.0;	
	
	
	
	
	//Start Loop propagation on time    
	for (int ktime=0; ktime<fpulse.Nmaxt; ktime++)
	{
		
		efieldz=imag(fpulse.efield.f[ktime]);
		
		wlarge.Zprop_REG(dt, efieldz );		
		wlarge.absorber1D( 0.1, 0.1, 1./6.);	   
		
		
		
		
		
		
		//***************************************
	    //     Write the error in the norm      
	    //*************************************** 
	    if (ktime%snap==0)
		{
			
			
			
			
	     	Energy1 = wlarge.kinetic_energy_finite() + wlarge.potential_energy();			
			wlarge1.mask_function1D( wlarge1, masc0, masc1, sigma_masc,0.);	
			
			
			cout << "\nLoop number = "         << ktime;
			cout << "  kont= "                 << kont;
			cout << "  Energy = "              << Energy1;
			cout << "  Norm = "                << wlarge.norm();
			cout << "  FNorm = "               << wlarge1.norm();
			cout << "  Error in the norm = "   << 1.-wlarge.norm();
			
			fprintf(out[1],"%d %e %e \n", ktime, wlarge.norm(), wlarge1.norm());
			
			
			
			
			
			
			//wlarge.snapshot(out2, snaper1); //Snapshot full
			//wlarge1.snapshot(out3, snaper1); //Snapshot mask
			
			
			
			
			
			
			
			/*
			//******************************************************
			//     Building projection on scattering wave      
			//*************************************** 				
			for (int kmom=0; kmom<Nmom; kmom++)
			{
				
				
				wlarge2.kmomenta= kmin+dk*kmom;
				wlarge2.Continuum_WF();
				
				
				
				complex proj=wlarge2.projection(wlarge2, wlarge);				
				fprintf(out4,"%e %e %e \n", real(proj), imag(proj), abs(proj));
				
				
				
				
				if (kmom%500==0) {
					cout << "\nkmom= "                << kmom;
					cout << "  k= "                 << wlarge2.kmomenta;	
					cout << "  'ScatteringNorm'= "  << wlarge2.norm();
					cout << "  Adip= "              << abs(proj);				
				}
				
				
				
			  }//End propagation*/
			kont+=1;
		
			
			
		}//End snapshot
		
		
		
	}//End loop propagation on time
	
	
	
	
	//*/
	
	
	
	
	
	
	Energy1 = wlarge.kinetic_energy_finite() + wlarge.potential_energy();			
	wlarge1.mask_function1D( wlarge1, masc0, masc1, sigma_masc,0.);	
	
	
	cout << "\nLoop number = "         << fpulse.Nmaxt;
	cout << "  Energy = "              << Energy1;
	cout << "  Norm = "                << wlarge.norm();
	cout << "  FNorm = "               << wlarge1.norm();
	cout << "  Error in the norm = "   << 1.-wlarge.norm();
	
	
	
	
	fprintf(out[1],"%d %e %e \n", fpulse.Nmaxt, wlarge.norm(), wlarge1.norm());
	
	
	
	
	wlarge.snapshot(  out[2], snaper1); //Snapshot full
	wlarge1.snapshot( out[3], snaper1); //Snapshot mask
	
	
	
	
	
	
	
	
	double right_dist=0;
	double left_dist=0;	
	complex *proj=(complex*)mkl_malloc(Nz*sizeof(complex),16);		
	
	int tid;
	int chunk =10;
	
#pragma omp parallel shared(chunk,dk) private(tid)
	{
        
		waveUniformZ1D wlarge2;		
		wlarge2.initialize( lNz, dz);	
		

		
		
		//Potential well
		wlarge2.set_potential_hlike1D( charge_nuclei, soft_core );  		
		
		
		
#pragma omp for schedule(dynamic,chunk)	
	
	//******************************************************
	//     Building projection on scattering wave      
	//*************************************** 				
	for (int kmom=0; kmom<Nmom; kmom++)
	{
		wlarge2.kmomenta= kmin+dk*kmom;
		wlarge2.Continuum_WF();
		
		
		
		proj[kmom] = wlarge2.projection(wlarge2, wlarge);
		
		
		
		
		if (kmom%500==0) {
			cout << "\nkmom= "                << kmom;
			cout << "  k= "                 << wlarge2.kmomenta;	
			cout << "  'ScatteringNorm'= "  << wlarge2.norm();
			cout << "  Adip= "              << abs(proj[kmom]);
		}
		
		
	
		

	}//End propagation
	
	}	
	
	
	
	
	
	
	
	//Save momentum distribution and calculing asymmetry
	for (int kmom=0; kmom<Nmom; kmom++) 
	{
		double k= kmin+dk*kmom;
		if (k<0) {
			left_dist+= abs(proj[kmom])*dk*k*k;
		}
		else {
			right_dist+= abs(proj[kmom])*dk*k*k;
		}
		
		fprintf(out[4],"%e %e %e\n", real(proj[kmom]), imag(proj[kmom]), abs(proj[kmom]));
		
		fprintf(out[5],"%e \n",k);
	
	}
	
	
	
	
	fprintf(out[6],"%d %e %e %e %e\n",
			kont,fpulse.cep0[0]-pi/2.,
			left_dist, right_dist,
			(left_dist-right_dist)/(left_dist+right_dist));	
	
	
	
	
	
	
	
	
	mkl_free(proj);	
	for(int i=0;i<Nfile;i++)
		fclose(out[i]);
/*	fclose(out1);
	fclose(out2);
	fclose(out3);	
	fclose(out4);
	fclose(out5);
	fclose(out6);	*/
	
	cout << "\nEnd of the program"<<endl;
}
