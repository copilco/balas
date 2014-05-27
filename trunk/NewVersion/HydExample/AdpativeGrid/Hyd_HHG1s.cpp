//
//  3D Simulation of Multiphoton Regimen on cylindrical approximation
//					Ground state Hydrogen atom
//
//  Created by Alexis Chacon on 26/12/11.
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
#include <time.h>

double dipole_acceleration( waveUniform2D &w);

int iparam=1;
int main (int argc, char * const argv[])
{   
	
  iparam = atoi(argv[1]);
  cout << "iparam= " << iparam;

   cout << "\n\n//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "TDSE Simulation . ..." << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
   	
	
	clock_t ct0,ct1;
	
	
	fstream axes("axes.txt",ios::out);
    fstream qaxes("qaxes.txt",ios::out);
	fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    fstream out3("out3.txt",ios::out);    
    fstream out4("out4.txt",ios::out);
    fstream out5("out5.txt",ios::out);	
    fstream out6("out6.txt",ios::out);	
	
	
	fstream B1axes("B1axes.txt",ios::out);
    fstream B1qaxes("B1qaxes.txt",ios::out);
	fstream B1out0("B1out0.txt",ios::out);
    fstream B1out1("B1out1.txt",ios::out);
    fstream B1out2("B1out2.txt",ios::out);
    fstream B1out3("B1out3.txt",ios::out);    
    fstream B1out4("B1out4.txt",ios::out);
    fstream B1out5("B1out5.txt",ios::out);	
    fstream B1out6("B1out6.txt",ios::out);	
	
	
	fstream B2axes("B2axes.txt",ios::out);
    fstream B2qaxes("B2qaxes.txt",ios::out);
	fstream B2out0("B2out0.txt",ios::out);
    fstream B2out1("B2out1.txt",ios::out);
    fstream B2out2("B2out2.txt",ios::out);
    fstream B2out3("B2out3.txt",ios::out);    
    fstream B2out4("B2out4.txt",ios::out);
    fstream B2out5("B2out5.txt",ios::out);	
    fstream B2out6("B2out6.txt",ios::out);	
	
	
	
	FILE *laserout, *out7, *out8, *out9, *out10, *out11;
	laserout=fopen("outa.txt","w");
	
	out7	= fopen("out7.txt","w");
	out8	= fopen("out8.txt","w");
	out9	= fopen("out9.txt","w");
	out10	= fopen("out10.txt","w");
	out11	= fopen("out11.txt","w");	

	//=================================//	
	// LASER OBJECT DEFINITION 
	//.................................//	
	
	//Laser parameters
	double E0		= 0.01 + .01*iparam ;//0.04;
	double I0		= E0*E0*3.5e16;			//Intensity w/cm^2
	double wfreq	= 0.057;//0.55185513-0.38881693 -0.003;//1.634e-3;//0.057;				//fequency a.u.
	
	
	int ncycle		= 4;				//Number Cycles / be equivalent to time-bandwidth
	double cep		= 0*pi/2. ;			//CEP
	double period0  = dospi/wfreq;       // Period
	double ir_period = period0;
	
	
	
	//Creating LASER PULSE	
	int npulses      = 1;	                       // Number of pulses		
	string env_name ="rsin2";              
	
										/**** 
											IMPORTANT The Name of the envelop, may be: 
											rect, sin2, gauss or rsin2.
											*****/ 
	
	double tstart     = 0;					// Start time of the first pulse
	double realdt     = 0.05;				// Temporary increase
	
	double offset     = 55.;					// Time before the pulse atomic unit
	double outset     = 55.;				// Time after the pulse atomic unit
	
	
		
	laser fpulse(npulses);					// Constructor of Laser Pulses
	
	
	
	//First pulse Laser
	fpulse.I0[0]		= I0;					// Intensity W/cm^2 
	fpulse.e[0]			= 0.00;				// Elliptical of the pulse
	fpulse.w0[0]		= wfreq;              // Central frequency
   	fpulse.cycles0[0]	= ncycle;				// Cycles number
   	fpulse.cep0[0]		= cep;				// Carrier Envelop Phase
	fpulse.phi_rel[0]	= 0.;					// Relative phase between the polarization Ex and Ey
	
		
	fpulse.envelope[0]	= env_name;            // Envelop name
	
	
   	fpulse.laser_pulses(realdt, tstart,  offset,  outset);  // Making the linear polarization pulse Ey
	
	
	//Printing laser parameters
	cout << "\n//****************************************//"<<endl;
	cout << "           LASER PARAMETERS     "<<endl;
	cout << "//****************************************//"<<endl;	
	cout << "\nTimeStep= "		<<realdt ;
	cout << "    Ntime=  "		<<fpulse.Nmaxt ;	
	cout << "\nIntensity=  "	<<fpulse.I0[0] ;
	cout << "\nFrequency=  "	<<fpulse.w0[0] ;	
	cout << "\nCycles=  "		<<fpulse.cycles0[0] ;	
	cout << "\nCEP= "			<<fpulse.cep0[0];		
	cout << "\nElliptical= "	<<fpulse.e[0] ;	
	cout << "\nRelative_Elliptical_Phase= "   <<fpulse.phi_rel[0]<<endl;		
	
	
	
	
	// Save the laser pulse
	for (int ktime=0; ktime<fpulse.Nmaxt; ktime++)
		fprintf(laserout,"%e %e %e %e %e %e\n",
				fpulse.g.t[ktime],
				real(fpulse.efield.f[ktime]), real(fpulse.avector.f[ktime]),
				imag(fpulse.efield.f[ktime]), imag(fpulse.avector.f[ktime]),
				imag(fpulse.env[0].f[ktime]));
    //End save the laser pulse	
	
	
	
	
	//Laser parameters // 
	fprintf(out9,"%e %e %e %e %e %e ", 
			fpulse.I0[0], fpulse.w0[0], 
			fpulse.cycles0[0], fpulse.cep0[0], 
			fpulse.e[0], fpulse.phi_rel[0]);
	
	
	
	
	
	
	//=================================//	
	// SPATIAL PARAMETERs 
	//.................................//	
    
	
	int smallNr		  = 300;//1000;//
	int smallNz		  = 600;//3000;//	
	
	
	
	
	int largeNr       = 800;
	int largeNz       = 1600;
    
	
	
	
    double dz		  = 0.3;
    double dr		  = 0.3;
    complex complexdt = realdt;//complex(dr*dr,0.);
    
	
	
	
	
	// Hydrogen parameters     
	double V0		= 7.565;
    double alpha	= 0.160;	
	
	
	
	
	//Rmax parameter and snaper
	double smallRmax  = smallNr*dr;	
	double largeRmax  = largeNr*dr;
	
	int snap          = floor(fpulse.Nmaxt/40);	
	int snaper1       = 4;
	int snaper2       = 4;
	
	
	
	
    //Print out the information on the screen 
	
	//cout.setf(ios::scientific);
	cout.precision(8);
	

	cout << "\n//****************************************//"<<endl;
	cout << "    Grid Parameters     "<<endl;
	cout << "//****************************************//"<<endl;					
	
	cout << "\nsmallNr= "     <<   smallNr      << endl;
	cout << "largeNr= "     <<   largeNr      << endl;	
    cout << "smallNz= "     <<   smallNz      << endl;
    cout << "largeNz= "     <<   largeNz      << endl;	
	cout << "dz= "          <<   dz           << endl;
    cout << "dt= "          <<   abs(complexdt)  << endl;
	cout << "smallRmax= "   <<   smallRmax    << endl;
	cout << "lageRmax= "    <<   largeRmax    << endl;	
    
	
	fprintf(out8,"%e %e %d %d %d %d \n",  dr, dz, smallNr, smallNz, largeNr, largeNz);
	

	
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
	HankelMatrix HGround(smallNr, smallRmax);	
    HankelMatrix HH( largeNr, largeRmax );

	
	
    waveUniform2D w;
	waveUniform2D w0;
    waveUniform2D wg0;
	waveUniform2D wg1;	
    waveUniform2D wg2;
	waveUniform2D wg3;	
	
	waveUniform2D wGround;	
	waveUniform2D wGround0;		
	waveUniform2D wGround1;		
	waveUniform2D wGround2;		
	waveUniform2D wGround3;		
	
	
	
	//Declare Object Hankel wave
	waveH2D w1;
	
	
	
	wGround.initialize(  HGround, smallNz, dz);
	wGround0.initialize( HGround, smallNz, dz);	
	wGround1.initialize( HGround, smallNz, dz);
	wGround2.initialize( HGround, smallNz, dz);	
	wGround3.initialize( HGround, smallNz, dz);	
	
    w.initialize(  HH, largeNz, dz);
    w0.initialize( HH, largeNz, dz);
	w1.initialize( HH, largeNz, dz);

    wg0.initialize( HH, largeNz, dz);
	wg1.initialize( HH, largeNz, dz);	
    wg2.initialize( HH, largeNz, dz);
	wg3.initialize( HH, largeNz, dz);		
	
	//return 0;	
	
	//Load ground state with binary function
	FILE *bwave,*bwave0,*bwave1,*bwave2,*bwave3;
	//bwave=fopen("/Users/alexisagustinchaconsalazar/Dropbox/CODE/BALAS/TestPoshTeller/Bound/Grid300x600/binwave0.bin","rb");
	//bwave=fopen("/Users/alexisagustinchaconsalazar/Dropbox/CODE/BALAS/TestPoshTeller/Bound/Grid1000x3000/Bigbinwave0.bin","rb") ;				
	///Users/alexis/Google Drive/LaserInteraction/BALAS/TestPoshTeller/Bound
	
	
	/*
	bwave=fopen("/Users/alexis/Google\ Drive/LaserInteraction/BALAS/TestPoshTeller/Bound/dr03/Grid300X600/binwave0.bin","rb") ;
	bwave0=fopen("/Users/alexis/Google\ Drive/LaserInteraction/BALAS/TestPoshTeller/Bound/dr03/Grid300X600/binwave0.bin","rb") ;
	bwave1=fopen("/Users/alexis/Google\ Drive/LaserInteraction/BALAS/TestPoshTeller/Bound/dr03/Grid300X600/binwave1.bin","rb") ;
	bwave2=fopen("/Users/alexis/Google\ Drive/LaserInteraction/BALAS/TestPoshTeller/Bound/dr03/Grid300X600/binwave2.bin","rb") ;
	bwave3=fopen("/Users/alexis/Google\ Drive/LaserInteraction/BALAS/TestPoshTeller/Bound/dr03/Grid300X600/binwave2.bin","rb") ;
	*/
	
	bwave=fopen("../Bound/binwave0.bin","rb") ;
	bwave0=fopen("../Bound/binwave0.bin","rb") ;
	bwave1=fopen("../Bound/binwave1.bin","rb") ;
	bwave2=fopen("../Bound/binwave2.bin","rb") ;
	bwave3=fopen("../Bound/binwave2.bin","rb") ;
	
	
	
	
	wGround.binread( bwave);
	wGround0.binread(bwave0);	
	wGround1.binread(bwave1);
	wGround2.binread(bwave2);	
	wGround3.binread(bwave3);	
	
	//wGround.rho_gaussian( 0., 0., 2., 2. );
	
	
    // Set hydrogen potential	
    wGround.set_potential_hlike2D(   1., 0.);
    wGround0.set_potential_hlike2D(  1., 0.);
    wGround1.set_potential_hlike2D(  1., 0.);	
    wGround2.set_potential_hlike2D(  1., 0.);
    wGround3.set_potential_hlike2D(  1., 0.);	
	
	
    w.set_potential_hlike2D(    1., 0.);
    w0.set_potential_hlike2D(   1., 0.); 
	w1.set_potential_hlike2D(	1., 0.);
    
    wg0.set_potential_hlike2D(  1., 0.);
	wg1.set_potential_hlike2D(  1., 0.); 	
    wg2.set_potential_hlike2D(  1., 0.); 
	wg3.set_potential_hlike2D(  1., 0.);
	
	
	
	// Save axes	
    w.saveAxes(axes,snaper1,snaper2);
	w1.saveQAxes(qaxes,snaper1,snaper2);
    
	placeWF( w,   wGround);
	placeWF( wg0, wGround0);
	placeWF( wg1, wGround1);
	placeWF( wg2, wGround2);
	placeWF( wg3, wGround3);	

	
	
	
	
	cout << "\n//****************************************//"<<endl;
	cout << "    Check norm and ground-state eigen energy     "<<endl;
	cout << "//****************************************//"<<endl;					
	
    cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm();    
	
	
	
	double EnergyTemp	= wGround0.kinetic_energy_finite()+wGround0.potential_energy();
	double EnergyTemp1	= wGround1.kinetic_energy_finite()+wGround1.potential_energy();
	double EnergyTemp2	= wGround2.kinetic_energy_finite()+wGround2.potential_energy();
	double EnergyTemp3	= wGround3.kinetic_energy_finite()+wGround3.potential_energy();
	
	
	
    cout <<"\nE0LargeG= "<<EnergyTemp<<endl;
	cout << "E0SmallG = " << w.kinetic_energy_finite()+w.potential_energy();
	
    cout <<"\nE1LargeG= "<<EnergyTemp1<<endl;
	cout << "E1SmallG = " << wg1.kinetic_energy_finite()+wg1.potential_energy();	

    cout <<"\nE2LargeG= "<<EnergyTemp2<<endl;
	cout << "E2SmallG = " << wg2.kinetic_energy_finite()+wg2.potential_energy();	

    cout <<"\nE3LargeG= "<<EnergyTemp3<<endl;
	cout << "E3SmallG = " << wg3.kinetic_energy_finite()+wg3.potential_energy();	
	
	
	fprintf(out9," %e  %e  %e  %e ",abs(EnergyTemp), EnergyTemp1,EnergyTemp2,EnergyTemp3);
	
    
	
	

	
	//=================================//	
	//		Temporary propagation 
	//.................................//		
	
	// Mask parameters	
    double masc0       = 20.;
    double masc1       = 35.;
	double sigma_masc  = 10. ;
	double Energy1     = 0.0;
	double Energy2     = 0.0;
	
	
	
	
	
	
	cout << "\n//****************************************//"<<endl;
	cout << "    Check norm and ground-state eigen energy     "<<endl;
	cout << "//****************************************//"<<endl;				
	

	cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm() << endl;
    cout <<"Energy= "<<w.kinetic_energy_finite()+w.potential_energy();	
		
	
	
	
	
	cout << "\n//****************************************//"<<endl;
	cout << "           Real Temporary evolution     "<<endl;
	cout << "//****************************************//"<<endl;			

	
	Energy1 = 0.0;
	complexdt = realdt;//complex(.05,0.);
	cout << "\ndt ="<<complexdt<< "  Ntime= "<<  fpulse.Nmaxt << endl;
	
	
	// ****************************** 
	// Time Flags //
	
	cout << "\nTime Flags \n";
	double t0		= fpulse.g.t[0];
	
	double ir_t0	= fpulse.clock0[1]+fpulse.cycles0[1]*ir_period*.3;
	double ir_tf	= fpulse.clock0[1]+fpulse.cycles0[1]*ir_period*.5;//fpulse.clock2[1];	
	double ir_tfs	= ir_tf + fpulse.cycles0[1]*ir_period;
	

	
	int NtimeEvol   = floor( ( ir_tfs - t0 )/fpulse.dt)+1;
	
	int ir_t0_index = floor( ( ir_t0 - t0 )/fpulse.dt )+1;
	int ir_tf_index = floor( ( ir_tf - t0 )/fpulse.dt )+1;
	
	
	cout << "\nIR-t0= "		<< ir_t0;	
	cout << "    IR-tf= "	<< ir_tf << "  Long-t= " << ir_tf-ir_t0 << endl<<endl;
	cout << "NtimeEvol= "	<< NtimeEvol <<endl;	
	
	
	complex P0, P1, P2, P3;
	
	//Preparing Crank-Nicholson arrrays    
	w.PrepareCrankArraysOnRho(complexdt    );
	w.PrepareCrankArraysOnZ(  complexdt/2. );
	//w.PrepareCrankArraysOnZLaserAV(  complexdt/2. );
	
	double adip		= 0.;
	double TotNorm  = 0.;

	//Start Loop propagation on time    
	for (int ktime=0; ktime<fpulse.Nmaxt; ktime++)
	{
		
		
		//w.parallel_evolution( complexdt);
		
		//Velocity Gauge
		/*
	     w.Zprop_PAG(complexdt/2., imag(fpulse.avector.f[ktime]) );
		 w.Rprop(complexdt);
		 w.Zprop_PAG(complexdt/2., imag(fpulse.avector.f[ktime]) );
		 //*/
		
		
		/*//Length Gauge
		w.Zprop(complexdt/2. );
		w.Rprop_REG( complexdt, imag(fpulse.efield.f[ktime]) );
		w.Zprop(complexdt/2. );
		
		//				w.Rprop(complexdt);*/
		
		
		
		//Velocity Gauge		
		//w.parallel_evolution_laser_PA( complexdt, imag(fpulse.avector.f[ktime]) );


		
		if(ktime==ir_t0_index || ktime==ir_tf_index)
		{			
			snaper1		= 4;
			snaper2		= 4;			
			
			largeNr		= 1500;
			largeNz		= 3000;	
			
			if(ktime==ir_tf_index)
			{
				largeNr=2000;
				largeNz=4000;
			};
			
			largeRmax  = largeNr*dr;			
			
			
			HH.resize_HankelMatrix(largeNr, largeRmax);
			cout << "\nNew size of the function \n ";
			cout << "\nHH-Nr= " << HH.Nr ;
			
			w.resize_waveUniform2D(HH, largeNz, dz );	
			
			cout << "\nNew dr= " << w.dr;
			
			w.set_potential_hlike2D(   1.,  0.);	
			
			
			w.PrepareCrankArraysOnRho(complexdt    );
			w.PrepareCrankArraysOnZ(  complexdt/2. );
			
			
	     	Energy1		= w.kinetic_energy_finite() + w.potential_energy();						
			cout << "\n\nCheking-Resize at ktime = " << ktime << "  Energy = " << Energy1 << "  Norm = " << w.norm() << "  Error in the norm = " << 1.-w.norm() << endl;
			
			
			w0.resize_waveUniform2D(HH, largeNz, dz );						
			w0.set_potential_hlike2D(    1.,  0.);
			
			wg0.resize_waveUniform2D(HH, largeNz, dz );						
			wg0.set_potential_hlike2D(    1.,  0.);	
			
			wg1.resize_waveUniform2D(HH, largeNz, dz );						
			wg1.set_potential_hlike2D(    1.,  0.);
			
			wg2.resize_waveUniform2D(HH, largeNz, dz );						
			wg2.set_potential_hlike2D(    1.,  0.);	
			
			wg3.resize_waveUniform2D(HH, largeNz, dz );						
			wg3.set_potential_hlike2D(    1.,  0.);	
			
			w1.resize_waveH2D(HH, largeNz, dz);
			w1.set_potential_hlike2D(	  1.,  0.);
			
			
			if(ktime==ir_t0_index)
			{
				w.saveAxes(B1axes,snaper1,snaper2);
				w1.saveQAxes(B1qaxes,snaper1,snaper2);	
			};
			
			if(ktime==ir_tf_index)
			{
				w.saveAxes(B2axes,snaper1,snaper2);
				w1.saveQAxes(B2qaxes,snaper1,snaper2);	
			};			
			
			cout << "\n........................\nThe redimension of the wavefunctions have been done succefully ....\n\n";
			
		};
		
		
		
		//Length Gauge
		w.parallel_evolution_laser_RE( complexdt, imag(fpulse.efield.f[ktime]) );
		w.absorber(0.0, 0.05, 0.05, 0.05, 1./6.);	   
		
		if(ktime%100==0)
		{
			
			P0= projection( wg0,  w );			
			P1= projection( wg1,  w );
			P2= projection( wg2,  w );			
			P3= projection( wg3,  w );
			
			TotNorm = w.norm();
			
			fprintf(out10,"%e %e %e %e %e %e %e\n", 
					 fpulse.g.t[ktime], imag(fpulse.efield.f[ktime]) 
					 ,norm(P0),norm(P1),norm(P2),norm(P3),TotNorm);
			
			//cout << "\n\nNorm= " << norm(P0) << "\n";
			
		};
		
		
		//***************************************
	    //     Write the error in the norm      
	    //*************************************** 		
	    if (ktime%snap==0 && ktime<ir_t0_index)
		{						
	     	Energy1 = w.kinetic_energy_finite()+w.potential_energy();						
			cout << "Loop number = " << ktime << "  Energy = " << Energy1 << "  Norm = " << w.norm() << "  Error in the norm = " << 1.-w.norm() << endl;
			
			
			
			w.mask_function2D( w0, masc0, masc1, sigma_masc);				
	       	out1 << ktime << "  " << 1.-w.norm() << " " << Energy1 << endl;			
			out2 << ktime << "  " <<  w.norm() <<"   "<<w0.norm()<<"   "<< 1.-abs(projection(w, wGround ))  << endl;
			
			
			
			//Snapshot full
			w.snapshot(out3, snaper1, snaper2);
			
			
			
			//Snapshot mask
			w0.snapshot(out4, snaper1, snaper2);
			
			

			//Making interpolation to creat Momentum distribution by Hankel Transform
			//interpU2H(w0,w1);
			//w1.qsnapshot(out5,HH, snaper1, snaper2);
			
		
			out6 << ktime <<"  " << fpulse.g.t[ktime] << " " << imag(fpulse.efield.f[ktime]) <<endl;			
			
			
			
		};
				
		
		if (ktime%500==0 && ktime>ir_t0_index &&  ktime<ir_tf_index )
			//if (ktime%5==0 && ktime>ir_t0_index &&  ktime<ir_tf_index )		
		{	
			
	     	Energy1		= w.kinetic_energy_finite() + w.potential_energy();						
			cout << "Loop number = " << ktime << "  Energy = " << Energy1 << "  Norm = " << w.norm() << "  Error in the norm = " << 1.-w.norm() << endl;
			
			
			w.saveAxes(B1axes,	  snaper1,snaper2);
			w1.saveQAxes(B1qaxes, snaper1,snaper2);	
			
			
			w.mask_function2D( w0, masc0, masc1, sigma_masc);				
	       	B1out1 << ktime << "  " << 1.-w.norm()	<< " "	<<	Energy1 << endl;			
			B1out2 << ktime << "  " <<  w.norm()	<<"   "	<<	w0.norm()<<"   "<< 1.-abs(projection(w, wGround ))  << endl;
			
			
			
			//Snapshot full
			w.snapshot(B1out3, snaper1, snaper2);
			
			
			
			//Snapshot mask
			w0.snapshot(B1out4, snaper1, snaper2);
			
			
			
			//Making interpolation to creat Momentum distribution by Hankel Transform
			//interpU2H(w0,w1);
			//w1.qsnapshot(out5,HH, snaper1, snaper2);
			
			
			B1out6 << ktime <<"  " << fpulse.g.t[ktime] << " " << imag(fpulse.efield.f[ktime]) <<endl;			
			
		};		
		
		
		
		if (ktime%snap==0 && ktime>ir_tf_index)
			//if (ktime%10==0 && ktime>ir_tf_index)		
		{						
	     	Energy1		= w.kinetic_energy_finite() + w.potential_energy();						
			cout << "Finals Loop number = " << ktime << "  Energy = " << Energy1 << "  Norm = " << w.norm() << "  Error in the norm = " << 1.-w.norm() << endl;
			
			
			w.mask_function2D( w0, masc0, masc1, sigma_masc);				
	       	B2out1 << ktime << "  " << 1.-w.norm() << " " << Energy1 << endl;			
			B2out2 << ktime << "  " <<  w.norm() <<"   "<<w0.norm()<<"   "<< 1.-abs(projection(w, wGround ))  << endl;
			
			
			
			//Snapshot full
			w.snapshot(B2out3, snaper1, snaper2);
			
			
			
			//Snapshot mask
			w0.snapshot(B2out4, snaper1, snaper2);
			
			
			
			//Making interpolation to creat Momentum distribution by Hankel Transform
			//interpU2H(w0,w1);
			//w1.qsnapshot(out5,HH, snaper1, snaper2);
			
			
			B2out6 << ktime <<"  " << fpulse.g.t[ktime] << " " << imag(fpulse.efield.f[ktime]) <<endl;			
			
		};		
		
		
		
		if (ktime%2==0) 
		{
			adip =dipole_acceleration( w);
			
			fprintf(out7,"%e %e %e \n",fpulse.g.t[ktime], imag(fpulse.efield.f[ktime]), adip);
		};
		//*/
		
	};//End loop propagation on time
	

	

	
/*	
	//Applying mask functions
	w.mask_function2D( w0, masc0, masc1, sigma_masc);		
	out1 << fpulse.Nmaxt << "  " << 1.-w.norm() << " " << Energy1 << endl;	
	out2 << fpulse.Nmaxt << "  " << 1.-w.norm() <<" "<< w.norm()-w0.norm()<<" "<< 1.-abs(projection(w, wGround ))  << endl;
	
	
	
	
	
	//Snapshot full
	w.snapshot(out3,snaper1,snaper2);
	
	
	
	//Snapshot mask
	w0.snapshot(out4,snaper1,snaper2);
	
	
	
	//Making interpolation to creat Momentum distribution by Hankel Transform
	interpU2H(w0,w1);
	w1.qsnapshot(out5,HH,snaper1,snaper2);	
	
	
	//*/
	//Saving binary data whole wavefunction
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
	
	
	cout << "\nEND OF THE PROGRAM\n"<<endl;
	
}

double dipole_acceleration( waveUniform2D &w )
{
	
	double adip		= 0.;
	double norm		= 0.;
	double rprime	= 0.;
	double edensity = 0.;
	double ac		= 0.;
	double rtemp	= 0.;
	double temp     = 0.;
	double VolElem  = 0.;
	
	for(int j=0;j<w.Nr;j++)
	{
		rtemp= w.r[j];
		for(int i=0;i<w.Nz;i++)
		{
			VolElem     = w.dz*w.dr*rtemp;
			rprime		= sqrt( rtemp*rtemp + w.z[i]*w.z[i] );
			edensity	= pow( abs( w.phi[w.index(j,i)] ), 2. );
			
			ac			= w.z[i]/pow(rprime,3.);
			
			
			adip+= VolElem*ac*edensity;
			norm+= VolElem*edensity;
		};
	};
	
	return double(adip/norm);
};

