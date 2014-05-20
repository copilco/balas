//
//  mainTestImaginaryCrankUniProp2D.cpp
//  
//
//  Created by Alexis Chacón on 26/12/11.
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



int main(int argc, char *argv[])
{   

	
	
	/*fstream axis("axis.txt",ios::out);
    fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    fstream out3("out3.txt",ios::out);	
    //fstream outEne("outEne.txt",ios::out);
	*/
	
	
	fstream axes("axes.txt",ios::out);
    fstream qaxes("qaxes.txt",ios::out);
	fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    fstream out3("out3.txt",ios::out);    
    fstream out4("out4.txt",ios::out);
    fstream out5("out5.txt",ios::out);	
    fstream out6("out6.txt",ios::out);
	
    FILE *bwave0;
    bwave0 = fopen("binwave0.bin","wb") ;	
    
	FILE *bwave1;
    bwave1 = fopen("binwave1.bin","wb") ;		
    
	FILE *bwave2;
    bwave2 = fopen("binwave2.bin","wb") ;
    
	FILE *bwave3;
    bwave3 = fopen("binwave3.bin","wb") ;	
	
    FILE *bwave4;
    bwave4 = fopen("binwave4.bin","wb") ;
	
	
	
	
	
	//=================================//	
	// SPATIAL PARAMETERs 
	//.................................//  
    int Nr     = 300;
    int Nz     = 600;	
	
	
    double dz  = 0.3;
    double dr  = 0.3;	
	
	
	
    //Time step
	complex dt = -0.05*I;
    
	
	
	
    int Ntime=50000;
    int snap=Ntime/500;
    
    
	
	
	
	//cout.setf(ios::scientific);
	cout.precision(8);

	
	
	
    //Gaussian parameters
    double z0    = 1.5;	
	double r0     = 0.;
    double sigma = 4.;
    double Rmax  = Nr*dr;
	
	
	
    //Print out the information on the screen    
    cout << "\nNz= "       << Nz;
	cout << " Nr= "        << Nr;
	
    cout << "\ndz= "        << dz;
    cout << "  dr= "        << dr;
	
    cout << " dt= "        << dt;
    cout << " Ntime= "     << Ntime;
    cout << " Snapshots= " << snap;
    cout << " z0= "        << z0;
    cout << " Sigma= "     << sigma << endl;
    
    
	
	
	
    
    ////////////////////////
    //  Declare objects  //
    ///////////////////////
	HankelMatrix HGround(Nr, Rmax);	
	
	
	//Declare wavefunction 
    waveUniform2D w;
    waveUniform2D w1;
    waveUniform2D w2;
    waveUniform2D w3;	
    waveUniform2D w4;	
	
	
	//Declare Object Hankel wave
	waveH2D wq;
    
    w.initialize(  HGround, Nz, dz);
    w1.initialize(  HGround, Nz, dz);
    w2.initialize(  HGround, Nz, dz);	
    w3.initialize(  HGround, Nz, dz);
    w4.initialize(  HGround, Nz, dz);		
    wq.initialize(  HGround, Nz, dz);	
    
    
	
	
    //Initial test function	
    w.rho_gaussian(r0, z0, sigma, sigma );
    w1.rho_gaussian(r0, z0, sigma, sigma );	
    w2.rho_gaussian(r0, z0, sigma, sigma );	
    w3.rho_gaussian(r0, z0, sigma, sigma );	
    w4.rho_gaussian(r0, z0, sigma, sigma );		
    //w.binread(bwave0);
	
		
	
	//Saving potential well 
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			w3.phi[w3.index(j,i)] = w3.phi[w3.index(j,i)]*w3.z[i];	
	
	
	
    // Hydrogen parameters 
    double zmol2 = 0.7;
    double zmol1 = 3.-zmol2;
	
	
    double soft_core1     = 0.00001;
    double soft_core2     = 0.0001;	
    
	double R              = 2.2;
	double shift          = 0.0;
	
	
	// Hydrogen parameters     
	
//	set_potential_hlike2D(1.,0.);//
    w.set_potential_hlike2D(  1., 0.);
	
	w1.set_potential_hlike2D( 1., 0.);
	
	w2.set_potential_hlike2D( 1., 0.);
	w3.set_potential_hlike2D( 1., 0.);
	
	w4.set_potential_hlike2D( 1., 0.);	
	
	
	
	//Saving potential well 
	for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out2 << w.pot[w.index(j,i)] << endl;
	
	
	
	
    //Check the norm	
    cout << "\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "  Initial norm= " << w.norm() << endl;
    
   
	
	
	
	
	
    //*******************************************
    //  Save initial wavefunction and axis  //
    //*******************************************    
    w.saveAxes(axes,1,1);
    
      
    double ene1 = 0.;
    double ene2 = 0.;
    
	
	
	
	//w.PrepareCrankArrays( dt );
	w.PrepareCrankArraysOnRho( dt );
	w.PrepareCrankArraysOnZ( dt/2. );	
	
	
	
	//////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////	
	int kont		 = 0;
	int askiper		 = 1;
	double temp_flag = 5e-14;
	double ierror    = 0.;
	
    for (int ktime=0; ktime<5000; ktime++)
	{
		
		w.Zprop(  dt/2. );
		w.Rprop(  dt );
		w.Zprop(  dt/2. );
		
		
		w.normalize();	   
		
		
		/************************************
			Write the error in the norm  
		************************************/
		if(ktime%200 == 0)
		{
			  cout << "\nLoop number0= " << ktime;			
			  w.snapshot(out1, 1, 1);
			
			ene1   = w.kinetic_energy_finite()+w.potential_energy();
			ierror = abs(ene1-ene2);
			
			double static_dip = w.static_dipole2D();
			
			cout << "  Energy0=  " << ene1 << "  StaticDipole0= " << static_dip << "   Error= " << ierror << endl;	
			w.snapshot(out0,askiper, askiper);
			
			
			out5 << ene1 << endl;			
			
			if(ierror<=temp_flag)
			{
				kont+=1;
				if(kont==5)
					break;
			}
					
		}
		
		ene2=w.kinetic_energy_finite()+w.potential_energy();
		
		
	}//End propagation

	w.binwrite(bwave0 );
	
	out4 << ene2 << endl;	
	
	cout << "\n\nStating the Calculation of Bound w1 \n";

	//w1.PrepareCrankArrays( dt );
	w1.PrepareCrankArraysOnRho( dt );
	w1.PrepareCrankArraysOnZ(dt/2. );
	
	
	kont=0;
    for (int ktime=0; ktime<55000; ktime++)
	{
		
		w1.Zprop(  dt/2. );
		w1.Rprop(  dt );
		w1.Zprop(  dt/2. );
		
		

		project_out(w, w1 );		
		w1.normalize();	   
		
		
		/************************************
		 Write the error in the norm  
		 ************************************/
		if(ktime%1000 == 0){
			cout << "\nLoop number1= " << ktime;	
			
			ene1 = w1.kinetic_energy_finite()+w1.potential_energy();
			double static_dip = w1.static_dipole2D();
			
			ierror = abs(ene1-ene2);			
			
			cout << "  Energy1=  " << ene1 << "  StaticDipole1= " << static_dip << "   Error= " << ierror << endl;	
			w1.snapshot(out0,askiper, askiper );
			
			
			if(ierror<=temp_flag)
			{
				kont+=1;
				if(kont==5)
					break;
			}			
			
			
			
		}
		
		ene2=w1.kinetic_energy_finite()+w1.potential_energy();
		
		
	}//End propagation
	
	w1.binwrite(bwave1 );	

	cout << "\n\nStating the Calculation of Bound w2 \n ";
	out4 << ene2 << endl;		
	
	
	
	//w2.PrepareCrankArrays( dt );
	w2.PrepareCrankArraysOnRho( dt );
	w2.PrepareCrankArraysOnZ(dt/2. );	
	
	
	kont=0;
    for (int ktime=0; ktime<85000; ktime++)
	{
		
		w2.Zprop(  dt/2. );
		w2.Rprop(  dt );
		w2.Zprop(  dt/2. );
		
		
		project_out(w,  w2 );				
		project_out(w1, w2 );		
		w2.normalize();	   
		
		
		/************************************
		 Write the error in the norm  
		 ************************************/
		if(ktime%1000 == 0){
			cout << "\nLoop number2= " << ktime;	
			
			ene1 = w2.kinetic_energy_finite() + w2.potential_energy();
			ierror = abs(ene1-ene2);			
			double static_dip = w2.static_dipole2D();			
			cout << "  Energy2=  " << ene1 << "  StaticDipole2= " << static_dip << "   Error= " << ierror << endl;	
			w2.snapshot(out0, askiper, askiper );
			
			out5 << ene1 << endl;			
			
			if(ierror<=temp_flag)
			{
				kont+=1;
				if(kont==5)
					break;
			}			
			
		}
		
		ene2=w2.kinetic_energy_finite()+w2.potential_energy();
		
		
	}//End propagation
	
	w2.binwrite(bwave2 );
	
	out4 << ene2 << endl;	
	
	cout << "\n\nStating the Calculation of Bound w3\n ";
	
	/*
	w3.PrepareCrankArrays( dt );
	kont=0;
    for (int ktime=0; ktime<95000; ktime++)
	{
		
		w3.Zprop(  dt/2. );
		w3.Rprop(  dt );
		w3.Zprop(  dt/2. );
		
		
		project_out(w,  w3 );				
		project_out(w1, w3 );
		project_out(w2, w3 );		
		w3.normalize();	   
		
		
		/************************************
		 Write the error in the norm  
		 ************************************/
	
	/*
		if(ktime%1000 == 0)
		{
			cout << "\nLoop number3= " << ktime;	
			
			ene1 = w3.kinetic_energy_finite() + w3.potential_energy();
			ierror = abs(ene1-ene2);			
			double static_dip = w3.static_dipole2D();			
			cout << "  Energy3=  " << ene1 << "  StaticDipole3= " << static_dip << "   Error= " << ierror << endl;	
			w3.snapshot(out0, askiper, askiper );
			
			
			out5 << ene1 << endl;			
			
			if(ierror<=temp_flag)
			{
				kont+=1;
				if(kont==5)
					break;
			}			
			
		}
		
		ene2=w3.kinetic_energy_finite()+w3.potential_energy();
		
		
	}//End propagation
	
	w3.binwrite(bwave3 );
	out4 << ene2 << endl;		
	
	cout << "\n\nStating the Calculation of Bound w4\n ";
	
	
	w4.PrepareCrankArrays( dt );
	kont=0;
    for (int ktime=0; ktime<90000; ktime++)
	{
		
		w4.Zprop(  dt/2. );
		w4.Rprop(  dt );
		w4.Zprop(  dt/2. );
		
		
		project_out(w,  w4 );				
		project_out(w1, w4 );
		project_out(w2, w4 );	
		project_out(w3, w4 );		
		w4.normalize();	   
		
		
		/************************************
		 Write the error in the norm  
		 ************************************/
	/*
		if(ktime%5000 == 0){
			cout << "\nLoop number4= " << ktime;	
			
			ene1 = w4.kinetic_energy_finite() + w4.potential_energy();
			double static_dip = w4.static_dipole2D();
			
			cout << "  Energy4=  " << ene1 << "  StaticDipole4= " << static_dip << "   Error= " << ierror << endl;	
			w4.snapshot(out0, askiper, askiper );
			
			out5 << ene1 << endl;
			
			if(ierror<=temp_flag)
			{
				kont+=1;
				if(kont==5)
					break;
			}			
			
		}
		
		ene2=w4.kinetic_energy_finite()+w4.potential_energy();
		
		
	}//End propagation
	
	out4 << ene2 << endl;		
	w4.binwrite(bwave4 );	*/
	
	
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
