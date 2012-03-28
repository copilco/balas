//
//  mainTestImaginaryCrankUniProp2D.cpp
//  
//
//  Created by Alexis Chacón on 26/12/11.
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



int main(int argc, char *argv[])
{   

	
	
	/*fstream axis("axis.txt",ios::out);
    fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    fstream out3("out3.txt",ios::out);	
    //fstream outEne("outEne.txt",ios::out);
	*/
	
	
    FILE *out0;
    FILE *out1;	
    FILE *out2;
    FILE *out3;
    FILE *out4;
    FILE *out5;	
	
    FILE *bwave0;
    FILE *bwave1;
    
	
	
	
    out0  = fopen("out0.txt","w");
    out1  = fopen("out1.txt","w");
    out2  = fopen("out2.txt","w");
    out3  = fopen("out3.txt","w");    
    out4  = fopen("out4.txt","w");
    out5  = fopen("out5.txt","w");    	
	
    bwave0 = fopen("/Users/alexis/Dropbox/BALAS/TestingUniformeZ1D/binwave.bin","rb") ;
    bwave1 = fopen("binwavel.bin","wb") ;	
	
	
	
	
    //////////////////
    //  Parameters  //
    //////////////////    
    int Nz     = 15000;
    double dz  = 0.05;

    //Time step
	complex dt = -0.0025*I;
    
    
    int Ntime=50000;
    int snap=Ntime/500;
    
    
	
	
	
    //Gaussian parameters
    double z0    = 0.;		
    double sigma = 10.;
    
   
    //Print out the information on the screen    
    cout << "\nNz= "       << Nz;
    cout << " dz= "        << dz;
    cout << " dt= "        << dt;
    cout << " Ntime= "     << Ntime;
    cout << " Snapshots= " << snap;
    cout << " z0= "        << z0;
    cout << " Sigma= "     << sigma << endl;
    
    
	
	
	
    
    ////////////////////////
    //  Declare objects  //
    ///////////////////////    
    waveUniformZ1D w;
    w.initialize(Nz,dz);
    
    
    
    
    //Initial test function	
    //w.gaussian1D(z0, sigma );
    w.binread(bwave0);
    
    
	
	
    // Hydrogen parameters 
    double charge_nuclei = 1.;
    double soft_core     = 2.;
    
    w.set_potential_hlike1D( charge_nuclei, soft_core );
    
	
	
	
	
	
    //Check the norm	
    cout << "\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "  Initial norm= " << w.norm() << endl;
    
   
	
	
	
	
	
    //*******************************************
    //  Save initial wavefunction and axis  //
    //*******************************************    
    for(int i=0;i<Nz;i++)
		fprintf(out0,"%e \n", w.z[i]);
    
      
    double ene1 = 0.;
    double ene2 = 0.;
    
	
	
	
	
	//////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////	
    for (int ktime=0; ktime<15000; ktime++)
	{
		
		w.Zprop( dt );
		w.normalize();	   
		
		
		/************************************
			Write the error in the norm  
		************************************/
		if(ktime%200 == 0){
			cout << "\nLoop number= " << ktime;			
			w.snapshot(out1, 1);
			
			ene1 = w.kinetic_energy_finite()+w.potential_energy();
			
			cout << "  Energy:  " << ene1 << "   Error: " << log10(abs(ene1-ene2)) << endl;	
					
			}
		
		ene2=w.kinetic_energy_finite()+w.potential_energy();
		
		
	}//End propagation
	
	
	
	
	
	
	
    waveUniformZ1D wlarge;
    wlarge.initialize(150000,dz);	
	
	wlarge.placeWF(wlarge,w);
    wlarge.set_potential_hlike1D( charge_nuclei, soft_core );	

    
	for(int i=0;i<wlarge.Nz;i++)
		fprintf(out3,"%e %e\n", wlarge.z[i],wlarge.pot[i]);		
	
	
	
	
	
	
	
	//////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////	
	int kont=1;
    for (int ktime=0; ktime<15000; ktime++)
	{
		
		wlarge.Zprop( dt );
		wlarge.normalize();	   
		
		
		/************************************
		 Write the error in the norm  
		 ************************************/
		if(ktime%200 == 0){
			cout << "\nLoop number= " << ktime;			
			wlarge.snapshot(out2, 1);
			
			ene1 = wlarge.kinetic_energy_finite()+wlarge.potential_energy();
			
			cout <<" kont= " <<kont << "  Energy:  " << ene1 << "   Error: " << log10(abs(ene1-ene2)) << endl;	
			
			
			if (abs(ene1-ene2)<5.0e-15)
				kont+=1;
			
			if (kont>=10) 
				break;
			
		}
		
		ene2=wlarge.kinetic_energy_finite()+wlarge.potential_energy();
		
		
	}//End propagation	
	
	
	
	
	
	
	fprintf(out4,"%12.5e %12.5e %12.5e",dz,ene1,log10(abs(ene1-ene2)));
	wlarge.binwrite(bwave1);	
	
	
	
	fclose(out0);
	fclose(out1);
	fclose(out2);
	fclose(out3);	
	fclose(out4);
	fclose(out5);	
	
}
