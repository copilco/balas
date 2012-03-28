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
	
    bwave0 = fopen("/Users/alexis/Dropbox/BALAS/TestingUniformeZ1D/binwavel.bin","rb") ;
	
	
	
	
    //////////////////
    //  Parameters  //
    //////////////////    
    int Nz     = 150000;
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
    waveUniformZ1D wlarge;
    waveUniformZ1D wlarge1;	
    waveUniformZ1D wlarge2;		
	
	
    w.initialize(Nz,dz);
    wlarge.initialize(150000,dz);	
    wlarge1.initialize(150000,dz);
    wlarge2.initialize(150000,dz);	
	
	
	
	
    // Hydrogen parameters 
    double charge_nuclei = 1.;
    double soft_core     = 2.;	
	
	
	
	//Potential
    wlarge.set_potential_hlike1D( charge_nuclei, soft_core );	    
    //wlarge1.set_potential_hlike1D( charge_nuclei, soft_core );    
    
	//return 0;    
    //Initial test function	
    w.binread(bwave0);
    w.gaussian1D(z0, sigma );

	//Place wavefunction
	wlarge.placeWF(wlarge,w);    
    

    
	
	
	
	
	
    //Check the norm	
    cout << "\nOriginal norm= " << wlarge.norm() << endl;
    wlarge.normalize();
    cout << "  Initial norm= " << wlarge.norm() << endl;
    
	
	
	
	
	
	
    //*******************************************
    //  Save initial wavefunction and axis  //
    //*******************************************    
    for(int i=0;i<wlarge.Nz;i++){
		fprintf(out0,"%e \n", wlarge.z[i]);
		wlarge.phi[i]=-wlarge.z[i]*wlarge.phi[i];
    }
    
	
	
	
	
	
    //*******************************************
    //  Building Momentum Wave function   //
    //*******************************************    
	int kont = 1;
	int Nmom = 10000;
	
	double k    =  0.;
	double dk   =  0.001;
	double kmin = -dk*(Nmom-1)/2.;
    for (int kmom=0; kmom<Nmom; kmom++)
	{
		wlarge1.kmomenta= kmin+dk*kmom;
		wlarge1.Continuum_WF();
		
		complex proj=wlarge.projection(wlarge1,wlarge);
		
		fprintf(out1,"%e %e %e %e\n",wlarge1.kmomenta,real(proj),imag(proj),abs(proj));
		
		if (kmom%500==0) {
			cout << "\nkmom= " << kmom;
			cout << "    k= "  << wlarge1.kmomenta;	
			cout << "    'ScatteringNorm'= "  << wlarge1.norm();
			cout << "    Adip= "  << abs(proj);				
		}
		
	}//End propagation	
	



	fclose(out0);
	fclose(out1);
	fclose(out2);
	fclose(out3);	
	fclose(out4);
	fclose(out5);	
	
}
