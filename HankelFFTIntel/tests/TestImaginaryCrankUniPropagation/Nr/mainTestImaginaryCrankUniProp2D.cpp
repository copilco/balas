//
//  mainTestImaginaryCrankUniProp2D.cpp
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



int main(int argc, char *argv[])
{   
   	
	cout << "\n\n//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "mainTestImaginaryCrankUniProp2D. Running example..." << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
	
	
	fstream axis("axis.txt",ios::out);
    fstream out0("out0.txt",ios::out);
    fstream out1("out1.txt",ios::out);
    fstream out2("out2.txt",ios::out);
    //fstream outEne("outEne.txt",ios::out);
	
    //FILE *out3;
    //FILE *out4;
    //FILE *bwave;
    
    //out3=fopen("out3.txt","w");		 
    //out4=fopen("out4.txt","w");    
    //bwave=fopen("binwave.bin","wb") ;
    

	FILE *outEne;
	outEne = fopen("outEne.txt","w+");
    
    //////////////////
    //  Parameters  //
    //////////////////    
    int Nr=50+20*atoi(argv[1]);//100;//250;
    int Nz=1500;//400;
    
    double dz  = 0.3;
    double dr  = 0.3;
    complex dt = complex(0.,-0.05);
    
    
    int Ntime=5000;
    int snap=Ntime/100;
    
    
    //Gaussian parameters
    double Rmax  = 150;//Nr*dr;
    double rho0  = 0.;
    double z0    = 0.;		
    double vr0   = 0.;
    double vz0   = 0.;	
    double sigma = 20.;
    
   
    //Print out the information on the screen    
    cout << "\nNr= " << Nr << endl;
    cout << "Nz= " << Nz << endl;
    cout << "dz= " << dz << endl;
    cout << "dt= " << dt << endl;
    cout << "Ntime= " << Ntime << endl;
    cout << "Snapshots= " << snap << endl;
    cout << "Rmax= " << Rmax << endl;
    cout << "rho0= " << rho0 << endl;
    cout << "z0= " << z0 << endl;
    cout << "v0r= " << vr0 << endl;
    cout << "v0z= " << vz0 << endl;
    cout << "Sigma= " << sigma << endl;
    
    
    
    ////////////////////////
    //  Declare objects  //
    ///////////////////////    
    HankelMatrix HH(Nr,Rmax);
    waveUniform2D w;
    w.initialize(HH,Nz,dz);
    
    
    
    
    //Initial test function	
    w.rho_gaussian(rho0, z0, sigma, sigma );
    
    
    
    // Hydrogen parameters 
    double charge_nuclei = 1.;
    double soft_core = 0.0073685;
    
    w.set_potential_hlike2D( charge_nuclei, soft_core );
    
    //Check the norm	
    cout << "\n\nOriginal norm= " << w.norm() << endl;
    w.normalize();
    cout << "Initial norm= " << w.norm() << endl;
    
   
    //*******************************************
    //  Save initial wavefunction and axis  //
    //*******************************************
    for(int j=0;j<Nr;j++)
      	axis << w.r[j] << endl;
    
    for(int j=0;j<Nz;j++)
       	axis << w.z[j] << endl;
    
    for(int j=0;j<Nr;j++)
		for(int i=0;i<Nz;i++)
			out0 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;
    
    
    
    ////////////////////////////
    //  Prepare Crank Arrays  //
    ////////////////////////////
    w.PrepareCrankArrays(dt);
      
    double ene1 = 0.;
    double ene2 = 0.;
    
	
	//////////////////////////////////////
    //  Start imaginary temporal loop  //
    /////////////////////////////////////

	
    for (int ktime=0; ktime<Ntime; ktime++)
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
		
		ene1 = w.kinetic_energy_finite()+w.potential_energy();
		
		cout << "Energy (Expected value):  " << ene1 << "   Error: " << log10(abs(ene1-ene2)) << endl;
		//cout << "Norm after Transform  " << w.norm() << endl;
		
		//outEne << ene1 << "  " << log10(abs(ene1-ene2)) << endl;
		
		ene2=ene1;
		
		/*
		if(ktime%(Ntime/(snap-1)) == 0)
			for (int j=0; j<Nr; j++)
				for(int i=0; i<Nz; i++)
					out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl		  
		*/			
		
	}//End propagation
	
	fprintf(outEne,"%d %12.5e %12.5e",Nr,ene1,log10(abs(ene1-ene2)));
    
	//w.binwrite(bwave);
	
    out0.close();
    out1.close();
    out2.close();
    //out3.close();
    //out4.close();
	//outEne.close();
	//bwave.close();
	fclose(outEne);
	
	
}
