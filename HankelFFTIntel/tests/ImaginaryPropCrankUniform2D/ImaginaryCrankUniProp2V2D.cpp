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
    
    
    //////////////////
    //  Parameters  //
    //////////////////    
    int Nr=250;
    int Nz=400;
    
    double dz  = 0.3;
    double dr  = 0.3;
    complex dt = complex(0.,-0.05);
    
    
    int Ntime=4000;
    int snap=Ntime/100;
    
    
    
    
    //Gaussian parameters
    double Rmax  = Nr*dr;
    double rho0  = 0.;
    double z0    = 0.;		
    double vr0   = 3.;
    double vz0   = 0.;	
    double sigma = 2.;
    
    
    
    
    
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
      	out0 << w.r[j] << endl;
    
    for(int j=0;j<Nz;j++)
       	out0 << w.z[j] << endl;
    
    for(int j=0;j<Nr;j++)
	for(int i=0;i<Nz;i++)
       	  out1 << abs(conj(w.phi[w.index(j,i)])*w.phi[w.index(j,i)])*w.r[j] << endl;
    
    
    
    ////////////////////////////
    //  Prepare Crank Arrays  //
    ////////////////////////////
    w.PrepareCrankArrays(dt);
    
    
    
    FILE *out3;
    FILE *out4;
    FILE *bwave;
    
    out3=fopen("out3.txt","w");		 
    out4=fopen("out4.txt","w");    
    bwave=fopen("binwave.bin","wb") ;
    
    
    double Energy      =-0.5;
    double ErrorConvergence = 0.;
    
    
    
    //start propagation
    for (int ktime=0; ktime<Ntime; ktime++)
      {
		  w.Zprop(  dt/2. );
		  w.Rprop(  dt );
		  w.Zprop(  dt/2. );
		  w.normalize();	   
	   
		  //************************************
		  //  Write the error in the norm  
		  //************************************
		  if (ktime%snap==0)
		  {
			  double Energy1 = w.kinetic_energy_finite()+w.potential_energy();
			  cout << "Both error  norm= "           << abs(1.-w.norm());
			  cout << "  Energy= "                   << Energy1;			
			  cout << "  ErrorTheoricalEnergy= "     << abs(Energy-Energy1);
			  cout << "  ErrorTheoricalEnergyConv= " << abs(Energy1-ErrorConvergence)<<endl;	
			  
			  out1 << ktime << "  " << abs(1.-w.norm()) << endl;
			  
			  for (int j=0; j<Nr; j++)
				  for(int i=0; i<Nz; i++)
					  fprintf(out3,"%e \n", w.r[j]*abs( w.phi[w.index(j,i)] )*abs( w.phi[w.index(j,i)] ) );
			  
			  fprintf(out4, "%e %e %e\n ", Energy1, abs(Energy1 - ErrorConvergence), abs(1.-w.norm()) );
		  }//
		  ErrorConvergence = w.kinetic_energy_finite()+w.potential_energy();

      }//End propagation
	
    w.binwrite(bwave);
	
    out0.close();
    out1.close();
    out2.close();
    //out3.close();
    //out4.close();
}
