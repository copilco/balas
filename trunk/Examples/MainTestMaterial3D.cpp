#include <iostream>
//#include "vector.h"
#include "grid.h"
#include "wavef3.h"
#include "material.h"
#include "propagator.h"

using namespace std;

int main()
{
	
	/****************************/
	//These are grid parameters
	
	int nx=100;
	int ny=100;
	int nt=100;
	
	double dx=0.3;
	double dy=0.3;
	double dt=0.3;
	/****************************/
	
	grid g1;                                    //Define a grid object
	g1.set_grid(nx,ny,nt,dx,dy,dt);             //Inittialize  
	
	//Sone info
	
	cout << "Size in n1  "<< g1.n1 << "\n";     
	cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
	cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";
	
	/****************************/  
	//A wavefunction can be initialized using a grid as parameter.
	
	field w1; 
	w1.put_on_grid(g1);
	
	material h;
	h.put_on_grid(g1);
	h.set_v_one_over_rho();
	
	/*****************************************************/ 
	// Start the gaussian
	
	gaussian(w1);
	w1.normalize(); //Normalize the wavefunction
	/*****************************************************/
	//Do the FFT in 3D
	
	FILE *out4,*out5,*out6;
	
	out4=fopen("out4.txt","w");
	out5=fopen("out5.txt","w");
	out6=fopen("out6.txt","w");
	
	
	double factor=1./w1.n1/w1.n2/w1.n3;     
	double kinetic;
	double dz=0.02;
	
	int skiper1=5;
	int skiper2=5;
	int skiper3=5;
	
	for(int ktime=1;ktime<=25;ktime++)
	{
		//prop_dispersion_1(w1,h,dt);
		snapshot(out4,w1, skiper1, skiper2, skiper3);
		prop_dispersion_diffraction(w1,h,dt);
		
		fprintf(out6,"%d %e\n",ktime,w1.expected_x1());
		
		absorber(w1,0.1,0.1,0.1,0.1,0.1,0.,1./6.);
		cout << "\n Norm after the way on and way back "<<w1.norm();// << " norm with factor " << w1.norm()*factor*factor;  
		//w1.normalize();
	}
	
	
}
