//Constants
#ifndef CONSTANT_H
#define CONSTANT_H
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <new>
#include <vector>
#include <string>
#include <complex>
//#define complex complex<double>
//#define MKL_Complex16 std::complex<double>
#define complex complex<double>

using namespace std;

//complex *RotToMol(complex *d,double phi, double tetha, double psi);
//complex *RotToLab(complex *d,double phi, double theta, double psi);

const complex I=complex(0.,1.);
const double lightC_au = 137.036;
const double one_by_lightC_au = 1./lightC_au;

const double pi = 3.141592653589793238463383;
const double dospi = 6.2831853071795862;
const double charge_e1ectron_au = -1.;
const double mass_proton_au=1836.;
const double BohrRadius = 5.29177211e-11;
const double time_SI = 2.418884326505e-17;   //Equivalencia de tiempo en segundos de una unidad atómica
const double lightC_SI = 2.99792458e8;

enum gauge_t { lengthgauge, velocitygauge, othergauge };

/*
complex *RotToMol(complex *d,double phi, double theta, double psi)
{
    // Rotate with direct Euler angles matrix R
    
    complex *dr=new complex[3];
    
    dr[0] = (cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi))*d[0] + (cos(psi)*sin(phi)-cos(theta)*cos(phi)*sin(psi))*d[1] + (sin(psi)*sin(theta))*d[2]; 
    dr[1] = (-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi))*d[0] + (-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi))*d[1] + (cos(psi)*sin(theta))*d[2]; 
    dr[2] = (sin(theta)*sin(phi))*d[0] - (sin(theta)*cos(phi))*d[1] + cos(theta)*d[2];
	
    return dr;
    
}


complex *RotToLab(complex *d,double phi, double theta, double psi)
{
    // Rotate with inverse Euler angles matrix R^-1
    
    complex *dr=new complex[3];
    
    dr[0] = (cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi))*d[0] + (sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi))*d[1] + (sin(phi)*sin(theta))*d[2]; 
    dr[1] = (cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi))*d[0] + (sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi))*d[1] + (cos(phi)*sin(theta))*d[2]; 
    dr[2] = (sin(theta)*sin(psi))*d[0] + (sin(theta)*cos(phi))*d[1] + cos(theta)*d[2];
	
    return dr;
    
}
*/
#endif  /* CONSTANT_H */


