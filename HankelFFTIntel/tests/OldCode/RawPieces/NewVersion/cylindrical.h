#include <iostream>
#include <math.h>
#include <complex.h>
#include <vector.h>
#define complex complex<double>
using namespace std;

double cnorm( vector<complex> &phi, vector<double> &rho, int _nr, int _nz, double _dr, double _dz );
void cnormalize( vector<complex> &phi, vector<double> &rho, int _nr, int _nz, double _dr, double _dz ); 
void Zprop( vector<complex> &phi, vector<double> &v, complex _dt, int _nz, int _nr, double _dz );
void Rprop( vector<complex> &phi, vector<double> &v, vector<double> &rho, complex _dt, int _nz, int _nr, double _dr );
void tridag( vector<complex> &a, vector<complex> &b, vector<complex> &c, vector<complex> &r, vector<complex> &u, int n);
void trid_simple( complex a, vector<complex> &b, complex c, vector<complex> &r, vector<complex> &u, int n);
int index(int j, int i, int nr);
void zeros(vector<complex> &phi);