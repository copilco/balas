#include <iostream>
#include "grid.h"

int main()
{
  
  int nx=300;
  int ny=1;
  int nz=3;
  
  double dx=0.1;
  double dy=1.;
  double dz=0.2;

  grid g1;
  
  g1.set_grid(nx,ny,nz,dx,dy,dz);
  
  cout << "Size in n1  "<< g1.n1 << "\n";
  cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
  cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";
  cout << "kfactor     "<< g1.kfactor2 << "\n";
  
}
