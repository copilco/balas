#include "arrai.h"
#include <stdio.h>

int main()
{
  int N=100;
  arrai a(N,2.);

  for (int i=0;i<a.N;i++)
    printf("%e \n",a.v[i]);

}
