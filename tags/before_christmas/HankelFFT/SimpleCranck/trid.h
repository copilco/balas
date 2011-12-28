#include <iostream.h>
#include <complex.h>
#include <fftw.h>
#include <alloc.h>
//#include <malloc.h>  //Para el ibook
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define complex complex<double>

void tridag(complex a,complex b[],complex c,complex r[],complex u[], int n)
{
    int j;
    complex bet, *gam;
    gam=(complex *)malloc((n+1)*sizeof(complex));
    if(!gam) printf("\nerror de colocacion en la asignacion complex");
    if(b[1]==0.0) printf("\nerror de colocacion en la asignacion complex"); //nrerror ("error 1 en tridag");
    u[1]=r[1]/(bet=b[1]);
    for(j=2;j<=n;j++)
    {
        gam[j]=c/bet;
        bet=b[j]-a*gam[j];
        if (bet==0.0) printf("\nerror de colocacion en la asignacion complex");//nrerror("error 2 en tridag");
        u[j]=(r[j]-a*u[j-1])/bet;
    }
    for(j=(n-1);j>=1;j--)
        u[j]-=gam[j+1]*u[j+1];
    free(gam);
}
