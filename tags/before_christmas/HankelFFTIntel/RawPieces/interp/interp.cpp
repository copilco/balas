/*
 *  interp.cpp
 *  
 *
 *  Created by Alejandro de la Calle on 18/11/11.
 *  
 *
 */
#include <iostream>
#include <complex>
#include <new>
#include <vector>
#include "HankelMatrix.h"
#include "constant.h"

#define complex complex<double>

using namespace std;

void tridager(complex *a00, complex *b, complex *c00, complex *r, complex *u,int n);
void nrerror (char error_text[]);

int main()
{
	FILE *in0, *out0, *out1;
	int Nr=1000;
	double *x, *xint;
	complex *y,*yint, *s, *u;
	
	in0=fopen("in0.txt","w+");
	out0=fopen("out0.txt","w+");
	out1=fopen("out1.txt","w+");
	
	
	x=new double[Nr];
	y=new complex[Nr];
	s=new complex[Nr];
	u=new complex[Nr-1];
	xint=new double[Nr];
	yint=new complex[Nr];
	
	double Rmax=100.;
	double dx=Rmax/((double)Nr);
	double dt=0.01;

	
	complex yp1=complex(0.0,0.0);
	complex ypn=complex(0.0,0.0);
	
	HankelMatrix HH(Nr,Rmax);
	
	for(int i=0;i<Nr;i++)
	{
		//x[i]=(-Nr/2+i)*dx;
		x[i]=HH.r[i];
		xint[i]=i*dx;
		y[i]=exp(-(x[i]-Rmax/2.)*(x[i]-Rmax/2.)/2./2.);
	}
	

	double norm=0.0;
	for(int i=0;i<Nr;i++)
	{
		norm+=HH.dr[i]*x[i]*real(conj(y[i])*y[i]);
	}
	
	for(int i=0;i<Nr;i++)
	{
		y[i]=y[i]/sqrt(norm);
	}
	
	norm=0.0;
	for(int i=0;i<Nr;i++)
	{
		norm+=HH.dr[i]*x[i]*real(conj(y[i])*y[i]);
	}
	
	printf("Original normalization = %e\n",norm);
	
	for(int i=0;i<Nr;i++)
	{
		fprintf(in0,"%e %e\n",x[i],abs(y[i]));
	}

	
	
	double sig;
	complex p,qn,un;
	
	
	if (abs(yp1) > 0.99e30)
		
		s[1]=u[1]=0.0;
	
	else {
		s[1] =complex(-0.5,-0.5);
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	
	for (int i=2;i<=Nr-1;i++)
	{	
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*s[i-1]+2.0;
		s[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	
	if (abs(ypn) > 0.99e30)
		qn=un=0.0;
	else 
	{
		qn=complex(0.5,0.5);
		un=(3.0/(x[Nr]-x[Nr-1]))*(ypn-(y[Nr]-y[Nr-1])/(x[Nr]-x[Nr-1]));
	}
	
	s[Nr]=(un-qn*u[Nr-1])/(qn*s[Nr-1]+1.0);
	
	for (int k=Nr-1;k>=1;k--)	
	{
		s[k]=s[k]*s[k+1]+u[k];
	}
	
	int klo,khi,k; 
	double h,b,a;
	
	
	for(int i=0;i<Nr;i++)
	{
		klo=1;
		khi=Nr; 
		
		while (khi-klo > 1)
		{
			k=(khi+klo) >> 1;
			if (x[k] > xint[i])
				khi=k;
			else klo=k;
		}
		
		h=x[khi]-x[klo];
		if (h == 0.0) 
			printf("Bad x input to routine, bro!");
		
		a=(x[khi]-xint[i])/h;
		b=(xint[i]-x[klo])/h;
		yint[i]=a*y[klo]+b*y[khi]+((a*a*a-a)*s[klo]+(b*b*b-b)*s[khi])*(h*h)/6.0;
		
	}
	
	for(int i=0;i<Nr;i++)
	{
		fprintf(out0,"%e %e\n",xint[i],abs(yint[i]));	
	}
	
	norm=0.0;
	for(int i=0;i<Nr;i++)
	{
		norm+=dx*xint[i]*real(conj(yint[i])*yint[i]);
	}
	
	printf("Interpolate normalization = %e\n",norm);
	
	// Now propagate
	complex *yyint, *az, *cz, *vl, *sl;
	yyint=new complex[Nr];
	az=new complex[Nr];
	cz=new complex[Nr];
	vl=new complex[Nr];
	sl=new complex[Nr];
	
	
	for(int ktime=0;ktime<5000;ktime++)
	{
		for(int i=1;i<Nr-1;i++)
        {	
			
            az[i]=complex(0.,-dt/4./dx/dx+dt/8./dx/xint[i]);
            cz[i]=complex(0.,-dt/4./dx/dx-dt/8./dx/xint[i]);
            
            vl[i]=complex(1.,dt/2./dx/dx);//+complex(0.,v[index(k,j,i)])*dt/6.+complex(0.,vr[k*(ndimzr+2)+j])*dt/4.;
            sl[i]=-az[i]*yint[i-1]+(conj(vl[i]))*yint[i]-cz[i]*yint[i+1];
            
            
        }
		
		az[0]=complex(0.,-dt/4./dx/dx+dt/8./dx/xint[0]);
		cz[0]=complex(0.,-dt/4./dx/dx-dt/8./dx/xint[0]);
		
		vl[0]=complex(1.,dt/2./dx/dx);//+complex(0.,v[index(k,j,i)])*dt/6.+complex(0.,vr[k*(ndimzr+2)+j])*dt/4.;
		sl[0]=-az[0]*yint[0]+(conj(vl[0]))*yint[0]-cz[0]*yint[1];
		
		
		az[Nr-1]=complex(0.,-dt/4./dx/dx+dt/8./dx/xint[Nr-1]);
		cz[Nr-1]=complex(0.,-dt/4./dx/dx-dt/8./dx/xint[Nr-1]);
		
		vl[Nr-1]=complex(1.,dt/2./dx/dx);//+complex(0.,v[index(k,j,i)])*dt/6.+complex(0.,vr[k*(ndimzr+2)+j])*dt/4.;
		sl[Nr-1]=-az[Nr-1]*yint[Nr-2]+(conj(vl[Nr-1]))*yint[Nr-1]-cz[Nr-1]*yint[Nr-1];
		
		
		tridager(az,vl,cz,sl,yyint,Nr);
		
		for(int i=0;i<Nr;i++)
        {
            fprintf(out1,"%d %e %e\n",ktime,x[i],abs(yyint[i]));
            yint[i]=yyint[i];
        }
		
		norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=dx*xint[i]*real(conj(yint[i])*yint[i]);
		}
		
		printf("Normalization = %e\n",norm);
		
		
	}
	
	
	
	
	
	delete[] x,y,xint,yint, s,u, yyint;
	delete[] az, cz, vl, sl;
	
	fclose(in0);
	fclose(out0);
	fclose(out1);
	
	
}



void tridager(complex *a00, complex *b, complex *c00, complex *r, complex *u,int n)
/*
 Solves for a vector u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n],
 b[1..n], c[1..n], and r[1..n] are input vectors and are not modi?ed. */
{
	//int j;
	complex bet,*gam2;
	gam2=(complex *)malloc((n+2)*sizeof(complex));
	if(!gam2) cout<<"\nerror de colocacion en la asignacion complex";
	//gam=vector(1,n); One vector of workspace, gam is needed.
	if (b[1] == 0.0) nrerror("Error 1 in tridager **");
	
	u[1]=r[1]/(bet=b[1]);
	for (int j=2;j<n;j++)
    {
		gam2[j]=c00[j-1]/bet;
		bet=b[j]-a00[j]*gam2[j];
		//printf("%d %e %e\n",j,bet.real(),bet.imag());
		if (bet == 0.0) nrerror("Error 2 in tridager **"); //Algorithm fails; see below.
		u[j]=(r[j]-a00[j]*u[j-1])/bet;
    }
	for (int j=(n-2);j>=0;j--)
		u[j] -= gam2[j+1]*u[j+1];
	free(gam2);
}


void nrerror (char error_text[])
{
	fprintf(stderr, "numerical recipes runtime error..\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr,"... now exiting system..\n");
}