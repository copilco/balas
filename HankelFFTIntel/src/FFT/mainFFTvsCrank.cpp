/*
 *  mainFFTvsCrank.cpp
 *  
 *
 *  Created by Alejandro de la Calle on 02/12/11.
 *  
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "constant.h"

void tridager(complex *a00, complex *b, complex *c00, complex *r, complex *u,int n);
void nrerror (char error_text[]);


int main()
{
	
	cout << "\n\nmain FFT vs. Crank. Runnig example..." << endl;
	
	fstream out1("out1.txt",ios::out);
	fstream outFFT("outFFT.txt",ios::out);
	fstream outUniform("outUniform.txt",ios::out);
	
	int Nz=1000;

	/////////////////////////
	// Set axis parameters //
	/////////////////////////
	
	double dz=0.05;
	double sigma=10;
	double z0=0.;
	
	double dw=dospi/Nz/dz;
	double wmax=pi/dz;
	
	complex *phiFFT=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	complex *phiCrank=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	complex *phiCrank2=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	double *z=new double[Nz];
	double *w=new double[Nz];
	double *fase=new double[Nz];
	
	cout << "\nNumber of points: " << Nz << endl;
	cout << "dz: " << dz << endl;
	cout << "zmax: " << Nz/2.*dz << endl;
	cout << "Sigma: " << sigma << endl;
	cout << "z0: " << z0 << endl;
	cout << "dw: " << dw << endl;
	cout << "wmax: " << wmax << endl;
	
	//////////////////////////
	//  Set FFT parameters  //
	//////////////////////////
	
	MKL_LONG status;
	
	DFTI_DESCRIPTOR_HANDLE hand=0;
	double fscale = dz/sqrt(dospi);// 1.;//dz/sqrt(dospi);
	double bscale = dw/sqrt(dospi);//1./Nz;//dw/sqrt(dospi);
	
	
	status = DftiCreateDescriptor(&hand,DFTI_DOUBLE, DFTI_COMPLEX,1,(MKL_LONG)Nz);
	status = DftiSetValue(hand, DFTI_FORWARD_SCALE, fscale);
    status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, bscale);
	status = DftiCommitDescriptor(hand);
		
	
	
	///////////////////
	//  Create axes  //
	///////////////////
	
	for(int j=0;j<Nz;j++)
	{
		z[j]=(-Nz/2.+j)*dz;
	}
	
	for(int k=0;k<Nz/2;k++)
		w[k]=(k)*dw;
	
	for(int k=Nz/2;k<Nz;k++)
		w[k]=-wmax+(k-Nz/2)*dw;
	
	/////////////////////////
	// Create wavefunction //
	/////////////////////////
	
	for(int j=0;j<Nz;j++)
	{
		phiFFT[j]=exp(-(z[j]-z0)*(z[j]-z0)/sigma/sigma);
		phiCrank[j]=phiFFT[j];
		out1 << z[j] << "  " << phiFFT[j] << endl; 
		
	}
	
	double norm1=0.0;
	double norm2=0.0;
	for(int i=0;i<Nz;i++)
	{
		norm1+=dz*real(conj(phiFFT[i])*phiFFT[i]);
		norm2+=dz*real(conj(phiCrank[i])*phiCrank[i]);
	}
	
	for(int i=0;i<Nz;i++)
	{
		phiFFT[i]=phiFFT[i]/sqrt(norm1);
		phiCrank[i]=phiCrank[i]/sqrt(norm2);
	}
	
	norm1=0.0;
	norm2=0.0;
	for(int i=0;i<Nz;i++)
	{
		norm1+=dz*real(conj(phiFFT[i])*phiFFT[i]);
		norm2+=dz*real(conj(phiCrank[i])*phiCrank[i]);
	}
	
	cout << "\nInitial FFT normalization: " << norm1 << endl;
	cout << "Initial Crank normalization: " << norm2 << endl;
	
	
	//////////////////////////
	// Prepare Crank arrays //
	//////////////////////////
	
	complex *az=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	complex *cz=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	complex *vl=(complex*)mkl_malloc(Nz*sizeof(complex),16);
	complex *sl=(complex*)mkl_malloc(Nz*sizeof(complex),16);

	
	/////////////////////////
	// Start temporal loop //
	/////////////////////////
	
	double dt=0.01;
	
	cout << "\nStart temporal loop.\n";
	cout << "dt: " << dt << endl;
	
	for(int ktime=0;ktime<1;ktime++)
	{
		cout << "\n\nNumber of loop: " << ktime << endl;
		
		/////////////////////////
		//    FFT evolution    //
		/////////////////////////
		
		
		status=DftiComputeForward(hand,phiFFT);
		
		double qnorm=0;
		for (int i=0;i<Nz;i++)
		{
			qnorm+=dw*real(conj(phiFFT[i])*phiFFT[i]);//* (dz/sqrt(dospi))*(dz/sqrt(dospi)) ;
			//phiFFT[i]=phiFFT[i]*dw/sqrt(dospi);
		}
			
		cout << "\nError FFT normalization z: " << 1.-qnorm << endl;
		
		//for(int j=0;j<Nz;j++)
		//	cout << abs(phiFFT[j]) << endl;
		
		for (int i=0;i<Nz;i++)
		{
			fase[i]=dospi*dospi*w[i]*w[i]*dt/2.;
			phiFFT[i]=phiFFT[i]*exp(I*fase[i]);
		}
		
		
	
		status=DftiComputeBackward(hand,phiFFT);
		
		/*
		for (int i=0;i<Nz;i++)
		{
			phiFFT[i]=phiFFT[i]*dz/sqrt(dospi);
		}
		*/
		
		/////////////////////////
		//   Crank evolution   //
		/////////////////////////
		
		
		for(int i=1;i<Nz-1;i++)
        {	
			
            az[i]=complex(0.,-dt/4./dz/dz);
            cz[i]=complex(0.,-dt/4./dz/dz);
            
            vl[i]=complex(1.,dt/2./dz/dz);
            sl[i]=-az[i]*phiCrank[i-1]+(conj(vl[i]))*phiCrank[i]-cz[i]*phiCrank[i+1];
            
            
        }
		
		az[0]=complex(0.,-dt/4./dz/dz);
		cz[0]=complex(0.,-dt/4./dz/dz);
		
		vl[0]=complex(1.,dt/2./dz/dz);
		sl[0]=-az[0]*phiCrank[0]+(conj(vl[0]))*phiCrank[0]-cz[0]*phiCrank[1];
		
		
		az[Nz-1]=complex(0.,-dt/4./dz/dz);
		cz[Nz-1]=complex(0.,-dt/4./dz/dz);
		
		vl[Nz-1]=complex(1.,dt/2./dz/dz);
		sl[Nz-1]=-az[Nz-1]*phiCrank[Nz-2]+(conj(vl[Nz-1]))*phiCrank[Nz-1]-cz[Nz-1]*phiCrank[Nz-1];
		
		
		tridager(az,vl,cz,sl,phiCrank2,Nz);

		for(int i=0;i<Nz;i++)
        {
            phiCrank[i]=phiCrank2[i];
        }
		
		
		//////////////////////////////////
		// Load the data into the files //
		//////////////////////////////////
		
	
		if((ktime%10)==0)
		{
			for (int i=0; i<Nz; i++)
			{
				outFFT << z[i] << "  " << phiFFT[i] << endl; 
				outUniform << z[i] << "  " << phiCrank[i] << endl; 
			}
		}	
		
		
		//////////////////////////////////////
		// Check the error in normalization //
		//////////////////////////////////////
		
		
		norm1=0.0;
		norm2=0.0;
		
		//phiFFT[3]=complex(1.,0.);
		
		for(int i=0;i<Nz;i++)
		{
			norm1+=dz*real(conj(phiFFT[i])*phiFFT[i]);
			norm2+=dz*real(conj(phiCrank[i])*phiCrank[i]);
		}
		
		
		cout << "\nError FFT Después de volver normalization z: " << 1.-norm1 << endl;
		cout << "Error Crank normalization z: " << 1.-norm2 << endl;
		
		
	}
	
	
	
	
	status=DftiFreeDescriptor(&hand);
	
	out1.close();
	outFFT.close();
	outUniform.close();
	
	mkl_free(phiFFT);
	mkl_free(phiCrank);
	delete[] fase;
	mkl_free(az);
	mkl_free(cz);
	mkl_free(vl);
	mkl_free(sl);
	
}


void tridager(complex *a00, complex *b, complex *c00, complex *r, complex *u,int n)
/*
 Solves for a vector u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n],
 b[1..n], c[1..n], and r[1..n] are input vectors and are not modi?ed. */
{
	
	complex bet,*gam;
	gam=(complex *)mkl_malloc((n+2)*sizeof(complex),16);
	if(!gam) cout << "\nError de colocacion en la asignacion complex\n";
	
	if (b[1] == 0.0) nrerror("Error 1 in tridager **");
	
	u[1]=r[1]/(bet=b[1]);
	for (int j=2;j<n;j++)
    {
		gam[j]=c00[j-1]/bet;
		bet=b[j]-a00[j]*gam[j];
		//printf("%d %e %e\n",j,bet.real(),bet.imag());
		if (bet == 0.0) nrerror("Error 2 in tridager **"); //Algorithm fails; see below.
		u[j]=(r[j]-a00[j]*u[j-1])/bet;
    }
	for (int j=(n-2);j>=0;j--)
		u[j] -= gam[j+1]*u[j+1];
	mkl_free(gam);
}


void nrerror (char error_text[])
{
	
//	stderr << "Runtime error...\n";
//	stderr << error_text << endl;
//	stderr << "... now exting the system...\n";
	
	
	
	fprintf(stderr, "numerical recipes runtime error..\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr,"... now exiting system..\n");
	
 
}
