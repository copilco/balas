#include "cylindrical.h"
vector<complex> gam1;
vector<complex> bz;
vector<complex> rv_z;
vector<complex> phi_z;

vector<complex> ar;
vector<complex> br;
vector<complex> cr;
vector<complex> rv_r;
vector<complex> phi_r;


int index(int j, int i, int nr)
{	            //k         j    i
	//int index=k*nz*nr + j*nr + i;
	int index = j*nr + i;	
	return index;
}


double cnorm( vector<complex> &phi, vector<double> &rho, int _nr, int _nz, double _dr, double _dz ){
	
	int nr = _nr;
	int nz = _nz;	
	
	double dr = _dr;
	double dz = _dz;	
	
	double norm0 = 0.;
	
	for(int i=0; i< nr; i++)	
		for(int j=0; j<nz; j++)
			norm0+= dr*dz*rho[i]*abs(phi[index(i,j,nz)])*abs(phi[index(i,j,nz)]);
	
	return norm0;
}


void cnormalize( vector<complex> &phi, vector<double> &rho, int _nr, int _nz, double _dr, double _dz )
{
	int nr =_nr;
	int nz =_nz;	
	
	double dr = _dr;
	double dz = _dz;
	
	double norm0 = 0.;
	
	norm0 = cnorm( phi, rho, nr, nz, dr, dz );
	
	for(int i=0; i< nr*nz; i++ )
		phi[i] = phi[i]/sqrt(norm0);
	
} 



void trid_simple( complex a, vector<complex> &b, complex c, vector<complex> &r, vector<complex> &u, int n)
{
	int j;
	gam1.resize( n, 0. );
	
	complex bet=b[0];
	u[0]=r[0]/bet;
	
	for(j=1;j<n;j++)
	{
		gam1[j] = c/bet;
		bet     = b[j] - a*gam1[j];
		
		
		u[j]  = (r[j]-a*u[j-1])/bet;
	}
	
	for(j=(n-2);j>=0;j--)
		u[j]-= gam1[j+1]*u[j+1];
	
}


void tridag( vector<complex> &a, vector<complex> &b, vector<complex> &c, vector<complex> &r, vector<complex> &u, int n)
{
	int j;
	gam1.resize( n, 0. );
	
	complex bet = b[0];
	u[0]=r[0]/bet;
	
	for(j=1;j<n;j++)
	{
		gam1[j] = c[j-1]/bet;
		bet     = b[j] - a[j]*gam1[j];
		
		u[j]  = (r[j]-a[j]*u[j-1])/bet;
	}
	
	for(j=(n-2);j>=0;j--)
		u[j]-=gam1[j+1]*u[j+1];
	
}



//*******************************
//===============================//
//Z_operator 1/2*dt


void Zprop( vector<complex> &phi, vector<double> &v, complex _dt, int _nz, int _nr, double _dz )
{
	int nz     = _nz;
	int nr     = _nr;
	
	complex dt  = _dt;
	double  dz  = _dz;
	
	complex az = complex( 0., -1./dz/dz/4. )*dt;
	complex cz = complex( 0., -1./dz/dz/4. )*dt;

	complex a1 =complex(1.,0.);
	complex a2 =complex(0.,0.);		
	
	bz.resize(    nz, 0. );
	rv_z.resize(  nz, 0. );
	phi_z.resize( nz, 0. );
	
	
	//*******************************
	//===============================//
	//Z_operator 1/2*dt
	for (int i=0; i<nr; i++ ){
		
		//Left part Z
		for(int j=0; j<nz; j++)	{
			a2= complex(0.,1./2./dz/dz  + v[index(i,j,nz)]/4.) ;
			bz[j]	=	a1 + a2*dt;
		}
		
		
		//Right part Z
		rv_z[0]		=   conj(bz[0])*phi[index(i,0,nz)]   + 
						conj( cz  )*phi[index(i,1,nz)];	
		
		
		for (int j=1; j<nz-1; j++) 
			rv_z[j] =	conj( az    )*phi[index(i,j-1,nz)] + 
						conj( bz[j] )*phi[index(i,j,nz)]   + 
						conj( cz    )*phi[index(i,j+1,nz)];
		
		
		rv_z[nz-1]	=	conj( az  )*phi[index(i,nz-2,nz)]  +
						conj( bz[index(i,nz-1,nz)]  )*phi[index(i,nz-1,nz)];	
		//Finishing right part Z
		
		
		
		//Zeros on Z
		zeros(phi_z);		
		
		//Solving Triagonal Matrix	for Z
		trid_simple(az,bz,cz,rv_z,phi_z,nz);
		
		//Save function 
		for (int j=0; j<nz; j++)
			phi[index(i,j,nz)] =	phi_z[j];		//psi0[index(i,j,nz)] = psi0[index(i,j,nz)];//
	}	
}
//******************************************// End Z_propagator
//==========================================//


//**********************************
//==================================//
//RHO_Operator dt	Complete Version 

void Rprop( vector<complex> &phi, vector<double> &v, vector<double> &rho, complex _dt, int _nz, int _nr, double _dr )
{	
	int nz      = _nz;
	int nr      = _nr;
	
	complex dt  = _dt;
	double dr   = _dr;
	
	ar.resize(    nr, 0. );
	br.resize(    nr, 0. );
	cr.resize(    nr, 0. );
	
	rv_r.resize(  nr, 0. );
	phi_r.resize( nr, 0. );
	
	complex a1 = complex(1.,0.);
	complex a2 = complex(0.,0.);	
	
	
	//**********************************
	//**********************************
	//==================================//
	//RHO_Operator dt	
	for ( int j=0; j<nz; j++ ){
		
		//Left part Rho
		for( int i=0; i<nr; i++ ){						
			ar[i]     =   (complex( 0., - 1./4./dr/dr  ) +
						  complex( 0., + 1./8./dr/rho[i]  ))*dt;

			a2		  =   complex(0., 1./2./dr/dr  + v[index(i,j,nz)]/4.)*dt ;
			br[i]	  =   a1 + a2;
			
			cr[i]     =   (complex( 0., - 1./4./dr/dr  ) +
						  complex( 0., - 1./8./dr/rho[i]  ))*dt;
		}//End left part
		
		
		
		//Right part Rho 
		complex phi0 = phi[index(1,j,nz)];//complex(0.,0.);//
		
		rv_r[0]		 =	conj( ar[0]  )*phi0 +
						conj( br[0]  )*phi[index(0,j,nz)]   + 
						conj( cr[0]  )*phi[index(1,j,nz)];			
		
		for (int i=1; i<nr-1; i++) 
			rv_r[i]	=	conj( ar[i]  )*phi[index(i-1,j,nz)] + 
						conj( br[i]  )*phi[index(i,j,nz)]   + 
						conj( cr[i]  )*phi[index(i+1,j,nz)];
		
		rv_r[nr-1]	 =	conj( ar[nr-1]  )*phi[index(nr-2,j,nz)]  +
						conj( br[index(nr-1, j, nz)]  )*phi[index(nr-1, j, nz)];
		//Finishing right
		
		
		
		//Zeros
		zeros(phi_r);
		
		
		//Solving Triagonal Matrix			
		tridag(ar, br, cr, rv_r, phi_r, nr);
		
		
		for (int i=0; i<nr; i++)
			phi[index(i,j,nz)] =	phi_r[i];//phi[index(i,j,nz)];
	}
}
//******************************************// End Rho_propagator
//==========================================//


void zeros(vector<complex> &phi){
	int n=phi.size();
	for (int i=0; i<n; i++)
		phi[i]=complex(0.,0.);
	
}