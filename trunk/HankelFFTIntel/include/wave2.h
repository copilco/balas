/*
 *  wave.h
 *  
 *
 *  Created by Camilo Ruiz Méndez on 29/11/11.
 *  Copyright 2011 USAL. All rights reserved.
 *
 */

class waveH
{
	
public:
	
	complex *phiHank;
	complex *F,*f;
	complex *F2,*f2;
	complex *Fretrieved;
	
	int Nr;
	int space_indicator;  //Defines in which space you are.
	double *r,*v,*dr,*dv;
	
	
	void initialize(Hankel HH)
	{
		Nr=HH.Nr;
		space_indicator=0;//We are in position rho space
		
		phiHank=new complex[Nr];
		
		F=new complex[Nr];
		f=new complex[Nr];
		
		F2=new complex[Nr];
		f2=new complex[Nr];
		
		Fretieved=new complex[Nr];
		
		r=new double[Nr];
		v=new double[Nr];
		dr=new double[Nr];
		dv=new double[Nr];
		
		
		for(int i=0;i<Nr;i++)
		{
			phiHank[i]=complex(0.,0.);
			
			F[i]=complex(0.,0.);
			f[i]=complex(0.,0.);
			
			F2[i]=complex(0.,0.);
			f2[i]=complex(0.,0.);
			
			Fretrieved[i]=complex(0.,0.);
		}
		
		for(int i=0;i<Nr;i++)
		{
			r[i]=HH.r[i];
			v[i]=HH.v[i];
			dr[i]=HH.dr[i];
			dv[i]=HH.dv[i];
		}
		
	}
	
	double norm()
	{
		double norm=0.0;
		for(int i=0;i<Nr;i++)
		{
			norm+=dr[i]*r[i]*real(conj(phiHank[i])*phiHank[i]);
		}
		return norm;
	}
	
	void normalize()
	{
		double norm1=norm();
		for(int i=0;i<Nr;i++)
		{
			phiHank[i]=phiHank[i]/sqrt(norm1);
		}
	}
	

}
