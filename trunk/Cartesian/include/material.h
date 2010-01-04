#include <iostream>
#include <vector.h>

class material: public grid
	{
		
	public:
		
		double kpp; //kprime 
		double k0;		
		
		double rho_c;
		double q2;		
		
		double index_n2; //index
		double x_kd;		
		
		double t_kd;
		double rho_n;
		
		double w_n;
		
		double a1;
		double b1;
		
		double a2;
		double b2;
		
		double a3;
		double b3;
		
		vector<double> v;  
		int something;
		
		
		/***************************************/
		// Methods.
		
		void initialize()
		{
			cout << "size 1 "<< n1;
			cout << "size 2 "<< n2;
			cout << "size 3 "<< n3;
			
			v.resize(n1*n2*n3, 0.);
			
			cout << " size v "<< n1*n2*n3;
			kpp=0.2; 
			k0=-1; 
			
			rho_c=1.;  
			index_n2=1;  
			
			x_kd=1.;  
			t_kd=1;  
			
			if(n1==1)
			{
				a1=0.;     //If the dimension is 1
				b1=0.;     //just a round zero for the operator
				
			}
			else
			{
				a1=kpp/2.;     //factor for the momentum
				b1=0.;     //factor for the momentum
			}
			
			if(n2==1)
			{
				a2=0.;     //If the dimension is 1
				b2=0.;     //just a round zero for the operator
			}
			else
			{
				a2=1./2./k0;     //factor for the momentum
				b2=0.;     //factor for the momentum
			}
			
			if(n3==1)
			{
				a3=0.;     //If the dimension is 1
				b3=0.;     //just a round zero for the operator
			}
			else
			{
				a3=1./2./k0;     //factor for the momentum
				b3=0.;     //factor for the momentum
			}
			
		}
		
		void put_on_grid(grid g)
		{
			
			set_grid(g.n1,g.n2,g.n3,g.dx1,g.dx2,g.dx3 );
			
			symmetry_n1=g.symmetry_n1;
			symmetry_n2=g.symmetry_n2;
			symmetry_n3=g.symmetry_n3;
			
			initialize();
		}
		
		/******************************************************/
				
		
		/******************************************************/
		// 
		/******************************************************/
		void set_v_one_over_rho()
		{
			cout << "V. size: " << (int) v.size() << endl;
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{ 
						//index(k,j,i);
						v[index(k,j,i)]= 1./(x2[j]*x2[j]+x3[k]*x3[k] ); //  1/(x2^2+x3^2)
					}
			
		}//End 3D pot
		
		/******************************************************/
		
		
		
	};
