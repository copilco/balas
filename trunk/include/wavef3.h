#include <iostream>
//#include <grid.h>
#include "fftw3.h"



class field: public grid
	{
	public:
		
		int space_indicator;
		
		fftw_complex *w;
		
		fftw_complex *w12;
		fftw_complex *w13;
		fftw_complex *w23;
		
		fftw_complex *w1;
		fftw_complex *w2;
		fftw_complex *w3;
		
		fftw_plan p3DF, p3DB;
		
		fftw_plan p2D_12F, p2D_12B;
		fftw_plan p2D_13F, p2D_13B;
		fftw_plan p2D_23F, p2D_23B;
		
		fftw_plan p1D_1F, p1D_1B;
		fftw_plan p1D_2F, p1D_2B;
		fftw_plan p1D_3F, p1D_3B;
		
		double qnorm;
		double expected_kin;
		
		int something;
		
		void initialize()
		{
			w = (fftw_complex*) fftw_malloc(n1*n2*n3 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			w1 = (fftw_complex*) fftw_malloc(n1 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			w2 = (fftw_complex*) fftw_malloc(n2 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			w3 = (fftw_complex*) fftw_malloc(n3 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			w12 = (fftw_complex*) fftw_malloc(n1*n2 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			w13 = (fftw_complex*) fftw_malloc(n1*n3 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			w23 = (fftw_complex*) fftw_malloc(n2*n3 * sizeof(fftw_complex));
			if(!w) cout<<"\nError allocating array, please check ...";
			
			p3DF  =fftw_plan_dft_3d(n3,n2,n1,w,w,-1, FFTW_ESTIMATE);
			p3DB  =fftw_plan_dft_3d(n3,n2,n1,w,w, 1, FFTW_ESTIMATE);
			
			p2D_12F  =fftw_plan_dft_2d(n2,n1,w12,w12,-1, FFTW_ESTIMATE);
			p2D_12B  =fftw_plan_dft_2d(n2,n1,w12,w12, 1, FFTW_ESTIMATE);
			
			p2D_13F  =fftw_plan_dft_2d(n3,n1,w13,w13,-1, FFTW_ESTIMATE);
			p2D_13B  =fftw_plan_dft_2d(n3,n1,w13,w13, 1, FFTW_ESTIMATE);
			
			p2D_23F  =fftw_plan_dft_2d(n3,n2,w23,w23,-1, FFTW_ESTIMATE);
			p2D_23B  =fftw_plan_dft_2d(n3,n2,w23,w23, 1, FFTW_ESTIMATE);
			
			p1D_1F  =fftw_plan_dft_1d(n1,w1,w1,-1, FFTW_ESTIMATE);
			p1D_1B  =fftw_plan_dft_1d(n1,w1,w1, 1, FFTW_ESTIMATE);
			
			p1D_2F  =fftw_plan_dft_1d(n2,w2,w2,-1, FFTW_ESTIMATE);
			p1D_2B  =fftw_plan_dft_1d(n2,w2,w2, 1, FFTW_ESTIMATE);
			
			p1D_3F  =fftw_plan_dft_1d(n3,w3,w3,-1, FFTW_ESTIMATE);
			p1D_3B  =fftw_plan_dft_1d(n3,w3,w3, 1, FFTW_ESTIMATE);
			
			
			qnorm=0.;
			expected_kin=0.;
			
			//p1DF_1 = fftw_plan_dft_1d(n1,w[0], w[0], FFTW_FORWARD, FFTW_ESTIMATE);
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
		// Functions of the wavefunction itself
		double norm()
		{
			/* Norm */
			double norm0=0.0;
			
			for(int i=0;i<n1*n2*n3;i++)
				norm0+=dx1*dx2*dx3*(w[i][0]*w[i][0]+w[i][1]*w[i][1]);
			
			return norm0;  
		}
		
		/******************************************************/
		//Normalize.
		
		void normalize()
		{
			double norm0=norm();
			//cout << "\n norm imput "<<norm0;
			for(int i=0;i<n1*n2*n3;i++)
			{
				w[i][0]/=sqrt(norm0);
				w[i][1]/=sqrt(norm0);
			}
			
			norm0=norm();
			//cout << "\n norm output "<<norm0;
			
		}
		
		/******************************************************/
		//Expected value x1
		double expected_x1()
		{
			double expec_x1=0.;
			
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{ 
						expec_x1+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x1[i]*w[index(k,j,i)][1]);
					}
			
			return expec_x1;
		}
		
		//Expected value x1
		double expected_x2()
		{
			double expec_x2=0.;
			
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{ 
						expec_x2+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x2[j]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x2[j]*w[index(k,j,i)][1]);
					}
			
			return expec_x2;
		}
		
		//Expected value x3
		double expected_x3()
		{
			double expec_x3=0.;
			
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{ 
						expec_x3+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x3[k]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x3[k]*w[index(k,j,i)][1]);
					}
			
			return expec_x3;
		}
		
		double q_norm()
		{
			
			qnorm=0.;
			double factor=1./ n1/ n2/ n3;
			
			/**********************/
			fftw_execute( p3DF);
			/**********************/
			
			for(int k=0;k< n3;k++)
				for(int j=0;j< n2;j++)
					for(int i=0;i< n1;i++)
					{
						qnorm+= dq1* dq2* dq3*
						( w[ index(k,j,i)][0]* w[ index(k,j,i)][0]+ w[ index(k,j,i)][1]* w[ index(k,j,i)][1])
						* kfactor1* kfactor1* kfactor2* kfactor2* kfactor3* kfactor3;
						
						w[ index(k,j,i)][0]*=factor;
						w[ index(k,j,i)][1]*=factor;
						
					}
			
			/**********************/
			fftw_execute( p3DB);
			/**********************/
			return qnorm;
		}
		
		void fft_over1_F()
		{
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
				{
					//Copy into the 1d array
					for(int i=0;i<n1;i++)
					{
						w1[i][0]=w[index(k,j,i)][0];
						w1[i][1]=w[index(k,j,i)][1];
					}
					
					//Do the fft1d with the factor
					fftw_execute(p1D_1F);
					
					for(int i=0;i<n1;i++)
					{
						w[index(k,j,i)][0]=w1[i][0];//*kfactor1;
						w[index(k,j,i)][1]=w1[i][1];//*kfactor1;
					}					
				}//end double loop
			
		}//end of fft_over1_F()
		
		void fft_over1_B()
		{
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
				{
					//Copy into the 1d array
					for(int i=0;i<n1;i++)
					{
						w1[i][0]=w[index(k,j,i)][0];
						w1[i][1]=w[index(k,j,i)][1];
					}
					
					//Do the fft1d with the factor
					fftw_execute(p1D_1B);
					
					for(int i=0;i<n1;i++)
					{
						w[index(k,j,i)][0]=w1[i][0];//*kfactor1;
						w[index(k,j,i)][1]=w1[i][1];//*kfactor1;
					}					
				}//end double loop
			
		}//end of fft_over1_F()
		
		
		void fft_over2_F()
		{
			for(int k=0;k<n3;k++)
				for(int i=0;i<n1;i++)
				{
					//Copy into the 1d array
					for(int j=0;j<n2;j++)
					{
						w2[j][0]=w[index(k,j,i)][0];
						w2[j][1]=w[index(k,j,i)][1];
					}
					
					//Do the fft1d with the factor
					fftw_execute(p1D_2F);
					
					for(int j=0;j<n2;j++)
					{
						w[index(k,j,i)][0]=w2[j][0];//*kfactor1;
						w[index(k,j,i)][1]=w2[j][1];//*kfactor1;
					}					
				}//end double loop
			
		}//end of fft_over1_F()
		
		void fft_over2_B()
		{
			for(int k=0;k<n3;k++)
				for(int i=0;i<n1;i++)
				{
					//Copy into the 1d array
					for(int j=0;j<n2;j++)
					{
						w2[j][0]=w[index(k,j,i)][0];
						w2[j][1]=w[index(k,j,i)][1];
					}
					
					//Do the fft1d without the factor
					fftw_execute(p1D_2B);
					
					for(int j=0;j<n2;j++)
					{
						w[index(k,j,i)][0]=w2[j][0];//*kfactor1;
						w[index(k,j,i)][1]=w2[j][1];//*kfactor1;
					}					
				}//end double loop
			
		}//end of fft_over2_B()
		
		
		void fft_over3_F()
		{
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{
					//Copy into the 1d array
					for(int k=0;k<n3;k++)
					{
						w3[k][0]=w[index(k,j,i)][0];
						w3[k][1]=w[index(k,j,i)][1];
					}
					
					//Do the fft1d without the factor
					fftw_execute(p1D_3F);
					
					for(int k=0;k<n3;k++)
					{
						w[index(k,j,i)][0]=w3[k][0];//*kfactor1;
						w[index(k,j,i)][1]=w3[k][1];//*kfactor1;
					}					
				}//end double loop
			
		}//end of fft_over3_F()
		
		void fft_over3_B()
		{
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{
					//Copy into the 1d array
					for(int k=0;k<n3;k++)
					{
						w3[k][0]=w[index(k,j,i)][0];
						w3[k][1]=w[index(k,j,i)][1];
					}
					
					//Do the fft1d without the factor
					fftw_execute(p1D_3B);
					
					for(int k=0;k<n3;k++)
					{
						w[index(k,j,i)][0]=w3[k][0];//*kfactor1;
						w[index(k,j,i)][1]=w3[k][1];//*kfactor1;
					}					
				}//end double loop
			
		}//end of fft_over3_B()
		
		
		void fft_over12_F()
		{

			for(int k=0;k<n3;k++)
			{
				
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{
						w12[in12(j,i)][0]=w[index(k,j,i)][0];
						w12[in12(j,i)][1]=w[index(k,j,i)][1];
					}
					
				//Do the fft1d without the factor
				fftw_execute(p2D_12F);
				
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{
						w[index(k,j,i)][0]=w12[in12(j,i)][0];
						w[index(k,j,i)][1]=w12[in12(j,i)][1];
					}					
				}//end double loop
			
		}//end of fft_over12_F()
		
		void fft_over12_B()
		{
			
			for(int k=0;k<n3;k++)
			{
				
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{
						w12[in12(j,i)][0]=w[index(k,j,i)][0];
						w12[in12(j,i)][1]=w[index(k,j,i)][1];
					}
				
				//Do the fft1d without the factor
				fftw_execute(p2D_12B);
				
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{
						w[index(k,j,i)][0]=w12[in12(j,i)][0];
						w[index(k,j,i)][1]=w12[in12(j,i)][1];
					}					
			}//end double loop
			
		}//end of fft_over12_F()
		
		
		void fft_over13_F()
		{
			
			for(int j=0;j<n2;j++)
			{
				
				for(int k=0;k<n3;k++)
					for(int i=0;i<n1;i++)
					{
						w13[in13(k,i)][0]=w[index(k,j,i)][0];
						w13[in13(k,i)][1]=w[index(k,j,i)][1];
					}
				
				//Do the fft2d without the factor
				fftw_execute(p2D_13F);
				
				for(int k=0;k<n3;k++)
					for(int i=0;i<n1;i++)
					{
						w[index(k,j,i)][0]=w13[in13(k,i)][0];
						w[index(k,j,i)][1]=w13[in13(k,i)][1];
					}					
			}//end double loop
			
		}//end of fft_over13_F()
		
		void fft_over13_B()
		{
			
			for(int j=0;j<n2;j++)
			{
				
				for(int k=0;k<n3;k++)
					for(int i=0;i<n1;i++)
					{
						w13[in13(k,i)][0]=w[index(k,j,i)][0];
						w13[in13(k,i)][1]=w[index(k,j,i)][1];
					}
				
				//Do the fft2d without the factor
				fftw_execute(p2D_13B);
				
				for(int k=0;k<n3;k++)
					for(int i=0;i<n1;i++)
					{
						w[index(k,j,i)][0]=w13[in13(k,i)][0];
						w[index(k,j,i)][1]=w13[in13(k,i)][1];
					}					
			}//end double loop
			
		}//end of fft_over13_F()
		
		void fft_over23_F()
		{
			
			for(int i=0;i<n1;i++)
			{
				
				for(int k=0;k<n3;k++)
					for(int j=0;j<n2;j++)
					{
						w23[in23(k,j)][0]=w[index(k,j,i)][0];
						w23[in23(k,j)][1]=w[index(k,j,i)][1];
					}
				
				//Do the fft2d without the factor
				fftw_execute(p2D_23F);
				
				for(int k=0;k<n3;k++)
					for(int j=0;j<n2;j++)
					{
						w[index(k,j,i)][0]=w23[in23(k,j)][0];
						w[index(k,j,i)][1]=w23[in23(k,j)][1];
					}					
			}//end double loop
			
		}//end of fft_over23_F()
		
		void fft_over23_B()
		{
			
			for(int i=0;i<n1;i++)
			{
				
				for(int k=0;k<n3;k++)
					for(int j=0;j<n2;j++)
					{
						w23[in23(k,j)][0]=w[index(k,j,i)][0];
						w23[in23(k,j)][1]=w[index(k,j,i)][1];
					}
				
				//Do the fft2d without the factor
				fftw_execute(p2D_23B);
				
				for(int k=0;k<n3;k++)
					for(int j=0;j<n2;j++)
					{
						w[index(k,j,i)][0]=w23[in23(k,j)][0];
						w[index(k,j,i)][1]=w23[in23(k,j)][1];
					}					
			}//end double loop
			
		}//end of fft_over23_F()
		
		
		
	};
