/*
 *  laser.cpp laser pulses train differents characteristics by pulse
 *
 *  Created by Alexis Chacón and Camilo Ruiz
 *
 */

#ifndef LASER_H
#define LASER_H
#include <stdlib.h>
#include <string>
#include "timegrid.h"
#include "timeobject.h"
#include "constant.h"
#include<vector>
using namespace std;

class laser 
{

public:

	int Npulses;				//Total pulses number
	int Nmaxt;					//Total iteration number

	double t01;					//Initial time to first pulse
	double dt;					//Time step
	double blaser;				//Time before of the pulses
	double alaser;				//Time after of the pulses

	double major0;				//Major value
	double minus0;				//Minus value

	vector<double> clock0;      //Initial time per pulse
	vector<double> clock1;      //The time until half per pulse
	vector<double> clock2;		//The time until final per pulse
	vector<double> twidth;		//Time bandwidth per pulse
	
	vector<int> kclock0;		//k-th value to the initial time per pulse 
	vector<int> kclock1;		//k-th value to the time until half per pulse
	vector<int> kclock2;		//k-th value to the time until final per pulse
	
	
	vector<double> I0;			//Intensity per pulse
	vector<double> e;			//Ellipticity per pulse
	vector<double> E0x;			//Electric Field Component in the x-direction
	vector<double> E0y;			//Electric Field Component in the y-direction
	
	vector<double> w0;			//Frequency per pulse
	vector<double> period0;		//Period per pulse
	vector<double> cycles0;		//Cycles number per pulse
	
	vector<double> cep0;		//Carrier envelope phase per pulse
	vector<double> phi_rel;		//Relative phase between Ex and Ey components of the electric field of the laser pulse
	vector<double> delay0;		//Delay between consecutive pulses

	vector<string> envelope;	//Envelope name used
	
	timegrid g;				//Object of the time grid
	timeobject efield;		//Object to the electric field for total pulses
	timeobject avector;		//Object to the vector potential for total pulses

	timeobject *ef;			//Object to the electric field per pulse
	timeobject *env;			//Object to the envelope per pulse
	timeobject *av;			//Object to the vector potential per pulse
	
	timeobject *av_int;		//Object to the integral of the vector potential per pulse
	timeobject *avsq_int;	//Object to the square integral of the vector potential per pulse	
	
	/*==========================*/
	/*      MAIN FUNCTIONS
	/*==========================*/
	laser(int _Npulses);      							                                // Creator Object
	//~laser();																			//Destructor
	void laser_pulses(double _dt, double _t01, double _blaser, double _alaser);	        // LASER PULSES TRAIN
	
	/*=========================*/
	/*   SECUNDARY FUNCTIONS
	/*=========================*/
	inline void Initialize_Amplitude_Period();			// Initialize max amplitude and period
	inline void Set_StartTime_EndTime();		        // "Start" and "end" time per pulse
	void Laser_Grid();									// Object time axis and field_set
	void Evaluation_Laser_Pulse();						// Evaluation pulses train
	void Sum_Pulses();									// Pulses sum
	void Set_Vector_Potential();						// Vector potential per pulse
	void Vector_Potential();							// Total vector potential
	void put_on_envelope();								// Generator of Envelope of the pulses
	void Set_Av_Integral();
	
	double qmajor(vector<double>& v);					//Finding the maximum value of a vector
	double qminus(vector<double>& v);					//Finding the minimum value ot a vector
};


//==============================================================================//
			/*=== MAIN FUNCTIONS ===*/

/*===========================================================          		
		/*=== OBJECT'S CONSTRUCTOR LASER  ===*/
laser::laser(int _Npulses)
  {

	  
	  Npulses = _Npulses;
	  clock0.resize(Npulses,0.0);
	  clock1.resize(Npulses,0.0);	  
	  clock2.resize(Npulses,0.0);
	  
	  kclock0.resize(Npulses,0.0);
	  kclock1.resize(Npulses,0.0);
	  kclock2.resize(Npulses,0.0);
	  
	  twidth.resize(Npulses,0.0);
	  
	  I0.resize(Npulses,0.0);
	  e.resize(Npulses,0.0);	  
	  
	  envelope.resize(Npulses);
	  
	  for (int kpulse=0; kpulse<Npulses; kpulse++) 
		 envelope[kpulse] ="gauss";
	  
	  
	  
	  E0x.resize(Npulses,0.0);
	  E0y.resize(Npulses,0.0);	  
     
	  
	  
	  w0.resize(Npulses,0.0);
	  period0.resize(Npulses,0.0);
     
	  
	  
	  cycles0.resize(Npulses,0.0); 
	  cep0.resize(Npulses,0.0);	  
	  phi_rel.resize(Npulses,0.0);
	  
	  
	  
	  ef  = new timeobject[Npulses];     // reserve memory
	  env = new timeobject[Npulses];     // reserve memory
	  av  = new timeobject[Npulses];
     
	  
	  
	  av_int    = new timeobject[Npulses];		// reserve memory
	  avsq_int  = new timeobject[Npulses];		// reserve memory	  	  
	  
	  if (Npulses > 1 )
		  delay0.resize(Npulses-1,0.0);	  		  
  }//End initialize variable  



/*laser::~laser(){
	delete avsq_int;
	delete av_int;
	delete ef;
	delete env;
	delete av;
}*/


/*===========================================================          		
 /*===  LASER PULSES FUNCTION  ===*/
void laser::laser_pulses(double _dt, double _t01, double _blaser, double _alaser)
{
	//return;
   	/*======= Pulse's Parameter ======*/
   	t01    = _t01;										//	Start time first pulse
   	dt     = abs(_dt);			                        //	Time step 	
   	blaser = abs(_blaser);	                        //	Time before the laser or the pulse train
   	alaser = abs(_alaser);	                        //	Time after the laser or the pulse train 	
	
   	//================================//
	Initialize_Amplitude_Period();					//	Initialize max amplitude and period        
   	Set_StartTime_EndTime();	 	     	        //	"Start" and "end" time by each pulse
	
   	major0=qmajor(clock2);	     	                //	Maximum time to all pulse
   	minus0=qminus(clock0);          	            //	Minimun time to all pulse     
	
   	Laser_Grid();       	     	                //	Objet time axis      
   	Evaluation_Laser_Pulse();          	            //	Evaluation pulses
	
   	Sum_Pulses();	     	     	                //	Pulses sum (build pulses train)
   	Set_Vector_Potential();		     	            //	Vector potential by each pulse
   	Vector_Potential();								//	Vector potential to pulses train 

}//End laser_pulses

/*============ END MAIN FUNCTIONS ===============*/



/*=================================================================*/
		/*========= SECUNDARY FUNCTIONS ===========*/

//== Function initialize max amplitude and period ==// 
inline void laser::Initialize_Amplitude_Period()
{
	
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{
		double factor      =  1.0/(1.0+e[kpulse]*e[kpulse]);           //The ellipticity is e = E0x/E0y, Other condition
		E0x[kpulse]        =  e[kpulse]*sqrt(I0[kpulse]*factor/3.5e16);
		E0y[kpulse]        =  sqrt(I0[kpulse]*factor/3.5e16); 		
		period0[kpulse]    =  dospi/w0[kpulse];
		twidth[kpulse]	   =  cycles0[kpulse]*period0[kpulse];
	}
}


//== Function "start" and "end" time by each pulse ==//
inline void laser::Set_StartTime_EndTime()
{
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{		
		if(kpulse==0)
		{
			clock0[kpulse]   = t01;
			kclock0[kpulse]  = floor( (abs(clock0[kpulse]) + blaser)/dt) + 1;
			
			clock1[kpulse]   = clock0[kpulse] + twidth[kpulse]/2.0;
			kclock1[kpulse]  = floor( (abs(clock1[kpulse]) + blaser)/dt ) + 1;
			
			clock2[kpulse]   = clock0[kpulse] + twidth[kpulse];
			kclock2[kpulse]  = floor( (abs(clock2[kpulse]) + blaser)/dt ) + 1;
		    
		}
		else
		{

			clock0[kpulse]  = clock2[kpulse-1] + delay0[kpulse-1]
			                    - (twidth[kpulse-1] + twidth[kpulse])/2.0;
			
			kclock0[kpulse] = floor((abs(clock0[kpulse])+ blaser)/dt) + 1;
			
			
			clock1[kpulse]  = clock0[kpulse] + twidth[kpulse]/2.0;
			kclock1[kpulse] = floor( (abs(clock1[kpulse]) + blaser)/dt ) + 1;
			
			clock2[kpulse]  = clock0[kpulse] + twidth[kpulse];
			kclock2[kpulse] = floor( (abs(clock2[kpulse]) + blaser)/dt ) + 1;
			
		} 
	}//End loop
}//End function


//== Function object time axis and field set  ==//
void laser::Laser_Grid()   
{	

	Nmaxt=floor((major0+alaser-(minus0-blaser) )/dt)+1;
	if (Nmaxt%2!=0) 
		Nmaxt=Nmaxt+1;	

	
	g.set_grid(Nmaxt, dt, minus0-blaser);      
	efield.put_on_grid(g);
	avector.put_on_grid(g);
	
	for (int kfield=0; kfield <Npulses;kfield++){
		env[kfield].put_on_grid(g);
		ef[kfield].put_on_grid(g);
		av[kfield].put_on_grid(g);
	}
	cout << "Ntime= "<<Nmaxt<<endl;
	cout << "Memory= "<<((2+3*Npulses)*24+8)*Nmaxt*1e-6 << "  Mb"<<endl;
}


//== Funtion that evaluate of the pulses train ==//
void laser::Evaluation_Laser_Pulse()
{
	
	put_on_envelope();
	
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
			
			double arg1=w0[kpulse]*(g.t[ktime]-clock0[kpulse])+cep0[kpulse];                 // Phase x
			double arg2=w0[kpulse]*(g.t[ktime]-clock0[kpulse])+cep0[kpulse]+phi_rel[kpulse]; // Phase y
			
				
			real(ef[kpulse].f[ktime])  = real(env[kpulse].f[ktime])*sin(arg1);			   
			imag(ef[kpulse].f[ktime]) = imag(env[kpulse].f[ktime])*sin(arg2);
			
			
			real(av[kpulse].f[ktime]) = real(ef[kpulse].f[ktime] );
			imag(av[kpulse].f[ktime]) = imag(ef[kpulse].f[ktime] );
	   }
	}
}



void laser::put_on_envelope(){
	
	//string _name_env=envelope;
	vector<string> envelopes (5);
	
	envelopes[0]="rect";
	envelopes[1]="sin2";		
	envelopes[2]="rsin2";
	envelopes[3]="gauss";	
	envelopes[4]="konst";
	
	
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{	
	//Rectangle Envelope 
	if (envelope[kpulse]==envelopes[0]) {		
		

			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				if (g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse]){					
					real(env[kpulse].f[ktime]) = 0.0;
					imag(env[kpulse].f[ktime]) = 0.0;
				} 
				else{
					real(env[kpulse].f[ktime]) = E0x[kpulse];			   
					imag(env[kpulse].f[ktime]) = E0y[kpulse];			   
					
				}
				
			}
		}//End Rectangle Envelope
		
		

		//sin2 Envelope 
		if (envelope[kpulse]==envelopes[1]) 			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				double arg0=w0[kpulse]*(g.t[ktime]-clock0[kpulse])/2./cycles0[kpulse];            //Envelope Argument
				real(env[kpulse].f[ktime]) = E0x[kpulse]*sin(arg0)*sin(arg0);			   
				imag(env[kpulse].f[ktime]) = E0y[kpulse]*sin(arg0)*sin(arg0);			   
				
			}
		//End sin2 Envelope
		
		
		//Rectangle Multiplied sin2 Envelope (rsin2) 
		if (envelope[kpulse]==envelopes[2]) {		
			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				if (g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse]){					
					real(env[kpulse].f[ktime]) = 0.0;
					imag(env[kpulse].f[ktime]) = 0.0;
				} 
				else{
					double arg0=w0[kpulse]*(g.t[ktime]-clock0[kpulse])/2./cycles0[kpulse];            //Envelope Argument
					
					real(env[kpulse].f[ktime]) = E0x[kpulse]*sin(arg0)*sin(arg0);			   
					imag(env[kpulse].f[ktime]) = E0y[kpulse]*sin(arg0)*sin(arg0); 
					
				}
				
			}			
		}//End Rectangle  Multiplied sin2 Envelop (rsin2)
		

		
		//Gauss Envelope 
		if (envelope[kpulse]==envelopes[3])
			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				real(env[kpulse].f[ktime]) = E0x[kpulse]*
				                             exp( -64.*log10(2.)*( g.t[ktime] - (clock0[kpulse]+twidth[kpulse]/2. ) )*
					                         ( g.t[ktime] - (clock0[kpulse]+twidth[kpulse]/2. ) ) /twidth[kpulse]/twidth[kpulse] );		   
				
				imag(env[kpulse].f[ktime]) = E0y[kpulse]*
				                             exp( -64.*log10(2.)*(g.t[ktime]- (clock0[kpulse]+twidth[kpulse]/2.) )*
					                         ( g.t[ktime]-(clock0[kpulse]+twidth[kpulse]/2.) )/twidth[kpulse]/twidth[kpulse] );			   
				
			}//End Gauss Envelope	
		
		
		
		//Constant Envelope 
		if (envelope[kpulse]==envelopes[4]) {		
			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				real(env[kpulse].f[ktime]) = E0x[kpulse];			   
				imag(env[kpulse].f[ktime]) = E0y[kpulse];			   
				
			}
						
		}//End Constant Envelope		
	}//End number of pulses
	
	
}



//== Function pulses sum (build pulses train) ==//
void laser::Sum_Pulses()
{
	//Start loop sum pulses
	for(int ktime=0;ktime<Nmaxt;ktime++)
	{
		real(efield.f[ktime])=0.0;
		imag(efield.f[ktime])=0.0;		
		for(int kpulse=0;kpulse<Npulses;kpulse++)
		{
			real(efield.f[ktime])+= real(ef[kpulse].f[ktime]);
			imag(efield.f[ktime])+= imag(ef[kpulse].f[ktime]);

		}	
	}//End sum set pulses       
}



//== Function vector potential by each pulse ==//
void laser::Set_Vector_Potential()
{
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{
		av[kpulse].integrateRK4();//av[kpulse].integrate();      //
		double aa= real(av[kpulse].f[0]);
		double bb= imag(av[kpulse].f[0]);
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
			real(av[kpulse].f[ktime])-= aa; 
			imag(av[kpulse].f[ktime])-= bb;
			
			real(av[kpulse].f[ktime])*=-1; 	//According to Landau & Lifshitz in The Classical Theory of Fields,
			imag(av[kpulse].f[ktime])*=-1; 	//and Remetter Attosecond Wave Packet Interference, vector potential has negative sign.
			
		}  
	}
}





//== Function Integral of the vector potential by each pulse ==//
void laser::Set_Av_Integral()
{
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{		
		av_int[kpulse].put_on_grid(g);
		avsq_int[kpulse].put_on_grid(g);
		
		
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
			
			real(av_int[kpulse].f[ktime])=real(av[kpulse].f[ktime]); 
			imag(av_int[kpulse].f[ktime])=imag(av[kpulse].f[ktime]);
			
			
			real(avsq_int[kpulse].f[ktime]) = real(av[kpulse].f[ktime])*real(av[kpulse].f[ktime]); 
			imag(avsq_int[kpulse].f[ktime]) = imag(av[kpulse].f[ktime])*imag(av[kpulse].f[ktime]);
			
		}
		
		av_int[kpulse].integrateRK4();//av[kpulse].integrate(); 
		avsq_int[kpulse].integrateRK4();//av[kpulse].integrate(); 				
		//*/
	}
	
} //End integral of the vector potential





//== Function: vector potential to pulses train ==//
void laser::Vector_Potential()
{
	for(int ktime=0;ktime<Nmaxt;ktime++)
	{
		real(avector.f[ktime]) = 0.0;
		imag(avector.f[ktime]) = 0.0;		
		for(int kpulse=0;kpulse<Npulses;kpulse++)
		{
			real(avector.f[ktime])+= real(av[kpulse].f[ktime]);
			imag(avector.f[ktime])+= imag(av[kpulse].f[ktime]);
		}	
	}//End sum set pulses       
}



//Start function major
double laser::qmajor(vector<double>& v)
{
	int N=v.size();
	double may = v[0];
	
	for (int k=1;k<N;k++)   
		if (v[k]>may) may=v[k];
	
	return may;
}//End major

//Start function minus
double laser::qminus(vector<double>& v)
{
	int N=v.size();
	double min=v[0];
	
	for (int k=1;k<N;k++)
		if (v[k]<min) min=v[k];	
	
	return min;
}//End minus

/*========================== END  SECUNDARY FUNCTIONS =============================*/

#endif
//END
