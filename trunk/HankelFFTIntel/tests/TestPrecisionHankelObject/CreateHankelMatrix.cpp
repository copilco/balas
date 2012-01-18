//
//  mainTestPrecisionHankelObject.cpp
//  
//
//  Created by Alejandro de la Calle on 28/09/11.
//  
//
//

#include <iostream>
#include <math.h>
#include <complex>
#include "mkl.h"
#include "mkl_dfti.h"
#include "HankelMatrix.h"
#include "constant.h"


int main()
{
	cout << "\n\n//////////////////////////////////////////////////" << endl;
    cout << "//////////////////////////////////////////////////" << endl;
    cout << "mainTestPrecisionHankelObject. Running example..." << endl;
    cout << "//////////////////////////////////////////////////" << endl;
	cout << "//////////////////////////////////////////////////" << endl;
	
    // Number of points
	int Nr=1000;
    
    // Invoke HankelMatrix
	HankelMatrix HH(Nr,200.);
	
    // Get files with the object's data.
    HH.getAxisBin();
    
    HH.getHankelMatrixBin();
    

}