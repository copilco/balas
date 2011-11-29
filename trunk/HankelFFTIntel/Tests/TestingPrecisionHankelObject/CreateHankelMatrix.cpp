/*
 *  Hankel1D.cpp
 *  
 *
 *  Created by Alejandro de la Calle on 28/09/11.
 *  
 *
 */

#include <iostream>
#include <math.h>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"
#include "constant.h"
#include "arrai.h"
#include "HankelMatrix.h"


int main()
{
	// Parameters!!
		
	int Nr=1000;
	HankelMatrix HH(Nr,200.);
	

}