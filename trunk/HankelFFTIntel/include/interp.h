//
//  interp.h
//  
//
//  Created by de la Calle Negro Alejandro on 26/12/11.
//  
//

#ifndef INTERP_H
#define INTERP_H

#include <iostream>
#include <math.h>
#include "mkl.h"
#include "arrai.h"
#include "HankelMatrix.h"
#include "constant.h"
#include "tools.h"
#include "waveH2D.h"
#include "waveUniform2D.h"

void interpH2U(waveH2D &wHank, waveUniform2D &wU);
void interpU2H(waveUniform2D &wU, waveH2D &wHank);
void CheckDfError(int num);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////   Two functions to interpolate arrays from wave and waveUniform class   ////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void interpH2U(waveH2D &wHank, waveUniform2D &wU)
{
    
    /////////////////////////////////////////////////
	// Declare parameters for interpolate routines //
	/////////////////////////////////////////////////
    
    int Nr = wHank.Nr;
    int Nz = wHank.Nz; 
    
    DFTaskPtr taskReal;                     // Data Fitting task descriptor
    DFTaskPtr taskImag; 
    DFTaskPtr celltask;
    double *auxInReal = new double[Nr*Nz];        // Auxiliary arrays
    double *auxInImag = new double[Nr*Nz];
    double *auxOutReal = new double[Nr*Nz];
    double *auxOutImag = new double[Nr*Nz];
    MKL_INT *cells = new MKL_INT[wU.Nr];
    double *sites = new double[2];
    double *scoeffReal;
	double *scoeffImag;
	
    
	int errcode = 0;
    
	////// Initializing parameters for Data Fitting task //////
    
    MKL_INT sorder = DF_PP_CUBIC;                     // spline order
    MKL_INT stype = DF_PP_NATURAL;                    // spline type
    
    ///// Parameters describing interpolation interval //////
    
    MKL_INT nx = wHank.Nr;                           // number of break points
    MKL_INT xhint = DF_NON_UNIFORM_PARTITION;                  // additional info about break points
    
    
    ////// Parameters describing function //////
    MKL_INT ny = Nz;                      // number of functions
    MKL_INT yhint = DF_NO_HINT;                      // additional info about function
    
    
    ////// Parameters describing spline coefficients storage //////
    
    MKL_INT nscoeff = ny*(nx-1)*sorder;                    // number of spline coefficients
    MKL_INT scoeffhint = DF_NO_HINT;                       // additional info about spline coefficients
    scoeffReal = new double[nscoeff];                                // array of spline coefficients
    scoeffImag = new double[nscoeff];    
	
    ////// Parameters describing boundary conditions type //////
    
    MKL_INT bc_type = DF_BC_FREE_END;                    // boundary conditions type
    double *bc;                                          // array of boundary conditions
    bc = NULL;
    
    ////// Parameters describing internal conditions //////
    
    MKL_INT ic_type = DF_NO_IC;          // internal conditions type
    double *ic;                          // array of internal conditions
    ic = NULL;
    
    ////// Parameters describing interpolation sites //////
    
    MKL_INT nsite = Nr;                      // number of interpolation sites
    MKL_INT sitehint = DF_UNIFORM_PARTITION;    // additional info about interpolation sites
    
    sites[0]=wU.r[0];
    sites[1]=wU.r[wU.Nr-1];
    
    
    ////// Additional info about structure of arrays x and y //////
	
    double *datahint;                   // additional info about structure
    datahint=0;
    
    ////// Parameter describing size of array for derivative orders //////
	
    MKL_INT ndorder = 1;                    // size of array describing derivative orders
    MKL_INT dorder[] = {1};                 // spline values will be computed
    
    ////// Parameter describing interpolation results storage //////
    
    MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;       // interpolation results storage format
    
    
    ////////////////////
	// Cell searching //
	////////////////////
    
    ////// Parameter describing array of cell indices //////
    
    MKL_INT *cell_idx;                  // indices of cells containing interpolation sites
    cell_idx = NULL;
    
    cell_idx=cells;
    
    ////// Create cell search task //////
    
    errcode = dfdNewTask1D( &celltask, nx, wHank.r, xhint, 0, 0, 0 );
    CheckDfError(errcode);
    
    errcode = dfdSearchCells1D( celltask, DF_METHOD_STD, nsite, sites,
                               sitehint, datahint, cell_idx);
    CheckDfError(errcode);
    
    ////// Delete cell search task //////
    
    errcode = dfDeleteTask( &celltask );
	CheckDfError(errcode);
	
    
    //for(int k=0;k<Nr;k++)
    //    printf("Number of cell %d: %d\n",k,cells[k]);
    
    ////////////////////////////
	// Create the interpolant //
	////////////////////////////
    
    for(int j=0;j<Nr;j++)
        for(int i=0;i<Nz;i++)
        {
            auxInReal[wHank.index(i,j)]=wHank.phi[wHank.index(i,j)].real();
            auxInImag[wHank.index(i,j)]=wHank.phi[wHank.index(i,j)].imag();
        }
    
	
    ////// Create Data Fitting task //////
    errcode = dfdNewTask1D(&taskReal, nx, wHank.r, xhint, ny, auxInReal, yhint );
    CheckDfError(errcode);
    
    errcode = dfdNewTask1D(&taskImag, nx, wHank.r, xhint, ny, auxInImag, yhint );
    CheckDfError(errcode);
    
    
    ////// Edit task parameters for Hermite cubic spline construction //////
    errcode = dfdEditPPSpline1D( taskReal, sorder, stype, bc_type, 0,
                                ic_type, 0, scoeffReal, scoeffhint );
    CheckDfError(errcode);
    
    errcode = dfdEditPPSpline1D( taskImag, sorder, stype, bc_type, 0,
								ic_type, 0, scoeffImag, scoeffhint );
    CheckDfError(errcode);
	
    ////// Construct Natural cubic spline using STD method //////
    errcode = dfdConstruct1D( taskReal, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);
    
    errcode = dfdConstruct1D( taskImag, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);
    
    
    
    ////////////////////////////
    // Evaluate the polynom   //
    ////////////////////////////
    
    ////// Interpolate using PP method //////
    errcode = dfdInterpolate1D( taskReal, DF_INTERP, DF_METHOD_PP,
                               nsite, sites, sitehint, ndorder,
                               dorder, datahint, auxOutReal, rhint, cell_idx );
    CheckDfError(errcode);
	
    
    errcode = dfdInterpolate1D( taskImag, DF_INTERP, DF_METHOD_PP,
                               nsite, sites, sitehint, ndorder,
                               dorder, datahint, auxOutImag, rhint, cell_idx );
    CheckDfError(errcode);
    
    for(int j=0;j<Nr;j++)
        for(int i=0;i<Nz;i++)
        {
            wU.phi[wU.index(i,j)].real()=auxOutReal[wU.index(i,j)];
            wU.phi[wU.index(i,j)].imag()=auxOutImag[wU.index(i,j)];
        }
    
    
    ////// Delete Data Fitting task //////
    errcode = dfDeleteTask( &taskReal );
    CheckDfError(errcode);
    
    errcode = dfDeleteTask( &taskImag );
    CheckDfError(errcode);
    
    
    
    ///////////////////////////
	// Free auxiliary arrays //
	///////////////////////////
    
    delete[] auxInReal, auxInImag, auxOutReal, auxOutImag;
    delete[] scoeffReal, scoeffImag;
    delete[] sites, cells;
    
}



///////////////////////////////////////////////////////////////////////////////////////////////


void interpU2H(waveUniform2D &wU, waveH2D &wHank)
{
    
    /////////////////////////////////////////////////
	// Declare parameters for interpolate routines //
	/////////////////////////////////////////////////
    
    int Nr = wU.Nr;
    int Nz = wU.Nz; 
    
    DFTaskPtr taskReal;                     // Data Fitting task descriptor
    DFTaskPtr taskImag; 
    DFTaskPtr celltask;
    double *auxInReal = new double[Nr*Nz];        // Auxiliary arrays
    double *auxInImag = new double[Nr*Nz];
    double *auxOutReal = new double[Nr*Nz];
    double *auxOutImag = new double[Nr*Nz];
    MKL_INT *cells = new MKL_INT[wU.Nr];
    double *x_idx = new double[2];
    double *scoeffReal;
	double *scoeffImag;
    
    
    int errcode = 0;
    
    ////// Initializing parameters for Data Fitting task //////
    
    MKL_INT sorder = DF_PP_CUBIC;                     // spline order
    MKL_INT stype = DF_PP_NATURAL;                    // spline type
    
    ///// Parameters describing interpolation interval //////
    
    MKL_INT nx = Nr;                           // number of break points
    MKL_INT xhint = DF_UNIFORM_PARTITION;                  // additional info about break points
    
    x_idx[0]=wU.r[0];
    x_idx[1]=wU.r[Nr-1];
    
    ////// Parameters describing function //////
    MKL_INT ny = Nz;                      // number of functions
    MKL_INT yhint = DF_NO_HINT;                      // additional info about function
    
    
    ////// Parameters describing spline coefficients storage //////
    
    MKL_INT nscoeff = ny*(nx-1)*sorder;                    // number of spline coefficients
    MKL_INT scoeffhint = DF_NO_HINT;                       // additional info about spline coefficients
    scoeffReal = new double[nscoeff];                                // array of spline coefficients
    scoeffImag = new double[nscoeff];
	
    ////// Parameters describing boundary conditions type //////
    
    MKL_INT bc_type = DF_BC_FREE_END;                    // boundary conditions type
    double *bc;                                          // array of boundary conditions
    bc = NULL;
    
    ////// Parameters describing internal conditions //////
    
    MKL_INT ic_type = DF_NO_IC;          // internal conditions type
    double *ic;                          // array of internal conditions
    ic = NULL;
    
    ////// Parameters describing interpolation sites //////
    
    MKL_INT nsite = wHank.Nr;                      // number of interpolation sites
    MKL_INT sitehint = DF_NON_UNIFORM_PARTITION;    // additional info about interpolation sites
    
    ////// Additional info about structure of arrays x and y //////
    
    double *datahint;                   // additional info about structure
    datahint=0;
    
    ////// Parameter describing size of array for derivative orders //////
    
    MKL_INT ndorder = 1;                    // size of array describing derivative orders
    MKL_INT dorder[] = {1};                 // spline values will be computed
    
    ////// Parameter describing interpolation results storage //////
    
    MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;       // interpolation results storage format
    
    
    ////////////////////
	// Cell searching //
	////////////////////
    
    ////// Parameter describing array of cell indices //////
    
    MKL_INT *cell_idx;                  // indices of cells containing interpolation sites
    cell_idx = 0;
    
    cell_idx=cells;
    
    ////// Create cell search task //////
    
    errcode = dfdNewTask1D( &celltask, nx, x_idx, xhint, 0, 0, 0 );
    CheckDfError(errcode);
    
    errcode = dfdSearchCells1D( celltask, DF_METHOD_STD, nsite, wHank.r,
                               sitehint, datahint, cell_idx);
    CheckDfError(errcode);
    
    ////// Delete cell search task //////
    
    errcode = dfDeleteTask( &celltask );
    CheckDfError(errcode);
    
    
    //for(int k=0;k<Nr;k++)
    //    printf("Number of cell %d: %d\n",k,cells[k]);
    
    ////////////////////////////
	// Create the interpolant //
	////////////////////////////
    
    for(int i=0;i<Nz;i++)
        for(int j=0;j<Nr;j++)
        {
            auxInReal[wU.index(i,j)]=wU.phi[wU.index(i,j)].real();
            auxInImag[wU.index(i,j)]=wU.phi[wU.index(i,j)].imag();
        }
	
    
    ////// Create Data Fitting task //////
    errcode = dfdNewTask1D(&taskReal, nx, x_idx, xhint, ny, auxInReal, yhint );
    CheckDfError(errcode);
    
    errcode = dfdNewTask1D(&taskImag, nx, x_idx, xhint, ny, auxInImag, yhint );
    CheckDfError(errcode);
    
    
    ////// Edit task parameters for Hermite cubic spline construction //////
    errcode = dfdEditPPSpline1D( taskReal, sorder, stype, bc_type, 0,
                                ic_type, 0, scoeffReal, scoeffhint );
    CheckDfError(errcode);
    
    errcode = dfdEditPPSpline1D( taskImag, sorder, stype, bc_type, 0,
                                ic_type, 0, scoeffImag, scoeffhint );
    CheckDfError(errcode);
    
    ////// Construct Natural cubic spline using STD method //////
    errcode = dfdConstruct1D( taskReal, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);
    
    errcode = dfdConstruct1D( taskImag, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);
    
    
    
    ////////////////////////////
    // Evaluate the polynom   //
    ////////////////////////////
    
    ////// Interpolate using PP method //////
    errcode = dfdInterpolate1D( taskReal, DF_INTERP, DF_METHOD_PP,
                               nsite, wHank.r, sitehint, ndorder,
                               dorder, datahint, auxOutReal, rhint, cell_idx );
    CheckDfError(errcode);
    
    
    errcode = dfdInterpolate1D( taskImag, DF_INTERP, DF_METHOD_PP,
                               nsite, wHank.r, sitehint, ndorder,
                               dorder, datahint, auxOutImag, rhint, cell_idx );
    CheckDfError(errcode);
    
    for(int i=0;i<Nz;i++)
        for(int j=0;j<Nr;j++)
        {
            wHank.phi[wHank.index(i,j)].real()=auxOutReal[wHank.index(i,j)];
            wHank.phi[wHank.index(i,j)].imag()=auxOutImag[wHank.index(i,j)];
        }   
	
    
    ////// Delete Data Fitting task //////
    errcode = dfDeleteTask( &taskReal );
    CheckDfError(errcode);
    
    errcode = dfDeleteTask( &taskImag );
    CheckDfError(errcode);
    
    
    
    ///////////////////////////
	// Free auxiliary arrays //
	///////////////////////////
    
    delete[] auxInReal, auxInImag, auxOutReal, auxOutImag;
    delete[] x_idx, cells;
	delete[] scoeffReal, scoeffImag;
	
    
}

///////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
/////////////////////   A function for check errors   ////////////////////////
//////////////////////////////////////////////////////////////////////////////


void CheckDfError(int num)
{
    switch(num)
    {
        case DF_ERROR_NULL_TASK_DESCRIPTOR:
        {
            printf( "Error: null task descriptor (code %d).\n", num );
            break;
        }
        case DF_ERROR_MEM_FAILURE:
        {
            printf( "Error: memory allocation failure in DF functionality (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_NX:
        {
            printf( "Error: the number of breakpoints is invalid (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_X:
        {
            printf( "Error: the array which contains the breakpoints is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_X_HINT:
        {
            printf( "Error: invalid flag describing structure of partition (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_NY:
        {
            printf( "Error: invalid dimension of vector-valued function y (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_Y:
        {
            printf( "Error: the array which contains function values is invalid (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_Y_HINT:
        {
            printf( "Error: invalid flag describing structure of function y (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_SPLINE_ORDER:
        {
            printf( "Error: invalid spline order (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_SPLINE_TYPE:
        {
            printf( "Error: invalid type of the spline (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_IC_TYPE:
        {
            printf( "Error: invalid type of internal conditions used in the spline construction (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_IC:
        {
            printf( "Error: array of internal conditions for spline construction is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_BC_TYPE:
        {
            printf( "Error: invalid type of boundary conditions used in the spline construction (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_BC:
        {
            printf( "Error: array which presents boundary conditions for spline construction is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_PP_COEFF:
        {
            printf( "Error: array of piece-wise polynomial spline coefficients is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_PP_COEFF_HINT:
        {
            printf( "Error: invalid flag describing structure of the piece-wise polynomial spline coefficients (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_PERIODIC_VAL:
        {
            printf( "Error: function values at the end points of the interpolation interval are not equal as required in periodic boundary conditions (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_DATA_ATTR:
        {
            printf( "Error: invalid attribute of the pointer to be set or modified in Data Fitting task descriptor with EditIdxPtr editor (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_DATA_IDX:
        {
            printf( "Error: index of pointer to be set or modified in Data Fitting task descriptor with EditIdxPtr editor is out of range (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_NSITE:
        {
            printf( "Error: invalid number of interpolation sites (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_SITE:
        {
            printf( "Error: array of interpolation sites is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_SITE_HINT:
        {
            printf( "Error: invalid flag describing structure of interpolation sites (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_NDORDER:
        {
            printf( "Error: invalid size of array that defines order of the derivatives to be computed at the interpolation sites (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_DORDER:
        {
            printf( "Error: array defining derivative orders to be computed at interpolation sites is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_DATA_HINT:
        {
            printf( "Error: invalid flag providing a-priori information about partition and/or interpolation sites (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_INTERP:
        {
            printf( "Error: array of spline based interpolation results is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_INTERP_HINT:
        {
            printf( "Error: invalid flag defining structure of spline based interpolation results (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_CELL_IDX:
        {
            printf( "Error: array of indices of partition cells containing interpolation sites is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_NLIM:
        {
            printf( "Error: invalid size of arrays containing integration limits (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_LLIM:
        {
            printf( "Error: array of left integration limits is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_RLIM:
        {
            printf( "Error: array of right integration limits is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_INTEGR:
        {
            printf( "Error: array of spline based integration results is not defined (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_INTEGR_HINT:
        {
            printf( "Error: invalid flag defining structure of spline based integration results (code %d).\n", num );
            break;
        }
        case DF_ERROR_BAD_LOOKUP_INTERP_SITE:
        {
            printf( "Error: bad site provided for interpolation with look-up interpolator (code %d).\n", num );
            break;
        }
        default: break;
    }
    
    if(num < 0) {
        exit(1);
    }
}


#endif   // INTERP_H
