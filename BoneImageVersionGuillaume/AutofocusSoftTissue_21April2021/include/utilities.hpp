#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "mex.h"
#include "structure.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>

/* function to point toward NULL before deallocate him
 * Necessary when pointer is used for MEX output.
 *  */
void*
nulifyAndFreePtr( double* pTr );
/* function to zero reset Input_st/Output_st 's tables
 * when looping over model parameters.
 *  */
void*
resetImage(const Input_st *pData, Output_st *pImage);
/* function to allocate a 2D Table used for our Image to Build
 *  */
double **
alloc2D (const int N, const int M);

/* function to desallocate a 2D table
 *  */
void
free2D (double **pTr, const int N);

/* function to calculate the square of a double x
 * */
double
sqr (const double& x);

/* function to calculate 2d cartesian distance between 2 point(x,z)
 * */
double
distance (const double& x1, const double& z1, const double& x2, const double& z2);

/* function to allocate and initialize The Output_st structure with 0 using calloc.
 * */
void*
initializeImage(const Input_st *pData, Output_st *pImage);

/* function to desalocate the Output_st structure using free.
 * */
void*
freeImage(Output_st *pImage, const int NX, const int NZ);

/* A function to get prhs (matlab mxarray pointer) pointing to 1D mxArray and fill The Expsetup_st structure.
 * It contains all constants from the experience, See include/Structure.h
 * ShowInput = 1 means display the structure content.
 *  */
void*
fillSetup(ExpSetup_st &rSetup, const mxArray *pMxArray, const int ShowInput);

/* Function to fill 1D Table from prhs pointing to 1D mxArray. Example: X/Z Position.
 *  */
int
fillPointerFromMxArray(double **pData, const mxArray *pMxArray);

/* Function to fill 1D Table from prhs pointing to 3D mxArray (Raw Signal).
 * Signalpart = 0/1 means we fill Inphase/Quadratic part of the Signal.
 * */
double*
fillSignalFromMxArray(Input_st *pData, const mxArray *pMxArray, const int& SignalPart);

#endif
