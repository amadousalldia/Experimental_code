#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "mex.h"
#include "structure.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>


void*
nulifyAndFreePtr( double* pTr );
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

/* Calculates maximum value (&Maximum) and its position (&Index) in array
 * */
void
maxInArray(double *pTr, const int N, double &Maximum, int &Index);

/* Function to allocate a int **pTr pointing to a 2D Table made of int.
 * */
int**
newTable2D (const int N, const int M);

/* Function to desallocate a int **pTr pointing to a 2D Table.
 * */
void
deleteTable2D (int **pTr, const int N);

/* function to allocate and initialize The Output_st structure with 0 using calloc.
 * */
void*
initializeImage(const Input_st *pData, Output_st *pImage);

/* function to desalocate the Output_st structure using free.
 * */
void*
freeImage(Output_st *pImage, const int NX, const int NZ);

/* function to copy a 2D table using calloc
 * */
double**
copyTable2D(double **pTrToCopy, const int N, const int M);

/* function to desalocate a double **pTr pointing to 2D table made of double using free.
 * */
void**
freeTable2D(double **pTrToFree, const int N);

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
/* Function to write out to a nlhs pointing to 1D mxArray containg updated Parabolic
 * Parameters
 *  */
mxArray*
fillOutParabola(Parabola_st &rParabola);

/* Function to save to binary file an Expesetup_st */
void*
saveSetup_stToBinary(ExpSetup_st *Setup);

/* Function to load to binary file an Expesetup_st */
void*
loadSetup_stFromBinary(ExpSetup_st *Setup);

/* Function to save to binary file an Input_st  */
void*
saveInput_stToBinary(Input_st *Data);

/* Function to load to binary file an Input_st  */
void*
loadInput_stFromBinary(Input_st *Data);

/* Function to save to binary file an Output_st  */
void*
saveImageToBinary(double **pImage, const int nX, const int nZ, const char* FileName);

/* Function to save to binary file an Output_st  */
void*
loadImageFromBinary(double **pImage, const int nX, const int nZ, const char* FileName);

void*
saveParabola_stToBinary(Parabola_st &rParabola, const char* FileName);

Parabola_st
loadParabola(const char* FileName);

/* Open file FileName with specialified mode Mode to check if exist
 */
void*
isFileExist(const std::string& FileName, const char *Mode);
#endif
