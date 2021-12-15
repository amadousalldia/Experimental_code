#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "mex.h"
#include "structure.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>

/* function to calculate the square of a double x
 * */
double
sqr (const double& x);

/* function to calculate 2d cartesian distance between 2 point(x,z)
 * */
double
distance (const double& x1, const double& z1, const double& x2, const double& z2);

/* function to set to zero an array of double.
 * */
void*
setArrayToZero(double* pArray, const int N);

/* Calculates maximum value (&Maximum) and its position (&Index) in array
 * */
void
maxInArray(double *pTr, int N, double &Maximum, int &Index);

/* Function to allocate a int **pTr pointing to a 2D Table made of int.
 * */
int**
newTable2D (const int N, const int M);

/* Function to desallocate a int **pTr pointing to a 2D Table.
 * */
void
deleteTable2D (int **pTr, const int N);
/* function to allocate a 2D Table used for our Image to Build
 *  */
double **
alloc2D (const int N, const int M);

/* function to desallocate a 2D table
 *  */
void
free2D (double **pTr, const int N);

/* A function to get prhs (matlab mxarray pointer) pointing to 1D mxArray and fill The Expsetup_st structure.
 * It contains all constants from the experience, See include/Structure.h
 *  */
void*
initializeImage(const Input_st *pData, Output_st *pImage);

void*
freeImage(Output_st *pImage, const int NX);

/* A function to get prhs (matlab mxarray pointer) pointing to 1D mxArray and fill The Expsetup_st structure.
 * It contains all constants from the experience, See include/Structure.h
 * Set Showinput to 1 if you want to cout values for debugging.
 *  */
void*
fillSetup(ExpSetup_st &rSetup, const mxArray *pMxArray, const int ShowInput);

/* Function to get prhs pointing to 1D mxArray and fill A Parabola_st Structure containing parabolic
 * parameters for discontinuities (Perios, Endeos etc...)
 * */
int
fillParabola(Parabola_st &rParabola, const mxArray *pMxArray);

/* Function to fill 1D Table from prhs pointing to 1D mxArray. Example: X/Z Position.
 *  */
int
fillPointerFromMxArray(double **pData, const mxArray *pMxArray);

/* Function to fill Parabola_st from prhs pointing to 2D mxArray.
 *  */
void*
fillSubAperture(SubAperture_st *pData, const mxArray *pMxArray, const double Pitch);

/* function to free Parabola_st content
 * */
void*
freeSubAperture(SubAperture_st *pData);

/* Function to fill 1D Table from prhs pointing to 3D mxArray (Raw Signal).
 * */
double*
fillSignalFromMxArray (Input_st *pData, const mxArray *pMxArray, const int& SignalPart);

/* Function to write out to a nlhs pointing to 1D mxArray containg updated Parabolic 
 * Parameters
 *  */
mxArray*
fillOutParabola(Parabola_st &rParabola);

/* function to save Angle at Pixel into a binary file
 *  */
void*
saveAngleAtPix(const Input_st *Data, Output_st *Image,
                const char*  Reception, const char* Transmission);

/* function to load Angle at Pixel into a binary file
 *  */
void*
loadAngleAtPix(Input_st *Data, Output_st *Image,
                const char* Reception, const char* Transmission);

/* function to save TimeArrival and Emerging Angles into a binary file to
 * avoid calculations at each doppler image.
 *  */
void*
saveTimeAndAngle(Input_st *Data, Output_st *Image,
                const char* Reception, const char* Transmission);

/* function to load TimeArrival and Emerging Angles from a binary file to
 * avoid calculations at each doppler image.
 *  */
void*
loadTimeAndAngle(Input_st *data, Output_st *Image,
                const char* Reception, const char* Transmission);

/* function to save Parabola parameters into a binary file when avoiding
 * calculation of time arrivals.
 * */
void*
saveParabola( Parabola_st *Para, const char* Parabola);

/* function to load Parabola parameters from a binary file when avoiding
 * calculation of time arrivals.
 * */
void*
loadParabola( Parabola_st *Para, const char* Parabola);

/* function to clean the matrix of selected AngleAtPixR (Data->pListenedAngleAtPixR)
 * before it is given back to the Matlab. Avoid to sort NaN and non NaN when doing VFI.
 * */
void*
cleanAngleAtPixel(const ExpSetup_st *Setup, const Input_st *pData,
                Output_st *pImage, const int IdxLayer, const int ReconTo);
#endif
