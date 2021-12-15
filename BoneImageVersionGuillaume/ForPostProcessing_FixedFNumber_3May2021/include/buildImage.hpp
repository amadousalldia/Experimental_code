#ifndef DEF_IMAGERECONSTR
#define DEF_IMAGERECONSTR

#include <math.h>

#include "structure.h"
#include "omp.h"
#include "iostream"
#include "mex.h"
#include "utilities.hpp"

/* Delay and Sum algorthim for B-Mode image.
* */
void*
buildImage(const ExpSetup_st *Setup, const Input_st *Data, Output_st *Image);

/* Delay and Sum algorithm for Vector Flow analysis.
* */
void*
buildImageForVFI(const ExpSetup_st *Setup, const Input_st *Data,
                     Output_st *Image);

/* At a given pixel and proposed listening angle, the function searches which receptor is
 * aligned with the refracted beam.
 * Since we take into account refraction, it cannot be perfect and some pixel cannot be oserved
 * regarding an error of acceptance (We chose 5deg). If no receptor is found, return 0.
 * Then, with a given F-number we seek for the width of such listening subaperture centered around the just previously calculated receptor.
 * return the number of elements for this active subaperture (if zero, no construction is done.)
 * and update through reference its range (minR, maxR).
 * IdX, IdZ, IdL are the index of the (X,Z)pixel and (L)ListenAngle.
 * *minR, maxR are the first, last element of the suabterture updated along the function.
 * FNumber is the F-number.
 * Maxdeviationfromlistenangle is the max deviation from the ListenAngle we allow.
 * */
int
getActiveSubAperture(const int IdX, const int IdZ, const int IdL,
                const Output_st* Image, const Input_st* Data,
                int *minR, int *maxR, const double FNumber,
                const double MaxDeviationFromListenAngle);

/* Calculate F-number for a given active subaperture regarding given illuminated pixel.
 * Elemxmin and Elemxmax are the first and last X positions of a subaperture.
 * Pixelx, Pixelz gives the position of the illuminated pixel.
 * Assumes Element depth is 0.
 * */
double
getFnumber(const double ElemXMin, const double ElemXMax, const double PixelX, const double PixelZ);

/* Calculate the Hamming function from a number of elements.
 * N is the Number of elements per subaperture.
 * Flag = 0,1,2, we use no/hann/hamming filters;
* */
void*
getWindowFunction(double* Window, const int N, const int Flag);

#endif
