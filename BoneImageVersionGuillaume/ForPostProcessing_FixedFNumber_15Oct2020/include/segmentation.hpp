#ifndef DEF_SEGMENTATION
#define DEF_SEGMENTATION

#include <iostream>
#include <math.h>
#include <string.h>
#include "utilities.hpp"


/* Segmentation is used to find the interfaces, which are the lines with highest amplitudes from
 * left to right of the image
 * */
int*
createContour(double **Im, const int Nx, const int Nz);

/* Segmentation is used to find the interfaces,
 * which are the lines with highest amplitudes from left to right of the image
 * */
void
segmentation(Parabola_st *Parabola, double **Im, const int NX, const int NZ, double *X, double *Z);

#endif
