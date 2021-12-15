#ifndef MINIMIZATION_FUNCTION_HPP
#define MINIMIZATION_FUNCTION_HPP

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "vecteur.hpp"
#include "structure.h"

#define _RECIPES_INFINITE 1.0E+100
#define SHFT(a,b,c,d)   (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)       ((b) > 0.0 ? fabs(a) : -fabs(a))
#define CGOLD           0.3819660
#define ZEPS            1.0e-10
#define TEST            fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)
#define GOLD            1.618034
#define GLIMIT          100.0
#define TINY            1.0e-20
#define MAX(a,b)        ((a) > (b) ? (a) : (b))

/* Brent Algorithm: Finds the minimum (xmin) of the function f(x) with a given tol accuracy 
 * through a ITMAX number of iteration. It returns the found minimum and the updated xmin by reference.
 * (Numerical Recipes in C, William H. Press, page 403)
 * given 3 points ax, bx, cx : parabolic fit */
double
brent (const double ax, const double bx, const double cx, double (*f)(double, void *),
                void * localBuffer, const double tol, double * xmin, const int ITMAX);

void*
testMetric(const Output_st *pImage, const Input_st *pData,
        double *rIntensity, double *rVarIm, double *rGradBrenner);
#endif
