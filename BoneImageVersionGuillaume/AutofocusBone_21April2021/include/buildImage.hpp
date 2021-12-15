#ifndef DEF_IMAGERECONSTR
#define DEF_IMAGERECONSTR

#include <math.h>

#include "structure.h"
#include "omp.h"
#include "iostream"
#include "utilities.hpp"

/* Function that calculated the image amplitude at each point but applying Delay and Sum algorithm.
* */
void*
buildImage(ExpSetup_st *Setup, Input_st *Data, Output_st *Image);
#endif
