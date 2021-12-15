#ifndef CALCULATE_TIME_ARR_HPP
#define CALCULATE_TIME_ARR_HPP

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>

#include "minimization.hpp"
#include "rayTheory.hpp"
#include "omp.h"

#define __min(a,b) ((a)>(b) ? (b) : (a)) // function to find the minimum between two numbers
#define __max(a,b) ((a)<(b) ? (b) : (a)) // function to find the maximum between two numbers
#define IterMaxBrent  100
#define AccuracyBrent 1e-8

/* Function that find the minimum travel time to cross the Lens/Tissue medium.
 * Calls Brent algo and Traveltimeisotropic.
 */
void*
getTimeArrivalTissue(const ExpSetup_st GlobalSetup, const Input_st *Data, Output_st *Image);

/* Function to refine the range of research before calling Brent. When a source is very far from the
 * probe, we calculate the angle needed to a virtual source to reach the left [s][0] and right [s][2] of the probe.
 * The starting point for Brent is then the middle between these 2 angles.
 * */
double**
refineCriticalAngleTable(const Input_st *Data, const double CriticalMin, const double CriticalMax);
#endif
