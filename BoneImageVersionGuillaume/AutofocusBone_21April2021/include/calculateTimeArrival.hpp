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
#define nTestPriorTx  100  // Number of test to infer Prior Angle to help Brent converging properly.
#define nTestPriorRx  100  // Number of test to infer Prior Angle to help Brent converging properly.
#define AccuracyBrent 1e-8

/* Function that find the minimum travel time to cross the Lens/Tissue medium.
 * Calls Brent algo and Traveltimeisotropic.
 */
void*
getTimeArrivalTissue (const ExpSetup_st GlobalSetup, 
					    const Input_st * Data, Output_st * Image);

/* Function that find the minimum travel time to cross the Tissue/Bone medium.
 * Calls Brent algo and TraveltimeAnisotropic.
 */
void*
getTimeArrivalBone(const ExpSetup_st GlobalSetup, 
                    const Input_st * Data, Output_st * Image);

/* Function to refine the range of research before calling Brent. When a source is very far from the
 * probe, we calculate the angle needed to a virtual source to reach the left [s][0] and right [s][2] of the probe.
 * The starting point for Brent is then the middle between these 2 angles.
 * */
double**
refineCriticalAngleTable(const Input_st *Data, const double CriticalMin, const double CriticalMax);

/* A function that calls brent algorithm after grid searching the Prior Angle
 * by testing functionToTest (can be travelTimeThroughTissue[Bone,BoneMarrow]).
 * When estimating the path of a ray crossing Tissue->Bone or Tissue->Bone->Marrow,
 * Brent can fall in a local minima leading to wrong results.
 * With a given a priori information, brent converges properly.
 * If AnglePrior is provided (!=0), grid search is avoided.
 * */
double
runBrentWithPriorAngleEstimation(ExpSetup_st &Setup, double& AngleMin,
                double AnglePrior, double (*functionToTest)(double, void *));
#endif
