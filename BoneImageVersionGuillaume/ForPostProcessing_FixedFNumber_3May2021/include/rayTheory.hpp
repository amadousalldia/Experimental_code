#ifndef RAYTHEORY_HPP
#define RAYTHEORY_HPP

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "vecteur.hpp"
#include "utilities.hpp"
#include "minimization.hpp"
#include "structure.h"

/* Apply Snell Descartes law for an Isotropic medium
 * Return 1 of refraction occured anf update the refracted angle
 * */
int
refractionIsotropic(const vecteur& Ui, const vecteur& N, const double& RatioC1C2, vecteur * Ur);

/* Calculate the distance between a parabola (a,b,c) and a point.
 * */
double
rayLength(const double& a, const double& b, const double& c, const vecteur& Pt, const vecteur& Ur);

double
snellEq(const double& SinTheta, const vecteur& Ui, const vecteur& N,
                const double& AnisoCoef, const double& AnisoShape,
                const double& CBone, const double& CIso);

/* Apply Snell Descartes law for an VTI Anisotropic medium
 * Return 1 of refraction occured anf update the refracted angle
 * */
int
refractionAnisotropic (const vecteur& Ui, const vecteur& N, vecteur *Ur,
                double *SinSqrGrAngle, ExpSetup_st *Setup);

/* Calculate the travel time of a ray between the source element and the pixel
 * assuming the medium is isotropic
 */
double
travelTimeInTissue(double Theta, void *pBuffer);

/* Extract the angle of refraction at pixel
 * when bouncing point is within the Tissue.
 */
double
getAngleAtScatterWithinTissue(const double& Theta, const double& RatioClCt);

/* Calculate the travel time of a ray between the source element and the pixel
 * assuming the medium is anisotropic
 */
double
travelTimeThroughTissueBone(double Theta, void *pBuffer);

/* Extract the angle of refraction at pixel
 * when bouncing point is within the bone
 */
double
getAngleAtScatterWithinBone (const double& Theta, void *pBuffer);

/* Calculate the travel time of a ray between the source element and the pixel
 * assuming we are crossing Lens/Tissue/Bone/Marrow.
 */
double
travelTimeThroughTissueBoneMarrow( double Theta, void *pBuffer);

/* Extract the Angle of refraction at Pixel
 * when bouncing point is within the Marrow.
 */
double
getAngleAtScatterWithinMarrow (const double& Theta, void *pBuffer);
#endif
