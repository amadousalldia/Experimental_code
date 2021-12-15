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

double
travelTimeInTissue(double Theta, void *pBuffer);

double
RayLengthUntilTissue(double Theta, void *pBuffer);

#endif
