#include "rayTheory.hpp"
using namespace std;

double
travelTimeInTissue(double Theta, void *pBuffer)
{
	ExpSetup_st *pSetup     = (ExpSetup_st *)pBuffer;
    	// first path, from sources to the end of the lens
    	double tanTheta 	= tan(Theta);
        // point reached at the interface lens/tissue
    	double x                = (pSetup->XStart + (pSetup->LensThick - pSetup->ZStart) * tanTheta);
    	double z                = pSetup->LensThick;

        // time for the ray to reach the interface (travel in the lens)
    	double TravelTime = distance(pSetup->XStart, pSetup->ZStart, x, z) / pSetup->CLens;

    	// time from the end (x) of the lens to the image point in tissue
    	TravelTime += distance(x, z, pSetup->XEnd, pSetup->ZEnd) / pSetup->CTissue;

    	// Traveltime is the total time from the element to the point
	return TravelTime;
}

double
RayLengthUntilTissue(double Theta, void *pBuffer)
{
	ExpSetup_st *pSetup     = (ExpSetup_st *)pBuffer;
    // first path, from sources to the end of the lens
    double tanTheta 	= tan(Theta);
    // point reached at the interface lens/tissue
    double x                = (pSetup->XStart + (pSetup->LensThick - pSetup->ZStart) * tanTheta);
    double z                = pSetup->LensThick;

    double RayL = distance(pSetup->XStart, pSetup->ZStart, x, z);

    RayL += distance(x, z, pSetup->XEnd, pSetup->ZEnd);

	return RayL;
}
