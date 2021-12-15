#include "buildImage.hpp"
using namespace std;

void*
buildImage(const ExpSetup_st *Setup,const Input_st *Data, Output_st *Image)
{
    int i, j, iTx, Tx, r, Index, OffsetT, OffsetSig, OffsetR;
    double SumI, SumQ, TotalTime, Weight, TransmitTime;

    omp_set_num_threads(Setup->NCore);
    #pragma omp parallel for schedule(static) \
    shared(Setup, Data, Image) \
    private( i,j, iTx, Tx, OffsetT, TransmitTime, \
                    r, OffsetR, TotalTime, Index, \
                    OffsetSig, Weight, SumI, SumQ )

    for (i = 0; i < Data->NX; ++i)
	{
        for (j = Data->ZMin; j < Data->ZMax; ++j)
        {
            for (iTx = 0; iTx < Data->NTX; ++iTx)
            {
                SumI = 0., SumQ = 0.; // total intensity for pixel

                Tx 	= Data->pTxToUse[iTx] - 1;
                OffsetT = i*Data->NS*Data->NZ + j*Data->NS + Tx;
                
//                 if (fabs(Image->pAngleT[OffsetT]) < Setup->HalfOpenAngle)
//                 { // if angle at transmit element lies within acceptance angle
                
                // Total transmitted time.
                TransmitTime = Image->pTimeT[OffsetT]  - Data->pAddDelay[Tx];
                for (r = 0; r < Data->NR; ++r)
                {
                    OffsetR = i*Data->NR*Data->NZ + j*Data->NR + r;
                    if (fabs(Image->pAngleR[OffsetR]) < Setup->HalfOpenAngle)
                    { // if angle at receive element lies within acceptance angle
                        TotalTime = (Image->pTimeR[OffsetR] + TransmitTime) * Setup->SamplingFreq;
                        Index     = int(TotalTime);

                        if (Index >= 1 && Index < Data->TimeSig -1)
                        { // get data
                            OffsetSig = Tx*Data->TimeSig*Data->NR + r*Data->TimeSig + Index;
                            Weight    = 1. - (TotalTime - floor(TotalTime));
                            SumI     += Data->pSignalI[OffsetSig] * Weight + Data->pSignalI[OffsetSig + 1] * (1. - Weight); // linear interpolation
                            SumQ     += Data->pSignalQ[OffsetSig] * Weight + Data->pSignalQ[OffsetSig + 1] * (1. - Weight); // linear interpolation

                        }
                    }
                }
//                 }
                Image->pImI[i][j] += SumI;
                Image->pImQ[i][j] += SumQ;
            }
        }
    }
    return NULL;
}

void*
buildImageForVFI(const ExpSetup_st *Setup, const Input_st *Data, Output_st *Image)
{
    int i, j, iTx, Tx, iRx, r, Index, OffsetT, OffsetSig, OffsetR, OffsetBeam;
    double sumApI, sumApQ, TotalTime, Weight, TempI, TempQ, TransmitTime;
    int minR, maxR, nElemInSubR;

    omp_set_num_threads(Setup->NCore);
    #pragma omp parallel for schedule(static) \
    shared(Setup, Data, Image) \
    private(i, j, iTx, Tx, OffsetT, TransmitTime, iRx, sumApI, sumApQ, r, OffsetR, \
            TotalTime, Index, OffsetSig, Weight, TempI, TempQ, OffsetBeam, minR, maxR, nElemInSubR )
    for (i = 0; i < Data->NX; ++i)
	{
        double* pWindow = new double [Data->NR];
        for (j = Data->ZMin; j < Data->ZMax; ++j)
        {
            for (iTx = 0; iTx < Data->NTX; ++iTx)
            { // Nb of Transmission.
                Tx 	= int(Data->pTxToUse[iTx]) - 1;
                OffsetT = i*Data->NS*Data->NZ + j*Data->NS + Tx;
                TransmitTime = Image->pTimeT[OffsetT] - Data->pAddDelay[Tx];

                for (iRx = 0; iRx < Data->NLA; ++iRx)
                { // Nb of listening angles
                    sumApI = 0.; sumApQ = 0.; minR = 0.; maxR = 0.; nElemInSubR = 0;
                    // Search for the active listening subaperture for a given listening angle anf F-number.
                    nElemInSubR = getActiveSubAperture(i, j, iRx, Image, Data, &minR, &maxR, Setup->Fnumber, Setup->MaxDevListenAngle);
                    // Initialize Window at each iRx
                    setArrayToZero(pWindow,Data->NR);
                    // If we have listening elements (size of subaperture!=0) we build.
                    if (nElemInSubR > 0) // test to check if we use at least 10 elements.
                    {
                        // Apply apodization on listening SubAperture.
                        getWindowFunction(pWindow, nElemInSubR, Setup->SubAperApodis);
                        for (r = minR; r <= maxR; r++)
                        { // receiver index within Rx aperture
                            OffsetR   = i*Data->NR*Data->NZ + j*Data->NR + r;
                            TotalTime = (Image->pTimeR[OffsetR] + TransmitTime) * Setup->SamplingFreq;
                            Index     = int(TotalTime);

                            if (Index >= 1 && Index < Data->TimeSig -1)
                            { // Delay and Sum
                                OffsetSig = Tx*Data->TimeSig*Data->NR + r*Data->TimeSig + Index;
                                Weight    = 1. - (TotalTime - floor(TotalTime));
                                TempI     = Data->pSignalI[OffsetSig] * Weight + Data->pSignalI[OffsetSig + 1] * (1. - Weight);
                                TempQ     = Data->pSignalQ[OffsetSig] * Weight + Data->pSignalQ[OffsetSig + 1] * (1. - Weight);
                                sumApI    += TempI*pWindow[r-minR];
                                sumApQ    += TempQ*pWindow[r-minR];
                            }
                        }
                        OffsetBeam = Tx*Data->NX*Data->NZ*Data->NLA +
                                        iRx*Data->NX*Data->NZ +
                                        i*Data->NZ +
                                        j;

                    Image->pBeamQ[OffsetBeam] += sumApQ;
                    Image->pBeamI[OffsetBeam] += sumApI;
                    }
                }
            }
        }
        delete [] pWindow;
    }
    return NULL;
}

int
getActiveSubAperture(const int IdX, const int IdZ, const int IdL, const Output_st* pImage, const Input_st* pData,
                     int *rMinR, int *rMaxR, const double FNumber, const double MaxDeviationFromListenAngle)
{
    double diffAngle, minDiffAngle = 0.5236; // Difference between LisenAngle/AngleAtPix, Initial value for minimum seeking. (30deg)
    int IdxClosest                 = 0;                           // Idx of the listening element.

    // We first get the element which is aligned with the given listenAngle from a given Pixel(x,z).
    for (int r = 0; r < pData->NR; r++)
    {// Looping over elements to find r which observe a pixel with the closest value to ListenAngle.

        // Difference between the angle of refraction of the wave and the ListenAngle.
        diffAngle = fabs(pData->pListenAngle[IdL] - pImage->pAngleAtPixR[IdX*pData->NR*pData->NZ + IdZ*pData->NR + r]);
        if ( diffAngle < minDiffAngle )
        {
                minDiffAngle = diffAngle; // keep minimum difference.
                IdxClosest   = r;         // IdxCloset = index of r that minimzes.
        }
    }
    // If there are no element listening with that angle (Or too far from that listening angle, i.e.  > tresholdAngle, we dont construct image.
    if (fabs(minDiffAngle) > MaxDeviationFromListenAngle) return 0; // Not very useful, but for safety.

    // Else IdxClosest is the center of an active listening subaperture.
    // Using the F-number definition and its value imposed by user, we get the width of this subaperture centered at IdxClosest.
    // Approximation: The subaperture is symmetrical around IdxClosest.
    else {
        // We keep angles at Pixel when the center is identified.
        pImage->pListenedAngleAtPixR[IdX*pData->NLA*pData->NZ + IdZ*pData->NLA + IdL] =
                pImage->pAngleAtPixR[IdX*pData->NR*pData->NZ + IdZ*pData->NR + IdxClosest];

        // Euclidian distance between center of sub and a given pixel.
        double FocalDistance       = distance(pData->pXR[IdxClosest], 0, pData->pX[IdX], pData->pZ[IdZ]);
        double ActiveApertureWidth = FocalDistance/FNumber;

        // Get the 1st element within Activeaperturewidth left to Idxclosest;
        for (int i = 0; i < IdxClosest; i++)
        {
            if (sqrt(sqr(pData->pXR[i]-pData->pXR[IdxClosest])) >= ActiveApertureWidth/2) *rMinR = i;
        }

        // Get the 1st element within Activeaperturewidth right to Idxclosest;
        for (int i = IdxClosest; i < pData->NR; i++)
        {
            if (sqrt(sqr(pData->pXR[i]-pData->pXR[IdxClosest])) <= ActiveApertureWidth/2) *rMaxR = i;
        }

        // Check if the subaperture is symmetrical. If yes we move.
        if (*rMaxR-IdxClosest==IdxClosest-*rMinR) return abs(*rMaxR - *rMinR) + 1;

        // if not we check who is the longest and scale down to the other side value.
        if (*rMaxR-IdxClosest > IdxClosest-*rMinR) *rMaxR = 2*IdxClosest - *rMinR;
        else *rMinR = 2*IdxClosest - *rMaxR;

        // Number of element within our Sub.
        return abs(*rMaxR - *rMinR) + 1;
    }
}

double
getFnumber(const double ElemXMin, const double ElemXMax, const double PixelX, const double PixelZ)
{
    return distance((ElemXMax+ElemXMin)/2, 0, PixelX, PixelZ) / (ElemXMax-ElemXMin);
}

void*
getWindowFunction(double* pWindow,  const int N, const int Flag)
{
    // N = L -1;
    // Definition from Discrete time signal processing. Oppenheim et al., 1999
    const double PI = 3.14159265358979323846;

    switch (Flag)
    {
        case 0 :
                // No window. Apply 1
                // N = L -1; so we just do < instead of <=.
                for (int i = 0; i < N; i++)
                {
                    pWindow[i] = 1.;
                }
                break;
        case 1 :
                // Hamming window.
                for (int i = 0; i < N; i++)
                {
                    pWindow[i] = 0.53836 - 0.46164*cos(2*PI*i/N);
                }
                break;
        case 2 :
                // Hann window.
                for (int i = 0; i < N; i++)
                {
                    pWindow[i] = 0.5 - 0.5*cos(2*PI*i/N);
                }
                break;
    }
    return NULL;
}
