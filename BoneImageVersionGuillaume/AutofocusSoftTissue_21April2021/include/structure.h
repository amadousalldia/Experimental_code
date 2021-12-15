#ifndef STRUCTURE_H
#define STRUCTURE_H

/* Structure containing Signal, Pixel position,
 *  Receivers and Transmetters and there dimsnesions
 * */
struct Input_st {
        double *pXR, *pZR, *pXS, *pZS, *pX, *pZ;
        double *pSignalI, *pSignalQ;
        double *pAddDelay, *pTxToUse;
        int NTX, NR, NS, NX, NZ, ZMax, ZMin;
        int TimeSig, SrcSig, RcptSig;
};

/* Structure containing Model parameters, Local Pixel positions for different functions,
 * Experience settings, and Parameters for MEX behavior.
 * */
struct ExpSetup_st {

        double XStart = 0, XEnd = 0, ZStart = 0, ZEnd = 0;      // Local Pixel poititions.
        double LensThick, LensPitch;             		// Model/Media parameters.
        double AngleCriticalMin = 0., AngleCriticalMax = 0.;
        double HalfOpenAngle, SamplingFreq, CLens, CTissue;
        int NCore;

        // Number of Members within this Structure.
        const int NSetupMember = 7;
};

/* Structure containing the built image and outputs (angle, time arrival, beamform)
 * */
struct Output_st {
        double *pTimeT, *pTimeR;                // Time Arrival Trans/Receiv
        double *pAngleT, *pAngleR;              // Emerging angles
        double *pRayLengthT, *pRayLengthR;      // Ray Length
        double **pIm, **pImI, **pImQ;           // Enveloppe/I/Q Image
};
#endif
