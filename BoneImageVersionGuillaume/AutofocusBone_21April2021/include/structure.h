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

/* Structure containing Parabolic parameters
 * and its Minimum
 * */
struct Parabola_st {
    // y = ax2 + bx + c, MinIdx is the index in Data->pZ of the Parabola Min.
    double a = 0, b = 0, c = 0;
    int IdxMin = 0, IdxMax = 0;
    // Number of Members within the structure.
    const int NParameters = 5;
};

/* Structure containing Model parameters, Local Pixel positions for different functions,
 * Experience settings, and Parameters for MEX behavior.
 * */
struct ExpSetup_st {

    Parabola_st *pPeri; // Pointing to Parabola.
    double XStart = 0, XEnd = 0, ZStart = 0, ZEnd = 0;
    double RatioLT = 0., RatioTB = 0., AngleCriticalMin =0., AngleCriticalMax =0.;
    double AnisoCoef, AnisoShape, LensThick, LensPitch;
    double DProbePerios, CorticalThickMin, CorticalThickMax;
    double HalfOpenAngle, SamplingFreq, CLens, CTissue, CBoneAxial, CBoneRadial;   		// Media Parameters
    int NCore, WhichModel;

    // Number of Members within this Structure.
    const int NSetupMember = 14;
};

/* Structure containing the built image and outputs (angle, time arrival, beamform)
 * */
struct Output_st {
    double *pTimeT, *pTimeR;                // Time Arrival Trans/Receiv
    double *pAngleT, *pAngleR;              // Emerging angles
    double *pRayLengthT, *pRayLengthR;      // Ray Length
    double **pIm, **pImI, **pImQ;                           // Enveloppe Image
};
#endif
