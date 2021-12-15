#ifndef STRUCTURE_H
#define STRUCTURE_H

/* Structure containing info defining A SubAperture in a sense
 * of fixed subaperture defined by the user.
 * Not used when the code uses F-number constrain.
 * */
struct SubAperture_st {
    int **pElements;                // Idx of elements of the subs.
    int *pCenterIdx;                // Idx of the center of the sub.
    int nSub, nElemPerSub;          // Number of subs, elements per subs.
    double CenterCorrection = 0.;   // Indicates if a sub is even or odd (Correction = PITCH/2).
};
/* Structure containing Signal, Pixel position,
 *  Receivers and Transmetters and there dimsnesions
 * */
struct Input_st {
    double *pXR, *pZR, *pXS, *pZS, *pX, *pZ;
    double *pSignalI, *pSignalQ, *pAddDelay;
    double *pTxToUse, *pListenAngle;
    int NTX, NR, NS, NX, NZ, ZMax = 0, ZMin = 0;
    int TimeSig, SrcSig, RcptSig, NLA;
};

/* Structure containing Parabolic parameters
 * and its Minimum/Maximum.
 * */
struct Parabola_st {
    // y = ax2 + bx + c
    double a = 0, b = 0, c = 0;
    int IdxMin = 0, IdxMax = 0;
    // Number of Members within the structure.
    const int NParameters = 5;
};

/* Structure containing Model parameters, Local Pixel positions for different functions,
 * Experience settings, and Parameters for MEX behavious.
 * */
struct ExpSetup_st {

    const Parabola_st *pPeri; // Pointing to Parabola.
    const Parabola_st *pEndo;
    double XStart = 0, XEnd = 0, ZStart = 0, ZEnd = 0;              // Pixel poititions.
    double RatioLT = 0., RatioTB = 0., AngleCriticalMin = 0., AngleCriticalMax =0.;
    double AnisoCoef, AnisoShape, LensThick, LensPitch;             // Model/Media parameters.
    double DProbePerios, CorticalThickMin, CorticalThickMax;
    double HalfOpenAngle, Fnumber, MaxDevListenAngle;
    double CLens, CTissue, CMarrow, CBoneAxial, CBoneRadial;

    int SamplingFreq;
    int NCore, ReconTo, NeedTravelTime, SubAperApodis;

    // Number of Members to load from Matlab. The rest (5 members) are inferred.
    const int NSetupMember = 19;
};

/* Structure containing the built image and outputs (angle, time arrival, beamform)
 * */
struct Output_st {
        double *pTimeT, *pTimeR;                // Time Arrival Trans/Receiv
        double *pAngleT, *pAngleR;              // Emerging angles T,R
        double *pAngleAtPixT, *pAngleAtPixR;    // Angle at Pixel
        double **pImI, **pImQ;                  // I/Q Image
        double *pBeamI, *pBeamQ;                // I/Q Beamform
        double *pListenedAngleAtPixR;           // Effective angles at pixel in reception regarding Angle of Listenting.
};
#endif
