#include <string.h>
#include <fstream>
#include <iostream>
#include "buildImage.hpp"
#include "calculateTimeArrival.hpp"
#include "minimization.hpp"
#include "segmentation.hpp"
#include "utilities.hpp"
#include "structure.h"

#define __min(a,b) ((a)>(b) ? (b) : (a)) // function to find the minimum between two numbers
#define __max(a,b) ((a)<(b) ? (b) : (a)) // function to find the maximum between two numbers
#define NPLHS 6                          // Nb of Matlab outputs
#define NPRHS 12                         // Nb of Matlab inputs

using namespace std;

/* Main code */
void
mexFunction (int nlhs, mxArray * plhs[],
             int nrhs, const mxArray * prhs[])
{
    if (nrhs!=NPRHS) mexErrMsgIdAndTxt("Main:Input","Number of inputs must be %u", NPRHS);
    if (nlhs!=NPLHS) mexErrMsgIdAndTxt("Main:Output","Number of outputs must be %u", NPLHS);
    // Count time elapsed at each pointed steps.

    // Declaring Structures needed.
    ExpSetup_st Setup; // Contains Media Informations and Run Setup.
    Input_st    *Data; // Contains Pixel info and Raw Signal.
    Output_st   *ImageTissue, *ImageBone; // Contains Image to build (Time, Angle, Image, BeamForm)
    Parabola_st PeriParab; // Contains Parabola parameters and its Min.
    // Allocating dynamically. The rest to allocate is done on the flow after each segmentation.
    Data        = new Input_st;
    
    // Parameter we want to test
    double      *ModelToTest;
    
    /* Getting arguments and store them ############################################
     * #############################################################################
     * */

    int Counter = 0;
    fillSetup(Setup, prhs[Counter], 0); Counter++;

    Data->NR   = fillPointerFromMxArray(&(Data->pXR), prhs[Counter]); Counter++;
    int TestNR = fillPointerFromMxArray(&(Data->pZR), prhs[Counter]); Counter++;
    if (TestNR!=Data->NR) mexErrMsgIdAndTxt("Main:Input","XR and ZR have different dimensions");

    Data->NS   = fillPointerFromMxArray(&(Data->pXS), prhs[Counter]); Counter++;
    int TestNS = fillPointerFromMxArray(&(Data->pZS), prhs[Counter]); Counter++;
    if (TestNS!=Data->NS) mexErrMsgIdAndTxt("Main:Input","XS and ZS have different dimensions");

    Data->NX   = fillPointerFromMxArray(&(Data->pX), prhs[Counter]); Counter++;
    Data->NZ   = fillPointerFromMxArray(&(Data->pZ), prhs[Counter]); Counter++;
    if (Data->NX<=0 || Data->NZ<=0)    mexErrMsgIdAndTxt("Main:Input","X or Z are empty");
    if (Data->pZ[0] < Setup.LensThick) mexErrMsgIdAndTxt("Main:Input", "Z[0] < LensThickness");
    if (Data->pX[1] > 1. || Data->pZ[1] > 1. ||
                    //Data->pXS[1] > 1. || Data->pZS[1] > 1. ||
                    Data->pXR[1] > 1. || Data->pZR[1] > 1.)
            mexErrMsgIdAndTxt("Main:Input", "X, Z, XS, ZS, XR and ZR must be in meters");

    double *SignalI = fillSignalFromMxArray(Data, prhs[Counter], 0); Counter++;
    double *SignalQ = fillSignalFromMxArray(Data, prhs[Counter], 1); Counter++;

    Data->NTX  = fillPointerFromMxArray(&(Data->pTxToUse), prhs[Counter]); Counter++;
    int NDelay = fillPointerFromMxArray(&(Data->pAddDelay),prhs[Counter]); Counter++;
    if (NDelay!=Data->NS) mexErrMsgIdAndTxt("Main:Input","Nb of Sources and delays for virtual sources are different");

    // Refering Parabola to Setup structure.
    Setup.pPeri = &PeriParab;
    
    // Model to test
    int nModel  = fillPointerFromMxArray(&(ModelToTest),prhs[Counter]); Counter++;
    if (nModel<=0) mexErrMsgIdAndTxt("Main:Input","Nb of model to test is empty");

    if (Setup.WhichModel == 1) cout << "Testing for Axial Velocity" << endl;
    if (Setup.WhichModel == 2) cout << "Testing for Radial Velocity" << endl;
    if (Setup.WhichModel == 3) cout << "Testing for Aniso Shape Coef" << endl;
    if (Setup.WhichModel == 4) cout << "Testing for Isotropic Velocity model" << endl;
    cout << "Testing for " << nModel << " model values" << endl;
    
    // Resolution in both axes
    double 	ResolutionZ = Data->pZ[2]-Data->pZ[1];
    double 	ResolutionX = Data->pX[2]-Data->pX[1];

    /* Main code ################################################################
     * Allocation and Initialization  ###########################################
     * */

    // Dynamic allcation.
    double          *OutTissue, *OutBone, *OutInt, *OutVar, *OutGrad;
    ImageTissue     = new Output_st;
    ImageBone       = new Output_st;

    // Dimension of the image.
    const mwSize DImage[] = { (mwSize)Data->NZ, (mwSize)Data->NX };
    const mwSize DMetric[] = { (mwSize)nModel};

    // Allocation of plhs pointing to Matlab mxarray.
    Counter = 0;
    plhs[Counter] = mxCreateNumericArray(2, DImage, mxDOUBLE_CLASS, mxREAL);
    OutTissue     = mxGetPr(plhs[Counter]);
    initializeImage(Data, ImageTissue); Counter++;

    plhs[Counter] = mxCreateNumericArray(2, DImage, mxDOUBLE_CLASS, mxREAL);
    OutBone       = mxGetPr(plhs[Counter]);
    initializeImage(Data, ImageBone); Counter++;
    
    plhs[Counter] = mxCreateNumericArray(1, DMetric, mxDOUBLE_CLASS, mxREAL);
    OutInt        = mxGetPr(plhs[Counter]); Counter++;
    
    plhs[Counter] = mxCreateNumericArray(1, DMetric, mxDOUBLE_CLASS, mxREAL);
    OutVar        = mxGetPr(plhs[Counter]); Counter++;
    
    plhs[Counter] = mxCreateNumericArray(1, DMetric, mxDOUBLE_CLASS, mxREAL);
    OutGrad       = mxGetPr(plhs[Counter]); Counter++;
    /* Main code ################################################################
     * Tissue Imaging first #####################################################
     * */

    // Min and Max depths to Image, We do DAS only between ZMin and Zmax.
    Data->ZMin = 0;
	Data->ZMax = __min(round((Setup.DProbePerios - Data->pZ[0])/ResolutionZ), Data->NZ);
    // Time Arrival Calculation within Tissue.
    getTimeArrivalTissue(Setup, Data, ImageTissue);

    // Build Image within Tissue.
	buildImage(&Setup, Data, ImageTissue);

    
    double XDist, ZPerios, ZDist;
    // muting region deeper than assumed max depth of periosteum to help segmentation.
    for (int i = 0; i < Data->NX; i++)
    {
        for (int j = 0; j < Data->NZ; j++)
        {
            ZDist = ResolutionZ*j + Data->pZ[0];                // depth in m
            if (ZDist > Setup.DProbePerios) ImageTissue->pIm[j][i] = 0;
        }
    }
        
    
    // Segmentation of image Im to find the Perios (Dijkstra filter)
    segmentation(&PeriParab, ImageTissue->pIm, Data->NX, Data->NZ, Data->pX, Data->pZ);

    // Filling mxArray for Image displaying.
	for (int i = 0; i < Data->NX; i++)
	{
        for (int j = 0; j < Data->NZ; j++)
		{
            *OutTissue = ImageTissue->pIm[j][i];
			OutTissue++;
		}
	}
    // Output the Perios.
    plhs[Counter] = fillOutParabola(PeriParab);

    // Deallocating Image part.
    freeImage(ImageTissue, Data->NX, Data->NZ);
    delete [] ImageTissue;

    /* Main code ################################################################
     * Now we calculated the tissue part, we can test the bone model ############
     */
    
    Data->ZMin = __max(PeriParab.IdxMin + Setup.CorticalThickMin / ResolutionZ,0);
    if (PeriParab.IdxMax >= 0) Data->ZMax = __min(
        PeriParab.IdxMax + Setup.CorticalThickMax / ResolutionZ, Data->NZ );
    else Data->ZMax = Data->NZ;

    // Metric allocation
    double MaxInt = 0, MaxVarNorm = 0, MaxGrad = 0;
    double *pIntensity = (double*)calloc(nModel, sizeof(double));
    double *pVarNormIm = (double*)calloc(nModel, sizeof(double));
    double *pGradBrenn = (double*)calloc(nModel, sizeof(double));
    
    
    // Looping over models
    for (int v = 0; v < nModel; v++)
    {
        // model to test
        if (Setup.WhichModel==1)
        {
            Setup.CBoneAxial = ModelToTest[v];
            Setup.AnisoCoef  = (Setup.CBoneAxial - Setup.CBoneRadial)/Setup.CBoneAxial;
        }        
        else if (Setup.WhichModel==2)
        {
            Setup.CBoneRadial = ModelToTest[v];
            Setup.AnisoCoef  = (Setup.CBoneAxial - Setup.CBoneRadial)/Setup.CBoneAxial;
        }
        else if (Setup.WhichModel==3) Setup.AnisoShape  = ModelToTest[v];
        else if (Setup.WhichModel==4)
        {
            Setup.CBoneAxial = ModelToTest[v];
            Setup.CBoneRadial = ModelToTest[v];
            Setup.AnisoCoef  = 0;
        }
        
        // Set to zero travel time and Image tables before testing
        resetImage (Data, ImageBone);

        // Time Arrival Calculation within Tissue.
        getTimeArrivalBone(Setup, Data, ImageBone);

        // Build Image within Tissue.
        buildImage(&Setup, Data, ImageBone);

        
        // muting of pixels in ImageBone that are not included in region defined by CorticalThickMin CorticalThickMax.
        for (int i = 0; i < Data->NX; i++)
        {
            XDist   = ResolutionX*i + Data->pX[0]; // lateral distance in m (x-axis)
            ZPerios = (PeriParab.a*sqr(XDist) + PeriParab.b*XDist + PeriParab.c); // depth of the periosteum at lateral distance in m
            for (int j = 0; j < Data->NZ; j++)
            {
                ZDist = ResolutionZ*j + Data->pZ[0];
                if (ZDist < ZPerios+Setup.CorticalThickMin) ImageBone->pIm[j][i] = 0;
                if (ZDist > ZPerios+Setup.CorticalThickMax) ImageBone->pIm[j][i] = 0;
            }
        }
        
        
        
        // testing model 
        testMetric(ImageBone, Data,
        &pIntensity[v], &pVarNormIm[v], &pGradBrenn[v]);

        // Keep max value
        if (pIntensity[v] > MaxInt) MaxInt = pIntensity[v];
        if (pVarNormIm[v] > MaxVarNorm) MaxVarNorm = pVarNormIm[v];
        if (pGradBrenn[v] > MaxGrad) MaxGrad = pGradBrenn[v];
    }
    
    // Now we calculated everything, lets find the max of all tests.
    double MaxMetric = 0., SumMetric = 0.;
    int    idxBest   = -1;
    for (int v = 0; v < nModel; v++)
    {
        // Normalization
        pIntensity[v] /= MaxInt;
        pVarNormIm[v] /= MaxVarNorm;
        pGradBrenn[v] /= MaxGrad;
        
        // Give back the 3 metrics
        *OutInt  = pIntensity[v]; OutInt++;
        *OutVar  = pVarNormIm[v]; OutVar++;
        *OutGrad = pGradBrenn[v]; OutGrad++;
        
        // Sum of all normalized metrics
        SumMetric = pIntensity[v] + pVarNormIm[v] + pGradBrenn[v];
        
        // Keep only the max and its index
        if (SumMetric > MaxMetric)
        {
            MaxMetric = SumMetric;
            idxBest = v;
        }
    }

    // Now we need to calculate again the image for the output purpose
    if (Setup.WhichModel==1)
    {
        Setup.CBoneAxial = ModelToTest[idxBest];
        Setup.AnisoCoef  = (Setup.CBoneAxial - Setup.CBoneRadial)/Setup.CBoneAxial;
    }        
    else if (Setup.WhichModel==2)
    {
        Setup.CBoneRadial = ModelToTest[idxBest];
        Setup.AnisoCoef  = (Setup.CBoneAxial - Setup.CBoneRadial)/Setup.CBoneAxial;
    }
    else if (Setup.WhichModel==3) Setup.AnisoShape  = ModelToTest[idxBest];
    else if (Setup.WhichModel==4)
    {
        Setup.CBoneAxial = ModelToTest[idxBest];
        Setup.CBoneRadial = ModelToTest[idxBest];
    }

    cout << "Best model found = " << ModelToTest[idxBest] << endl;
    
    resetImage (Data, ImageBone);

    // Time Arrival calculation within the Bone.
    getTimeArrivalBone(Setup, Data, ImageBone);

    // Build of the Image within the Bone between the Perios and ZMax.
    buildImage(&Setup, Data, ImageBone);

    
    // muting of pixels in ImageBone that are not included in region defined by CorticalThickMin CorticalThickMax.
    for (int i = 0; i < Data->NX; i++)
    {
        XDist   = ResolutionX*i + Data->pX[0]; // lateral distance in m (x-axis)
        ZPerios = (PeriParab.a*sqr(XDist) + PeriParab.b*XDist + PeriParab.c); // depth of the periosteum at lateral distance in m
        for (int j = 0; j < Data->NZ; j++)
        {
            ZDist = ResolutionZ*j + Data->pZ[0];
            if (ZDist < ZPerios+Setup.CorticalThickMin) ImageBone->pIm[j][i] = 0;
            if (ZDist > ZPerios+Setup.CorticalThickMax) ImageBone->pIm[j][i] = 0;
        }
    }
    
    
    // Filling mxArray for Image displaying.
    for (int i = 0; i < Data->NX; i++)
    {
        for (int j = 0; j < Data->NZ; j++)
        {
            *OutBone = ImageBone->pIm[j][i];
            OutBone++;
        }
    }
    // Deallocating Image part.
    nulifyAndFreePtr(pIntensity);
    nulifyAndFreePtr(pVarNormIm);
    nulifyAndFreePtr(pGradBrenn);

    // Deallocating Image part.
    freeImage(ImageBone, Data->NX, Data->NZ);
    delete [] ImageBone;
    delete [] Data;
    delete [] SignalI;
    delete [] SignalQ;

} // End
