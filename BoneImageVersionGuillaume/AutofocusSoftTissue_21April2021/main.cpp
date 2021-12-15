#include <string.h>
#include <fstream>
#include <iostream>
#include "buildImage.hpp"
#include "calculateTimeArrival.hpp"
#include "minimization.hpp"
#include "utilities.hpp"
#include "structure.h"

#define __min(a,b) ((a)>(b) ? (b) : (a)) // function to find the minimum between two numbers
#define __max(a,b) ((a)<(b) ? (b) : (a)) // function to find the maximum between two numbers

#define NPLHS 4                         // Nb of Matlab outputs
#define NPRHS 12                        // Nb of Matlab inputs

using namespace std;

/* Main code */
void
mexFunction (int nlhs, mxArray * plhs[],
             int nrhs, const mxArray * prhs[])
{
    if (nrhs!=NPRHS) mexErrMsgIdAndTxt("Main:Input","Number of inputs must be %u", NPRHS);
    if (nlhs!=NPLHS) mexErrMsgIdAndTxt("Main:Output","Number of outputs must be %u", NPLHS);

    // Declaring Structures needed.
    ExpSetup_st Setup;                                  // Contains Media Informations and Run Setup.
    Input_st    *Data;                                  // Contains Pixel info and Raw Signal.
    Output_st   *ImageTissue; 			    // Contains Image to build (Time, Angle, Image, BeamForm)
    // Allocating dynamically. The rest to allocate is done on the flow after each segmentation.
    Data        = new Input_st;

    double      *VelocityModel;
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

    // Reading velocity model to test
    int nVelocity = fillPointerFromMxArray(&(VelocityModel),prhs[Counter]); Counter++;
    if (nVelocity<=0) mexErrMsgIdAndTxt("Main:Input","Nb of model to test is empty");
    cout << "Testing for " << nVelocity << " model values" << endl;
    // Resolution in both axes
    double 	ResolutionZ = Data->pZ[2]-Data->pZ[1];
    double 	ResolutionX = Data->pX[2]-Data->pX[1];

    /* Main code ################################################################
     * Allocation and Initialization  ###########################################
     * */

    // Dynamic allcation.
    double          *OutTissue, *OutInt, *OutVar, *OutGrad;
    ImageTissue     = new Output_st;

    // Dimension of the outputs.
    const mwSize DImage[]  = { (mwSize)Data->NZ, (mwSize)Data->NX };
    const mwSize DMetric[] = { (mwSize)nVelocity};
    
    // Allocation of plhs pointing to Matlab mxarray.
    Counter = 0;
    plhs[Counter] = mxCreateNumericArray(2, DImage, mxDOUBLE_CLASS, mxREAL);
    OutTissue     = mxGetPr(plhs[Counter]);
    initializeImage(Data, ImageTissue);

    Counter++;
    plhs[Counter] = mxCreateNumericArray(1, DMetric, mxDOUBLE_CLASS, mxREAL);
    OutInt        = mxGetPr(plhs[Counter]);
    
    Counter++;
    plhs[Counter] = mxCreateNumericArray(1, DMetric, mxDOUBLE_CLASS, mxREAL);
    OutVar        = mxGetPr(plhs[Counter]);
    
    Counter++;
    plhs[Counter] = mxCreateNumericArray(1, DMetric, mxDOUBLE_CLASS, mxREAL);
    OutGrad       = mxGetPr(plhs[Counter]);
    
    // Metric allocation
    double MaxInt = 0, MaxVarNorm = 0, MaxGrad = 0;
    
    Data->ZMin = 0;
    Data->ZMax = Data->NZ;
    
    double *pIntensity = (double*)calloc(nVelocity, sizeof(double));
    double *pVarNormIm = (double*)calloc(nVelocity, sizeof(double));
    double *pGradBrenn = (double*)calloc(nVelocity, sizeof(double));
    for (int v = 0; v < nVelocity; v++)
    {
        // model to test
        Setup.CTissue = VelocityModel[v];
        // Set to zero travel time and Image tables before testing
        resetImage (Data, ImageTissue);

        // Time Arrival Calculation within Tissue.
        getTimeArrivalTissue(Setup, Data, ImageTissue);

        // Build Image within Tissue.
        buildImage(&Setup, Data, ImageTissue);

        // testing model 
        testMetric(ImageTissue, Data,
        &pIntensity[v], &pVarNormIm[v], &pGradBrenn[v]);

        // Keep max value
        if (pIntensity[v] > MaxInt) MaxInt = pIntensity[v];
        if (pVarNormIm[v] > MaxVarNorm) MaxVarNorm = pVarNormIm[v];
        if (pGradBrenn[v] > MaxGrad) MaxGrad = pGradBrenn[v];
    }
    
    // Now we calculated everything, lets find the max of all tests.
    double MaxMetric = 0., SumMetric = 0.;
    int    idxBest   = -1;
    for (int v = 0; v < nVelocity; v++)
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
    Setup.CTissue = VelocityModel[idxBest];
    cout << "Best model found = " << Setup.CTissue << endl;
    
    resetImage (Data, ImageTissue);

    // Time Arrival Calculation within Tissue.
    getTimeArrivalTissue(Setup, Data, ImageTissue);

    // Build Image within Tissue.
    buildImage(&Setup, Data, ImageTissue);

    // Filling mxArray for Image displaying.
    for (int i = 0; i < Data->NX; i++)
    {
        for (int j = 0; j < Data->NZ; j++)
        {
            *OutTissue = ImageTissue->pIm[j][i];
            OutTissue++;
        }
    }
    // Deallocating Image part.
    nulifyAndFreePtr(pIntensity);
    nulifyAndFreePtr(pVarNormIm);
    nulifyAndFreePtr(pGradBrenn);
    freeImage(ImageTissue, Data->NX, Data->NZ);
    delete [] ImageTissue;
    delete [] Data;
    delete [] SignalI;
    delete [] SignalQ;

} // End
