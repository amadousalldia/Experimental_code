#include "utilities.hpp"
using namespace std;

void*
nulifyAndFreePtr( double* pTr )
{
    pTr = NULL;
    free(pTr);
    return NULL;
}
void*
resetImage(const Input_st *pData, Output_st *pImage)
{
    int DimXZS = pData->NX*pData->NZ*pData->NS;
    int DimXZR = pData->NX*pData->NZ*pData->NR;
    for (int x = 0; x < pData->NX; x++)
    {
        for (int z = 0; z < pData->NZ; z++)
        {
            pImage->pImI[x][z] = 0.;
            pImage->pImQ[x][z] = 0.;
            pImage->pIm[z][x] = 0.;
        }
    }
    for (int xzs = 0; xzs < DimXZS; xzs++)
    {
        pImage->pTimeT[xzs]  = 0.;
        pImage->pAngleT[xzs] = 0.;
        pImage->pRayLengthT[xzs] = 1.;
    }
    for (int xzr = 0; xzr < DimXZR; xzr++)
    {
        pImage->pTimeR[xzr]  = 0.;
        pImage->pAngleR[xzr] = 0.;
        pImage->pRayLengthR[xzr] = 1.;
    }
    return NULL;
}

double **
alloc2D (const int N, const int M)
{
	double **pTr = new double * [N];
	for (int i = 0; i < N; i++)
	{
		pTr[i] = new double [M];
		memset(pTr[i], 0, M*sizeof(double));
	}
	return pTr;
}

void
free2D (double **pTr, const int N)
{
	for (int i = 0; i < N; i++)
	{
		delete [] pTr[i];
	}
	delete [] pTr;
}

double
sqr (const double& x)
{
	return x * x;
}

double
distance (const double& x1, const double& z1, const double& x2, const double& z2)
{
    return sqrt( sqr(x2 - x1) + sqr(z2 - z1) );
}

void*
initializeImage(const Input_st *pData, Output_st *pImage)
{
    int DimXZS      = pData->NX*pData->NZ*pData->NS;
    int DimXZR      = pData->NX*pData->NZ*pData->NR;
    pImage->pTimeT  = (double*)calloc(DimXZS,sizeof(double));
    pImage->pTimeR  = (double*)calloc(DimXZR,sizeof(double));
    pImage->pAngleT = (double*)calloc(DimXZS,sizeof(double));
    pImage->pAngleR = (double*)calloc(DimXZR,sizeof(double));
    pImage->pRayLengthT = (double*)calloc(DimXZS,sizeof(double));
    pImage->pRayLengthR = (double*)calloc(DimXZR,sizeof(double));
    pImage->pImI    = (double**)calloc(pData->NX,sizeof(double*));
    pImage->pImQ    = (double**)calloc(pData->NX,sizeof(double*));
    for (int x = 0; x < pData->NX; x++)
    {
            pImage->pImI[x] = (double*)calloc(pData->NZ,sizeof(double));
            pImage->pImQ[x] = (double*)calloc(pData->NZ,sizeof(double));
    }
    // The ordering of pIm is inverted compared to I/Q. Faster for output (no reshape)
    // Thus we must follow the Matlab convention.
    pImage->pIm     = (double**)calloc(pData->NZ,sizeof(double*));
    for (int z = 0; z < pData->NZ; z++)
    {
            pImage->pIm[z]  = (double*)calloc(pData->NX,sizeof(double));
    }
    return NULL;
}

void*
freeImage(Output_st *pImage, const int NX, const int NZ)
{
    free(pImage->pTimeT);
    free(pImage->pTimeR);
    free(pImage->pAngleT);
    free(pImage->pAngleR);
    free(pImage->pRayLengthT);
    free(pImage->pRayLengthR);
    for (int x = 0; x < NX; x++)
	{
		free(pImage->pImI[x]);
		free(pImage->pImQ[x]);
	}
        free(pImage->pImI);
        free(pImage->pImQ);

    for (int z = 0; z < NZ; z++)
	{
		free(pImage->pIm[z]);
	}
    free(pImage->pIm);

    return NULL;
}

void*
fillSetup(ExpSetup_st &rSetup, const mxArray *pMxArray, const int ShowInput)
{
    const mwSize *pMxArrayDims;
    double       *pMxData;
    pMxArrayDims = mxGetDimensions(pMxArray);
    if (pMxArrayDims[0]!=1 || pMxArrayDims[1]!=rSetup.NSetupMember)
    {
        mexErrMsgIdAndTxt("FillSetup:Input",
            "Argument #1 must be %u mxArray, See ExpSetup_st in include/Structures.h",
            rSetup.NSetupMember );
    }

    // Reading pointer to the first block of mxArray..
    pMxData = mxGetPr(pMxArray);

    // Filling Setup structure.
    int Counter = 0;
    rSetup.LensThick         = (double)pMxData[Counter]; Counter++;
    rSetup.LensPitch         = (double)pMxData[Counter]; Counter++;
    rSetup.HalfOpenAngle     = (double)pMxData[Counter]; Counter++;
    rSetup.SamplingFreq      = (double)pMxData[Counter]; Counter++;
    rSetup.CLens             = (double)pMxData[Counter]; Counter++;
    rSetup.CTissue           = (double)pMxData[Counter]; Counter++;
    rSetup.NCore             = (int)pMxData[Counter]; Counter++;

    rSetup.AngleCriticalMax  = asin(rSetup.CLens/rSetup.CTissue);
    rSetup.AngleCriticalMin  = -1*rSetup.AngleCriticalMax;
    // Checking if all Parameters are in meters.
    if (rSetup.LensThick > 1.) mexErrMsgIdAndTxt("FillSetup:Input","LensThick must be in meters");
    if (rSetup.LensPitch > 1.) mexErrMsgIdAndTxt("FillSetup:Input","LensPitch must be in meters");

    // Display parameters for checking
    if (ShowInput==1)
    {
            cout << "LensThick="        << rSetup.LensThick << endl;
            cout << "LensPitch="        << rSetup.LensPitch << endl;
            cout << "AngleCriticalMax=" << rSetup.AngleCriticalMax << endl;
            cout << "AngleCriticalMin=" << rSetup.AngleCriticalMin << endl;
            cout << "HalfOpenAngle="    << rSetup.HalfOpenAngle << endl;
            cout << "SamplingFreq="     << rSetup.SamplingFreq << endl;
            cout << "CLens="            << rSetup.CLens << endl;
            cout << "CTissue="          << rSetup.CTissue << endl;
            cout << "NCore="            << rSetup.NCore << endl;
    }
    return NULL;
}

int
fillPointerFromMxArray(double **pData, const mxArray *pMxArray)
{
    const mwSize *pMxArrayDims;
    pMxArrayDims = mxGetDimensions(pMxArray);
    if (pMxArrayDims[0]!=1 || pMxArrayDims[1] <=0)
    {
            mexErrMsgIdAndTxt("FillPointerFromMxArray:Input",
                            "Argument with Array must be 1*N mxArray");
    }
    *pData = mxGetPr(pMxArray);
    return (int)pMxArrayDims[1];
}

double*
fillSignalFromMxArray (Input_st *pData, const mxArray *pMxArray, const int& SignalPart)
{
    const mwSize *pMxArrayDims;
    mwSize       NDimArray;
    double       *pMxData;
    pMxArrayDims = mxGetDimensions(pMxArray);
    NDimArray    = mxGetNumberOfDimensions(pMxArray);
    if (NDimArray!=3)
    {
            mexErrMsgIdAndTxt("FillSignalFromMxArray:Input",
                            "Argument with I/Q raw signal must be N*M*M");
    }
    pData->TimeSig = pMxArrayDims[0];
    pData->RcptSig = pMxArrayDims[1];
    pData->SrcSig  = pMxArrayDims[2];
    if (pData->TimeSig<=0 || pData->RcptSig<=0 || pData->SrcSig<=0)
    {
            mexErrMsgIdAndTxt("FillSignalFromMxArray:Input",
                            "Signal I/Q has a wrong dimension");
    }

    pMxData = mxGetPr(pMxArray);

    double *Signal = new double [pData->TimeSig*pData->RcptSig*pData->SrcSig](); // Allocate and Initialize
	for (int s = 0; s < pData->SrcSig; ++s)
	{
		for (int r = 0; r < pData->RcptSig; ++r)
		{
			for (int t = 0; t < pData->TimeSig; ++t)
			{
                int offset      = s*pData->TimeSig*pData->RcptSig + r*pData->TimeSig + t;
                Signal[offset]  = *pMxData;
				++pMxData;
			}
		}
	}
    // Choosing between I/Q output.
    if (SignalPart==0) pData->pSignalI = Signal;
    else pData->pSignalQ = Signal;
    return Signal;
}
