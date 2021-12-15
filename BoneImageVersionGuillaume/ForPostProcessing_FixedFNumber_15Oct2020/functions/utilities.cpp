#include "utilities.hpp"
using namespace std;

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
setArrayToZero(double* pArray, const int N)
{
    for (int i = 0; i < N; i++)
    {
        pArray[i] = 0.;
    }
    return NULL;
}

void
maxInArray(double *pTr, int N, double &Maximum, int &Index)
{
   Maximum = pTr[0];
   Index   = 0;
   for (int i = 1; i < N; i++)
   {
       if (pTr[i] > Maximum)
       {
           Maximum = pTr[i];
           Index   = i;
       }
   }
}

int**
newTable2D (const int N, const int M)
{
	int **pTr = new int * [N];
	for (int i = 0; i < N; i++)
	{
		pTr[i] = new int [M];
		memset(pTr[i], 0, M * sizeof(int));
	}
	return pTr;
}

void
deleteTable2D (int **pTr, const int N)
{
	for (int i = 0; i < N; i++)
	{
		delete [] pTr[i];
	}
	delete [] pTr;
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

void*
initializeImage(const Input_st *pData, Output_st *pImage)
{
    const int DimXZS   = pData->NX*pData->NZ*pData->NS;
    const int DimXZR   = pData->NX*pData->NZ*pData->NR;
    const int DimXZLA  = pData->NX*pData->NZ*pData->NLA;
    const int DimXZLAS = pData->NX*pData->NZ*pData->NLA*pData->NS;
    // Initiate by 0. the mxarray already allocated.
    // Allocation of NON MEX outputs.
    pImage->pAngleAtPixR = (double*)calloc(DimXZR, sizeof(double));
    // Allocation of MEX outputs
    for (int s = 0; s < DimXZS; s++)
    {
        pImage->pTimeT[s]       = 0.;
        pImage->pAngleT[s]      = 0.;
        pImage->pAngleAtPixT[s] = 0.;
    }
    for (int r = 0; r < DimXZR; r++)
    {
        pImage->pTimeR[r]        = 0.;
        pImage->pAngleR[r]       = 0.;
    }
    for (int b = 0; b < DimXZLAS; b++)
    {
        pImage->pBeamI[b] = 0.;
        pImage->pBeamQ[b] = 0.;
    }
    for (int rx = 0; rx < DimXZLA; rx++)
    {
        pImage->pListenedAngleAtPixR[rx]  = NAN; // For straighforward selection when doing VFI.
    }

    // Allocating for Non-mxarray type.
    pImage->pImI     = (double**)calloc(pData->NX,sizeof(double));
    pImage->pImQ     = (double**)calloc(pData->NX,sizeof(double));
    for (int x = 0; x < pData->NX; x++)
    {
        pImage->pImI[x] = (double*)calloc(pData->NZ,sizeof(double));
        pImage->pImQ[x] = (double*)calloc(pData->NZ,sizeof(double));
        for (int z = 0; z < pData->NZ; z++)
        {
            pImage->pImI[x][z] = 0.;
            pImage->pImQ[x][z] = 0.;
        }
    }
    return NULL;
}
void*
freeImage(Output_st *pImage, const int NX)
{
    // Free non MEX outputs. The rest has to stay for Matlab.
    free(pImage->pAngleAtPixR);
    for (int x = 0; x < NX; x++)
	{
		free(pImage->pImI[x]);
		free(pImage->pImQ[x]);
	}
    free(pImage->pImI);
    free(pImage->pImQ);
    return NULL;
}
void*
fillSetup(ExpSetup_st &rSetup, const mxArray *pMxArray, const int ShowInput)
{
    const mwSize *pMxArrayDims;
    double       *pMxData;
    pMxArrayDims = mxGetDimensions(pMxArray);
    if (pMxArrayDims[0]!=1 || pMxArrayDims[1]!=(mwSize)rSetup.NSetupMember)
    {
        mexErrMsgIdAndTxt("FillSetup:Input",
            "Argument #1 must be %u mxArray, See ExpSetup_st in include/Structures.h",
            rSetup.NSetupMember );
    }

    // Reading pointer to the first block of mxArray..
    pMxData = mxGetPr(pMxArray);

    // Filling Setup structure.
    int Counter = 0;
    rSetup.AnisoShape        = (double)pMxData[Counter]; Counter++;
    rSetup.LensThick         = (double)pMxData[Counter]; Counter++;
    rSetup.LensPitch         = (double)pMxData[Counter]; Counter++;
    rSetup.DProbePerios      = (double)pMxData[Counter]; Counter++;
    rSetup.CorticalThickMin  = (double)pMxData[Counter]; Counter++;
    rSetup.CorticalThickMax  = (double)pMxData[Counter]; Counter++;
    rSetup.HalfOpenAngle     = (double)pMxData[Counter]; Counter++;
    rSetup.Fnumber           = (double)pMxData[Counter]; Counter++;
    rSetup.MaxDevListenAngle = (double)pMxData[Counter]; Counter++;
    rSetup.SamplingFreq      = (double)pMxData[Counter]; Counter++;
    rSetup.CLens             = (double)pMxData[Counter]; Counter++;
    rSetup.CTissue           = (double)pMxData[Counter]; Counter++;
    rSetup.CMarrow           = (double)pMxData[Counter]; Counter++;
    rSetup.CBoneAxial        = (double)pMxData[Counter]; Counter++;
    rSetup.CBoneRadial       = (double)pMxData[Counter]; Counter++;
    rSetup.NCore             = (int)pMxData[Counter]; Counter++;
    rSetup.ReconTo           = (int)pMxData[Counter]; Counter++;
    rSetup.NeedTravelTime    = (int)pMxData[Counter]; Counter++;
    rSetup.SubAperApodis     = (int)pMxData[Counter]; Counter++;

    rSetup.AnisoCoef         = (rSetup.CBoneAxial - rSetup.CBoneRadial)/rSetup.CBoneAxial;
    rSetup.RatioLT           = rSetup.CLens/rSetup.CTissue;
    rSetup.RatioTB           = rSetup.CTissue/rSetup.CBoneAxial;
    rSetup.AngleCriticalMax  = asin(rSetup.RatioLT);
    rSetup.AngleCriticalMin  = -1*rSetup.AngleCriticalMax;
    
    // Testing values
    if (rSetup.ReconTo < 1 || rSetup.ReconTo > 3)
            mexErrMsgIdAndTxt("FillSetup:Input", "ReconTo must be =1,2,3");
    if (rSetup.SubAperApodis < 0 || rSetup.SubAperApodis > 2)
            mexErrMsgIdAndTxt("FillSetup:Input", "SubAperture Apodization must be = 0,1,2");
    if (rSetup.NeedTravelTime < 0 || rSetup.NeedTravelTime > 2)
            mexErrMsgIdAndTxt("FillSetup:Input", "NeedTravelTime must be = 0,1,2");
    if (rSetup.NCore < 1)
            mexErrMsgIdAndTxt("FillSetup:Input", "Number of CPU  must be > 0");
    if (rSetup.LensPitch < 0 || rSetup.LensPitch > 1e-3)
            mexErrMsgIdAndTxt("FillSetup:Input", "Pitch must be > 0 or in [m]");
    if (rSetup.DProbePerios < 0 || rSetup.DProbePerios > 1)
            mexErrMsgIdAndTxt("FillSetup:Input", "A Priori DProbePerios must be >= 0 or in [m]");
    if (rSetup.CorticalThickMin < 0 || rSetup.CorticalThickMin > 100e-3)
            mexErrMsgIdAndTxt("FillSetup:Input", "A Priori CorticalThickMin must be >= 0 or in [m]");
    if (rSetup.CorticalThickMax < 0 || rSetup.CorticalThickMax > 100e-3)
            mexErrMsgIdAndTxt("FillSetup:Input", "A Priori CorticalThickMax must be >= 0 or in [m]");
    if (rSetup.CorticalThickMax < rSetup.CorticalThickMin)
            mexErrMsgIdAndTxt("FillSetup:Input", "CorticalThickMin is > CorticalThickMax ");
    if (rSetup.AngleCriticalMax < rSetup.AngleCriticalMin)
            mexErrMsgIdAndTxt("FillSetup:Input", "AngleCriticalMin is > AngleCriticalMin ");
    if (rSetup.Fnumber <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "Fnumber must be > 0");
    if (rSetup.MaxDevListenAngle <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "MaxDeviationFromListenAngle must be > 0");
    if (rSetup.CLens <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "CLens must be > 0");
    if (rSetup.CTissue <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "CTissue must be > 0");
    if (rSetup.CMarrow <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "CMarrow must be > 0");
    if (rSetup.CBoneRadial <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "CBoneRadial must be > 0");
    if (rSetup.CBoneAxial <= 0)
            mexErrMsgIdAndTxt("FillSetup:Input", "CBoneAxial must be > 0");
    if (rSetup.SamplingFreq < 1e6 )
            mexErrMsgIdAndTxt("FillSetup:Input", "SamplingFreq must be in [Hz]");

    // See Inputs for debugging
    if (ShowInput==1)
    {
        cout << "AnisoCoef="        << rSetup.AnisoCoef << endl;
        cout << "AnisoShape="       << rSetup.AnisoShape << endl;
        cout << "LensThick="        << rSetup.LensThick << endl;
        cout << "LensPitch="        << rSetup.LensPitch << endl;
        cout << "RatioLT="          << rSetup.RatioLT << endl;
        cout << "RatioTB="          << rSetup.RatioTB << endl;
        cout << "AngleCriticalMax=" << rSetup.AngleCriticalMax << endl;
        cout << "AngleCriticalMin=" << rSetup.AngleCriticalMin << endl;
        cout << "DProbePerios="     << rSetup.DProbePerios << endl;
        cout << "CorticalThickMin=" << rSetup.CorticalThickMin << endl;
        cout << "CorticalThickMax=" << rSetup.CorticalThickMax << endl;
        cout << "HalfOpenAngle="    << rSetup.HalfOpenAngle << endl;
        cout << "SamplingFreq="     << rSetup.SamplingFreq << endl;
        cout << "CLens="            << rSetup.CLens << endl;
        cout << "CTissue="          << rSetup.CTissue << endl;
        cout << "CMarrow="          << rSetup.CMarrow << endl;
        cout << "CBoneAxial="       << rSetup.CBoneAxial << endl;
        cout << "CBoneRadial="      << rSetup.CBoneRadial << endl;
        cout << "NCore="            << rSetup.NCore << endl;
        cout << "ReconTo="          << rSetup.ReconTo << endl;
        cout << "NeedTravelTime="   << rSetup.NeedTravelTime << endl;
        cout << "SubAp. Apod.="     << rSetup.SubAperApodis << endl;
        cout << "F-number ="        << rSetup.Fnumber << endl;
        cout << "Max error Angle =" << rSetup.MaxDevListenAngle << endl;
    }
    return NULL;
}

int
fillParabola(Parabola_st &rParabola, const mxArray *pMxArray)
{
    int          TestParab = 0;
    const mwSize *pMxArrayDims;
    double       *pMxData;
    pMxArrayDims = mxGetDimensions(pMxArray);
    if (pMxArrayDims[0]!=1 || pMxArrayDims[1]!=(mwSize)rParabola.NParameters)
    {
        mexErrMsgIdAndTxt("FillParabola:Input",
                        "Argument # must be 1*5 mxArray, See Parabola_st in include/structures.h");
    }
    pMxData = mxGetPr(pMxArray);

    int Counter = 0;
    rParabola.a      = (double)pMxData[Counter]; Counter++;
    rParabola.b      = (double)pMxData[Counter]; Counter++;
    rParabola.c      = (double)pMxData[Counter]; Counter++;
    rParabola.IdxMin = (double)pMxData[Counter]; Counter++;
    rParabola.IdxMax = (double)pMxData[Counter]; Counter++;

    if (fabs(rParabola.a) + fabs(rParabola.b) + fabs(rParabola.c)
                    + fabs(rParabola.IdxMax) + fabs(rParabola.IdxMin) !=0) TestParab++;
    return TestParab;
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
void*
fillSubAperture(SubAperture_st *pData, const mxArray *pMxArray, const double Pitch)
{
    const mwSize *pMxArrayDims;
    mwSize       NDimArray;
    double       *pMxData;

    pMxArrayDims = mxGetDimensions(pMxArray);
    NDimArray    = mxGetNumberOfDimensions(pMxArray);
    if (NDimArray!=2)
    {
        mexErrMsgIdAndTxt("FillSubAperture:Input",
            "Argument for SubAerture must be nSub*nElemPerSub");
    }
    // Infer number of subs.
    pData->nSub       = pMxArrayDims[0];
    // Infer number of elements per subaperture (Assume they all have the same dimensions.
    pData->nElemPerSub = pMxArrayDims[1];
    // Allocation for elements in Subs.
    pData->pElements   = newTable2D(pData->nSub, pData->nElemPerSub);
    // Allocation for center.
    pData->pCenterIdx  = (int*)calloc(pData->nSub,sizeof(int));

    // fetching data.
    pMxData = mxGetPr(pMxArray);
    for (int j = 0; j < pData->nElemPerSub; j++)
    {
        // Getting elements.
        for (int i = 0; i < pData->nSub; i++)
        {
            pData->pElements[i][j] = (int)*pMxData -1; // From Matlab to C.
            ++pMxData;
            // Infer median of each subaperture.
            // Even number of elements in subarrays.
            if (pData->nElemPerSub%2 == 0)
            {
                pData->pCenterIdx[i]    = ( pData->pElements[i][ (pData->nElemPerSub-1)/2 ] +
                    pData->pElements[i][ pData->nElemPerSub/2 ] )/2;
                pData->CenterCorrection = Pitch/2; // Correction of center position is a half pitch.
            }
                // Odd number of elements in subarrays. No center Correction.
                else pData->pCenterIdx[i] = pData->pElements[i][ pData->nElemPerSub/2];
        }
    }
    return NULL;
}

void*
freeSubAperture(SubAperture_st *pData)
{
    free(pData->pCenterIdx);
    deleteTable2D(pData->pElements, pData->nSub);
    return NULL;
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
            "Argument with I/Q raw signal must be Time*nR*nTx");
    }
    pData->TimeSig = pMxArrayDims[0];
    pData->RcptSig = pMxArrayDims[1];
    pData->SrcSig  = pMxArrayDims[2];

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

mxArray*
fillOutParabola(Parabola_st &rParabola)
{
    const mwSize pMxArrayDim[] = {1,(mwSize)rParabola.NParameters};
    mxArray *pMxArray          = mxCreateNumericArray(2, pMxArrayDim, mxDOUBLE_CLASS, mxREAL);
    double  *pMxData           = mxGetPr(pMxArray);
    pMxData[0] = rParabola.a;
    pMxData[1] = rParabola.b;
    pMxData[2] = rParabola.c;
    pMxData[3] = (double)rParabola.IdxMin;
    pMxData[4] = (double)rParabola.IdxMax;

    return pMxArray;
}

void*
saveAngleAtPix(const Input_st *Data, Output_st *Image,
                const char* Reception, const char* Transmission)
{
    ofstream FileTransmiss(Transmission, ios::out | ios::binary | ios::trunc);
    ofstream FileReception(Reception,    ios::out | ios::binary | ios::trunc);
    double   *AngleR = NULL, *AngleT = NULL;

    AngleR = Image->pAngleAtPixR;
    AngleT = Image->pAngleAtPixT;

    if (FileTransmiss.is_open() && FileTransmiss.is_open())
    {
        for (int i = 0; i < Data->NX; ++i)
        {
            for (int j = Data->ZMin; j < Data->ZMax; ++j)
            {
                for (int r = 0; r < Data->NR; ++r)
                {
                    int OffsetR = i*Data->NR*Data->NZ + j*Data->NR + r;
                    FileReception.write (reinterpret_cast<char*>(&AngleR[OffsetR]), sizeof(double));
                }
                for (int s = 0; s < Data->NS; ++s)
                {
                    int OffsetT = i*Data->NS*Data->NZ + j*Data->NS + s;
                    FileTransmiss.write (reinterpret_cast<char*>(&AngleT[OffsetT]), sizeof(double));
                }
            }
        }
    }
    else mexErrMsgIdAndTxt("Utilitie:saveAngleAtPix", "Could not open file");
    FileReception.close();
    FileTransmiss.close();
    return NULL;
}
void*
loadAngleAtPix(Input_st *Data, Output_st *Image,
                const char* Reception, const char* Transmission)
{
    ifstream FileTransmiss(Transmission, ios::in | ios::binary);
    ifstream FileReception(Reception,    ios::in | ios::binary);
    double   *AngleR = NULL, *AngleT = NULL;

    AngleR = Image->pAngleAtPixR;
    AngleT = Image->pAngleAtPixT;

    if (FileTransmiss.is_open() && FileTransmiss.is_open())
    {
        for (int i = 0; i < Data->NX; ++i)
        {
            for (int j = Data->ZMin; j < Data->ZMax; ++j)
            {
                for (int r = 0; r < Data->NR; ++r)
                {
                    int OffsetR = i*Data->NR*Data->NZ + j*Data->NR + r;
                    FileReception.read (reinterpret_cast<char*>(&AngleR[OffsetR]), sizeof(double));
                }
                for (int s = 0; s < Data->NS; ++s)
                {
                    int OffsetT = i*Data->NS*Data->NZ + j*Data->NS + s;
                    FileTransmiss.read (reinterpret_cast<char*>(&AngleT[OffsetT]), sizeof(double));
                }
            }
        }
    }
    else mexErrMsgIdAndTxt("Utilitie:loadAngleAtPix", "Could not open file");
    FileReception.close();
    FileTransmiss.close();
    return NULL;
}

void*
saveTimeAndAngle(Input_st *Data, Output_st *Image,
                const char*  Reception, const char* Transmission)
{
    ofstream FileTransmiss(Transmission, ios::out | ios::binary | ios::trunc);
    ofstream FileReception(Reception,    ios::out | ios::binary | ios::trunc);
    double   *TimeT = NULL, *TimeR = NULL, *AngleR = NULL, *AngleT = NULL;

    TimeT  = Image->pTimeT; TimeR  = Image->pTimeR;
    AngleT = Image->pAngleT; AngleR = Image->pAngleR;
    if (FileTransmiss.is_open() && FileTransmiss.is_open())
    {
        for (int i = 0; i < Data->NX; ++i)
        {
            for (int j = Data->ZMin; j < Data->ZMax; ++j)
            {
                for (int r = 0; r < Data->NR; ++r)
                {
                    int OffsetR = i*Data->NR*Data->NZ + j*Data->NR + r;
                    FileReception.write (reinterpret_cast<char*>(&TimeR[OffsetR]),  sizeof(double));
                    FileReception.write (reinterpret_cast<char*>(&AngleR[OffsetR]), sizeof(double));
                }
                for (int s = 0; s < Data->NS; ++s)
                {
                    int offsetT = i*Data->NS*Data->NZ + j*Data->NS + s;
                    FileTransmiss.write (reinterpret_cast<char*>(&TimeT[offsetT]), sizeof(double));
                    FileTransmiss.write (reinterpret_cast<char*>(&AngleT[offsetT]), sizeof(double));
                }
            }
        }
    }
    else mexErrMsgIdAndTxt("Utilitie:saveTimeAndAngle", "Could not open file");
    FileReception.close();
    FileTransmiss.close();
    return NULL;
}

void*
loadTimeAndAngle(Input_st *Data, Output_st *Image,
                const char*  Reception, const char* Transmission)
{
    ifstream FileTransmiss(Transmission, ios::in | ios::binary);
    ifstream FileReception(Reception,    ios::in | ios::binary);
    double   *TimeT = NULL, *TimeR = NULL, *AngleR = NULL, *AngleT = NULL;

    TimeT  = Image->pTimeT; TimeR  = Image->pTimeR;
    AngleT = Image->pAngleT; AngleR = Image->pAngleR;
    if (FileTransmiss.is_open() && FileReception.is_open())
    {
        for (int i = 0; i < Data->NX; ++i)
        {
            for (int j = Data->ZMin; j < Data->ZMax; ++j)
            {
                for (int r = 0; r < Data->NR; ++r)
                {
                    int OffsetR = i*Data->NR*Data->NZ + j*Data->NR + r;
                    FileReception.read (reinterpret_cast<char*>(&TimeR[OffsetR]), sizeof(double));
                    FileReception.read (reinterpret_cast<char*>(&AngleR[OffsetR]), sizeof(double));
                }
                for (int s = 0; s < Data->NS; s++)
                {
                    int OffsetT = i*Data->NS*Data->NZ + j*Data->NS + s;
                    FileTransmiss.read (reinterpret_cast<char*>(&TimeT[OffsetT]), sizeof(double));
                    FileTransmiss.read (reinterpret_cast<char*>(&AngleT[OffsetT]), sizeof(double));
                }
            }
        }
    }
    else mexErrMsgIdAndTxt("Utilitie:loadTimeAndAngle", "Could not open file");
    FileReception.close();
    FileTransmiss.close();
    return NULL;
}

void*
saveParabola( Parabola_st *Para, const char* Parabola)
{
    ofstream FileParabola(Parabola, ios::out | ios::binary | ios::trunc);
    double A, B, C;
    int Min, Max;

    A = Para->a, B = Para->b, C = Para->c;
    Min = Para->IdxMin, Max = Para->IdxMax;
    if (FileParabola.is_open())
    {
        FileParabola.write (reinterpret_cast<char*>(&A),   sizeof(double));
        FileParabola.write (reinterpret_cast<char*>(&B),   sizeof(double));
        FileParabola.write (reinterpret_cast<char*>(&C),   sizeof(double));
        FileParabola.write (reinterpret_cast<char*>(&Min), sizeof(int));
        FileParabola.write (reinterpret_cast<char*>(&Max), sizeof(int));
    }
    else mexErrMsgIdAndTxt("Utilitie:saveParabola", "Could not open file");
    FileParabola.close();
    return NULL;
}

void*
loadParabola( Parabola_st *Para, const char* Parabola)
{
    ifstream FileParabola(Parabola, ios::in | ios::binary);
    double A, B, C;
    int Min, Max;

    if (FileParabola.is_open())
    {
        FileParabola.read (reinterpret_cast<char*>(&A),   sizeof(double));
        FileParabola.read (reinterpret_cast<char*>(&B),   sizeof(double));
        FileParabola.read (reinterpret_cast<char*>(&C),   sizeof(double));
        FileParabola.read (reinterpret_cast<char*>(&Min), sizeof(int));
        FileParabola.read (reinterpret_cast<char*>(&Max), sizeof(int));
    }
    else mexErrMsgIdAndTxt("Utilitie:loadParabola", "Could not open file");
    Para->a = A;
    Para->b = B;
    Para->c = C;
    Para->IdxMin = Min;
    Para->IdxMax = Max;
    FileParabola.close();
    return NULL;
}

void*
cleanAngleAtPixel(const ExpSetup_st *Setup, const Input_st *pData,
                Output_st *pImage, const int IdxLayer, const int ReconTo)
{
    if (IdxLayer!=0 && IdxLayer!=1 && IdxLayer!=2) mexErrMsgIdAndTxt("Cleaning Matrix", "Layer must be 1 or 2");
    double aPeri = Setup->pPeri->a;
    double bPeri = Setup->pPeri->b;
    double cPeri = Setup->pPeri->c;
    double aEndo = Setup->pEndo->a;
    double bEndo = Setup->pEndo->b;
    double cEndo = Setup->pEndo->c;

    switch (IdxLayer)
    {
        case 0 :
            // Tissue case Looping over pixels and put NaN when deeper than Periosteum.
            if (ReconTo>0) // If we dont segmentate=no parabola so we dont need cleaning.
            {
                for (int i = 0; i < pData->NX; i++)
                {

                    for(int j = 0; j < pData->NZ; j++)
                    {
                        // We clean if shallower than interface.
                        if (pData->pZ[j] > aPeri * sqr(pData->pX[i]) + bPeri * pData->pX[i] + cPeri)
                        {
                            for(int k = 0; k < pData->NLA; k++)
                            {
                                    pImage->pListenedAngleAtPixR[i*pData->NLA*pData->NZ + j*pData->NLA + k] = NAN;
                            }
                        }

                    }
                }
            }
            break;
        case 1 :
        // Bone case Looping over pixels and put NaN when shallowe than Periosteum and deeper than Endeosteum.
        if (ReconTo>1) // If we dont segmentate down to Marrow=no parabola so dont need cleaning.
        {
            for (int i = 0; i < pData->NX; i++)
            {
                for(int j = 0; j < pData->NZ; j++)
                {
                    // We clean if shallower than interface.
                    if (pData->pZ[j] < aPeri * sqr(pData->pX[i]) + bPeri * pData->pX[i] + cPeri ||
                    pData->pZ[j] > aEndo * sqr(pData->pX[i]) + bEndo * pData->pX[i] + cEndo)
                    {
                        for(int k = 0; k < pData->NLA; k++)
                        {
                            pImage->pListenedAngleAtPixR[i*pData->NLA*pData->NZ + j*pData->NLA + k] = NAN;
                        }
                    }

                }
            }
        }
        break;
        case 2 :
        // Marrow case Looping over pixels and put NaN when shallower than Endeosteum.
        for (int i = 0; i < pData->NX; i++)
        {
            for(int j = 0; j < pData->NZ; j++)
            {
                // We clean if shallower than interface.
                if (pData->pZ[j] < aEndo * sqr(pData->pX[i]) + bEndo * pData->pX[i] + cEndo)
                {
                    for(int k = 0; k < pData->NLA; k++)
                    {
                        pImage->pListenedAngleAtPixR[i*pData->NLA*pData->NZ + j*pData->NLA + k] = NAN;
                    }
                }
            }
        }
    }
    return NULL;
}