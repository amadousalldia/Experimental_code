#include "buildImage.hpp"
using namespace std;

/* sinc function */
double sinc(double x)
{
    x *= M_PI;
    if (x == 0.0) return 1.0;
    return sin(x)/x;
} /* sinc */

void*
buildImage(ExpSetup_st *Setup, Input_st *Data, Output_st *Image)
{
    int i, j, iTx, Tx, r, Index, OffsetT, OffsetSig, OffsetR;
    double SumI, SumQ, TotalTime, Weight, TimeTransmit, Tx_weigthing, Rx_weigthing, el_width_over_wavelength;

    el_width_over_wavelength = 1;
            
    omp_set_num_threads(Setup->NCore);
    #pragma omp parallel for schedule(static) \
    shared(Setup, Data, Image) \
    private( i, j, iTx, Tx, r, SumI, SumQ, \
                    TimeTransmit, OffsetR, OffsetT, \
                    TotalTime, Index, OffsetSig, Weight, Tx_weigthing, Rx_weigthing, el_width_over_wavelength)
                    
    for (i = 0; i < Data->NX; ++i)
	{
        for (j = Data->ZMin; j < Data->ZMax; ++j)
        {
            for (iTx = 0; iTx < Data->NTX; ++iTx)
            {
                Tx           = Data->pTxToUse[iTx] - 1;
                OffsetT      = i*Data->NS*Data->NZ + j*Data->NS +Tx;
                TimeTransmit = Image->pTimeT[OffsetT] - Data->pAddDelay[Tx];
                SumI = 0., SumQ = 0.;
                Tx_weigthing = cos(Image->pAngleT[OffsetT])*sinc(el_width_over_wavelength*sin(Image->pAngleT[OffsetT]))/Image->pRayLengthT[OffsetT];
                        
                for (r = 0; r < Data->NR; ++r)
                {
                    OffsetR = i*Data->NR*Data->NZ + j*Data->NR +r;
                    if (fabs(Image->pAngleR[OffsetR]) < Setup->HalfOpenAngle)
                    { // if angle lies within acceptance angle for this element
                        TotalTime = (Image->pTimeR[OffsetR] + TimeTransmit) * Setup->SamplingFreq;
                        Index     = int(TotalTime);
                        Rx_weigthing = cos(Image->pAngleR[OffsetR])*sinc(el_width_over_wavelength*sin(Image->pAngleR[OffsetR]))/Image->pRayLengthR[OffsetR];

                        if (Index >= 1 && Index < Data->TimeSig -1)
                        { // get data
                            OffsetSig = Tx*Data->TimeSig*Data->NR + r*Data->TimeSig + Index;
                            Weight    = 1. - (TotalTime - floor(TotalTime));
                            SumI     += (Data->pSignalI[OffsetSig] * Weight +
                                        Data->pSignalI[OffsetSig + 1] * (1. - Weight))*Rx_weigthing; // linear interpolation
                            SumQ     += (Data->pSignalQ[OffsetSig] * Weight +
                                        Data->pSignalQ[OffsetSig + 1] * (1. - Weight))*Rx_weigthing; // linear interpolation
                        }
                    }
                }
                Image->pImI[i][j] += SumI*Tx_weigthing;
                Image->pImQ[i][j] += SumQ*Tx_weigthing;
                //Image->pIm[j][i] += sqrt(sqr(SumI) + sqr(SumQ))*Tx_weigthing; // incoherent compounding
            }
            // For Matlab/C++ convention, i/j and j/i are inverted for the total image.
            Image->pIm[j][i] = sqrt(sqr(Image->pImI[i][j]) + sqr(Image->pImQ[i][j]));
        }
    }
    return NULL;
}
