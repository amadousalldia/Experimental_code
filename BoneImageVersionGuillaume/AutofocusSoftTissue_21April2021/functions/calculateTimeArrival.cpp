#include "calculateTimeArrival.hpp"

using namespace std;

void*
getTimeArrivalTissue(const ExpSetup_st GlobalSetup, const Input_st *Data, Output_st *Image)
{
        int i, j, r, s;                                 // private idx for loops
        double AngleMin;                                // Angle that minimizes TravelTimes

        double **BrentRange = refineCriticalAngleTable(Data, GlobalSetup.AngleCriticalMin, GlobalSetup.AngleCriticalMax);

        omp_set_num_threads(GlobalSetup.NCore);
        #pragma omp parallel for schedule(static) \
        shared(Data,Image) private(i,j,r,s,AngleMin)
        for (i = 0; i < Data->NX; i++)
        {
                ExpSetup_st Setup; // Need to be copied to each thread since we update some structure's elements.
                memcpy(&Setup, &GlobalSetup, sizeof(ExpSetup_st));
                double *pTimeR  = Image->pTimeR  + i*Data->NR*Data->NZ;
                double *pAngleR = Image->pAngleR + i*Data->NR*Data->NZ;
                double *pRayLengthR = Image->pRayLengthR + i*Data->NR*Data->NZ;
                double *pTimeT  = Image->pTimeT  + i*Data->NS*Data->NZ;
                double *pAngleT = Image->pAngleT + i*Data->NS*Data->NZ;
                double *pRayLengthT = Image->pRayLengthT + i*Data->NS*Data->NZ;

                Setup.XEnd = Data->pX[i];
                for (j = Data->ZMin; j < Data->ZMax; j++)
                {
                        Setup.ZEnd = Data->pZ[j];
                        for (r = 0; r < Data->NR; r++)
                        {
                                Setup.XStart = Data->pXR[r];
                                Setup.ZStart = Data->pZR[r];
                                *pTimeR      = brent(Setup.AngleCriticalMin, 0, Setup.AngleCriticalMax,
                                                        travelTimeInTissue, &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);
                                *pAngleR     = AngleMin;
                                *pRayLengthR = RayLengthUntilTissue(AngleMin, &Setup);
                                pTimeR++; pAngleR++; pRayLengthR++;
                        }
                        for (s = 0; s < Data->NS; s++)
                        {
                                Setup.XStart = Data->pXS[s];
                                Setup.ZStart = Data->pZS[s];
                                *pTimeT      = brent(BrentRange[s][0], BrentRange[s][1], BrentRange[s][2],
                                                travelTimeInTissue, &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);
                                *pAngleT     = AngleMin;
                                *pRayLengthT = RayLengthUntilTissue(AngleMin, &Setup);
                                pTimeT++; pAngleT++; pRayLengthT++;
                        }
                }
        }
        free2D(BrentRange, Data->NS);
        return NULL;
}

double**
refineCriticalAngleTable(const Input_st *Data, const double CriticalMin, const double CriticalMax)
{
 
        double **BrentRangeTable = alloc2D(Data->NS, 3);
        double margin_x = (Data->pXR[Data->NR-1]-Data->pXR[0])/2; // added by GR, May 12, 2020
        for (int i = 0; i < Data->NS; i++)
        {// Range of research is widened.
            if ( (fabs(Data->pZR[0] - Data->pZS[i]) > 0.) & (fabs(Data->pZR[Data->NR-1] - Data->pZS[i]) > 0.))
            { // virtual sources
                // Angle Between Source and the 1st element of the probe
                BrentRangeTable[i][0] = __max(
                        atan( (Data->pXR[0] - margin_x - Data->pXS[i]) / (Data->pZR[0] - Data->pZS[i]) ),
                        CriticalMin );  // modified by GR, May 12, 2020, absolute min angle
                // Angle Between Source and the last element of the probe
                BrentRangeTable[i][2] = __min(
                        atan( (Data->pXR[Data->NR-1] + margin_x - Data->pXS[i]) / (Data->pZR[Data->NR-1] - Data->pZS[i]) ),
                        CriticalMax );  // modified by GR, May 12, 2020, absolute max angle
            }
            else
            { // physical sources
                BrentRangeTable[i][0] = CriticalMin;  // modified by GR, May 12, 2020, absolute min angle
                BrentRangeTable[i][2] = CriticalMax;  // modified by GR, May 12, 2020, absolute max angle
            }
            // Mean of 2 angles to get the starting point for Brent.
            BrentRangeTable[i][1] = (BrentRangeTable[i][2] + BrentRangeTable[i][0])/2;
        }
        return BrentRangeTable;
}
