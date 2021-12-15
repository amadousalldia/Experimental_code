#include "calculateTimeArrival.hpp"
using namespace std;

void*
getTimeArrivalTissue(const ExpSetup_st GlobalSetup, const Input_st *Data,
                Output_st *Image)
{
    int i, j, r, s, OffsetR, OffsetT; //idx;
    double AngleMin;

    double **BrentRange = refineCriticalAngleTable(Data, GlobalSetup.AngleCriticalMin, GlobalSetup.AngleCriticalMax);

    omp_set_num_threads(GlobalSetup.NCore);
    #pragma omp parallel for schedule(static) shared(Data,Image, BrentRange) \
    private(i, j, r, s, OffsetR, OffsetT, AngleMin)
    for (i = 0; i < Data->NX; i++)
    {
        ExpSetup_st Setup; // Need to copy for each thread or change all function arguments...
        memcpy(&Setup, &GlobalSetup, sizeof(ExpSetup_st));

        Setup.XEnd = Data->pX[i];
        for (j = Data->ZMin; j < Data->ZMax; j++)
        {
            // Reception
            Setup.ZEnd = Data->pZ[j];
            for (r = 0; r < Data->NR; r++)
            {
                Setup.XStart = Data->pXR[r];
                Setup.ZStart = Data->pZR[r];
                OffsetR      = i*Data->NR*Data->NZ + j*Data->NR + r;
                Image->pTimeR[OffsetR]   = brent(Setup.AngleCriticalMin, 0, Setup.AngleCriticalMax,
                                                 travelTimeInTissue, &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);
                Image->pAngleR[OffsetR]  = AngleMin;

                // Angle at pixel for all elements.
                // Selection for a given angle of listening is done during beaforming
                // in function/buildImage.cpp with getActiveSubAperture function.
                Image->pAngleAtPixR[OffsetR] = getAngleAtScatterWithinTissue(AngleMin, Setup.RatioLT);

            }
            for (s = 0; s < Data->NS; s++)
            {
                Setup.XStart = Data->pXS[s];
                Setup.ZStart = Data->pZS[s];
                OffsetT      = i*Data->NS*Data->NZ + j*Data->NS + s;
                Image->pTimeT[OffsetT]  = brent(BrentRange[s][0], BrentRange[s][1], BrentRange[s][2],
                                                travelTimeInTissue, &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);
                Image->pAngleT[OffsetT] = AngleMin;

                // Angle at pixel for all elements.
                // Selection for a given angle of listening is done during beaforming
                // in function/buildImage.cpp with getActiveSubAperture function.
                Image->pAngleAtPixT[OffsetT] = getAngleAtScatterWithinTissue(AngleMin, Setup.RatioLT);
            }
        }
    }
    free2D(BrentRange, Data->NS);
    return NULL;
}

void*
getTimeArrivalBone(const ExpSetup_st GlobalSetup, const Input_st *Data, Output_st *Image)
{
    int i, j, r, s, OffsetR, OffsetT;
    double AngleMin, ZParabola;

    double **BrentRange = refineCriticalAngleTable(Data, GlobalSetup.AngleCriticalMin, GlobalSetup.AngleCriticalMax);

    omp_set_num_threads(GlobalSetup.NCore);
    #pragma omp parallel for schedule(static) shared(Data,Image, BrentRange) \
    private(i, j, r, s, OffsetR, OffsetT, ZParabola, AngleMin)
    for (i = 0; i < Data->NX; i++)
    {
        ExpSetup_st Setup; // Need to copy for each thread or change all function arguments...
        memcpy(&Setup, &GlobalSetup, sizeof(ExpSetup_st));

        Setup.XEnd = Data->pX[i];
        // If we are shallower than the parabola, we dont need to calculate Times and Angles.
        ZParabola = (Setup.pPeri->a * sqr(Setup.XEnd) + Setup.pPeri->b * Setup.XEnd + Setup.pPeri->c);

        for (j =  Data->ZMin; j < Data->ZMax; j++)
        {
            if (Data->pZ[j] > ZParabola)
            {
                Setup.ZEnd = Data->pZ[j];
                // Reception
                for (r = 0; r < Data->NR; r++)
                {
                    OffsetR       = i*Data->NR*Data->NZ + j*Data->NR + r;
                    Setup.XStart  = Data->pXR[r];
                    Setup.ZStart  = Data->pZR[r];
                    Image->pTimeR[OffsetR]       = brent( Setup.AngleCriticalMin, 0., Setup.AngleCriticalMax,
                                                   travelTimeThroughTissueBone, &Setup, AccuracyBrent, &AngleMin,
                                                   IterMaxBrent);
                    // Angle at pixel for all elements.
                    // Selection for a given angle of listening is done during beaforming
                    // in function/buildImage.cpp with getActiveSubAperture function.
                    Image->pAngleAtPixR[OffsetR] = getAngleAtScatterWithinBone(AngleMin, &Setup);
                    Image->pAngleR[OffsetR]      = AngleMin;
                }
                for (s = 0; s < Data->NS; s++)
                {
                    OffsetT      = i*Data->NS*Data->NZ + j*Data->NS + s;
                    Setup.XStart = Data->pXS[s];
                    Setup.ZStart = Data->pZS[s];
                    Image->pTimeT[OffsetT]  = brent(BrentRange[s][0], BrentRange[s][1], BrentRange[s][2],
                                                    travelTimeThroughTissueBone, &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);

                    // Angle at pixel for all elements.
                    // Selection for a given angle of listening is done during beaforming
                    // in function/buildImage.cpp with getActiveSubAperture function.
                    Image->pAngleAtPixT[OffsetT] = getAngleAtScatterWithinBone(AngleMin, &Setup);
                    Image->pAngleT[OffsetT]      = AngleMin;
                }
            }
        }
    }
    free2D(BrentRange, Data->NS);
	return NULL;
}

void*
getTimeArrivalMarrow(const ExpSetup_st GlobalSetup, const Input_st *Data, Output_st *Image)
{
    int i, j, r, s, OffsetR, OffsetT, OffsetPriorR, OffsetPriorT;
    double AngleMin, AnglePrior = 0., ZParabola;

    double **BrentRange = refineCriticalAngleTable(Data, GlobalSetup.AngleCriticalMin, GlobalSetup.AngleCriticalMax);

    omp_set_num_threads(GlobalSetup.NCore);
    #pragma omp parallel for schedule(static) shared(Data,Image) \
    private(i, j, r, s, OffsetR, OffsetT, OffsetPriorR, OffsetPriorT, ZParabola, AngleMin, AnglePrior)
	for (i = 0; i < Data->NX; i++)
	{
        // Reception
        ExpSetup_st Setup; // Need to copy for each thread or change all function arguments...
        memcpy(&Setup, &GlobalSetup, sizeof(ExpSetup_st));

        Setup.XEnd = Data->pX[i];
        ZParabola  = (Setup.pEndo->a*sqr(Data->pX[i]) + Setup.pEndo->b * Data->pX[i] + Setup.pEndo->c);
        for (j = Data->ZMin; j < Data->ZMax; j++)
		{
            if (Data->pZ[j] > ZParabola)
            {
                Setup.ZEnd = Data->pZ[j];
                for (r = 0; r < Data->NR; r++)
                {
                    Setup.XStart = Data->pXR[r];
                    Setup.ZStart = Data->pZR[r];
                    OffsetR      = i*Data->NR*Data->NZ + j*Data->NR + r;
                    OffsetPriorR = i*Data->NR*Data->NZ + (j-1)*Data->NR + r;
                    AnglePrior   = Image->pAngleR[OffsetPriorR];

                    Image->pTimeR[OffsetR]  = runBrentWithPriorAngleEstimation( Setup, AngleMin, AnglePrior,
                                              travelTimeThroughTissueBoneMarrow);
                    Image->pAngleR[OffsetR] = AngleMin;

                    // Angle at pixel for all elements.
                    // Selection for a given angle of listening is done during beaforming
                    // in function/buildImage.cpp with getActiveSubAperture function.
                    Image->pAngleAtPixR[OffsetR] = getAngleAtScatterWithinMarrow(AngleMin, &Setup);
                }
                for (s = 0; s < Data->NS; s++)
                {
                    Setup.XStart = Data->pXS[s];
                    Setup.ZStart = Data->pZS[s];
                    OffsetT      = i*Data->NS*Data->NZ + j*Data->NS + s;
                    OffsetPriorT = i*Data->NS*Data->NZ + (j-1)*Data->NS + s;
                    AnglePrior   = Image->pAngleT[OffsetPriorT];

                    Image->pTimeT[OffsetT]       = runBrentWithPriorAngleEstimationTransmission( Setup, AngleMin, AnglePrior,
                                                   travelTimeThroughTissueBoneMarrow, BrentRange, s);
                    Image->pAngleT[OffsetT]      = AngleMin;

                    // Angle at pixel for all elements.
                    // Selection for a given angle of listening is done during beaforming
                    // in function/buildImage.cpp with getActiveSubAperture function.
                    Image->pAngleAtPixT[OffsetT] = getAngleAtScatterWithinMarrow(AngleMin, &Setup);
                }
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

double
runBrentWithPriorAngleEstimationTransmission(ExpSetup_st &Setup, double& AngleMin,
                double AnglePrior, double (*functionToTest)(double, void *),
                double **BrentRange, int IdxS)
{
    // If we provide Angleprior, we dont do the grid search
    if (AnglePrior==0.)
    {
        double Theta, dTheta, TravelTime, TTmin = 1.;
        // step of Angle through iteration.
        dTheta = (BrentRange[IdxS][2] - BrentRange[IdxS][0]) / (nTestPriorTx - 1);
        // testing nTestPrior values, keep the AnglePrior minimzing the functionToTest.
        // helping Brent to converge.
        for (int i = 0; i < nTestPriorRx; i++)
        {
            Theta       = BrentRange[IdxS][0] + i * dTheta;
            TravelTime  = (*functionToTest)(Theta, &Setup);

            if (TravelTime < TTmin)
            {
                    AnglePrior 	= Theta;
                    TTmin           = TravelTime;
            }
        }
    }
    // We keep using AnglePrior, always safer.
    return  brent(BrentRange[IdxS][0], AnglePrior, BrentRange[IdxS][2],
            (*functionToTest),
            &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);
}
double
runBrentWithPriorAngleEstimation(ExpSetup_st &Setup, double& AngleMin,
                double AnglePrior, double (*functionToTest)(double, void *))
{
    // If we provide Angleprior, we dont do the grid search
    if (AnglePrior==0.)
    {
        double Theta, dTheta, TravelTime, TTmin = 1.;
        // step of Angle through iteration.
        dTheta = (Setup.AngleCriticalMax - Setup.AngleCriticalMin) / (nTestPriorRx - 1);
        // testing nTestPrior values, keep the AnglePrior minimzing the functionToTest.
        // helping Brent to converge.
        for (int i = 0; i < nTestPriorRx; i++)
        {
            Theta       = Setup.AngleCriticalMin + i * dTheta;
            TravelTime  = (*functionToTest)(Theta, &Setup);

            if (TravelTime < TTmin)
            {
                AnglePrior 	= Theta;
                TTmin           = TravelTime;
            }
        }
    }

    return  brent(Setup.AngleCriticalMin, AnglePrior, Setup.AngleCriticalMax,
            (*functionToTest),
            &Setup, AccuracyBrent, &AngleMin, IterMaxBrent);
}
