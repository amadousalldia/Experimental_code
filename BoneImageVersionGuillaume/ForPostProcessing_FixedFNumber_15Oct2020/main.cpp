// Main code for Image reconstruction //
#include <string.h>
#include <fstream>
#include <iostream>
#include <chrono>

// Global macro.
#define __min(a,b) ((a)>(b) ? (b) : (a)) // function to find the minimum between two numbers
#define __max(a,b) ((a)<(b) ? (b) : (a)) // function to find the maximum between two numbers

#include "buildImage.hpp"
#include "calculateTimeArrival.hpp"
#include "minimization.hpp"
#include "segmentation.hpp"
#include "structure.h"
#include "utilities.hpp"

#define Debug false
#define NLHSTissue 9
#define NLHSBone   20 // modified by GR, October 15 2020
#define NLHSMarrow 29
#define NPRHS      14

using namespace std;

void
mexFunction (int nlhs, mxArray * plhs[],
             int nrhs, const mxArray * prhs[])
{
    if (nrhs!=NPRHS) mexErrMsgIdAndTxt("Main:Input","Number of inputs must be %u", NPRHS);

    // Declaring Strctures needed.
    ExpSetup_st    Setup;                                  // Contains Media Informations and Run Setup.
    Input_st       *Data;                                  // Contains Pixel info and Raw Signal.
    Output_st      *ImageTissue, *ImageBone, *ImageMarrow; // Contains Image to build (Time, Angle, Image, BeamForm)
    Parabola_st    PeriParab, EndoParab;                   // Contains Parabola parameters and Min.

    // Allocating dynamically.
    Data        = new Input_st;
    ImageTissue = new Output_st;
    ImageBone   = new Output_st;
    ImageMarrow = new Output_st;

    const double PI = 3.14159265358979323846;
    /* Getting arguments and store them ############################################
     * #############################################################################
     * */

    // Counter for input reading.
    int Counter = 0;

    // Read Setup of experience.
    fillSetup(Setup, prhs[Counter], 0); Counter++;

    // Read Receptor positions.
    Data->NR   = fillPointerFromMxArray(&(Data->pXR), prhs[Counter]); Counter++;
    int TestNR = fillPointerFromMxArray(&(Data->pZR), prhs[Counter]); Counter++;
    if (TestNR!=Data->NR) mexErrMsgIdAndTxt("Main:Input","XR and ZR have different dimensions");

    // Read Sources/Tx positions.
    Data->NS   = fillPointerFromMxArray(&(Data->pXS), prhs[Counter]); Counter++;
    int TestNS = fillPointerFromMxArray(&(Data->pZS), prhs[Counter]); Counter++;
    if (TestNS!=Data->NS) mexErrMsgIdAndTxt("Main:Input","XS and ZS have different dimensions");

    // Read pixel positions.
    Data->NX   = fillPointerFromMxArray(&(Data->pX), prhs[Counter]); Counter++;
    Data->NZ   = fillPointerFromMxArray(&(Data->pZ), prhs[Counter]); Counter++;
    if (Data->NX<=0 || Data->NZ<=0) mexErrMsgIdAndTxt("Main:Input","X or Z are empty");

    // Read I/Q signal parts.
    double *SignalI = fillSignalFromMxArray(Data, prhs[Counter], 0); Counter++;
    double *SignalQ = fillSignalFromMxArray(Data, prhs[Counter], 1); Counter++;

    // Read Tx to use.
    Data->NTX = fillPointerFromMxArray(&(Data->pTxToUse), prhs[Counter]); Counter++;
    if (Data->NTX > Data->NS) mexErrMsgIdAndTxt("Main:Input", "TxtoUse cannot be higher than NS");

    // Read delays of virtual sources.
    int testAddDel = fillPointerFromMxArray(&(Data->pAddDelay),prhs[Counter]); Counter++;
    if (Data->NTX!=testAddDel) mexErrMsgIdAndTxt("Main:Input","AddDelayToFireing and TStart must have same length");

    // Reading listening angles.
    Data->NLA = fillPointerFromMxArray(&(Data->pListenAngle), prhs[Counter]); Counter++;
    if (Data->NLA <= 0) mexErrMsgIdAndTxt("Main:Input", "Listening angles input is empty");
    for (int i = 0; i < Data->NLA; i++)
    {
        if (fabs(Data->pListenAngle[i]) > PI) mexErrMsgIdAndTxt("Main:Input", "Listening angles must be in radian");
    }

    // Read Parabolae, returns 1 if we have A Priori, Otherwise 0.
    int TestPeri = fillParabola(PeriParab, prhs[Counter]); Counter++;
    int TestEndo = fillParabola(EndoParab, prhs[Counter]); Counter++;
    // Refering Parabola to Setup structre.
    Setup.pEndo = &EndoParab;
    Setup.pPeri = &PeriParab;

    if (Setup.NeedTravelTime==2 && (TestPeri==1 || TestEndo==1)) mexErrMsgIdAndTxt("Main:Inputs",
            "When we dont need to calculate Travel times (NeedTravelTime=2), provide a zero vector for parabolae! \n"
            "The parabolae parameters are actually loaded from .bin files to keep consitency. \n");

    // Testing how many nlhs we provided.
    if (Setup.ReconTo==1 && nlhs!=NLHSTissue) mexErrMsgIdAndTxt("Main:Output","Number of outputs must be %u", NLHSTissue);
    if (Setup.ReconTo==2 && nlhs!=NLHSBone)   mexErrMsgIdAndTxt("Main:Output","Number of outputs must be %u", NLHSBone);
    if (Setup.ReconTo==3 && nlhs!=NLHSMarrow) mexErrMsgIdAndTxt("Main:Output","Number of outputs must be %u", NLHSMarrow);

    /* Allocation of Output arguments for Matlab ################################################################
     * ##########################################################################################################
     * ##########################################################################################################
     * ##########################################################################################################
     * */
    if (Debug) cout << "Inputs read" << endl;
    // Declare output images and parabola parameters
    // The row/column convention is inverted from C++ to Matlab.
    const mwSize DimZX[]    = { (mwSize)Data->NZ, (mwSize)Data->NX };                                               // Image Dimension
    const mwSize DimSZX[] 	= { (mwSize)Data->NS, (mwSize)Data->NZ, (mwSize)Data->NX};                              // Transmission angle/time matrix
    const mwSize DimRZX[] 	= { (mwSize)Data->NR, (mwSize)Data->NZ, (mwSize)Data->NX };                             // Reception
    const mwSize DimLAZX[] 	= { (mwSize)Data->NLA, (mwSize)Data->NZ, (mwSize)Data->NX };                             // Reception
    const mwSize DimZXLAS[] = { (mwSize)Data->NZ, (mwSize)Data->NX, (mwSize)Data->NLA, (mwSize)Data->NS };

    // Pointers for plhs.
    double **ImTissue = NULL, **ImBone = NULL, **ImMarrow = NULL;
    double *outTissue = NULL, *outBone = NULL, *outMarrow = NULL, *pBuffer = NULL;
    if (Setup.ReconTo >= 1)
	{
        // Image in soft tissue
        Counter = 0;
        plhs[Counter]   = mxCreateNumericArray(2, DimZX, mxDOUBLE_CLASS, mxREAL);
        outTissue       = mxGetPr(plhs[Counter]); Counter++;
        ImTissue        = alloc2D(Data->NZ, Data->NX);

        // transmit times
        plhs[Counter]         	= mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer                 = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pTimeT     = pBuffer;

        // receive times
        plhs[Counter]           = mxCreateNumericArray(3, DimRZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer                 = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pTimeR     = pBuffer;

        // transmit angles
        plhs[Counter] 	        = mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer	                = mxGetPr(plhs[Counter]); Counter ++;
        ImageTissue->pAngleT    = pBuffer;

        // receive angles
        plhs[Counter]           = mxCreateNumericArray(3, DimRZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer                 = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pAngleR    = pBuffer;

        // beamformed image for each Tx/Rx pair
        plhs[Counter]           = mxCreateNumericArray(4, DimZXLAS, mxDOUBLE_CLASS, mxREAL);
        pBuffer                 = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pBeamI     = pBuffer;

        plhs[Counter]           = mxCreateNumericArray(4, DimZXLAS, mxDOUBLE_CLASS, mxREAL);
        pBuffer                 = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pBeamQ     = pBuffer;

        // Output the angle at scatters
        plhs[Counter] 		  = mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer               = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pAngleAtPixT = pBuffer;

        plhs[Counter] 		  = mxCreateNumericArray(3, DimLAZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer               = mxGetPr(plhs[Counter]); Counter++;
        ImageTissue->pListenedAngleAtPixR = pBuffer;

        initializeImage(Data, ImageTissue);
    }
    if (Setup.ReconTo >=2)
    {
        plhs[Counter] = mxCreateNumericArray(2, DimZX, mxDOUBLE_CLASS, mxREAL);
        outBone	      = mxGetPr(plhs[Counter]); Counter++;
        ImBone 	      = alloc2D(Data->NZ, Data->NX);

        // transmit times
        plhs[Counter] 	   = mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer            = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pTimeT  = pBuffer;

        // receive times
        plhs[Counter] 	   = mxCreateNumericArray(3, DimRZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer            = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pTimeR  = pBuffer;

        // transmit angles
        plhs[Counter] 	   = mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer            = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pAngleT = pBuffer;

        // receive angles
        plhs[Counter] 	   = mxCreateNumericArray(3, DimRZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer	           = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pAngleR = pBuffer;

        // beamformed image for each Tx/Rx pair
        plhs[Counter]      = mxCreateNumericArray(4, DimZXLAS, mxDOUBLE_CLASS, mxREAL);
        pBuffer            = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pBeamI  = pBuffer;

        plhs[Counter]      = mxCreateNumericArray(4, DimZXLAS, mxDOUBLE_CLASS, mxREAL);
        pBuffer            = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pBeamQ  = pBuffer;

        // Output the angle at scatters
        plhs[Counter] 		= mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer            	= mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pAngleAtPixT = pBuffer;

        plhs[Counter] 		= mxCreateNumericArray(3, DimLAZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer             = mxGetPr(plhs[Counter]); Counter++;
        ImageBone->pListenedAngleAtPixR = pBuffer;

        initializeImage(Data, ImageBone);
    }
    if (Setup.ReconTo >=3)
	{
        plhs[Counter]  = mxCreateNumericArray(2, DimZX, mxDOUBLE_CLASS, mxREAL);
        outMarrow      = mxGetPr(plhs[Counter]); Counter++;
        ImMarrow       = alloc2D(Data->NZ, Data->NX);

        // transmit times
        plhs[Counter] 		= mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer 	        = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pTimeT = pBuffer;

        // receive times
        plhs[Counter] 		= mxCreateNumericArray(3, DimRZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer	        	= mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pTimeR = pBuffer;

        // transmit angles
        plhs[Counter] 		 = mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer              = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pAngleT = pBuffer;

        // receive angles
        plhs[Counter] 		 = mxCreateNumericArray(3, DimRZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer 		     = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pAngleR = pBuffer;

        // beamformed image feor each Tx/Rx pair
        plhs[Counter]        = mxCreateNumericArray(4, DimZXLAS, mxDOUBLE_CLASS, mxREAL);
        pBuffer              = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pBeamI  = pBuffer;

        plhs[Counter]        = mxCreateNumericArray(4, DimZXLAS, mxDOUBLE_CLASS, mxREAL);
        pBuffer              = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pBeamQ  = pBuffer;

        // Output the angle at scatters
        plhs[Counter] = mxCreateNumericArray(3, DimSZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer       = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pAngleAtPixT = pBuffer;

        plhs[Counter] = mxCreateNumericArray(3, DimLAZX, mxDOUBLE_CLASS, mxREAL);
        pBuffer       = mxGetPr(plhs[Counter]); Counter++;
        ImageMarrow->pListenedAngleAtPixR = pBuffer;

        initializeImage(Data, ImageMarrow);

	}
    if (Debug) cout << "Pointing to outputs done" << endl;
    /* #######################################################################################
     * ################################# End of Setup ########################################
     * #######################################################################################

     * #######################################################################################
     * ############################# Beginning Calculation ###################################
     * #######################################################################################
    * */

    // Resolution in both axes
    double 	ResolutionZ = Data->pZ[2]-Data->pZ[1];
    double 	ResolutionX = Data->pX[2]-Data->pX[1];
    // Min and Max depths to Image
    Data->ZMin = 0;
    /* // modified by GR, Sept 11 2020
    // If we dont have the parabola we calculate down to a priori bone depth.
    if (TestPeri == 0) Data->ZMax = __min(round((Setup.DProbePerios - Data->pZ[0])/ResolutionZ), Data->NZ);
    // Else we cacluate to the deepest point of the Periosteum.
    else Data->ZMax = PeriParab.IdxMax; */
    Data->ZMax = __min(round((Setup.DProbePerios - Data->pZ[0])/ResolutionZ), Data->NZ); // modified by GR, Sept 11 2020

    // Time Arrival Calculation within Tissue.
    switch (Setup.NeedTravelTime)
    {
        case 0 :
            getTimeArrivalTissue(Setup, Data, ImageTissue);
            break;

        case 1 :
            getTimeArrivalTissue(Setup, Data, ImageTissue);
            // Save pTime, pAngle and Parabola structure.
            saveTimeAndAngle(Data, ImageTissue,
                            "TimeAngleTissueReception.bin",
                            "TimeAngleTissueTransmission.bin");
            saveAngleAtPix(Data, ImageTissue,
                            "AngleScattererTissueReception.bin",
                            "AngleScattererTissueTransmission.bin");
            saveParabola(&PeriParab, "PeriParab.bin");
            break;
        case 2 :
            loadTimeAndAngle(Data, ImageTissue,
                            "TimeAngleTissueReception.bin",
                            "TimeAngleTissueTransmission.bin");
            loadAngleAtPix(Data, ImageTissue,
                            "AngleScattererTissueReception.bin",
                            "AngleScattererTissueTransmission.bin");
            if (Setup.ReconTo!=1)
            {
                loadParabola(&PeriParab, "PeriParab.bin");
                // Data->ZMax = PeriParab.IdxMax; // modified by GR, Sept 11 2020
                if ( (Setup.DProbePerios-Data->pZ[0])/ResolutionZ < PeriParab.IdxMax ) mexErrMsgIdAndTxt(
                                "Main:Input","The deepest point of the Periosteum is shallower than a priori info DProbePerios.\n"
                                "Error of indexing in sight.\n"
                                "Rerun the code with NeedTravelTime=1 and a larger value for DProbePerios.");
            }
            break;
    }
    if (Debug) cout << "Arrival Time done" << endl;
    // Build Image within Tissue.
    buildImage(&Setup, Data, ImageTissue);
    // Build Image for VFI within Tissue.
    buildImageForVFI(&Setup, Data, ImageTissue);
    if (Debug) cout << "Reconstruction done" << endl;
    // Filling mxArray for Image displaying.
	for (int i = 0; i < Data->NX; i++)
	{
		for (int j = 0; j < Data->NZ; j++)
		{
            // j/i and i/j are inverted from the Image_st and the pointer pointing to the Matlab Image matrix. (Convention)
            ImTissue[j][i]  = sqrt( sqr( ImageTissue->pImI[i][j]) + sqr(ImageTissue->pImQ[i][j]) );
            *outTissue      = ImTissue[j][i];
			outTissue++;
		}
	}
    // Clean AngleAtPixR before doing VFI.
    if (Debug) cout << "Imaging Tissue done" << endl;
    // Tissue part finished, Now the Bone.
    if (Setup.ReconTo >1)
    {
        // Segmentation of image Im to find the Perios (Dijkstra)
        if (Setup.NeedTravelTime !=2 && TestPeri == 0) // modified by GR on September 30, 2020
        {
            segmentation(&PeriParab, ImTissue, Data->NX, Data->NZ, Data->pX, Data->pZ);
            // Beginning of the bone image in depth : higher point on the periosteum
            Data->ZMin = __max(PeriParab.IdxMin,0);
            if (PeriParab.IdxMin >= 0) Data->ZMax = __min(
                            PeriParab.IdxMax + Setup.CorticalThickMax / ResolutionZ, Data->NZ ); // modified by GR on September 30, 2020
            else Data->ZMax = Data->NZ;

            // Refering Parabola in Setup.
            Setup.pPeri = &PeriParab;
        }
        // After segmenting clean AngleAtPixR before doing VFI.
        cleanAngleAtPixel(&Setup, Data, ImageTissue, 0, Setup.ReconTo);
        // Time Arrival dow to the Bone.
        switch (Setup.NeedTravelTime)
        {
            case 0 :
                // If we provide Parabola we define ZMin, ZMax.
                if (TestPeri==1)
                {
                        Data->ZMin = PeriParab.IdxMin;
                        Data->ZMax = PeriParab.IdxMax + Setup.CorticalThickMax / ResolutionZ; // modified by GR on September 30, 2020
                }
                // Calculation of arrival times and angles.
                getTimeArrivalBone(Setup, Data, ImageBone);
                break;

            case 1 :
                // If we provide Parabola we define ZMin, ZMax.
                if (TestPeri==1)
                {
                        Data->ZMin = PeriParab.IdxMin;
                        Data->ZMax = PeriParab.IdxMax + Setup.CorticalThickMax / ResolutionZ; // modified by GR on September 30, 2020
                }
                // Calculation of arrival times and angles.
                getTimeArrivalBone(Setup, Data, ImageBone);
                // Save pTime, pAngle and Parabola structure.
                saveTimeAndAngle(Data, ImageBone,
                        "TimeAngleBoneReception.bin",
                        "TimeAngleBoneTransmission.bin");
                saveAngleAtPix(Data, ImageBone,
                        "AngleScattererBoneReception.bin",
                        "AngleScattererBoneTransmission.bin");
                saveParabola(&PeriParab, "PeriParab.bin");
                break;
            case 2 :
                // When loading we must get ZMin and Zmax from previously (case 1) segmentated images.
                loadParabola(&PeriParab, "PeriParab.bin");
                Data->ZMin = __max(PeriParab.IdxMin,0);
                if (PeriParab.IdxMin >= 0) Data->ZMax = __min(
                                PeriParab.IdxMax + Setup.CorticalThickMax / ResolutionZ, Data->NZ ); // modified by GR on September 30, 2020
                else Data->ZMax = Data->NZ;
                // Loading pTime and pAngles.
                loadTimeAndAngle(Data, ImageBone,
                        "TimeAngleBoneReception.bin",
                        "TimeAngleBoneTransmission.bin");
                loadAngleAtPix(Data, ImageBone,
                        "AngleScattererBoneReception.bin",
                        "AngleScattererBoneTransmission.bin");
                break;
                }
        if (Debug) cout << "Arrival Time Bone done" << endl;
        //debugTransmissionBone(Setup, Data);
        // Output the Perios.
		plhs[Counter] = fillOutParabola(PeriParab); Counter++;

        // Build of the Image within the Bone between the Perios and ZMax.
        buildImage(&Setup, Data, ImageBone);
        // Build image for VFI
        buildImageForVFI(&Setup, Data, ImageBone);
        if (Debug) cout << "Reconstruction Bone done" << endl;

        // Filling mxArray for Image displaying.
        for (int i = 0; i < Data->NX; i++)
        {
            for (int j = 0; j < Data->NZ; j++)
            {
                ImBone[j][i]  = sqrt( sqr( ImageBone->pImI[i][j]) + sqr(ImageBone->pImQ[i][j]) );
                *outBone      = ImBone[j][i];
                outBone++;
            }
        }


        // modified by GR, October 15 2020:        
        if (TestEndo == 0)
        {
            for (int i = 0; i < Data->NX; i++)
            {
                double XDist   = ResolutionX*i + Data->pX[0];                                   // lateral distance in mm (x-axis)
                double ZPerios = (PeriParab.a*sqr(XDist) + PeriParab.b*XDist + PeriParab.c); 	// depth of the periosteum at lateral distance iMm
                for (int j = 0; j < Data->NZ; j++)
                {
                    double ZDist = ResolutionZ*j + Data->pZ[0]; // depth in mm
                    if (ZDist < ZPerios + Setup.CorticalThickMin) ImBone[j][i] = 0;
                    // We ensure that the algorithm finds the endosteum in the good region of space.
                }
            }
            segmentation(&EndoParab, ImBone, Data->NX, Data->NZ, Data->pX, Data->pZ);
        }
        // Endosteum  outputs
        plhs[Counter]  = fillOutParabola(EndoParab); // modified by GR, October 15 2020



        if (Debug) cout << "Imaging Bone done" << endl;
        // Bone part finished, Now the Marrow.
        if (Setup.ReconTo > 2)
        {
            // In all case we calculate down to NZ.
            Data->ZMax   = Data->NZ;
            // If no Endeos provided, we set to zero the the bone Image shallower than the perios to find the Endeos by Segmentation.
            if (Setup.NeedTravelTime !=2 && TestEndo == 0) // modified by GR on September 30, 2020
            {
                for (int i = 0; i < Data->NX; i++)
                {
                    double XDist   = ResolutionX*i + Data->pX[0];                                   // lateral distance in mm (x-axis)
                    double ZPerios = (PeriParab.a*sqr(XDist) + PeriParab.b*XDist + PeriParab.c); 	// depth of the periosteum at lateral distance iMm
                    for (int j = 0; j < Data->NZ; j++)
                    {
                        double ZDist = ResolutionZ*j + Data->pZ[0]; // depth in mm
                        if (ZDist < ZPerios + Setup.CorticalThickMin) ImBone[j][i] = 0;
                        // We ensure that the algorithm finds the endosteum in the good region of space.
                    }
                }
                segmentation(&EndoParab, ImBone, Data->NX, Data->NZ, Data->pX, Data->pZ);
                // Min/Max of the Image in the Marrow.
                Data->ZMin   = __max(EndoParab.IdxMin,1);
                // Refering Parabola in Setup.
                Setup.pEndo = &EndoParab;
            }
            // Clean AngleAtPixR before doing VFI.
            cleanAngleAtPixel(&Setup, Data, ImageBone, 1, Setup.ReconTo);

            // Time Arrival calculation down to the Marrow.
            switch (Setup.NeedTravelTime)
            {
                case 0 :
                    if (TestEndo==1) Data->ZMin = EndoParab.IdxMin;
                    getTimeArrivalMarrow(Setup, Data, ImageMarrow);
                    break;

                case 1 :
                    if (TestEndo==1) Data->ZMin = EndoParab.IdxMin;
                    getTimeArrivalMarrow(Setup, Data, ImageMarrow);
                    saveTimeAndAngle(Data, ImageMarrow,
                            "TimeAngleMarrowReception.bin",
                            "TimeAngleMarrowTransmission.bin");
                    saveAngleAtPix(Data, ImageMarrow,
                            "AngleScattererMarrowReception.bin",
                            "AngleScattererMarrowTransmission.bin");
                    saveParabola(&EndoParab, "EndoParab.bin");
                    break;
                case 2 :
                    loadParabola(&EndoParab, "EndoParab.bin");
                    Data->ZMin   = EndoParab.IdxMin;
                    loadTimeAndAngle(Data, ImageMarrow,
                            "TimeAngleMarrowReception.bin",
                            "TimeAngleMarrowTransmission.bin");
                    loadAngleAtPix(Data, ImageMarrow,
                            "AngleScattererMarrowReception.bin",
                            "AngleScattererMarrowTransmission.bin");
                    break;
            }
            if (Debug) cout << "Arrival Time Marrow done" << endl;
            // Endosteum  outputs
            plhs[Counter]  = fillOutParabola(EndoParab);

            // Building of the marrow Image between the Endeos and ZMax.
            buildImage(&Setup, Data, ImageMarrow);
            // Build image for VFI
            buildImageForVFI(&Setup, Data, ImageMarrow);
            if (Debug) cout << "Reconstruction Marrow done" << endl;

            // Filling mxArray for Image displaying.
            for (int i = 0; i < Data->NX; i++)
            {
                for (int j = 0; j < Data->NZ; j++)
                {
                    ImMarrow[j][i]  = sqrt( sqr( ImageMarrow->pImQ[i][j]) + sqr(ImageMarrow->pImI[i][j]) );
                    *outMarrow      = ImMarrow[j][i];
                    outMarrow++;
                }
            }
            // Clean AngleAtPixR before doing VFI.
            cleanAngleAtPixel(&Setup, Data, ImageMarrow, 2, Setup.ReconTo);
            if (Debug) cout << "Imaging Marrow done" << endl;
            // Free Image of the Marrow
            free2D(ImMarrow, Data->NZ);
            freeImage(ImageMarrow, Data->NX);

        } // End of Marrow reconstruction.

        // Free memory for Bone.
        free2D(ImBone, Data->NZ);
        freeImage(ImageBone, Data->NX);

    } // End of Bone Imaging.

    // Free Memory for Tissue and all arguments
    free2D(ImTissue, Data->NZ);
    freeImage(ImageTissue,Data->NX);
    delete [] Data;
    delete [] ImageTissue;
    delete [] ImageBone;
    delete [] ImageMarrow;
    delete [] SignalI;
    delete [] SignalQ;
    if (Debug) cout << "Deallocation done" << endl;
} // End
