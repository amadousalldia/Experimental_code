#include "segmentation.hpp"
using namespace std;

int *
createContour(double **Im, const int Nx, const int Nz)

{
    double A, *CostFun1 = new double [2 * Nz], *CostFun2 = CostFun1 + Nz;
    int **CostPath = newTable2D(Nz, Nx);
    int x, z, n;
    int *ZIndexPath = new int [Nx];
    // from left to right
    memset(CostFun1, 0, 2 * Nz * sizeof(double));
    for (x = 0; x < Nx; x++)
    {
        for (z = 1; z < Nz - 1; z++)
        {
            maxInArray(CostFun1 + z - 1, 3, A, n);
            CostPath[z][x] = z + (n - 1);
            CostFun2[z]    = A + Im[z][x];
        }
        memcpy(CostFun1, CostFun2, Nz * sizeof(double));
    }

    maxInArray(CostFun2, Nz, A, n);
    ZIndexPath[Nx - 1] = n;
    for (x = Nx - 2; x >= 0; x--)
    {
        ZIndexPath[x] = CostPath[ZIndexPath[x + 1]][x+1];
    }
    // Desalocate.
    delete [] CostFun1;
    deleteTable2D(CostPath, Nz);
    return ZIndexPath;
}

void
segmentation (Parabola_st *Parabola, double ** Im, const int NX, const int NZ, double *X, double *Z)
{
    double resolutionZ = Z[2] - Z[1];
    double resolutionX = X[2] - X[1];
	int x, *Contour    = createContour(Im, NX, NZ);

    // Fit of Contour by a parabola based on Cramer's rule (https://en.wikipedia.org/wiki/Cramer's_rule)
    double M11 = 0, M12 = 0, M13 = 0, M21 = 0, M22 = 0, M23 = 0;
	double M31 = 0, M32 = 0, M33 = 0, Y1  = 0, Y2  = 0, Y3  = 0;
	double d, d2, dres;

    M11 = NX;
    for (x = 0; x < NX; x++)
    {
        d    = x;
        dres = d * resolutionX + X[0];
        d2   = dres * dres;
        M12  = M12 + dres;
        M13  = M13 + d2;
        M23  = M23 + dres * d2;
        M33  = M33 + d2 * d2;
        Y1   = Y1 + Contour[x] * resolutionZ + Z[0];
        Y2   = Y2 + dres * (Contour[x] * resolutionZ + Z[0]);
        Y3   = Y3 + d2 * (Contour[x] * resolutionZ + Z[0]);
    }
    M21 = M12;
    M22 = M13;
    M31 = M13;
    M32 = M23;

    double det_M = 0, det_M0 = 0, det_M1 = 0, det_M2 = 0;
    det_M  = M11 * (M22 * M33 - M23 * M32) - M21 * (M12 * M33 - M32 * M13) + M31 * (M12 * M23 - M22 * M13);
    det_M0 = Y1  * (M22 * M33 - M23 * M32) - Y2  * (M12 * M33 - M32 * M13) + Y3  * (M12 * M23 - M22 * M13);
    det_M1 = M11 * (Y2 * M33 - M23 * Y3)   - M21 * (Y1 * M33 - Y3 * M13)   + M31 * (Y1 * M23 - Y2 * M13);
    det_M2 = M11 * (M22 * Y3 - Y2 * M32)   - M21 * (M12 * Y3 - M32 * Y1)   + M31 * (M12 * Y2 - M22 * Y1);

    // now store parabola parameters
    Parabola->a = det_M2 / det_M;
    Parabola->b = det_M1 / det_M;
    Parabola->c = det_M0 / det_M;

    // Now we retrieve the Min and Max of the Parabola that
    // will define the beginning and end of the following image.
    int ParabMin = Contour[0], ParabMax = 0;
    for (x = 1; x < NX; x++)
    {
        if (Contour[x] < ParabMin) ParabMin = Contour[x];
        if (Contour[x] > ParabMax) ParabMax = Contour[x];
    }
    Parabola->IdxMin = ParabMin; Parabola->IdxMax = ParabMax;
    // Retrieve the raw (by pixel) segmentation before being defined as a Parabola.
    //memcpy(parabola_params.coords, Contour, NX * sizeof(int));

	delete [] Contour;
}
