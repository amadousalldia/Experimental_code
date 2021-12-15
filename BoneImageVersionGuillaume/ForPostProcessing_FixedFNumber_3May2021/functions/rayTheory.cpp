#include "rayTheory.hpp"
using namespace std;

int
refractionIsotropic (const vecteur& Ui, const vecteur& N, const double& RatioC1C2, vecteur * Ur)
{
    double 	Cos1, Sin1, Sin2, Cos2, Sin2Sqr;
    int 	IsRefrac = 0; 	        // if the function returns ret=0, no refracted ray has been found
    vecteur T(N.z, -N.x); 	        // tangent direction
    Cos1 	= dot(Ui, N); 	        // normal incident component
    Sin1 	= dot(Ui, T); 	        // tangent incident component
    Sin2 	= Sin1 / RatioC1C2; 	// Snell-Descartes' law
    Sin2Sqr = sqr(Sin2);
    if (Sin2Sqr <= 1.0) // if the refracted ray exists
	{
        Cos2 	= sign(Cos1) * sqrt(1.0 - Sin2Sqr);
        *Ur 	= Sin2 * T + Cos2 * N;
        IsRefrac= 1;
    }
    return IsRefrac;
}

double
rayLength(const double& a, const double& b, const double& c, const vecteur& Pt, const vecteur& Ur)
{
	double Lamda1, Lamda2;

	double Alpha = a * sqr(Ur.x);
	double Beta  = 2 * a * Pt.x * Ur.x + b * Ur.x - Ur.z;
	double Gamma = a * sqr(Pt.x) + b * Pt.x + c - Pt.z;
    //if (Alpha == 0)
    if (fabs(Alpha) < 1e-12) // modified by GR, 18 November 2020
    { // parabola is a straight line OR vertical propagation
        //if (Beta == 0.) return 0.; // ray is parallel to the interface
        if (fabs(Beta) < 1e-12) return 0.;  // modified by GR, 18 November 2020

        Gamma  = a * sqr(Pt.x) + b * Pt.x + c - Pt.z;
        Lamda1 = - Gamma / Beta; // Lamda1 is the length of the ray from Pt to interface.
	}
	else
    { // real parabola
        double Delta = sqr(Beta) - 4. * Alpha * Gamma;
        if (Delta < 0.) return 0.; // no solution

        Delta = sqrt(Delta);
        Lamda1 = 2. * Gamma / (- Beta - Delta);
        Lamda2 = 2. * Gamma / (- Beta + Delta);

        if (Lamda1 <= 0) Lamda1 = Lamda2;				// non acceptable solution
        else if (Lamda2 < Lamda1 && Lamda2 > 0) Lamda1 = Lamda2;	// min value expected
	}
	return Lamda1 >= 0 ? Lamda1 : 0;
}

double
snellEq(const double& SinTheta, const vecteur& Ui, const vecteur& N,
                const double& AnisoCoef, const double& AnisoShape,
                const double& CBone, const double& CIso)
{
    vecteur T(N.z, -N.x);
    double SinSqr 	= sqr(SinTheta);
    double CosSqr 	= 1. - SinSqr;
    // calculation of the phase velocity in the anisotropic medium
    double VWeakAniso 	= CBone * (1 - AnisoCoef * (AnisoShape * CosSqr * SinSqr + sqr(CosSqr)));
    // residual of Snell's law at the tissue/bone interface : must be 0 to satisfy exactly Snell's law
    double IsRefrac 	= fabs(dot(Ui, T)) / CIso - fabs(SinTheta) / VWeakAniso;

    return IsRefrac;

}

int
refractionAnisotropic(const vecteur& Ui, const vecteur& N, vecteur *Ur,
                         double *pSinSqrGrAngle, ExpSetup_st *Setup)
{
	vecteur T(N.z, -N.x); 	// tangent vector
	int     IsRefrac = 0; 	// if = 0, no refracted ray has been found
	double  Cos1     = dot(Ui, N);
	double  Sin1     = dot(Ui, T);

	if (fabs(Sin1) <= Setup->RatioTB)
	{ // if a refracted angle exists, we look for the phase angle given by Snell's law
		double xleft 	= 0.;   // minimum sin_t for search of phase angle
		double xright 	= 0.99; // maximum sin_t for search of phase angle
		double fxleft 	= snellEq(xleft,  Ui, N, Setup->AnisoCoef, Setup->AnisoShape, Setup->CBoneAxial, Setup->CTissue);
		double fxright 	= snellEq(xright, Ui, N, Setup->AnisoCoef, Setup->AnisoShape, Setup->CBoneAxial, Setup->CTissue);

		if (fxleft*fxright < 0.)
		{ // the 0 of the fx function is calculated with the false position method
			// sin2 is the sin squared of the phase angle in bone that satisfies Snell's law
			double Sin2 = falsePosition(snellEq, Setup, xleft, fxleft, xright, fxright, 1.e-8, Ui, N, Setup->CTissue);
			// We now want to find the group angle corresponding to sin2.
			double PhaseAngle = asin(Sin2);
			double Delta 	  = -Setup->AnisoCoef * Setup->AnisoShape;
			double Epsilon    = -Setup->AnisoCoef;
			double GrAngle 	  = PhaseAngle - ( Delta + 2*(Epsilon-Delta) * sqr(cos(PhaseAngle)) ) * sin(2*PhaseAngle);
			double SinGrAngle = sin(GrAngle);
			*pSinSqrGrAngle   = sqr(SinGrAngle);
			double CosGrAngle = cos(GrAngle);
			*Ur 		  = sign(Sin1)*SinGrAngle * T + sign(Cos1)*CosGrAngle * N; // unit vector describing the ray in bone

            IsRefrac++;
		}
	}
	return IsRefrac;
}

double
travelTimeInTissue(double Theta, void *pBuffer)
{
	ExpSetup_st *pSetup     = (ExpSetup_st *)pBuffer;
    // first path, from sources to the end of the lens
    double tanTheta 	= tan(Theta);
    // point reached at the interface lens/tissue
    double x                = (pSetup->XStart + (pSetup->LensThick - pSetup->ZStart) * tanTheta);
    double z                = pSetup->LensThick;

    // time for the ray to reach the interface (travel in the lens)
    double TravelTime = distance(pSetup->XStart, pSetup->ZStart, x, z) / pSetup->CLens;

    // time from the end (x) of the lens to the image point in tissue
    TravelTime += distance(x, z, pSetup->XEnd, pSetup->ZEnd) / pSetup->CTissue;

    // Traveltime is the total time from the element to the point
	return TravelTime;
}

double
getAngleAtScatterWithinTissue(const double& Theta, const double& RatioClCt)
{
    vecteur Urefrac, Uincident( sin(Theta), cos(Theta) );
    vecteur N(0.,1.);
    double Tmp = refractionIsotropic(Uincident, N, RatioClCt, &Urefrac);
    return asin(Urefrac.x);
}

double
travelTimeThroughTissueBone(double Theta, void *pBuffer)
{
	ExpSetup_st *pSetup = (ExpSetup_st *)pBuffer;

    // first path, from sources to the end of the lens
    vecteur Ur, Ui(sin(Theta), cos(Theta));
    vecteur N(0., 1.);

    double  tanTheta = tan(Theta);
    // point reached at the interface lens/tissue
    double  x        = pSetup->XStart + (pSetup->LensThick - pSetup->ZStart) * tanTheta;
    double  z        = pSetup->LensThick;
    vecteur Pe(pSetup->XStart, pSetup->ZStart);
    vecteur Pl = vecteur(x,z);

    // time for the ray to reach the interface (travel in the lens)
    double TimeTravel   = distance(Pe.x, Pe.z, Pl.x, Pl.z) / pSetup->CLens;
    // calculation of the refracted vector Ur at the tissue/bone interface
    int Ret             = refractionIsotropic(Ui, N, pSetup->RatioLT, &Ur);

    // if the ray is reflected and not refracted
    if (Ret == 0) return 1. + fabs(Theta);

    // At the end of the tissue, the ray intersects the periosteum (parabola with parameters a,b,c)
    double Lamda1 = rayLength(pSetup->pPeri->a, pSetup->pPeri->b, pSetup->pPeri->c, Pl, Ur); // Same as BEfore checked on 04/09/18.
	if (Lamda1==0.) return 1. + fabs(Theta);

    TimeTravel	+= Lamda1 / pSetup->CTissue; // time traveled in tissue
	vecteur Pb1 	= Pl + Lamda1*Ur; // total vector from the probe element to the periosteum
    // path in the bone
    Ui 	= Ur; // incident vector in the bone
    N 	= vecteur(2. * pSetup->pPeri->a * Pb1.x + pSetup->pPeri->b, -1.); // normal vector to the periosteum
    N 	= N / module(N); // normalization of N

    vecteur R(pSetup->XEnd, pSetup->ZEnd); 				// vector correponding to the point to reach

    // Calculation of the angle necessary to reach the targetted point
    double CosSqrGrAngle = 1.; // cos^2 group angle in bone
    if (module(R-Pb1) > 0) CosSqrGrAngle = sqr(dot(N,(R-Pb1)/module(R-Pb1))); // Correction by J.N 2019.
    // sin^2 group angle in bone
	double SinSqrGrAngle = 1. - CosSqrGrAngle;
    // Calculation of the group velocity at which the ray travels in the bone
    double GrVelocity = pSetup->CBoneAxial *
        (1 - pSetup->AnisoCoef * (pSetup->AnisoShape * CosSqrGrAngle * SinSqrGrAngle + sqr(CosSqrGrAngle)));
    // Time travelled in the bone
    TimeTravel += module(R - Pb1) /GrVelocity;
    // t is the total time from the element to the point
	return TimeTravel;
}

double
getAngleAtScatterWithinBone(const double& Theta, void *pBuffer)
{
	ExpSetup_st *pSetup = (ExpSetup_st *)pBuffer;

    vecteur Ur, Ui(sin(Theta), cos(Theta));
    vecteur N(0., 1.);

    double tanTheta 	= tan(Theta);
    double x                = (pSetup->XStart + (pSetup->LensThick - pSetup->ZStart) * tanTheta);
    double z                = pSetup->LensThick;

    //vecteur Pe(pSetup->XStart, pSetup->ZStart); // modified by GR, 18 November 2020
    vecteur Pl = vecteur(x,z);

    double Tmp      = refractionIsotropic(Ui, N, pSetup->RatioLT, &Ur);
    double Lamda1 	= rayLength(pSetup->pPeri->a, pSetup->pPeri->b, pSetup->pPeri->c, Pl, Ur);
    vecteur Pb1 	= Pl + Lamda1 * Ur;

    vecteur R(pSetup->XEnd, pSetup->ZEnd);
    vecteur Diff = (R-Pb1)/module(R-Pb1);

    return asin(Diff.x);
}

double
travelTimeThroughTissueBoneMarrow( double Theta, void *pBuffer)
{
	ExpSetup_st *Setup = (ExpSetup_st *)pBuffer;

    // first path, from sources to the end of the lens
    vecteur Ur, Ui(sin(Theta), cos(Theta)), N(0., 1.), Pe(Setup->XStart, Setup->ZStart);

    double TanTheta 	= tan(Theta);
    double x                = (Setup->XStart + (Setup->LensThick-Setup->ZStart) * TanTheta);
    double z                = Setup->LensThick;
    vecteur Pl(x,z);

    // time for the ray to reach the interface (travel in the lens)
    double TravelTime = distance(Pe.x, Pe.z, Pl.x, Pl.z) / Setup->CLens;
    // calculation of the refracted angle at the tissue/lens interface
    int Ret = refractionIsotropic(Ui, N, Setup->RatioLT, &Ur);
    if (Ret == 0) return 1. + fabs(Theta);

 	// At the end of the tissue, the ray intersects the periosteum (parabola with parameters a,b,c)
	double Lamda1 	= rayLength(Setup->pPeri->a, Setup->pPeri->b, Setup->pPeri->c, Pl, Ur);
	if (Lamda1==0.) return 1. + fabs(Theta);

	vecteur Pb1 	= Pl + Lamda1*Ur;
    // time traveled in the tissue
    TravelTime	+= Lamda1 / Setup->CTissue;

    // path in the bone:
	// must apply Snell's law at interface between soft tissue and anisotropic bone to find what is the phase angle	vecteur
    Ui 	= Ur;
    N 	= vecteur(2. * Setup->pPeri->a * Pb1.x + Setup->pPeri->b, -1.);
    N 	= N / module(N);
    double SinSqrGrAngle, GrVelocity; // phase angle in bone

    // Calculation of the refracted vector Ur in the bone
    Ret = refractionAnisotropic(Ui, N, &Ur, &SinSqrGrAngle, Setup);
	if (Ret == 0) return 1. + fabs(Theta);

   	GrVelocity = Setup->CBoneAxial*
                (1 - Setup->AnisoCoef*(Setup->AnisoShape*(1-SinSqrGrAngle) * SinSqrGrAngle + sqr(1-SinSqrGrAngle)));

    /* length of the ray path in bone
    * At the end of the bone, the ray intersects the endosteum (parabola with parameters A,B,C)
    */
    double Lamda2 = rayLength(Setup->pEndo->a, Setup->pEndo->b, Setup->pEndo->c, Pb1, Ur);
	if (Lamda2==0.) return 1. + fabs(Theta);

    TravelTime    += Lamda2/ GrVelocity; // time spent in the bone
    // Vector from the element to the targetted point (end of bone)
    vecteur Pb2 = Pb1 + Lamda2*Ur;
    // path in the marrow to reach (xr,zr)
    vecteur R(Setup->XEnd, Setup->ZEnd);
    // time spent in the marrow
    TravelTime += module(R - Pb2) / Setup->CMarrow;
    // Total time from the element to the point
	return TravelTime;
}

double
getAngleAtScatterWithinMarrow (const double& Theta, void *pBuffer)
{
	ExpSetup_st *Setup = (ExpSetup_st *)pBuffer;
    // first path, from sources to the end of the lens
    vecteur Ur, Ui(sin(Theta), cos(Theta)), N(0., 1.);
    double TanTheta 	= tan(Theta);
    double x                = Setup->XStart + (Setup->LensThick-Setup->ZStart)*TanTheta;
    double z                = Setup->LensThick;
    vecteur Pe(Setup->XStart, Setup->ZStart), Pl(x,z);

    // Path in tissue
    int Ret         = refractionIsotropic(Ui, N, Setup->RatioLT, &Ur);
	double Lamda1 	= rayLength(Setup->pPeri->a, Setup->pPeri->b, Setup->pPeri->c, Pl, Ur);
	vecteur Pb1 	= Pl + Lamda1*Ur;
	// path in the bone
    Ui 	= Ur;
    N 	= vecteur(2. * Setup->pPeri->a * Pb1.x + Setup->pPeri->b, -1.);
    N 	= N / module(N);
    // Phase angle within Bone.
    double SinSqrGrAngle, GrVelocity;

    // Calculation of the refracted vector Ur in the bone
    Ret             = refractionAnisotropic(Ui, N, &Ur, &SinSqrGrAngle, Setup);
    // Weak Anisotropy Assumption
   	GrVelocity = Setup->CBoneAxial*
                (1 - Setup->AnisoCoef*(Setup->AnisoShape*(1-SinSqrGrAngle) * SinSqrGrAngle + sqr(1-SinSqrGrAngle)));

    // Reaching Marrow interface
    double Lamda2 = rayLength(Setup->pEndo->a, Setup->pEndo->b, Setup->pEndo->c, Pb1, Ur);
    vecteur Pb2     = Pb1 + Lamda2*Ur;
    
    /*N 	        = vecteur(2. *Setup->pEndo->a * Pb2.x + Setup->pEndo->b, -1.);
    N 	        = N / module(N);
    vecteur R(Setup->XEnd, Setup->ZEnd);
    return asin(det(N,R-Pb2)/module(R-Pb2));*/ // modified by GR, 18 November 2020
    
    vecteur R(Setup->XEnd, Setup->ZEnd); // added by GR, 18 November 2020
    vecteur Diff = (R-Pb2)/module(R-Pb2); // added by GR, 18 November 2020

    return asin(Diff.x);
    
}
