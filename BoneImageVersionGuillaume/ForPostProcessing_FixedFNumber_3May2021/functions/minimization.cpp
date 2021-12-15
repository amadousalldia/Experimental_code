#include "minimization.hpp"
using namespace std;
double
brent (const double ax,
       const double bx,
       const double cx,
       double       (*f)(double, void *),
       void         *localBuffer,
       const double tol,
       double       *xmin,
       const int    ITMAX)

{
    int iter;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double tmp1, tmp2, e = 0.0;

    a = ax < cx ? ax : cx;
    b = ax > cx ? ax : cx;
    x = bx;
	w = bx;
	v = bx;
	fx = (*f)(x, localBuffer);
    fw = fx;
	fv = fx;
    if (fw == _RECIPES_INFINITE) {
    return fw;
	}

	d = 0.0;
        for (iter = 1; iter <= ITMAX; iter++) {
                xm = 0.5 * (a + b);
		tol1 = tol * fabs(x) + ZEPS;
                tol2 = 2.0 * tol1;

                if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                        *xmin = x;
                        return fx;
                }
                if (fabs(e) > tol1) {
                        tmp1 = x - w;
			tmp2 = x - v;
                        r = tmp1 * (fx - fv);
                        q = tmp2 * (fx - fw);
                        p = tmp2 * q - tmp1 * r;
                        q = 2.0 * (q - r);
                        if (q > 0.0) p = -p;
                        q = fabs(q);
                        etemp = e;
                        e = d;
                        if (TEST) {
				e = x >= xm ? a - x : b - x;
				d = CGOLD * e;
                        }
			else {
                                d = p / q;
                                u = x + d;
                                if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
                        }
                }
                else {
			e = x >= xm ? a - x : b - x;
			d = CGOLD * e;
                }
		u = fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d);
                fu = (*f)(u, localBuffer);
                if (fu == _RECIPES_INFINITE) {
			return fu;
		}
                if (fu <= fx) {
                        if (u >= x) a = x; 
			else        b = x;
                        SHFT(v, w, x, u)
                        SHFT(fv, fw, fx, fu)
                }
                else {
                        if (u < x) a = u; 
			else       b = u;
                        if (fu <= fw || w == x) {
                                v = w;
                                w = u;
                                fv = fw;
                                fw = fu;
                        }
                        else if (fu <= fv || v == x || v == w) {
                                v = u;
                                fv = fu;
                        }
                }
        }
        *xmin = x;
        return fx;
}

double
falsePosition (double (*f)(const double&, const vecteur&, const vecteur&,
                        const double&, const double&, const double&, const double&),
		ExpSetup_st *pExp,
            	double X, double FOfX, double X1, double FOfX1,
                const double& rAccuracy, vecteur Ui, vecteur N, const double& rSpeedIso)

{
	double X2, FOfX2;
    	do
	{
    		X2 	= (FOfX*X1 - FOfX1*X) / (FOfX - FOfX1);
    		FOfX2	= f(X2, Ui, N, pExp->AnisoCoef, pExp->AnisoShape, pExp->CBoneAxial, rSpeedIso);
        	if (FOfX * FOfX2 <= 0.)
		{
            		X1    = X2;
            		FOfX1 = FOfX2;
        	}
        	else
		{
            		X    = X2;
            		FOfX = FOfX2;
        	}
    	} while (fabs(FOfX2) > rAccuracy);
	return X2;
}

