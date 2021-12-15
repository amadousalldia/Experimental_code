#ifndef DEF_VECTEUR
#define DEF_VECTEUR

#include "utilities.hpp"
#include <iostream>
#include <cmath>
// class to manipluate vectors
class vecteur
{
public:
        vecteur (double lx = 0., double lz = 0.);
        ~vecteur();
        vecteur operator+(const vecteur &second);
        vecteur operator-(const vecteur &second);
        friend vecteur operator*(double l, const vecteur & second);
        vecteur operator/(double l);

        double x;
        double z;
};

double sign     ( double x );
double dot      ( const vecteur & first, const vecteur &second );
double det      ( const vecteur & first, const vecteur &second );
double module   ( const vecteur &v );

#endif
