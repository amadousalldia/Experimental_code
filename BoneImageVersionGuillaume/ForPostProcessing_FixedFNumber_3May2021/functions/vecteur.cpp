#include "vecteur.hpp"
// class to manipluate vectors with simple common mathematical laws
using namespace std;

// Constructeur ##############################
vecteur::vecteur (double lx, double lz)
{
        x = lx;
        z = lz;
}
//Destructeur ################################
vecteur::~vecteur()
{

}

// Operator #################################
vecteur
vecteur::operator+(const vecteur &second)
{
        return vecteur(x + second.x, z + second.z);
}

vecteur
vecteur::operator-(const vecteur & second)
{
        return vecteur(x - second.x, z - second.z);
}

// friendly one
vecteur
operator*(double l, const vecteur & second)
{
        return vecteur(l * second.x, l * second.z);
}

vecteur
vecteur::operator/(double l)
{
	return vecteur(x / l, z / l);
}

// Functions #####################################
// Function to find the sign of a number
double sign (double x)
{
        return x >= 0. ? 1 : -1;
}

// Function to make the scalar product between two vectors
double dot (const vecteur &first, const vecteur &second)
{
        return first.x * second.x + first.z * second.z;
}

double det (const vecteur &first, const vecteur &second)
{
        return first.x * second.z - second.x * first.z;
}
// Function to compute the module of a vector
double module (const vecteur & v)

{
        return sqrt( sqr(v.x) + sqr(v.z) );
}
