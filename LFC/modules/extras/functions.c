/***********************************************************
 * File functions.c
 *
 * Useful functions to compute observables for the harmonic
 * oscillator
 *
 * The externally accessible functions are
 * 	double my_arccosh(double x)
 * 		computes the arccosh
 *
 * 	double func_DE(double c0, double c1, double c2)
 * 		Computes DE for a given time
 *
 * Version 1.0
 * Author: Francesco Morelato
 *
 ************************************************************/

#define FUNCTIONS_C

#include<math.h>
#include<stdlib.h>

double my_arccosh(double x)
{
        return log( x + sqrt(x*x -1.0) );
}

double func_DE(double c0, double c1, double c2)
{
        return my_arccosh( (c2+c0)/(2.*c1) );
}

double func_M(double c, double DE, double t)
{
	return c/( exp(-t*DE)+exp((t-64)*DE) );
}
