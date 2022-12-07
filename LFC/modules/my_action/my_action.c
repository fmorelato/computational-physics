
/****************************************************************************
 *
 * File my_action.c
 *
 * Harmonic oscillator action
 *
 * The externally accessible functions are
 * 	double my_action()
 *  		Computes the action for a particle in a one-dimensional
 * 		harmonic oscillator which is located in xx[i] at time a*i
 *
 * 	double dS(int j, double new)
 * 		gives the difference of the action of the original vector
 * 		xx[N] and a new vector where the j-th entrance has been
 * 		changed into "newval"
 *
 * Version 1.0
 * Author: Francesco Morelato
 *
 ****************************************************************************/

#define MY_ACTION_C

#include<stdio.h>
#include<stdlib.h>
#include"global.h"
#include"my_action.h"

/*
double my_action(void){
	double S=0;
	int i;
	for(i=0; i<N; i++){
		S +=  (M/2.0) *( xx[i]*xx[i]*(2+W*W) - 2*xx[i]*xx[(i+1)%N] );
	}
	return S;
}
*/

double my_action(void){
	double S=0;
	int i;
	for(i=0; i<N; i++){
		S +=  (M/2.0) *( (xx[(i+1)%N] - xx[i])*(xx[(i+1)%N] - xx[i]) + W*W*xx[i]*xx[i] );
	}
	return S;
}

double dS(int j, double delta){
	return M*( (2+W*W)*delta*(2*xx[j]+delta)  - 2*(xx[(j+1)%N]+xx[(j+N-1)%N])*delta ) /2.;	/* new - old */
}
