
/****************************************************************************
 *
 * File Markov.c
 *
 * building of the Markov chain thanks to Metropolis method
 *
 * The externally accessible functions are
 * 	void sweep()
 * 		for a given physical system and a given path (hence a given
 * 		action), it gives the new path, as a Markov chain with the
 * 		Metropolis method
 *
 *
 * Version 1.0
 * Author: Francesco Morelato
 *
 ****************************************************************************/

#define MARKOV_C

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"random.h"
#include"global.h"
#include"Markov.h"
#include"my_action.h"

void sweep(int *acc){
	int i=0;
	double randv[2*N];

	ranlxd(randv, 2*N);
	for(i=0; i<N; i++){
		if(  randv[2*i+1] < exp( -dS(i , 2*Delta*(randv[2*i]-.5)) )  ){
			xx[i] += 2*Delta*(randv[2*i]-.5);
			*acc += 1;
		}
	}
}

	/*++++++++ Ineffective metropolis ++++++++*/
	/*
	 * use only to test
	 *
	 ******************************************/
/*
void sweep(int *acc){
	double randv[2];

	ranlxd(randv, 2);
	if(  randv[1] < exp( -dS(7 , 2*Delta*(randv[0]-.5)) )  ){
		xx[7] += 2*Delta*(randv[0]-.5);
		*acc += 1;
	}
}
*/
