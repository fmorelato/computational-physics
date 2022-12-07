/***************************************************************************
 *
 * File check_my_action.c
 *
 * This program checks that the module my_action.c works properly, i.e. 
 * it gives the right value of the action for a one-dimensional harmonic
 * oscillator
 *
 * Author: Francesco Morelato
 *
 **************************************************************************/

#define MAIN_PROGRAM
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"random.h"
#include"global.h"
#include"my_action.h"

int main(){
	int i;
	double c;
	double exact_action, dscr;
        double d_action;	/* delta action */
	double originalxx[N], randvals[N];

	int T1, T2;	/* boolean, 1 if any module works, 0 otherwise*/
	xx[0]=1.;	/* fill xx[N] with some easy values*/

	if(N<128){
		c = 1.03;
	} else {
		c = 1.003;
	}

	for(i=1; i<N; i++){	
		xx[i] = c*xx[i-1];
	}


		/*+++++ test my_action() +++++*/

	exact_action = (M/2) * ( (2.+W*W-2*c) * (xx[N-1]*xx[N-1]*xx[2] - 1 )/(c*c-1) + 2.*xx[N-1]*(xx[N-1]*c -1) );	/* exact value of the action computed by hand*/
	dscr = (exact_action - my_action())/exact_action;
	if( (dscr>-1.e-10)&&(dscr<1.e-10) ){
		T1=1;
		printf("\nTest passed => my_action() works correctly on this machine\n\n");
	} else {
		T1=0;
		printf("\nTest failed: my_action() gives incorrect results => do not use my_action on this machine\n\n");
	}
	
	/*
	printf("%f\n", my_action() );
	printf("%f\n", exact_action );
	*/

		/*+++++ test dS() +++++*/

	if(T1==1){	/* this check is executed only if my_action() works properly*/
		T2 = 1;
		rlxd_init(1,839460);		/*initialize with my badge number*/
		ranlxd(randvals, N);
		for(i=0; i<N; i++){
			originalxx[i] = xx[i];		/* save the old vector */
			randvals[i] = 2*Delta*(randvals[i] - .5);	/* N random numbers */ 
		}
		for(i=0; i<N; i++){
			d_action = dS(i, randvals[i]);
			xx[i] = xx[i]+randvals[i];
			dscr = (my_action()-exact_action) - d_action;
			if( dscr>.1e-8 || dscr<-1.e-8 ){
				T2=0;
				printf("\nTest failed: dS(int,double) gives incorrect results at i=%d => do not use my_action on this machine\n\n", i);
				printf("\nmy_action()-exact_action = %f\nd_action = %f\n", my_action()-exact_action, d_action);
				i = N;	/* if at least one control fails, everything is considered wrong and exit */
			}
			xx[i] = originalxx[i];	/* reset xx[i] */
		}
		if(T2==1){
			printf("\nTest passed => dS(int,double) works correctly on this machine\n\n");
		}
	}

	exit(EXIT_SUCCESS);
}
