#define MAIN_PROGRAM

#include<stdio.h>
#include<stdlib.h>
#include"random.h"
#include"Markov.h"
#include"my_action.h"
#include"global.h"

#define iter 600

double min(double a, double b){
	if(a>b){
		return b;
	} else {
		return a;
	}
}

double max(double a, double b){
	if(a<b){
		return b;
	} else {
		return a;
	}
}

int main(){
	int i=0, acc=0, partial_acc=0;
	FILE *doc_array;
	FILE *doc_c, *doc_h;
	FILE *doc_acc;
	double v_acc[iter];
	doc_c =   fopen("thermalization_c_data.txt", "wt" );	/* cold thermalization */
	doc_h =   fopen("thermalization_h_data.txt", "wt" );	/* hot  thermalization */
	doc_acc = fopen("acceptance_data.txt"      , "wt" );	/*      acceptance     */

		/*+++++++ COLD THERMALIZATION +++++++*/
	rlxd_init(1,42069);
	for(i=0; i<N; i++){
		xx[i]=0.;
	}
	for(i=0; i<iter; i++){
		fprintf( doc_c, "%d %f\n", i, my_action());
		sweep(&acc);
	}
	fprintf( doc_c, "%d %f\n", i, my_action());
	printf("the acceptance in the cold thermalization is %2.2f %%\n", 100*(float)acc/(N*iter));
	fclose(doc_c);

		/*+++++++ HOT  THERMALIZATION +++++++*/
	acc = 0;
	ranlxd( xx, N);
	for(i=0; i<N; i++){
		xx[i]=(xx[i]-0.5)*3;
	}
	for(i=0; i<iter; i++){
		sweep(&partial_acc);
		acc += partial_acc;
		v_acc[i] = 100*(float)partial_acc/N;
		partial_acc = 0;
		fprintf( doc_h, "%d %f\n", i, my_action());
	}
	fprintf( doc_h, "%d %f\n", i, my_action());
	printf("the acceptance in the hot thermalization is %2.2f %%\n", 100*(float)acc/(N*iter));
		/* moving average */
	for(i=0; i<iter; i++){
		int j, count=0;
		double media_mob=.0;
		for(j=max(i-4,0); j<min(i+4,iter); j++){
			media_mob += v_acc[j];
			count += 1;
		}
		v_acc[i] = media_mob/count;
		fprintf(doc_acc,"%d %2.2f\n", i, v_acc[i]);
	}
	fclose(doc_h);
	fclose(doc_acc);


		/*++++++++ FINAL ARRAY ++++++++*/
	doc_array = fopen("array.txt", "wt");
	for(i=0; i<N; i++){
		fprintf( doc_array, "%d %f\n", i, xx[i]);
	}
	fclose(doc_array);
	exit(EXIT_SUCCESS);
}
