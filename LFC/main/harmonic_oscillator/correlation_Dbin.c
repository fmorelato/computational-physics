/**************************************************
 *
 * File correlation.c
 *
 **************************************************/

#define MAIN_PROGRAM

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"global.h"
#include"random.h"
#include"Markov.h"

#define T_M   100000	/* number of sweeps after thermalization to compute gamma(t_M,t)  */
#define Dbin  30


void inttochar(int n, char *a){
        int i,j;
        char v[4];
        for(i=0; i<4; i++){
                a[i] = '0';
        }
        sprintf(v,"%d",n);
        for(i=1; i<4; i++){
                if(v[i]=='\0'){
                        break;
                }
        }
        for(j=4-i; j<4; j++){
                a[j]=v[j-4+i];
        }
        return;
}


int main(){
	int Nbin = T_M/Dbin;
	int i, j, k, acc=0;
	int max_tm = 10;
	FILE *gamma_doc;
	double c[N];		/* c[t] = correlation between a generic x_k and x_(k+t) */
	double **correlat;	/* correlat[N][Nbin] */
	char filename[29];
	strcpy(filename,"data_Dbin/gamma_data0000.txt");

	correlat = (double**) malloc(N*sizeof(double*));
	for(i=0; i<N; i++){
		correlat[i] = (double*) malloc(Nbin*sizeof(double));
		for(j=0; j<Nbin; j++){
			correlat[i][j] = 0;
		}
	}

	rlxd_init(1,251099);
	ranlxd(xx, N);
	for(i=0; i<N_therm; i++){	/* thermalization */
		sweep(&acc);
	}

		/*+++++ write the averaged corr within a bin in a matrix at each binned time  +++++*/
	for(j=0; j<Nbin; j++){
		acc=0;
		for(k=0; k<Dbin; k++){
			sweep(&acc);
			corr(c);
			for(i=0; i<N; i++){
				if(c[i]<0){
					c[i]=0;
				}
				correlat[i][j] += c[i];
			}
		}
		for(i=0; i<N; i++){
			correlat[i][j]/=(double)Dbin;
		}
	}


		/*++++++++ compute gamma(t_M,t)  ++++++++*/
	for(i=0; i<N; i++){	/* i=t */
		double prod, meansq, gamma0;
		char number[5];
		inttochar(i,number);
		for(j=0; j<4; j++){
			filename[20+j]=number[j];
		}
		gamma_doc = fopen(filename, "wt");

			/* gamma(0,t) */
		prod=.0;
		meansq=.0;
		gamma0=0;
		for(j=0; j<Nbin; j++){
			prod += correlat[i][j]*correlat[i][j];
			meansq += correlat[i][j];
		}
		prod/=Nbin;
		meansq/=Nbin;
		meansq*=meansq;
		gamma0 = prod - meansq;
		fprintf( gamma_doc, "%d %f\n", 0, 1. );
		
			/* gamma(j,t)/gamma(0,t) j>0 */
		for(j=1; j<max_tm; j++){	/* j=t_M */
			double prod=.0;
			double meansq=.0;
			for(k=0; k<Nbin-j; k++){
				prod += correlat[i][k]*correlat[i][k+j];
				meansq += correlat[i][k];
			}
			for(k=Nbin-j; k<Nbin; k++){
				meansq += correlat[i][k];
			}
			prod/=(Nbin-j);
			meansq/=Nbin;
			meansq*=meansq;
			fprintf( gamma_doc, "%d %f\n", j, (prod - meansq)/gamma0 );
			/*if((prod-meansq)>0){
				fprintf( gamma_doc, "%d %f\n", j, (prod - meansq)/gamma0 );
			} else {
				fprintf( gamma_doc, "%d %f\n", j, .0 );
			}*/
		}
		fclose( gamma_doc );
	}

	exit(EXIT_SUCCESS);	
}
