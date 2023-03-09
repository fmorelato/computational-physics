/**************************************************
 *
 * File correlation.c
 *
 **************************************************/

#define MAIN_PROGRAM

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"global.h"
#include"random.h"
#include"Markov.h"

#define Ns 300000		/* number of sweeps after thermalization to compute gamma(t_M,t)  */

void inttochar(int n, char *a)
{
        int i,j;
        char v[4];
        for(i=0; i<4; i++)
	{
                a[i] = '0';
        }
        sprintf(v,"%d",n);	/* converts the integer n into a string of chars v */
        for(i=1; i<4; i++)
	{
                if(v[i]=='\0')
		{
                        break;
                }
        }
        for(j=4-i; j<4; j++)
	{
                a[j]=v[j-4+i];
        }
        return;
}


int main()
{
	int i, j, k, l, acc=0;
	int max_tm = 50, Mt=60;
	FILE *gamma_doc;
	double c[N];		/* c[t] = correlation between a generic x_k and x_(k+t) */
	double **correlat;	/* correlat[N][Ns] */
	double **gamma;		/* gamma[max_tm][Nbin-max_tm] */
	char filename[24];
	strcpy(filename,"data/gamma_data0000.txt");

	correlat = (double**) malloc(N*sizeof(double*));
	gamma    = (double**) malloc(max_tm*sizeof(double*));
	for(i=0; i<N; i++)
	{
		correlat[i] = (double*) malloc(Ns*sizeof(double));
	}
	for(i=0; i<max_tm; i++){
		gamma[i] = (double*) malloc((Ns-max_tm)*sizeof(double));
	}

	rlxd_init(1,251099);
	ranlxd(xx, N);
	for(i=0; i<N_therm; i++)	/* thermalization */
	{
		sweep(&acc);
	}

		/*++++++++ write corr in a matrix at each Markov time  ++++++++*/
	printf("computing corr...\n");
	for(j=0; j<Ns; j++)
	{
		acc=0;
		sweep(&acc);
		corr(c);
		for(i=0; i<N; i++)
		{
			correlat[i][j] = c[i];
		}
	}


		/*++++++++ compute gamma(t_M,t)  ++++++++*/
	printf("\ncomputing gamma...\n");
	for(i=0; i<N; i++)	/* i=t */
	{
		double ym=0;		/* mean value of corr[i,*] */
		double gamma0=0;
		char number[5];
		inttochar(i,number);
		for(j=0; j<4; j++)
		{
			filename[15+j]=number[j];
		}
		gamma_doc = fopen(filename, "wt");

			/* computing ym */
		for(j=0; j<Ns; j++)
		{
			ym+=correlat[i][j];
		}
		ym/=Ns;

			/* filling of gamma[j][k] */
		for(j=0; j<max_tm; j++)
		{
			for(k=0; k< Ns-max_tm; k++)
			{
				gamma[j][k] = (correlat[i][k]-ym)*(correlat[i][k+j]-ym);
			}
		}

			/* at a fixed j compute gamma(j) and var_gamma(j) */
		/* gamma(0) */
		for(j=0; j<Ns-max_tm; j++)
		{
			gamma0 += gamma[0][j];
		}
		gamma0/=(Ns-max_tm);
		fprintf( gamma_doc, "0 1.0 0.0\n");
		/* gamma(j) */
		for(j=1; j<max_tm; j++)
		{
			double gamma_m = 0, var_gamma = 0;
			for(k=0; k<Ns-max_tm; k++)
			{
				gamma_m += gamma[j][k];
				var_gamma += gamma[j][k]*gamma[j][k];
			}
			gamma_m /= (Ns-max_tm);
			var_gamma = var_gamma/Ns - gamma_m*gamma_m*(1+2*Mt+2*j);


			for(l=1; l<Mt+j; l++)
			{
				double corr_gamma = 0;
				for(k=0; k<Ns-max_tm-l; k++)
				{
					corr_gamma += gamma[j][k]*gamma[j][k+l];
				}
				var_gamma += ( 2*corr_gamma/(Ns-max_tm-l) );
			}
			fprintf( gamma_doc, "%d %f %f\n", j, gamma_m/gamma0, sqrt(var_gamma/(Ns-max_tm))/gamma0 );
		}
		
		free(correlat[i]);
		fclose( gamma_doc );
	}
	free(correlat);

	exit(EXIT_SUCCESS);	
}
