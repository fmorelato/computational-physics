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

#define T_M 50000		/* number of sweeps after thermalization to compute gamma(t_M,t)  */

void inttochar(int n, char *a)
{
        int i,j;
        char v[4];
        for(i=0; i<4; i++)
	{
                a[i] = '0';
        }
        sprintf(v,"%d",n);
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
	int i, j, k, acc=0;
	int max_tm = 50;
	FILE *gamma_doc;
	double c[N];		/* c[t] = correlation between a generic x_k and x_(k+t) */
	double **correlat;	/* correlat[N][Nbin] */
	char filename[24];
	strcpy(filename,"data/gamma_data0000.txt");

	correlat = (double**) malloc(N*sizeof(double*));
	for(i=0; i<N; i++)
	{
		correlat[i] = (double*) malloc(T_M*sizeof(double));
	}

	rlxd_init(1,251099);
	ranlxd(xx, N);
	for(i=0; i<N_therm; i++)	/* thermalization */
	{
		sweep(&acc);
	}

		/*++++++++ write corr in a matrix at each Markov time  ++++++++*/
	printf("computing corr...\n");
	for(j=0; j<T_M; j++)
	{
		int ac=0;
		sweep(&ac);
		acc+=ac;
		corr(c);
		for(i=0; i<N; i++)
		{
			correlat[i][j] = c[i];
		}
		/*printf("\r %2.2f %%", 100*(double)j/(double)T_M );*/
	}
	printf("%f\n", (float)acc/(N*T_M));

	printf("\ncomputing gamma...\n");
		/*++++++++ compute gamma(t_M,t)  ++++++++*/
	for(i=0; i<N; i++)	/* i=t */
	{
		double prod, meansq, gamma0;
		char number[5];
		inttochar(i,number);
		for(j=0; j<4; j++)
		{
			filename[15+j]=number[j];
		}
		gamma_doc = fopen(filename, "wt");

			/* gamma(0,t) */
		prod=.0;
		meansq=.0;
		gamma0=0;
		for(j=0; j<T_M; j++)
		{
			prod += correlat[i][j]*correlat[i][j];
			meansq += correlat[i][j];
		}
		prod/=T_M;
		meansq/=T_M;
		meansq*=meansq;
		gamma0 = prod - meansq;
		fprintf( gamma_doc, "%d %f\n", 0, 1. );
		
			/* gamma(j,t)/gamma(0,t) j>0 */
		
		for(j=1; j<max_tm; j++)	/* j=t_M */
		{
			double prod=.0;
			double meansq=.0;
			for(k=0; k<T_M-j; k++)
			{
				prod += correlat[i][k]*correlat[i][k+j];
				meansq += correlat[i][k];
			}
			for(k=T_M-j; k<T_M; k++)
			{
				meansq += correlat[i][k];
			}
			prod/=(T_M-j);
			meansq/=T_M;
			meansq*=meansq;
			fprintf( gamma_doc, "%d %f\n", j, (prod - meansq)/gamma0 );
		}
		free(correlat[i]);
		fclose( gamma_doc );
	}
	free(correlat);

	exit(EXIT_SUCCESS);	
}
