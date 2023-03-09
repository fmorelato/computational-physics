/**********************************************************
 *
 * File observable.c
 *
 * compute the observable DE = E1-E0, namely the first energy
 * level of the harmonic oscillator system with 'M' mass
 * and 'W' pulsation using the two points correlation 
 * function
 *
 **********************************************************/

#define MAIN_PROGRAM

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"global.h"
#include"random.h"
#include"Markov.h"
#include"functions.h"

#define Dbin 50
#define Nbin 1000000

#define bar_lenght 50


void barra_avanzamento(int avanzamento,int totale)
{
	int l;
	printf("\r[");
	for(l=0; l<(avanzamento*bar_lenght)/totale; l++)
		printf("#");
	for(l=l; l<bar_lenght; l++)
		printf(" ");
	printf("] %2.2f %%", 100.0*avanzamento/(double)totale );
	fflush(stdout);
}


int main()
{
	int i, j, k, acc=0, maxnum=0;
	double c[N], cm[N], err_cm[N];
					/* c[t]  = correlation between a generic x_k and x_(k+t) 
					 * cm[t] = averaged c over all sweeps
					 */
	double var_DE[N], var_M[N];
	double *DEm_Jack, *Mm_Jack;
	double **correlat;	/* correlat[N][Nbin] */
	FILE *doc_c, *doc_E, *doc_M;
	double DE_tot=0, err_DE_tot=0, M_tot=0, err_M_tot=0, DE_J, M_J;
	double denom_DE=0, denom_M=0;

	char filename1[10], filename2[11], filename3[10];
	strcpy(filename1, "c(t)3.txt");
	strcpy(filename2, "DE(t)3.txt");
	strcpy(filename3, "M(t)3.txt");

	DEm_Jack = (double*) malloc(Nbin*sizeof(double));
	Mm_Jack = (double*) malloc(Nbin*sizeof(double));
	correlat = (double**) malloc(N*sizeof(double*));
	for(i=0; i<N; i++)
	{
		correlat[i] = (double*) malloc(Nbin*sizeof(double));
		for(j=0; j<Nbin; j++)
		{
			correlat[i][j] = 0;
		}
	}

	rlxd_init(1,271096);
	ranlxd(xx, N);
	for(i=0; i<N_therm; i++)	/* thermalization */
		sweep(&acc);

		/*+++++ write the averaged corr within a bin in a matrix at each binned time  +++++*/
	for(i=0; i<Nbin; i++)
	{
		acc=0;
		for(j=0; j<Dbin; j++)
		{
			sweep(&acc);
			corr(c);
			for(k=0; k<N; k++)
			{
				correlat[k][i]+=(c[k]/(double)Dbin);
			}
			barra_avanzamento(i*Dbin+j,Nbin*Dbin);
		}
	}
	printf("\r[");
	for(i=0; i<bar_lenght; i++)
		printf("#");
	printf("] 100.00 %%\n");

		/*++++++ evaluation of cm[i] and err_cm[i] at each time i ++++++*/
	for(i=0; i<N; i++)
	{
		cm[i] = 0;
		err_cm[i]=0;
		for(j=0; j<Nbin; j++)
		{
			cm[i]     += correlat[i][j];
			err_cm[i] += correlat[i][j]*correlat[i][j];
		}
		cm[i]/=(double)Nbin;
		err_cm[i] = err_cm[i]/(double)Nbin - cm[i]*cm[i];
	}

			/* Jacknife's correlation cluster */
	for(i=0; i<N; i++)
	{
		for(j=0; j<Nbin; j++)
		{
			correlat[i][j] = ( (double)Nbin*cm[i] - correlat[i][j] ) / (double)(Nbin-1);
		}
	}

	for(i=0; i<N; i++)
	{
		double meanDE = func_DE(cm[(i+N-1)%N], cm[i], cm[(i+1)%N]); 
		double meanM  = func_M(cm[i],meanDE,i);
		if( (maxnum==0)&&( (cm[(i+N-1)%N]+cm[(i+1)%N])/2 < cm[i] )&&(i>1) )
			maxnum = i;
		var_DE[i]=0;
		var_M[i]=0;
		for(j=0; j<Nbin; j++)
		{
			DE_J = func_DE(correlat[(i+N-1)%N][j], correlat[i][j], correlat[(i+1)%N][j]);
			M_J = func_M(correlat[i][j],DE_J,i);
			var_DE[i] += ( DE_J - meanDE )*( DE_J - meanDE );
			var_M[i] += ( M_J - meanM )*( M_J - meanM );
		}
		var_DE[i] =  (double)(Nbin-1)*var_DE[i]/(double)Nbin;
		var_M[i] =  (double)(Nbin-1)*var_M[i]/(double)Nbin;
	}

	
		/* overall estimation of DE and sigma DE */
	for(i=2; i<maxnum-1; i++)
	{
		denom_DE += 1.0/var_DE[i];
		denom_M += 1.0/var_M[i];
	}
	/*printf("%f %f \n", denom_DE, denom_M);*/
	for(j=0; j<Nbin; j++)
	{
		DEm_Jack[j]=0;
		Mm_Jack[j]=0;
		for(i=2; i<maxnum-1; i++)
		{
			DE_J = func_DE( correlat[(i+N-1)%N][j], correlat[i][j], correlat[(i+1)%N][j] );
			DEm_Jack[j] += DE_J/var_DE[i];
			Mm_Jack[j] += func_M(correlat[i][j],DE_J,i)/var_M[i];
		}
		DEm_Jack[j] /= denom_DE;
		Mm_Jack[j] /= denom_M;
	}

	for(j=0; j<Nbin; j++)
	{
		DE_tot += DEm_Jack[j];
		M_tot += Mm_Jack[j];
		err_DE_tot += DEm_Jack[j]*DEm_Jack[j];
		err_M_tot += Mm_Jack[j]*Mm_Jack[j];
		/*
		if( j%((int)Nbin/10)==0 )
			printf("DE = %f %f %f %f\n", DE_tot, err_DE_tot, M_tot, err_M_tot );
		*/
	}
	/*printf("DE = %f %f %f %f\n", DE_tot, err_DE_tot, M_tot, err_M_tot );*/
	DE_tot /= (double)Nbin;
	err_DE_tot = (double)(Nbin-1)*( err_DE_tot/(double)Nbin - DE_tot*DE_tot );
	err_DE_tot = sqrt(err_DE_tot);

	M_tot /= (double)Nbin;
	err_M_tot = (double)(Nbin-1)*( err_M_tot/(double)Nbin - M_tot*M_tot );
	err_M_tot = sqrt(err_M_tot);

	printf("DE = %f +/- %f\n", DE_tot, err_DE_tot);	
	printf("M  = %f +/- %f\n", M_tot, err_M_tot);	
		
		/*++++++ write cm[t] on file ++++++*/
	doc_c = fopen( filename1 ,"wt");
	doc_E = fopen( filename2,"wt");
	doc_M = fopen( filename3,"wt");
	for(i=0; i<N; i++){
		double V= func_DE( cm[(i+N-1)%N], cm[i], cm[(i+1)%N]);
		fprintf( doc_c, "%d %f %f\n", i, cm[i], sqrt(err_cm[i]/(double)Nbin) );
		fprintf( doc_E, "%d %f %f\n", i, V, sqrt( var_DE[i] ) );
		fprintf( doc_M, "%d %f %f\n", i, func_M(cm[i],V,i), sqrt( var_M[i] ) );
	}
	free(correlat);
	fclose(doc_c);
	fclose(doc_E);
	fclose(doc_M);

	exit(EXIT_SUCCESS);	
}
