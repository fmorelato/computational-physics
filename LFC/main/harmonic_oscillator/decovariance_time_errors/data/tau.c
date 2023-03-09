/**********************************************************
 * File tau.c
 *
 * this is used to  perform a fitting procedure to estimate
 * the parameter tau for the autocorrelation function.
 * This is a weighted fitting
 *
 * compile with the command
 * 	gcc -std=c99 -Wall -o tau tau.c -lm
 **********************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
	int max_tm=50;
	FILE* doc = fopen( "gamma_data0000.txt", "rt");
	FILE* pro = fopen( "prov.txt", "wt");
	double b, c, A=0, B=0;;
	double z[max_tm], sz[max_tm];
	int count=0;
	
		/*
		 * the first line contains an exact value so
		 * that I ignore it in the fitting procedure
		 */
	fscanf(doc, "%lf", &b);
	fscanf(doc, "%lf", &b);
	fscanf(doc, "%lf", &b);
	while(0==0)
	{
		fscanf(doc, "%lf", &b);
		fscanf(doc, "%lf", &b);
		if(b<0 || count > max_tm-1)
		{
			break;
		}
		count++;
		z[count] = log(b);	/*  log(gamma(i)) */
		fscanf(doc, "%lf", &c);
		sz[count] = c/b;	/*  sigma(i)/gamma(i)  */
		fprintf(pro, "%d %f\n", count, z[count] );
	}
	fclose(doc);
	fclose(pro);
	for(int i=1; i<count+1; i++)
	{
		A += i*z[i]/(sz[i]*sz[i]);
		B += i*i/(sz[i]*sz[i]);
	}

	printf( "tau = %f +/- %f\n", -B/A, pow(B,1.5)/(A*A) );
	exit(EXIT_SUCCESS);
}
