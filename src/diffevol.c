/*
	Differential Evolution Optimizer

    Copyright (C) 2023  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich
*/
/*****************************************************************              
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
 *  IEK-5 PHOTOVOLTAIK                                           *              
 *                                                               *              
 *        ########                _   _                          *              
 *     ##########                |_| |_|                         *              
 *    ##########     ##         _ _   _ _     ___ ____ _   _     *              
 *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *              
 *   #########     #####    _  | | | | | |    | | |   | |_| |    *              
 *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *              
 *    #          ######     \___/ \___/|_____|___\____|_| |_|    *              
 *     ##      #######      F o r s c h u n g s z e n t r u m    *              
 *       ##########                                              *              
 *                                                               *              
 *   http://www.fz-juelich.de/iek/iek-5/DE/Home/home_node.html   *              
 *****************************************************************
 *                                                               *
 *    Dr. Bart E. Pieters 2023                                   *
 *                                                               *             
 *****************************************************************/ 
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "emalloc.h"
#include "diffevol.h"

/* basic implmentation of the differential evolution algorithm */
#define RandFloat(a,b) ((a)+((b)-(a))*((double)rand())/RAND_MAX)
#define RandInt(a,b) ((a)+rand()%((b)-(a)+1))
// default settings
de_settings de_default={
	0.3,
	0.8,
	1e-3,
	5000,
	50,
	0
};
	
// compute standard deviation of the population	cost function
double StdDev(double *E, int M)
{
	int i;
	double std=0, m=0;
	for (i=0;i<M;i++)
		m+=E[i];
	m/=(double)M;
	
	for (i=0;i<M;i++)
		std+=(E[i]-m)*(E[i]-m);
		
	return std/(double)M;
}

// return 3 random indices from the population
void ThreeIndex(int N, int *i1, int *i2, int *i3)
{
	(*i1)=RandInt(0,N);
	(*i2)=(*i1);
	while ((*i2)==(*i1))
		(*i2)=RandInt(0,N);
	(*i3)=(*i1);
	while (((*i3)==(*i1))||((*i3)==(*i2)))
		(*i3)=RandInt(0,N);
}

// loop to evolve the population
int diffevol_loop(double (*fn)(double x[], void *args), void *args, int n, double **pop, double lb[], double ub[], double *E, int ibest, int *NE, const de_settings *S)
{
	int i, j, k, i1, i2, i3;
	int iter=S->Npop;
	double *x, Ebest, Ex;
	double D;
	D=StdDev(E, S->Npop);
	Ebest=E[ibest];
	if((x=Emalloc(n*sizeof(double)))==NULL)
		return ibest;
		
    while ((iter < S->max_iter)&&(D>S->eps))
    {
		for (i=0;i<S->Npop;i++)
		{
			// crossover
			ThreeIndex(S->Npop-1,&i1, &i2, &i3);
			k=RandInt(0,n-1);
			for (j=0;j<n;j++)
			{
				if (j==k)
				{
					x[j]=pop[i1][j]+S->F*(pop[i2][j]-pop[i3][j]);
					if (x[j]<lb[j])
						x[j]=ub[j]-fmod(lb[j]-x[j],ub[j]-lb[j]);
					if (x[j]>ub[j])
						x[j]=lb[j]+fmod(x[j]-ub[j],ub[j]-lb[j]);
				}
				else
				{
					if (RandFloat(0,1)<S->CR)
					{
						x[j]=pop[i1][j]+S->F*(pop[i2][j]-pop[i3][j]);
						if (x[j]<lb[j])
							x[j]=ub[j]-fmod(lb[j]-x[j],ub[j]-lb[j]);
						if (x[j]>ub[j])
							x[j]=lb[j]+fmod(x[j]-ub[j],ub[j]-lb[j]);
					}
					else
						x[j]=pop[i][j];
				}
			}
			Ex=fn(x,args);
			iter++;
			if (isfinite(E[i]))
			{
				if (Ex<E[i])
				{
					if (Ex<Ebest)
					{
						Ebest=Ex;
						ibest=i;
					}
					E[i]=Ex;
					for (j=0;j<n;j++)
						pop[i][j]=x[j];
				}
			}
		}
		D=StdDev(E, S->Npop);
    }
    free(x);
    if (NE)
		(*NE)=iter;
	return ibest;
}

// main routine, generates a population and evolves it
double diffevol(double (*fn)(double x[], void *args), void *args, int n, double xbest[], double lb[], double ub[], int *NE, const de_settings *S)
{
	double **pop;
	double *E, Ebest;
	int ibest, i, j;
	time_t t;
	srand((unsigned) time(&t));
	if ((pop=Emalloc(S->Npop*sizeof(double *)))==NULL)	// population parameters
		return INFINITY;
	if((E=Emalloc(S->Npop*sizeof(double)))==NULL)		// error vector
	{
		free(pop);
		return INFINITY;
	}
	
	// initialize the population
	if((pop[0]=Emalloc(n*sizeof(double)))==NULL)
	{
		free(pop);
		free(E);
		return INFINITY;
	}
	if (S->prime) // in case a good initial guess is known we can integrate it in the population
	{
		for (j=0;j<n;j++)
			pop[0][j]=xbest[j];
	}
	else
	{	
		for (j=0;j<n;j++)
			pop[0][j]=RandFloat(lb[j],ub[j]);
	}
	E[0]=fn(pop[0],args);
	if (!isfinite(E[0]))
		E[0]=INFINITY;
	Ebest=E[0];
	ibest=0;
    for (j=0;j<n;j++)
		xbest[j]=pop[ibest][j];	
	
	for (i=1;i<S->Npop;i++)
	{
		if((pop[i]=Emalloc(n*sizeof(double)))==NULL)
		{
			i--;
			for (;i>=0;i--)
				free(pop[i]);
			free(pop);
			free(E);
			return Ebest;
		}
		for (j=0;j<n;j++)
			pop[i][j]=RandFloat(lb[j],ub[j]);
		E[i]=fn(pop[i],args);
		if (!isfinite(E[i]))
			E[i]=INFINITY;
		if (E[i]<Ebest)
		{
			ibest=i;
			Ebest=E[i];
			for (j=0;j<n;j++)
				xbest[j]=pop[ibest][j];	
		}
	}
	
	ibest=diffevol_loop(fn, args, n, pop, lb, ub, E, ibest, NE, S);
	
    for (j=0;j<n;j++)
		xbest[j]=pop[ibest][j];
	Ebest=E[ibest];
	for (i=0;i<S->Npop;i++)
		free(pop[i]);
	free(pop);
	free(E);
	return Ebest;
}

// main routine: same as above but thos version takes a priming population to seed the initial population with
double diffevol_prime(double (*fn)(double x[], void *args), void *args, int n, double **prime, int Nprime, double xbest[], double lb[], double ub[], int *NE, const de_settings *S)
{
	double **pop;
	double *E, Ebest=INFINITY;
	int ibest=0, i, j;
	time_t t;
	srand((unsigned) time(&t));
	if((pop=Emalloc(S->Npop*sizeof(double *)))==NULL)	// population parameters
		return INFINITY;
	if((E=Emalloc(S->Npop*sizeof(double)))==NULL)		// error vector
	{
		free(pop);
		return INFINITY;
	}
	
	// initialize the population
	i=0;
	while ((i<Nprime)&&(i<S->Npop))
	{
		if((pop[i]=Emalloc(n*sizeof(double)))==NULL)
		{
			i--;
			for (;i>=0;i--)
				free(pop[i]);
			free(pop);
			free(E);
			return Ebest;
		}
		for (j=0;j<n;j++)
			pop[i][j]=prime[i][j];
		E[i]=fn(pop[i],args);
		if (!isfinite(E[i]))
			E[i]=INFINITY;
		if (i==0)
		{
			Ebest=E[0];
			ibest=0;	
		}
		else if (E[i]<Ebest)
		{
			ibest=i;
			Ebest=E[i];
		}
		i++;
	}
	
	for (;i<S->Npop;i++)
	{
		if((pop[i]=Emalloc(n*sizeof(double)))==NULL)
		{
			i--;
			for (;i>=0;i--)
				free(pop[i]);
			free(pop);
			free(E);
			return Ebest;
		}
		for (j=0;j<n;j++)
			pop[i][j]=RandFloat(lb[j],ub[j]);
		E[i]=fn(pop[i],args);
		if (!isfinite(E[i]))
			E[i]=INFINITY;
		if (i==0)
		{
			Ebest=E[0];
			ibest=0;	
		}
		else if (E[i]<Ebest)
		{
			ibest=i;
			Ebest=E[i];
		}
	}
	
	ibest=diffevol_loop(fn, args, n, pop, lb, ub, E, ibest, NE, S);
	
    for (j=0;j<n;j++)
		xbest[j]=pop[ibest][j];
	Ebest=E[ibest];
	for (i=0;i<S->Npop;i++)
		free(pop[i]);
	free(pop);
	free(E);
	return Ebest;
}
