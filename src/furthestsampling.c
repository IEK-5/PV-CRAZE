/*
    simple (naive?) implementation of furthest sampling
    Copyright (C) 2024  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich
*/
/*****************************************************************              
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
 +  IEK-5 PHOTOVOLTAIK                                           *              
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
 *   Copyright (C) Dr. Bart E. Pieters 2024                      *
 *                                                               *             
 *****************************************************************/ 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "emalloc.h"
#include "curve.h"

// compute mean
double FS_mean(double *d, int N)
{
	double S=0;
	int i;
	
	for (i=0;i<N;i++)
		S+=d[i];
	S/=(double)N;
	return S;
}

// compute variance
double FS_variance(double *d, int N, double mu)
{
	double S=0, s;
	int i;
	
	for (i=0;i<N;i++)
	{
		s=(d[i]-mu);
		S+=s*s;
	}
	S/=(double)N;
	return S;
}

// substract mean and scale to standard deviation
void FS_Scaler(double *d, int N) // in place scaler
{
	int i;
	double mu, s;
	mu=FS_mean(d, N);
	s=sqrt(FS_variance(d, N, mu));
	for (i=0;i<N;i++)
		d[i]=(d[i]-mu)/s;
}

// return distance²
double FS_dist(double x1, double y1, double x2, double y2)
{
	double dx,dy;
	dx=(x2-x1);
	dy=(y2-y1);
	return dx*dx+dy*dy;
}

// select element with index s (mark in sel array) and update the inverse distance array (di)
// used in routine FS_Select
int FS_SelectS(const double *x, const double *y, double *di, char *sel, int N, int s)
{
	double mind;
	int i=0, imin=-1;
	sel[s]=1; // mark s as selected
	while (i<N)
	{
		if (!sel[i])
			break;
		i++;
	}
	if (i<N)
	{
		di[i]+=1/FS_dist(x[s],y[s],x[i],y[i]);
		imin=i;
		mind=di[i];
		i++;
		for (;i<N;i++) // update inverse distances
		{
			if (!sel[i])
			{
				di[i]+=1/FS_dist(x[s],y[s],x[i],y[i]);
				if (di[i]<mind)
				{
					imin=i;
					mind=di[i];
				}
			}
		}
	}
	return imin;
}

// creates a character array which marks selected elements
// tries to select M elements, returns number of selected elements
int FS_Select(const double *x, const double *y, int N, char **SEL, int M, int off)
{
	char *sel;
	double *di;
	int s,i;
	
	if (M>N)
	{
		(*SEL)=NULL;
		return 0;
	}
	
	// alocate selection
	if((sel=Ecalloc(N,sizeof(char)))==NULL)
	{
		(*SEL)=NULL;
		return 0;
	}	
	(*SEL)=sel;
	// alocate inverse distance array
	if((di=Ecalloc(N,sizeof(double)))==NULL)
	{
		free(sel);
		(*SEL)=NULL;
		return 0;
	}	
	
	i=0;
	s=off; // initial point
	while (i<M)
	{
		s=FS_SelectS(x,y,di,sel,N,s);
		if (s<0)
			break;
		i++;
	}
	free(di);
	return i;
}

// select elements marked in the sel array from original data
double *FS_SelD(const double *d, int N, char *sel, int M)
{
	double *o;
	int i,j=0;
	if((o=Emalloc(M*sizeof(double)))==NULL)
		return NULL;
		
	for (i=0;i<N;i++)
	{
		 if (sel[i])
		 {
			 o[j]=d[i];
			 j++;
			 if (j==M)
				break;
		 }
	}
	if (j!=M)
	{
		free(o);
		o=NULL;
	}
	return o; 
}

// select elements marked in the sel array from original data, now for curvepars
curvepars *FS_SelC(const curvepars *d, int N, char *sel, int M)
{
	curvepars *o;
	int i,j=0;
	if((o=Emalloc(M*sizeof(curvepars)))==NULL)
		return NULL;
		
	for (i=0;i<N;i++)
	{
		 if (sel[i])
		 {
			 o[j]=d[i];
			 j++;
			 if (j==M)
				break;
		 }
	}
	if (j!=M)
	{
		free(o);
		o=NULL;
	}
	return o; 
}

// allocate and copy data
double * FS_CopyData(const double *in, int N)
{
	int i=0;
	double *copy;
	if((copy=Emalloc(N*sizeof(double)))==NULL)
		return NULL;
	for (i=0;i<N;i++)
		copy[i]=in[i];
	return copy;
}
// allocate and copy curvepars
curvepars * FS_CopyPars(const curvepars *in, int N)
{
	int i=0;
	curvepars *copy;
	if((copy=Emalloc(N*sizeof(curvepars)))==NULL)
		return NULL;
	for (i=0;i<N;i++)
		copy[i]=in[i];
	return copy;
}
/* Routine to select data far apart in the G,T plane
 * Input: irradiance (G) and temperature (T) and esp (C)
 * Output: selected dataset of M elements at irradiance (Gs), temperature (Ts), and with esp (Cs)
 * Elements selected such that the elements are far way from each other in the G-T plane
 * The offset value sets the first element to choose (i.e. 0 chooses the first in the list)
 */
int FurthestSampler(const double *G, const double *T, const curvepars *C, 
                    int N, double **Gs, double **Ts, curvepars **Cs, int M, int off)
{
	double *x, *y; // working copies of X and Y
	double *xs, *ys; // working copies of selection
	curvepars *cs; // working copies of selection
	char *sel;
	int m;
	
	(*Gs)=NULL;
	(*Ts)=NULL;
	(*Cs)=NULL;
	// sanity check
	if (!M)
	{
		fprintf(stderr, "Error: cannot select %d points in furtherst sampling\n", M);
		return 0;
	}
	
	// copy data
	x=FS_CopyData(G, N);
	if (x==NULL)
		return 0;	
			
	y=FS_CopyData(T, N);
	if (y==NULL)
	{
		free(x);
		return 0;
	}
	
	if (M>=N)
	{
		// just copy
		(*Cs)=FS_CopyPars(C, N);
		if ((*Cs)==NULL)
		{
			free(x);
			free(y);
			return 0;
		}
		(*Gs)=x;
		(*Ts)=y;
		return N;
	}
	
	off=abs(off)%N; // make sure offset is positive and smaller than N
	
	// scaling
	FS_Scaler(x, N);
	FS_Scaler(y, N);
	
	// select data
	m=FS_Select(x, y, N, &sel, M, off);
	if (m!=M)
	{
		free(x);
		free(y);
		if (sel)
			free(sel);
		fprintf(stderr, "Error: selection failed in furtherst sampling\n");
		return 0;
	}
	
	// make selected data arrays
	xs=FS_SelD(G, N, sel, m);
	if (xs==NULL)
	{
		free(x);
		free(y);
		free(sel);
		fprintf(stderr, "Error: sampling failed in furtherst sampling\n");
		return 0;
	}
	
	ys=FS_SelD(T, N, sel, m);
	if (ys==NULL)
	{
		free(x);
		free(y);
		free(sel);
		free(xs);
		fprintf(stderr, "Error: sampling failed in furtherst sampling\n");
		return 0;
	}
	
	cs=FS_SelC(C, N, sel, m);
	if (cs==NULL)
	{
		free(x);
		free(y);
		free(sel);
		free(xs);
		free(ys);
		fprintf(stderr, "Error: sampling failed in furtherst sampling\n");
		return 0;
	}
	
	// success!
	(*Gs)=xs;
	(*Ts)=ys;
	(*Cs)=cs;
	
	// clean up
	free(x);
	free(y);
	free(sel);
	return m;
}
