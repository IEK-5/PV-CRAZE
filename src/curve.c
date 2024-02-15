/*
    Curve analysis routines. Extract ESP and Flag faulty curves
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
#include <stdbool.h>
#include <math.h>
#include "emalloc.h"
#include "Error.h"
#include "curve.h"

#define MINPOINTS 7 // minimum number of points in an IV. I currently use a kernelsize of 7 for MPP
curvepars nullpars={0,0,0,0,0,0,0,0,0,0};

/* ******************************************************************
 * Preprocessing routines:
 * Routines to shape the input data:
 *   - Spike filter
 *   - remove duplicate values
 *   - filter out the active range
 *   - sort voltages in ascending order
 * ****************************************************************** */

/* third order median filter (for removing spikes) */
int MedianFilter3rd(curve *C)
{
	int i, index;
	double c;
	curve B;
	
	double x[3];
	double y[3];
	if (C->N<MINPOINTS+2)
		return 1; // insufficient data
	
	/* make a copy to work on */
	B.N=C->N;
	if((B.P=Emalloc(B.N*sizeof(point)))==NULL)
		return 2;
	for (i=0;i<B.N;i++)
		B.P[i]=C->P[i];
	for (index=1;index<B.N-1;index++)
	{
		/* construct arrays x and y at this position */
		for (i=0;i<3;i++)
		{
			x[i]=B.P[index+i-1].x;
			y[i]=B.P[index+i-1].y;
		}
		
		if (x[0]>x[1])
		{
			c=x[0];
			x[0]=x[1];
			x[1]=c;
		}
		if (x[1]>x[2])
		{
			c=x[1];
			x[1]=x[2];
			x[2]=c;
			if (x[0]>x[1])
			{
				c=x[0];
				x[0]=x[1];
				x[1]=c;
			}
		}
		if (y[0]>y[1])
		{
			c=y[0];
			y[0]=y[1];
			y[1]=c;
		}
		if (y[1]>y[2])
		{
			c=y[1];
			y[1]=y[2];
			y[2]=c;
			if (y[0]>y[1])
			{
				c=y[0];
				y[0]=y[1];
				y[1]=c;
			}
		}
		
		/* apply median filter */
		C->P[index].x=x[1]; 
		C->P[index].y=y[1]; 
	}
	free(B.P);
	/* remove the endpoints, we cannot filter spikes there */
	for (i=0;i<C->N-2;i++)
		C->P[i]=C->P[i+1];
	C->P[C->N-2].x=NAN;
	C->P[C->N-2].y=NAN;
	C->P[C->N-1].x=NAN;
	C->P[C->N-1].y=NAN;
	C->N-=2;
	return 0;
}

/* sometimes when measurements go bad we have identical values. At the very least this adds no information. */
int RemoveIdentical(curve *C, double EPS)
{
	int i, j;
	int r=0;
	int N0=C->N;
	for (i=0;i<C->N-1;i++)
	{
		if ((fabs(C->P[i].x-C->P[i+1].x)<EPS)||(fabs(C->P[i].y-C->P[i+1].y)<EPS))
		{
			if (fabs(C->P[i].y-C->P[i+1].y)>EPS)
				C->P[i].y=(C->P[i].y+C->P[i+1].y)/2;
			else if (fabs(C->P[i].x-C->P[i+1].x)>EPS)
				C->P[i].x=(C->P[i].x+C->P[i+1].x)/2;
			/* delete point i+1 */
			for (j=i+1;j<C->N;j++)
				C->P[j]=C->P[j+1];
			C->N--;
			r++;
		} 
	}
	/* pad remainder with NAN's */
	for (i=C->N;i<N0;i++)
	{
		C->P[i].x=NAN;
		C->P[i].y=NAN;
	}
	return r;
}
void IdenticalFilter(curve *C, double EPS)
{
	int r=1;
	while ((r>0)&&(C->N>1))
		r=RemoveIdentical(C, EPS);
}



// utility function for sorting
#define EPS 1e-12
int PointCompare (const void * a, const void * b)
{
	double A, B;
	A=(*(point*)a).x;
	B=(*(point*)b).x;
	A=(A-B)/(fabs(B)+fabs(A)+EPS);
	if (A>EPS)
		return 1;
	if (A<-EPS)
		return -1;
	return 0;
}


// find active range
void ActiveRangeIndex(curve C, int *start, int *end, int *Isign)
{
	int i=0, j;
	while((i<C.N-1)&&(C.P[i].x<0))
		i++;
	j=i+1;
	while((j<C.N)&&(C.P[j].y*C.P[i].y>=0))
		j++;
	(*Isign)=2*(C.P[i].y>0)-1;
	(*start)=i;
	(*end)=j;
}

// preprocessor
curve PreprocessCurve(curve C, int ext)
{
	curve N;
	int i=0,j=C.N,k,s=1;
	if (C.N==0)
	{
		N.P=NULL;
		N.N=0;
		return N;
	}
	
	qsort (C.P, C.N, sizeof(point), PointCompare);
	ActiveRangeIndex(C, &i, &j, &s);
	if (ext>0)
	{
		i-=ext;
		j+=ext;
		if (i<0)
			i=0;
		if (j>C.N)
			j=C.N;
	}
	
	N.N=(j-i);
	if (N.N==0)
	{
		N.P=NULL;
		return N;
	}
	
	if((N.P=Emalloc(N.N*sizeof(point)))==NULL)
	{
		N.P=NULL;
		N.N=0;
		return N;
	}
	
	for (k=i;k<j;k++)
	{
		N.P[k-i]=C.P[k];
		if (s<0)
			N.P[k-i].y*=-1.0;
	}
	return N;
}

// find the minimal and maximal values (bounding box)
void MinMax(curve C, point *min, point *max)
{
	int i;
	double xmin, ymin, xmax, ymax;
	if (!C.N)
	{
		min->x=0;
		min->y=0;
		max->x=0;
		max->y=0;
		return;
	}
		
	xmin=C.P[0].x;
	xmax=C.P[0].x;
	ymin=C.P[0].y;
	ymax=C.P[0].y;
	
	for (i=1;i<C.N;i++)
	{
		if (C.P[i].x<xmin)
			xmin=C.P[i].x;
		else if (C.P[i].x>xmax)
			xmax=C.P[i].x;
			
		if (C.P[i].y<ymin)
			ymin=C.P[i].y;
		else if (C.P[i].y>ymax)
			ymax=C.P[i].y;
	}
	min->x=xmin;
	min->y=ymin;
	max->x=xmax;
	max->y=ymax;

}
curve CheckAndFilterCurve(curve C, int *FLAGCODE, analyze_settings settings)
{
	curve Cn={0};
	point min, max;
	if (C.N==0)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFIVDATA);
		return Cn;
	}
	/* limit to active range, make sure currents are positive */
	Cn=PreprocessCurve(C,0);	
	if (Cn.N==0)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return Cn;
	}
	if (Cn.N<MINPOINTS)
	{
		//IVFLAGERR  FLAG_INSUFFIVDATA  "Insufficient data points in IV characteristics"
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFIVDATA);
		return Cn;		
	}
	MinMax(Cn, &min, &max);
	if (min.y/(max.y+1e-12)>1.0-settings.MinRelRange)
	{
		//IVFLAGERR  FLAG_INSUFF_I_RANGE "Relative Current range too small"
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_I_RANGE);
		return Cn;		
	}
		
	
	if (min.x/(max.x+1e-12)>1.0-settings.MinRelRange)
	{
		//IVFLAGERR  FLAG_INSUFF_V_RANGE "Relative Voltage range too small"
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_V_RANGE);
		return Cn;	
	}
		
	if (settings.Identicals)
		IdenticalFilter(&Cn, settings.IdentEps);
		
	if (settings.SpikeFilter)
	{
		int r;
		r=MedianFilter3rd(&Cn);
		if (r==2)	
		{
			(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
			return Cn;
		}
		else if (r==1)
		{
			(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFIVDATA);
			return Cn;		
		}		
		if (settings.Identicals)
			IdenticalFilter(&Cn, settings.IdentEps);
	}	
	
	if (Cn.N<MINPOINTS)
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFIVDATA);
	return Cn;	
}

/* symmetric 3rd order SG filters with kernelsize 7 */
/* determine a0 a1 a2 and a3 */ 
const double LSFIT_33[4][7] = {
				{-9.52380952380952328e-02,  1.42857142857142766e-01,  2.85714285714285532e-01,  3.33333333333333148e-01,  2.85714285714285532e-01,  1.42857142857142766e-01, -9.52380952380952328e-02},
				{ 8.73015873015873245e-02, -2.65873015873015428e-01, -2.30158730158729813e-01,  0.00000000000000000e+00,  2.30158730158729813e-01,  2.65873015873015428e-01, -8.73015873015873245e-02},
				{ 5.95238095238094997e-02,  0.00000000000000000e+00, -3.57142857142856984e-02, -4.76190476190475956e-02, -3.57142857142856984e-02,  0.00000000000000000e+00,  5.95238095238094997e-02},
				{-2.77777777777777485e-02,  2.77777777777777346e-02,  2.77777777777777415e-02,  0.00000000000000000e+00, -2.77777777777777415e-02, -2.77777777777777346e-02,  2.77777777777777485e-02}
			    };

// find the MPP using 3rd order Savitzky-Golay filter
void MppByLeastSsquaresFit(curve C, int i_mpp, double *Vmpp, double *Impp, int *FLAGCODE, analyze_settings settings)
{
	int i,j;
	int idx[7]={-3,-2,-1,0,1,2,3};
	double av[4], ai[4], ap[4];
	double p1=-1, p2=-1, t1, t2, D;
	(*Vmpp)=C.P[i_mpp].x;
	(*Impp)=C.P[i_mpp].y;
	
	if ((i_mpp<3)||(C.N-1-i_mpp<3))
	{
		//IVFLAGERR  FLAG_INSUFFMPP  "Insufficient data around mpp"
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFMPP);
		(*Vmpp)=C.P[i_mpp].x;
		(*Impp)=C.P[i_mpp].y;
		return;
	}
	i=i_mpp-4;
	// spread the points over a proper range
	while ((i>0)&&(C.P[i].x>(1.0-settings.MPP_f_V)*C.P[i_mpp].x))
	{
		if (C.P[i_mpp+idx[1]-1].x>(1.0-2.0*settings.MPP_f_V/3.0)*C.P[i_mpp].x)
		{
			if (C.P[i_mpp+idx[2]-1].x>(1.0-settings.MPP_f_V/3.0)*C.P[i_mpp].x)
				idx[2]--;
			idx[1]--;
		}
		idx[0]--;
		i--;
	}
	i=i_mpp+4;
	while ((i<C.N-1)&&(C.P[i].y>(1.0-settings.MPP_f_I)*C.P[i_mpp].y))
	{
		if (C.P[i_mpp+idx[5]+1].y>(1.0-2.0*settings.MPP_f_I/3.0)*C.P[i_mpp].y)
		{
			if (C.P[i_mpp+idx[4]+1].y>(1.0-settings.MPP_f_I/3.0)*C.P[i_mpp].y)
				idx[4]++;
			idx[5]++;
		}
		idx[6]++;
		i++;
	}
	
	
	for (i=0;i<4;i++)
	{
		av[i]=0;
		ai[i]=0;
		ap[i]=0;
		for (j=0;j<7;j++)
		{
			av[i]+=LSFIT_33[i][j]*C.P[i_mpp+idx[j]].x;
			ai[i]+=LSFIT_33[i][j]*C.P[i_mpp+idx[j]].y;
			ap[i]+=LSFIT_33[i][j]*C.P[i_mpp+idx[j]].x*C.P[i_mpp+idx[j]].y;
			
		}
	}
	D=ap[2]*ap[2]-3*ap[1]*ap[3];
	
	if (D<0)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_NOMPP);
		(*Vmpp)=C.P[i_mpp].x;
		(*Impp)=C.P[i_mpp].y;
		return;
	}
	
	t1=-(sqrt(D)+ap[2])/(3*ap[3]);
	
	if (fabs(t1)<3)
		p1=ap[0]+ap[1]*t1+2*ap[2]*t1*t1+3*ap[3]*t1*t1*t1;
	
	t2=(sqrt(D)-ap[2])/(3*ap[3]);
	
	if (fabs(t2)<3)
		p2=ap[0]+ap[1]*t2+2*ap[2]*t2*t2+3*ap[3]*t2*t2*t2;
		
	if (p1*p2>0)
	{
		/* either 0 or two solutions */
		//IVFLAGERR  FLAG_NOMPP  "No power maximum"
		//IVFLAGWARN  FLAG_MULTIPLEMPP  "more than one power maximum"
		if (p1<0)
			(*FLAGCODE)|=FLAG_CODE(FLAG_NOMPP);
		else
			(*FLAGCODE)|=FLAG_CODE(FLAG_MULTIPLEMPP);
		(*Vmpp)=C.P[i_mpp].x;
		(*Impp)=C.P[i_mpp].y;
		return;
	}
	
	if (p1>p2)
	{
		(*Vmpp)=av[0]+av[1]*t1+2*av[2]*t1*t1+3*av[3]*t1*t1*t1;
		(*Impp)=ai[0]+ai[1]*t1+2*ai[2]*t1*t1+3*ai[3]*t1*t1*t1;
	}
	else
	{
		(*Vmpp)=av[0]+av[1]*t2+2*av[2]*t2*t2+3*av[3]*t2*t2*t2;
		(*Impp)=ai[0]+ai[1]*t2+2*ai[2]*t2*t2+3*ai[3]*t2*t2*t2;	
	}
}

// find the MPP using SVD
extern void dgelsd_( int* m, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, double* s, double* rcond, int* rank,
                double* work, int* lwork, int* iwork, int* info );

void GetMPPLSQ(curve *C, int i_mpp, double *Vmpp, double *Impp, int *FLAGCODE, analyze_settings settings)
{
	double *M, *b, *s;
	int N, Nr, is, ie, i, j, k;
	
	int m, n, nrhs, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	double wkopt;
	double* work;
	double xs, xe;
    int *iwork;
    int order=3;
		
    
    // initial guess
	(*Vmpp)=C->P[i_mpp].x;
	(*Impp)=C->P[i_mpp].y;
	
	// select points for fit
	is=i_mpp;
	while ((is>0)&&(C->P[is].x>(1.0-settings.MPP_f_V)*C->P[i_mpp].x))
		is--;
	ie=i_mpp;
	while ((ie<C->N-1)&&(C->P[ie].y>(1.0-settings.MPP_f_I)*C->P[i_mpp].y))
		ie++;
	
	Nr=(ie-is);
	if (Nr<order+1)
	{
		/*IVFLAGERR  FLAG_INSUFFMPP  "Insufficient data around mpp"*/
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFMPP);
		return;
	}
	if (settings.am=='p')
	{
		xs=0.0;
		xe=(double)Nr;
		nrhs=3; // parametric fit
	}
	else
	{
		xs=	C->P[is].x;
		xe=	C->P[ie].x;
		nrhs=2; // function of voltage
	}
	N=(order+1)*Nr;
	if((M=Emalloc(N*sizeof(double)))==NULL)
	{
		//IVFLAGERR  FLAG_MALLOCFAIL  "Failed to analyze, malloc fail"
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	for (j=order;j>=0;j--)
		for (i=0;i<Nr;i++)
		{
			M[(order-j)*Nr+i]=1.0;
			for (k=0;k<j;k++)
			{
				if (settings.am=='p')
					M[(order-j)*Nr+i]*=(double)i; // parametric
				else
					M[(order-j)*Nr+i]*=C->P[is+i].x; // function of voltage
			}
		}
	if((b=Emalloc(nrhs*Nr*sizeof(double)))==NULL)
	{
		free(M);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	for (i=0;i<Nr;i++)
		b[i]=C->P[is+i].y;
	for (i=0;i<Nr;i++)
		b[Nr+i]=C->P[is+i].x*C->P[is+i].y;	
	if (settings.am=='p')
	{
		for (i=0;i<Nr;i++)
			b[i+2*Nr]=C->P[is+i].x;
	}
		
	if((s=Emalloc((order+1)*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	n=Nr;
	m=order+1;
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(14*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFMPP);
		return;
	}	
    lwork = -1;
    dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
                        iwork, &info );
	lwork = (int)wkopt;
	work = (double*)Emalloc( lwork*sizeof(double) );
	if (!work)
	{
		free(M);
		free(b);
		free(s);
		free(iwork);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info == 0 )
	{
#define A_i b[0]
#define B_i b[1]
#define C_i b[2]
#define D_i b[3]
#define A_p b[0+Nr]
#define B_p b[1+Nr]
#define C_p b[2+Nr]
#define D_p b[3+Nr]
#define A_v b[0+2*Nr]
#define B_v b[1+2*Nr]
#define C_v b[2+2*Nr]
#define D_v b[3+2*Nr]
		double t1,t2,t,D, P1=-1, P2=-1;
		int Ns=0;
		D=B_p*B_p-3*A_p*C_p;
		if (D>0)
		{
			t1=-(sqrt(D)+B_p)/(3*A_p);
			t2= (sqrt(D)-B_p)/(3*A_p);
			if ((t1>xs)&&(t1<xe))
			{
				Ns++;
				P1=((A_p*t1+B_p)*t1+C_p)*t1+D_p;
			}
			if ((t2>xs)&&(t2<xe))
			{
				Ns++;
				P1=((A_p*t2+B_p)*t2+C_p)*t2+D_p;
			}
			switch(Ns)
			{
				case 0:
					(*FLAGCODE)|=FLAG_CODE(FLAG_NOMPP); // cannot determine mpp
					goto FreeGetMPPLSQ;
				case 2:
					(*FLAGCODE)|=FLAG_CODE(FLAG_MULTIPLEMPP); // get the largest one
				case 1:
				default:
					if (P1>P2)
						t=t1;
					else
						t=t2;
					break;
			}
			if (settings.am=='p')
				(*Vmpp)=((A_v*t+B_v)*t+C_v)*t+D_v;
			else
				(*Vmpp)=t;
			(*Impp)=((A_i*t+B_i)*t+C_i)*t+D_i;
		}
		else	
		{
			(*FLAGCODE)|=FLAG_CODE(FLAG_NOMPP);
			goto FreeGetMPPLSQ;
		}	
#undef A_v
#undef B_v
#undef C_v
#undef D_v
#undef A_i
#undef B_i
#undef C_i
#undef D_i
#undef A_p
#undef B_p
#undef C_p
#undef D_p
	}
	else
	{
		//IVFLAGWARN  FLAG_SVDFAIL  "SVD failed, could not solve OLS mpp"
		(*FLAGCODE)|=FLAG_CODE(FLAG_SVDFAIL);
	}
FreeGetMPPLSQ:
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
}

#define LSQMPP   
void GetMpp(curve Cn, int *FLAGCODE, double *Vmpp, double *Impp, int *iim, analyze_settings settings)
{
	int i, im=0;
	double Pmax;
	Pmax=Cn.P[0].x*Cn.P[0].y;
	for (i=1;i<Cn.N;i++)
	{
		if (Cn.P[i].x*Cn.P[i].y>Pmax)
		{
			im=i;
			Pmax=Cn.P[i].x*Cn.P[i].y;
		}
	}
	switch (settings.am)
	{
		case 'f':
			// no noise accurate curve:
			
			/* MPP */
			if ((im==0)||(im==Cn.N-1))
			{
				(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFFMPP);
				(*Vmpp)=Cn.P[im].x;
				(*Impp)=Cn.P[im].y;
			}
			else
			{
				double t, a, b;
				/* solve quadratic function through 3 points around mpp and find the maximum */
				a=(Cn.P[im+1].x*Cn.P[im+1].y + Cn.P[im-1].x*Cn.P[im-1].y)/2-Pmax;
				b=(Cn.P[im+1].x*Cn.P[im+1].y - Cn.P[im-1].x*Cn.P[im-1].y)/2;
				/* solve derivative = 0 */
				t=-b/2/a; 
				if ((t<-1)||(t>1))
				{
					(*FLAGCODE)|=FLAG_CODE(FLAG_NOMPP);
					(*Vmpp)=Cn.P[im].x;
					(*Impp)=Cn.P[im].y;
				}
				else
				{
					a=(Cn.P[im+1].x + Cn.P[im-1].x)/2-Cn.P[im].x;
					b=(Cn.P[im+1].x - Cn.P[im-1].x)/2;
					(*Vmpp)=a*t*t+b*t+Cn.P[im].x;
					
					a=(Cn.P[im+1].y + Cn.P[im-1].y)/2-Cn.P[im].y;
					b=(Cn.P[im+1].y - Cn.P[im-1].y)/2;
					(*Impp)=a*t*t+b*t+Cn.P[im].y;
				}
			}
			break;
		case 's':
			MppByLeastSsquaresFit(Cn, im, Vmpp, Impp, FLAGCODE, settings);
			break;
		case 'r':
		case 'p':
		default:
			// normal
			GetMPPLSQ(&Cn, im, Vmpp, Impp, FLAGCODE, settings);
			break;
	}
	(*iim)=im;
}

void GetQMPPs(curve Cn, int *FLAGCODE, double Vmpp, double Impp, double *Vqmp1, double *Iqmp1, double *Vqmp2, double *Iqmp2, analyze_settings settings)
{
	int i, im=0;
	double Pmax;
	curve Cx={0};
	int f=0;
	
	if (Cn.N<=0)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_QMPM);
		(*Vqmp1)=0;
		(*Iqmp1)=0;
		(*Vqmp2)=0;
		(*Iqmp2)=0;
	}		
	Cx=PreprocessCurve(Cn,1);
	if (Cx.N==0)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		(*Vqmp1)=0;
		(*Iqmp1)=0;
		(*Vqmp2)=0;
		(*Iqmp2)=0;
	}		
		
	Cx.P[0].y-=Impp;
	Pmax=Cx.P[0].x*Cx.P[0].y;
	for (i=1;i<Cn.N;i++)
	{
		Cx.P[i].y-=Impp;
		if (Cx.P[i].x*Cx.P[i].y>Pmax)
		{
			im=i;
			Pmax=Cx.P[i].x*Cx.P[i].y;
		}
	}
	switch (settings.am)
	{
		case 'f':
			/* MPP */
			if ((im==0)||(im==Cn.N-1))
			{
				f|=FLAG_CODE(FLAG_INSUFFMPP);
				(*Vqmp1)=Cn.P[im].x;
				(*Iqmp1)=Cn.P[im].y;
			}
			else
			{
				double t, a, b;
				/* solve quadratic function through 3 points around mpp and find the maximum */
				a=(Cn.P[im+1].x*Cn.P[im+1].y + Cn.P[im-1].x*Cn.P[im-1].y)/2-Pmax;
				b=(Cn.P[im+1].x*Cn.P[im+1].y - Cn.P[im-1].x*Cn.P[im-1].y)/2;
				/* solve derivative = 0 */
				t=-b/2/a; 
				if ((t<-1)||(t>1)) // no solution
				{
					(*Vqmp1)=Cn.P[im].x;
					(*Iqmp1)=Cn.P[im].y;
				}
				else
				{
					a=(Cn.P[im+1].x + Cn.P[im-1].x)/2-Cn.P[im].x;
					b=(Cn.P[im+1].x - Cn.P[im-1].x)/2;
					(*Vqmp1)=a*t*t+b*t+Cn.P[im].x;
					
					a=(Cn.P[im+1].y + Cn.P[im-1].y)/2-Cn.P[im].y;
					b=(Cn.P[im+1].y - Cn.P[im-1].y)/2;
					(*Iqmp1)=a*t*t+b*t+Cn.P[im].y;
				}
			}
			break;
		case 's':
			MppByLeastSsquaresFit(Cx, im, Vqmp1, Iqmp1, &f, settings);
			break;
		case 'r':
		case 'p':
		default:
			// normal robust method
			GetMPPLSQ(&Cx, im, Vqmp1, Iqmp1, &f, settings);
	}
	(*Iqmp1)+=Impp;
	switch (f)
	{
		case 0:
			break;
		case FLAG_CODE(FLAG_INSUFFMPP):
				//IVFLAGERR  FLAG_INSUFF_QMPM  "Insufficient data to determine the qmp- point"
				(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_QMPM);
				break;
		case FLAG_CODE(FLAG_NOMPP):
				//IVFLAGERR  FLAG_NOQMPM  "the qmp- point does not exist"
				(*FLAGCODE)|=FLAG_CODE(FLAG_NOQMPM);
				break;
		case FLAG_CODE(FLAG_MULTIPLEMPP):
				//IVFLAGWARN  FLAG_MULTIPLEQMPM  "more than one qmp- points"
				(*FLAGCODE)|=FLAG_CODE(FLAG_MULTIPLEQMPM);
		default:
				break;
	} 
	f=0;				
	im=0;
	
	Cx.P[0].y+=Impp;
	Cx.P[0].x-=Vmpp;
	Pmax=Cx.P[0].x*Cx.P[0].y;
	for (i=1;i<Cn.N;i++)
	{
		Cx.P[i].y+=Impp;
		Cx.P[i].x-=Vmpp;
		if (Cx.P[i].x*Cx.P[i].y>Pmax)
		{
			im=i;
			Pmax=Cx.P[i].x*Cx.P[i].y;
		}
	}
	switch (settings.am)
	{
		case 'f':
			/* MPP */
			if ((im==0)||(im==Cn.N-1))
			{
				f|=FLAG_CODE(FLAG_INSUFFMPP);
				(*Vqmp2)=Cn.P[im].x;
				(*Iqmp2)=Cn.P[im].y;
			}
			else
			{
				double t, a, b;
				/* solve quadratic function through 3 points around mpp and find the maximum */
				a=(Cn.P[im+1].x*Cn.P[im+1].y + Cn.P[im-1].x*Cn.P[im-1].y)/2-Pmax;
				b=(Cn.P[im+1].x*Cn.P[im+1].y - Cn.P[im-1].x*Cn.P[im-1].y)/2;
				/* solve derivative = 0 */
				t=-b/2/a; 
				if ((t<-1)||(t>1))
				{
					(*Vqmp2)=Cn.P[im].x;
					(*Iqmp2)=Cn.P[im].y;
				}
				else
				{
					a=(Cn.P[im+1].x + Cn.P[im-1].x)/2-Cn.P[im].x;
					b=(Cn.P[im+1].x - Cn.P[im-1].x)/2;
					(*Vqmp2)=a*t*t+b*t+Cn.P[im].x;
					
					a=(Cn.P[im+1].y + Cn.P[im-1].y)/2-Cn.P[im].y;
					b=(Cn.P[im+1].y - Cn.P[im-1].y)/2;
					(*Iqmp2)=a*t*t+b*t+Cn.P[im].y;
				}
			}
			break;
		case 's':
			MppByLeastSsquaresFit(Cx, im, Vqmp2, Iqmp2, &f, settings);
			break;
		case 'r':
		case 'p':
		default:
			// normal
			GetMPPLSQ(&Cx, im, Vqmp2, Iqmp2, &f, settings);
	}
	
	(*Vqmp2)+=Vmpp;
	switch (f)
	{
		case 0:
			break;
		case FLAG_CODE(FLAG_INSUFFMPP):
				//IVFLAGERR  FLAG_INSUFF_QMPP  "Insufficient data to determine the qmp+ point"
				(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_QMPP);
				break;
		case FLAG_CODE(FLAG_NOMPP):
				//IVFLAGERR  FLAG_NOQMPP  "the qmp+ point does not exist"
				(*FLAGCODE)|=FLAG_CODE(FLAG_NOQMPP);
				break;
		case FLAG_CODE(FLAG_MULTIPLEMPP):
				//IVFLAGWARN  FLAG_MULTIPLEQMPP  "more than one qmp+ points"
				(*FLAGCODE)|=FLAG_CODE(FLAG_MULTIPLEQMPP);
		default:
				break;
	} 
	free(Cx.P);
}


// determine whether the extracted parameters can produce a monotone curve
int FlagMono(curvepars par)
{
	int section=1, FLAG=0;
	
	if (par.Gsc<0)
		FLAG|=section;
	section<<=1;	
	if (par.Isc<par.Iqmp1)
		FLAG|=section;
	section<<=1;	
	if (par.Iqmp1<par.Impp)
		FLAG|=section;
	section<<=1;		
	if (0>par.Vqmp1)
		FLAG|=section;
	section<<=1;	
	if (par.Vqmp1>par.Vmpp)
		FLAG|=section;
	section<<=1;				
	if (par.Impp<par.Iqmp2)
		FLAG|=section;
	section<<=1;	
	if (par.Iqmp2<0)
		FLAG|=section;
	section<<=1;			
	if (par.Vmpp>par.Vqmp2)
		FLAG|=section;
	section<<=1;	
	if (par.Vqmp2>par.Voc)
		FLAG|=section;	
	section<<=1;			
	if (par.Roc<0)
		FLAG|=section;
	section<<=1;	
	return FLAG;
}


int CheckPars(curvepars par, int *FLAGCODE)
{
	int R=1;
	if (!isfinite(par.Isc))
	{
		//IVFLAGERR  FLAG_ISC_NOT_FINITE  "Isc is not finite"
		(*FLAGCODE)|=FLAG_CODE(FLAG_ISC_NOT_FINITE);
		R=0;
	}
	if (!isfinite(par.Voc))
	{
		//IVFLAGERR  FLAG_VOC_NOT_FINITE  "Voc is not finite"
		(*FLAGCODE)|=FLAG_CODE(FLAG_VOC_NOT_FINITE);
		R=0;
	}
	if (!isfinite(par.Vmpp))
	{
		//IVFLAGERR  FLAG_VMPP_NOT_FINITE  "Vmpp is not finite"
		(*FLAGCODE)|=FLAG_CODE(FLAG_VMPP_NOT_FINITE);
		R=0;
	}
	if (!isfinite(par.Impp))
	{
		//IVFLAGERR  FLAG_IMPP_NOT_FINITE  "Impp is not finite"
		(*FLAGCODE)|=FLAG_CODE(FLAG_IMPP_NOT_FINITE);
		R=0;
	}
	
	if (par.Isc<0)
	{
		//IVFLAGERR  FLAG_ISC_NEGATIVE  "Isc is negative"
		(*FLAGCODE)|=FLAG_CODE(FLAG_ISC_NEGATIVE);
		R=0;
	}
	if (par.Gsc<0) // not an error
	{
		//IVFLAGWARN  FLAG_GSC_NEGATIVE  "Gsc is negative"
		(*FLAGCODE)|=FLAG_CODE(FLAG_GSC_NEGATIVE);
		R=1;// problem is second order at best so we warn and continue
	}
	if (par.Voc<0)
	{
		//IVFLAGERR  FLAG_VOC_NEGATIVE  "Voc is negative"
		(*FLAGCODE)|=FLAG_CODE(FLAG_VOC_NEGATIVE);
		R=0;
	}
	if (par.Roc<0) // this is, multivalued I vs V is not normal
	{
		//IVFLAGERR  FLAG_ROC_NEGATIVE  "Roc is negative"
		(*FLAGCODE)|=FLAG_CODE(FLAG_ROC_NEGATIVE);
		R=1;
	}
	if (par.Vmpp<0) 
	{
		//IVFLAGERR  FLAG_VMPP_NEG  "Vmpp is negative"
		(*FLAGCODE)|=FLAG_CODE(FLAG_VMPP_NEG);
		R=0;
	}
	if (par.Vmpp>par.Voc) // multivalued, i.e. error
	{
		//IVFLAGERR  FLAG_VMPP_GT_VOC  "Vmpp is larger than Voc"
		(*FLAGCODE)|=FLAG_CODE(FLAG_VMPP_GT_VOC);
		R=0;
	}
	if (par.Impp<0)
	{
		//IVFLAGERR  FLAG_IMPP_NEG  "Impp is negative"
		(*FLAGCODE)|=FLAG_CODE(FLAG_IMPP_NEG);
		R=0;
	}
	if (par.Impp>par.Isc)
	{
		//IVFLAGWARN  FLAG_IMPP_GT_ISC  "Impp is larger than Isc"
		(*FLAGCODE)|=FLAG_CODE(FLAG_IMPP_GT_ISC);
		R=0;
	}
	if (FlagMono(par))
	{
		//IVFLAGWARN  FLAG_IV_NOT_MONO  "IV parameters are not monotone"
		(*FLAGCODE)|=FLAG_CODE(FLAG_IV_NOT_MONO);
		R=0;
	}
	return R;
}



// compute the distance between two points with normalization factors for the x and y coordinates
double dist2(point p1, point p2, double norm_x, double norm_y)
{
	return (p1.x-p2.x)*(p1.x-p2.x)/(norm_x*norm_x)+(p1.y-p2.y)*(p1.y-p2.y)/(norm_y*norm_y);
}

point Interp(curve C, int *i0, double *l0, double l, double Voc, double Isc)
{
	double d=(*l0), dd;
	int i1=(*i0);
	point p;
	if ((*l0)>=l)
		return C.P[*i0];
	if ((*i0)>=C.N)
		return C.P[C.N-1];
	while (i1<C.N-1)
	{
		dd=sqrt(dist2(C.P[i1], C.P[i1+1], Voc, Isc));
		if (d+dd>l)
		{
			p.x=(C.P[i1].x*(d+dd-l)+C.P[i1+1].x*(l-d))/dd;
			p.y=(C.P[i1].y*(d+dd-l)+C.P[i1+1].y*(l-d))/dd;
			(*l0)=d;
			(*i0)=i1;
			return p;
		}
		d+=dd;
		i1++;
	}
	return C.P[C.N-1];
}

/* resample curve to have somehwat equidistant points in the normaliyed IV space */
curve Resample(curve C, int N, double Voc, double Isc)
{
	curve Cn={0};
	double d, dmean, dtot=0, l0;
	int i, i0;
	
	// integrate the curve length in normalized space
	
	for (i=0;i<C.N-1;i++)
	{
		d=sqrt(dist2(C.P[i], C.P[i+1], Voc, Isc));
		dtot+=d;
	}
	
	// mean distance
	dmean=dtot/(double)(N-1);
	// make a new curve atempting to make a uniform spacing
	i0=0;
	if((Cn.P=Emalloc(N*sizeof(point)))==NULL)
	{
		Cn.N=0;
		return Cn;
	}	
	Cn.N=N;
	// first and last poinst are as before
	Cn.P[0]=C.P[0];	
	Cn.P[N-1]=C.P[C.N-1];	
	d=dmean;
	l0=0;
	for (i=1;i<N-1;i++)
	{
		Cn.P[i]=Interp(C, &i0, &l0, d, Voc, Isc);
		d+=dmean;
	}
	return Cn;
}

/* Auxiliary routine: printing a matrix
void print_matrix(int n, int m, double* a)
{
	int i, j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
			printf("%e ", a[i+j*n]);
        printf( "\n" );
    }
}  */


/* fits a polynomal of order 1-4 to find the Voc
 * The current range is from the Iqmp2 to 0
 * Extrapolation using a polynomal is ill conditiond
 * so to avoid the wild oscillations at the end
 * we do a linear extrapolation from close to the end:
 *    Imin+settings.OC_f_I*Iqmp2
 */ 
#define RCEPS 1e-14
void GetOcpLSQ(curve *Cn, int *FLAGCODE, int order, double Iqmp2, double *Voc, double *Roc, analyze_settings settings)
{
	double *M, *b, *s;
	int N, Nr, i, j, k;
	
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	double wkopt;
	double* work;
	double Imin;
    int *iwork;
    
	// limit order to sensible range
	if (order>4)
		order=4;
	if (order<1)
		order=1;
		
	i=Cn->N-1;
	Imin=Iqmp2;
	while ((i>0)&&(Cn->P[i].y<Iqmp2))
	{
		if (Cn->P[i].y<Imin)
			Imin=Cn->P[i].y;
		i--;
	}
	Nr=(Cn->N-i);
	N=(order+1)*Nr;
	if ((Nr<order+1)||(Imin+settings.OC_f_I*Iqmp2>Iqmp2))
	{
		//IVFLAGERR  FLAG_INSUFF_VOC_RANGE  "Insufficient data to reliably determine Voc"
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_VOC_RANGE);
		(*Voc)=Cn->P[Cn->N-1].x;
		(*Roc)=1e-12;
		return;
	}
	rcond=Nr*RCEPS; // truncate svd
	if((M=Emalloc(N*sizeof(double)))==NULL)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	for (j=order;j>=0;j--)
		for (i=0;i<Nr;i++)
		{
			M[(order-j)*Nr+i]=1.0;
			for (k=0;k<j;k++)
				M[(order-j)*Nr+i]*=Cn->P[Cn->N-Nr+i].y;
		}
	if((b=Emalloc(Nr*sizeof(double)))==NULL)
	{
		free(M);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	for (i=0;i<Nr;i++)
		b[i]=Cn->P[Cn->N-Nr+i].x;
		
	if((s=Emalloc((order+1)*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	n=Nr;
	m=order+1;
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(14*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
    lwork = -1;
    dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
                        iwork, &info );
	lwork = (int)wkopt;
	work = (double*)Emalloc( lwork*sizeof(double) );
	if (work==NULL)
	{
		free(M);
		free(b);
		free(s);
		free(iwork);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}	
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*Voc)=Cn->P[Cn->N-1].x;
		(*Roc)=1e-12;
	}
	else
	{
		/* for a stable extrapolation we use the value and derivative at Imin+f*Iqmp2 */
		double V1=0, V2=0, I, In1, In2, dV;
		I=Imin+settings.OC_f_I*Iqmp2;
		V1=b[order];
		V2=b[order];
		In1=1.0;
		In2=1.0;
		for (i=order-1;i>=0;i--)
		{
			In1*=I;
			In2*=Imin;
			V1+=b[i]*In1;
			V2+=b[i]*In2;
		}
		dV=(V1-V2)/(I-Imin);
		(*Voc)=V1-I*dV;
		(*Roc)=-dV;
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
}
     
/* fits a polynomal of order 1-4 to find the Isc
 * The voltage range is from 0 to the Vqmp1
 * Extrapolation using a polynomal is ill conditiond
 * so to avoid the wild oscillations at the end
 * we do a linear extrapolation from close to the end:
 *    Vmin+settings.SC_f_V*Vqmp1
 */   
void GetScpLSQ(curve *Cn, int *FLAGCODE, int order, double Vqmp1, double *Isc, double *Gsc, analyze_settings settings)
{
	double *M, *b, *s;
	int N, Nr, i, j, k;
	
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	double Vmin;
	double wkopt;
	double* work;
    int *iwork;
    
	// limit order to sensible range
	if (order>4)
		order=4;
	if (order<1)
		order=1;
		
	i=0;
	Vmin=Cn->P[0].x;
	while ((i<Cn->N)&&(Cn->P[i].x<Vqmp1))
		i++;
	Nr=i;
	N=(order+1)*Nr;
	if ((Nr<order+1)||(Vmin+settings.SC_f_V*Vqmp1>Vqmp1))
	{
		//IVFLAGERR  FLAG_INSUFF_ISC_RANGE  "Insufficient data to reliably determine Isc"
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_ISC_RANGE);
		(*Isc)=Cn->P[0].y;
		(*Gsc)=1e12;
		return;
	}
	
	rcond=Nr*RCEPS; // truncate svd
	if((M=Emalloc(N*sizeof(double)))==NULL)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	for (j=order;j>=0;j--)
		for (i=0;i<Nr;i++)
		{
			M[(order-j)*Nr+i]=1;
			for (k=0;k<j;k++)
				M[(order-j)*Nr+i]*=Cn->P[i].x;
		}
	if((b=Emalloc(Nr*sizeof(double)))==NULL)
	{
		free(M);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	for (i=0;i<Nr;i++)
		b[i]=Cn->P[i].y;
		
	if((s=Emalloc((order+1)*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	n=Nr;
	m=order+1;
	//print_matrix(n, m, M);
	lda=n;
	ldb=n;
	if((iwork=Emalloc(14*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		free(iwork);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
    lwork = -1;
    dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
                        iwork, &info );
	lwork = (int)wkopt;
	work = (double*)Emalloc( lwork*sizeof(double) );
	if (!work)
	{
		free(M);
		free(b);
		free(s);
		free(iwork);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*Isc)=Cn->P[0].y;
		(*Gsc)=1e12;
	}
	else
	{
		/* for a stable extrapolation we use the value and derivative at Vmin+f*Vqmp1 */
		double I1=0, I2=0, V, Vn1, Vn2, dI;
		V=Vmin+settings.SC_f_V*Vqmp1;
		I1=b[order];
		I2=b[order];
		Vn1=1.0;
		Vn2=1.0;
		for (i=order-1;i>=0;i--)
		{
			Vn1*=V;
			Vn2*=Vmin;
			I1+=b[i]*Vn1;
			I2+=b[i]*Vn2;
		}
		dI=(I1-I2)/(V-Vmin);
		(*Isc)=I1-V*dI;
		(*Gsc)=-dI;
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
}


/* following routines OC SC extrapolation routines are not in use
 * The idea was to use the cubic hermite spline representation to find 
 * the SC/OC points and slope. Which in itself is simple and elegant
 * but as polynomal extrapolations are ill conditioned it got ugly and 
 * hacky to try to keep that under control.
 */

/* hermite cubic spline base equations */
#define H00(t) (2*(t)*(t)*(t)-3*(t)*(t)+1)
#define H10(t) ((t)*(t)*(t)-2*(t)*(t)+(t))
#define H01(t) (-2*(t)*(t)*(t)+3*(t)*(t))
#define H11(t) ((t)*(t)*(t)-(t)*(t))

/* hermite cubic spline base equations for the derivative */
#define H00dot(t) (6*(t)*(t)-6*(t))
#define H10dot(t) (3*(t)*(t)-4*(t)+1)
#define H01dot(t) (-6*(t)*(t)+6*(t))
#define H11dot(t) (3*(t)*(t)-2*(t))

void GetOcpLSQ_(curve *Cn, int *FLAGCODE, double Vmpp, double Vqmp2, double Iqmp2, double *Voc, double *Roc, analyze_settings settings)
{
	double *M, *b, *s;
	int N, Nr, i;
	
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond;
	double wkopt;
	double* work;
    int *iwork;
    
	i=Cn->N-1;
	
	while ((i>0)&&(Cn->P[i].y<Iqmp2))
		i--;
	Nr=(Cn->N-i);
	N=2*Nr;
	if (Nr<2)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_VOC_RANGE);
		(*Voc)=Cn->P[Cn->N-1].x;
		(*Roc)=1e-12;
		return;
	}
	
	if((M=Emalloc(N*sizeof(double)))==NULL)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	if((b=Emalloc(Nr*sizeof(double)))==NULL)
	{
		free(M);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	
	rcond=Nr*RCEPS; // truncated svd
	// hermite cubic spline based least square solution for Voc/Roc
	// least squares solution to:
	// V(t)-H00(t)*Vqmp2-H10(t)*(Vqmp2-Vmpp)/Iqmp2=H01(t)*Voc+H11(t)*Iqmp2*dV/dI	
	
	for (i=0;i<Nr;i++)
	{
		double t=(Iqmp2-Cn->P[Cn->N-Nr+i].y)/Iqmp2;
		b[i]=Cn->P[Cn->N-Nr+i].x-H00(t)*Vqmp2-H10(t)*(Vqmp2-Vmpp);
		M[i]=H01(t);
		M[i+Nr]=H11(t)*Iqmp2;
	}
		
	if((s=Emalloc(2*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	n=Nr;
	m=2;
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(14*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
    lwork = -1;
    dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
                        iwork, &info );
	lwork = (int)wkopt;
	work = (double*)Emalloc( lwork*sizeof(double) );
	if (work==NULL)
	{
		free(M);
		free(b);
		free(s);
		free(iwork);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*Voc)=Cn->P[Cn->N-1].x;
		(*Roc)=1e-12;
	}
	else
	{
		(*Voc)=b[0];
		(*Roc)=b[1];
		/*
		double t=(Iqmp2-Cn->P[Cn->N-1].y)/Iqmp2;
		double roc1, voc1, roc2, voc2;
		roc1=(Vqmp2*H00dot(t)+(Vqmp2-Vmpp)*H10dot(t)+b[0]*H01dot(t)+b[1]*H11dot(t)*Iqmp2)/Iqmp2;
		voc1=Cn->P[Cn->N-1].x+roc1*Cn->P[Cn->N-1].y;
		voc2=b[0];
		roc2=b[1];
		// here it gets ugly
		if (fabs((voc2-voc1)/(voc1-Cn->P[Cn->N-1].x))>settings.OC_f_I)
			(*Voc)=voc1;
		else
			(*Voc)=voc2;
		if (fabs((roc2-roc1)/roc1)>settings.OC_f_I)
			(*Roc)=roc1;
		else
			(*Roc)=roc2;
		*/
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
}

void GetScpLSQ_(curve *Cn, int *FLAGCODE, double Impp, double Vqmp1, double Iqmp1, double *Isc, double *Gsc, analyze_settings settings)
{
	double *M, *b, *s;
	int N, Nr, i;
	
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond;
	double wkopt;
	double* work;
    int *iwork;
    
    
	i=0;
	
	while ((i<Cn->N)&&(Cn->P[i].x<Vqmp1))
		i++;
	Nr=i;
	N=2*Nr;
	if (Nr<2)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_INSUFF_ISC_RANGE);
		(*Isc)=Cn->P[0].y;
		(*Gsc)=1e12;
		return;
	}
	
	rcond=Nr*RCEPS; // truncated svd
	
	if((M=Emalloc(N*sizeof(double)))==NULL)
	{
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	if((b=Emalloc(Nr*sizeof(double)))==NULL)
	{
		free(M);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	
	// hermite cubic spline based least square solution for Isc/Gsc
	// least squares solution to:
	// I(t)-H01(t)*Iqmp1-H11(t)*(Iqmp-Impp)=H00(t)*Isc+H10(t)*Vqmp1*dI/dV
	for (i=0;i<Nr;i++)
	{
		double t=Cn->P[i].x/Vqmp1;
		b[i]=Cn->P[i].y-H01(t)*Iqmp1+H11(t)*(Iqmp1-Impp);
		M[i]=H00(t);
		M[i+Nr]=-H10(t)*Vqmp1;
	}
		
	if((s=Emalloc(2*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	n=Nr;
	m=2;
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(14*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
    lwork = -1;
    dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
                        iwork, &info );
	lwork = (int)wkopt;
	work = (double*)Emalloc( lwork*sizeof(double) );	
	if (work==NULL)
	{
		free(M);
		free(b);
		free(s);
		free(iwork);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return;
	}
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*Isc)=Cn->P[0].y;
		(*Gsc)=1e12;
	}
	else
	{
		(*Isc)=b[0];
		(*Gsc)=b[1];
		/*
		double t=Cn->P[0].x/Vqmp1;
		double gsc1, isc1, gsc2, isc2;
		gsc1=(b[0]*H00dot(t)+b[1]*H10dot(t)*Vqmp1+Iqmp1*H01dot(t)-(Iqmp1-Impp)*H11dot(t))/Vqmp1;
		isc1=Cn->P[0].y+gsc1*Cn->P[0].x;
		isc2=b[0];
		gsc2=b[1];
		
		// here it gets ugly 
		if (fabs((isc2-isc1)/(isc1-Cn->P[0].y))>settings.SC_f_V)
			(*Isc)=isc1;
		else
			(*Isc)=isc2;
		if (fabs((gsc2-gsc1)/gsc1)>settings.SC_f_V)
			(*Gsc)=gsc1;
		else
			(*Gsc)=gsc2;
		*/
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
}
       

curvepars AnalyzeLSQ(curve C, curve *Cout, int *FLAGCODE, analyze_settings settings)
/* Analyze:
 * Input:
 * 	C		- IV data
 * 	am		- Analize method. either 'r' (default) for robust or 'f' for fragile. Only use fragile when your data is virtually noise free
 *  Cout	- Filtered IV out (optional, pass NULL to not export anything)
 *  FLAGCODE- Error flags
 */  
{
	curvepars par=nullpars;
	curve Cn={0};
	int im, f=0;
	Cn=CheckAndFilterCurve(C, &f, settings);
	if (f)
	{
		(*FLAGCODE)|=f;
		if (Cout)
		{
			Cout->N=Cn.N;
			Cout->P=Cn.P;
		}
		else
			free(Cn.P);
		return par;		
	}
	
	if (settings.Resample)
	{
		// replace Cn with a resampled one
		curve Cr;
		Cr=Resample(Cn,settings.Resample,Cn.P[Cn.N-1].x, Cn.P[0].y);			
		if (Cr.N)
		{
			free(Cn.P);
			Cn=Cr;
		}
		
	}	
	GetMpp(Cn, FLAGCODE, &par.Vmpp, &par.Impp, &im, settings);	
	if ((*FLAGCODE)&FLAG_CODE(FLAG_MALLOCFAIL))
	{
		if (Cout)
		{
			Cout->N=Cn.N;
			Cout->P=Cn.P;
		}
		else
			free(Cn.P);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return par;
	}
	GetQMPPs(Cn, FLAGCODE, par.Vmpp, par.Impp, &par.Vqmp1, &par.Iqmp1, &par.Vqmp2, &par.Iqmp2, settings);
	if ((*FLAGCODE)&FLAG_CODE(FLAG_MALLOCFAIL))	
	{
		if (Cout)
		{
			Cout->N=Cn.N;
			Cout->P=Cn.P;
		}
		else
			free(Cn.P);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return par;
	}	
	//GetScpLSQ_(&Cn, FLAGCODE, par.Impp, par.Vqmp1, par.Iqmp1, &par.Isc, &par.Gsc, settings);
	GetScpLSQ(&Cn, FLAGCODE, 3, par.Vqmp1, &par.Isc, &par.Gsc, settings);
	if ((*FLAGCODE)&FLAG_CODE(FLAG_MALLOCFAIL))	
	{
		if (Cout)
		{
			Cout->N=Cn.N;
			Cout->P=Cn.P;
		}
		else
			free(Cn.P);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return par;
	}
	//GetOcpLSQ_(&Cn, FLAGCODE, par.Vmpp, par.Vqmp2, par.Iqmp2, &par.Voc, &par.Roc, settings);
	GetOcpLSQ(&Cn, FLAGCODE, 3, par.Iqmp2, &par.Voc, &par.Roc, settings);	
	if ((*FLAGCODE)&FLAG_CODE(FLAG_MALLOCFAIL))
	{
		if (Cout)
		{
			Cout->N=Cn.N;
			Cout->P=Cn.P;
		}
		else
			free(Cn.P);
		(*FLAGCODE)|=FLAG_CODE(FLAG_MALLOCFAIL);
		return par;
	}
	
	/* do sanity check */
	CheckPars(par, FLAGCODE);
	
	if (Cout)
	{
		Cout->N=Cn.N;
		Cout->P=Cn.P;
	}
	else
		free(Cn.P);
	return par;
}

// this is broken, this routine checks whether a plain spline interpolation would be monotonous
// we should remove the monotone flags from the analyze routine or only chack actual monotone 
// properties w.r.t. the solar cell parameters 
/*int FlagMono(curvepars par)
{
	int section=1, FLAG=0;
	
	// I allow the sign of the derivative to be wrong
	// such massive errors already trip other flags anyway
	if (par.Gsc*(par.Isc-par.Iqmp1)<0)
		FLAG|=section;	
	if (par.Gsc*((par.Iqmp1-par.Impp)/par.Vqmp1)<0)
		FLAG|=section;		
	section=section<<1;
	
	if ((par.Iqmp1-par.Impp)*((par.Iqmp1-par.Impp)/par.Vqmp1)<0)
		FLAG|=section;	
	if (((par.Iqmp1-par.Impp)/par.Vqmp1)*(par.Impp/par.Vmpp)<0)
		FLAG|=section;		
	section=section<<1;
	
	if ((par.Impp-par.Iqmp2)*(par.Impp/par.Vmpp)<0)
		FLAG|=section;	
	if ((par.Impp/par.Vmpp)*(par.Iqmp2/(par.Vqmp2-par.Vmpp))<0)
		FLAG|=section;		
	section=section<<1;
	
	if (par.Iqmp2<0)
		FLAG|=section;	
	if ((par.Iqmp2/(par.Vqmp2-par.Vmpp))/par.Roc<0)
		FLAG|=section;		
	section=section<<1;
	
	return FLAG;
}*/
