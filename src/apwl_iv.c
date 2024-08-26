/*
    Advanbced Piece-Wise Linear IV model utilizing also the derivative information.
    Copyright (C) 2024  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2024                                   *
 *                                                               *             
 *****************************************************************/ 
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "emalloc.h"
#include "curve.h"
#include "fitmodels.h"
#include "curvemodel.h"

// internal struct
typedef struct {
	double V[13];
	double I[13];
} apwlpars;

// internal struct, some sort of fractal smoothing
typedef struct {
	double *V, *I;
	int N;
} spwlpars;

// Exported Function
// return curvepars according to the model fit.
curvepars APWLPars(generic_pars P, MODELSET *S)
{
	curvepars par;
// to en able SPWL:
//	if ((P.TYPE!=M_APWL)&&(P.TYPE!=M_SPWL))
	if (P.TYPE!=M_APWL)
	{
		fprintf(stderr,"Error: wrong parameter type in APWLPars\n");
		return nullpars;
	}
	par.Voc		= P.pars[0];
	par.Isc		= P.pars[1];
	par.Vmpp	= P.pars[2];
	par.Impp	= P.pars[3];
	par.Vqmp1	= P.pars[4];
	par.Iqmp1	= P.pars[5];
	par.Vqmp2	= P.pars[6];
	par.Iqmp2	= P.pars[7];
	par.Gsc		= P.pars[8];
	par.Roc		= P.pars[9];
	return par;
}	

// convert general to internal parameter struct
apwlpars Generic2APWLPars(generic_pars P, MODELSET *S)
{
	apwlpars Sd;
	int i;
	double dV, dv, dI[5];

//	if ((P.TYPE!=M_APWL)&&(P.TYPE!=M_SPWL))
	if (P.TYPE!=M_APWL)
	{
		int i;
		fprintf(stderr,"Error: wrong parameter type in Generic2APWLPars\n");
		for (i=0;i<13;i++)
		{
			Sd.V[i]=0;
			Sd.I[i]=0;
		}
		return Sd;
	}
	
	
	Sd.V[0] = 0;
	Sd.V[3] = P.pars[4];
	Sd.V[6] = P.pars[2];
	Sd.V[9] = P.pars[6];
	Sd.V[12]= P.pars[0];
	
	Sd.I[0] =P.pars[1];
	Sd.I[3] =P.pars[5];
	Sd.I[6] =P.pars[3];
	Sd.I[9] =P.pars[7];
	Sd.I[12]=0;
	
	dI[0]=-P.pars[8];
	dI[1]=-(P.pars[5]-P.pars[3])/P.pars[4];
	dI[2]=-P.pars[3]/P.pars[2];
	dI[3]=-P.pars[7]/(P.pars[6]-P.pars[2]);
	dI[4]=-1/P.pars[9];
	
	for (i=0;i<4;i++)
	{
		dV=Sd.V[3*(i+1)]-Sd.V[3*i];
		Sd.V[3*i+1]=Sd.V[3*i]+dV/3;
		Sd.V[3*i+2]=Sd.V[3*(i+1)]-dV/3;
		Sd.I[3*i+1]=Sd.I[3*i]+dI[i]*dV/3;
		Sd.I[3*i+2]=Sd.I[3*(i+1)]-dI[i+1]*dV/3;
		if ((Sd.I[3*i]>Sd.I[3*(i+1)])&&(dI[i]<=0)&&(dI[i+1]<=0))
		{
			// we should be monotone decreasing
			if (Sd.I[3*i+1]<Sd.I[3*i+2])
			{
				//need to fix
				dv=(Sd.I[3*(i+1)]-Sd.I[3*i])/(dI[i]+dI[i+1]);
				if (dv>dV/2)
					dv=0;
				Sd.V[3*i+1]=Sd.V[3*i]+dv;
				Sd.V[3*i+2]=Sd.V[3*(i+1)]-dv;
				Sd.I[3*i+1]=Sd.I[3*i]+dI[i]*dv;
				Sd.I[3*i+2]=Sd.I[3*(i+1)]-dI[i+1]*dv;				
			}
		}
	}
	return Sd;
}
double APWLC(double t, double points[])
{
	return (1-t)*points[0]+t*points[1];
} 

point APWL_i(double t, point *dp, void * pars)
{
	apwlpars P;
	point p={0,0};
	int i=0, j;
	
	P=*((apwlpars *)pars);
	// for t from 0-1 we go through 11 sections, each with a t' from 0 to 1
	t*=12.0;
	for (j=1;j<12;j++)
		i+=(t>(double)j);
	t=t-(double)i;
	p.x=APWLC(t, P.V+i);
	p.y=APWLC(t, P.I+i);
	if (dp)
	{
		dp->x=P.V[i+1]-P.V[i];
		dp->y=P.I[i+1]-P.I[i];
	}
	return p;
} 

spwlpars InitSmooth(apwlpars P)
{
	spwlpars S;
	int i,j;
	S.N=29;
	if((S.V=Emalloc(S.N*sizeof(double)))==NULL)
	{
		S.I=NULL;
		S.N=0;
		return S;
	}
	if((S.I=Emalloc(S.N*sizeof(double)))==NULL)
	{
		free(S.V);
		S.V=NULL;
		S.N=0;
		return S;
	}
	
	// split corners
	j=0;
	for (i=0;i<13;i++)
	{		
		if (i%3) // not a knot, will split corners only
		{
			double dV, dI;
			// new corner
			S.V[j]=(P.V[i]+P.V[i-1])/2;
			S.I[j]=(P.I[i]+P.I[i-1])/2;
			// new corner
			S.V[j+2]=(2*P.V[i]+P.V[i+1])/3;
			S.I[j+2]=(2*P.I[i]+P.I[i+1])/3;
			// new knot
			dV=S.V[j+2]-S.V[j];
			dI=S.I[j+2]-S.I[j];
			S.V[j+1]=P.V[i];
			S.I[j+1]=S.I[j]+(P.V[i]-S.V[j])*dI/dV;
			j+=3;
		}
		else
		{
			S.V[j]=P.V[i];
			S.I[j]=P.I[i];
			j+=1;
		}
		
	}
	return S;	
}
void FreeSmooth(spwlpars *S)
{
	if (S->V)
		free(S->V);
	if (S->I)
		free(S->I);
	S->N=0;
	S->V=NULL;
	S->I=NULL;
}
		
point SPWL_i(double t, point *dp, void * pars)
{
	spwlpars P;
	point p={0,0};
	int i=0, j;
	
	P=*((spwlpars *)pars);
	t*=(double)P.N-1;
	for (j=1;j<P.N-1;j++)
		i+=(t>(double)j);
	t=t-(double)i;
	p.x=APWLC(t, P.V+i);
	p.y=APWLC(t, P.I+i);
	if (dp)
	{
		dp->x=P.V[i+1]-P.V[i];
		dp->y=P.I[i+1]-P.I[i];
	}
	return p;
} 


// Exported Function
curve APWLToCurve(generic_pars P, int N, MODELSET *S)
{
	apwlpars par;
	curve C;
	par=Generic2APWLPars(P,S);
	C=ModelToCurve((void *) &par, &APWL_i, N,S);
	return C;
}

// Exported Function
curve SPWLToCurve(generic_pars P, int N, MODELSET *S)
{
	apwlpars par;
	spwlpars spar;
	curve C={0};
	par=Generic2APWLPars(P,S);
	spar=InitSmooth(par);
	if (spar.N!=0)
	{
		C=ModelToCurve((void *) &spar, &SPWL_i, N,S);
		FreeSmooth(&spar);
	}
	return C;
}
double APWLRMSErr_i(apwlpars P, curve C, MODELSET *S)
{
	return RMSErr(P.V[12], P.I[0], (void *) &P, &APWL_i, C,S);
}
// Exported Function
double APWLRMSErr(generic_pars P, curve C, MODELSET *S)
{
	apwlpars par;
	double E;
	par=Generic2APWLPars(P,S);
	E=APWLRMSErr_i(par, C,S);
	return E;
}

double SPWLRMSErr_i(spwlpars P, curve C, MODELSET *S)
{
	if ((P.V!=NULL)&&(P.I!=NULL))
		return RMSErr(P.V[P.N-1], P.I[0], (void *) &P, &SPWL_i, C,S);
	else
		return NAN;
}
// Exported Function
double SPWLRMSErr(generic_pars P, curve C, MODELSET *S)
{
	apwlpars par;
	spwlpars spar;
	double E;
	par=Generic2APWLPars(P,S);
	spar=InitSmooth(par);
	E=SPWLRMSErr_i(spar, C,S);
	FreeSmooth(&spar);
	return E;
}

curve APWLDistCurve_i(apwlpars P, curve C, MODELSET *S)
{
	return DistCurve(P.V[12], P.I[0], (void *) &P, &APWL_i, C,S);
}
// Exported Function
curve APWLDistCurve(generic_pars P, curve C, MODELSET *S)
{
	apwlpars par;
	curve d;
	par=Generic2APWLPars(P,S);
	d=APWLDistCurve_i(par, C,S);
	return d;
}

curve SPWLDistCurve_i(spwlpars P, curve C, MODELSET *S)
{
	return DistCurve(P.V[P.N-1], P.I[0], (void *) &P, &SPWL_i, C,S);
}
// Exported Function
curve SPWLDistCurve(generic_pars P, curve C, MODELSET *S)
{
	apwlpars par;
	spwlpars spar;
	curve d;
	par=Generic2APWLPars(P,S);
	spar=InitSmooth(par);
	d=SPWLDistCurve_i(spar, C,S);
	FreeSmooth(&spar);
	return d;
}

// Exported Function
double APWLIV(double V, generic_pars P, MODELSET *S)
{
	apwlpars par;
	par=Generic2APWLPars(P,S);
	return Solve_IV(V,P.pars[0],NULL,(void *) &par, &APWL_i, S);
}
// Exported Function
double APWLVI(double I, generic_pars P, MODELSET *S)
{
	apwlpars par;
	par=Generic2APWLPars(P,S);
	return Solve_IV(I,P.pars[1],NULL,(void *) &par, &APWL_i, S);
}
// Exported Function
double SPWLIV(double V, generic_pars P, MODELSET *S)
{
	apwlpars par;
	spwlpars spar;
	double I;
	par=Generic2APWLPars(P,S);
	spar=InitSmooth(par);
	I=Solve_IV(V,P.pars[0],NULL,(void *) &spar, &SPWL_i, S);
	FreeSmooth(&spar);
	return I;
}
// Exported Function
double SPWLVI(double I, generic_pars P, MODELSET *S)
{
	apwlpars par;
	spwlpars spar;
	double V;
	par=Generic2APWLPars(P,S);
	spar=InitSmooth(par);
	V=Solve_IV(I,P.pars[1],NULL,(void *) &spar, &SPWL_i, S);
	FreeSmooth(&spar);
	return V;
}
// Exported Function
generic_pars APWLFit(curvepars par, curve C, double *Efit, MODELSET *S)
{	
	apwlpars Sd;
	generic_pars P;
	
	P.TYPE=M_APWL;
	P.pars[0]=par.Voc;
	P.pars[1]=par.Isc;
	P.pars[2]=par.Vmpp;
	P.pars[3]=par.Impp;
	P.pars[4]=par.Vqmp1;
	P.pars[5]=par.Iqmp1;
	P.pars[6]=par.Vqmp2;
	P.pars[7]=par.Iqmp2;
	P.pars[8]=par.Gsc;
	P.pars[9]=par.Roc;
	
	
	Sd=Generic2APWLPars(P,S);
	(*Efit)=APWLRMSErr_i(Sd, C,S);
	return P;
}

// Exported Function
generic_pars SPWLFit(curvepars par, curve C, double *Efit, MODELSET *S)
{	
	apwlpars Sd;
	spwlpars spar;
	generic_pars P;
	
	P.TYPE=M_APWL;
	P.pars[0]=par.Voc;
	P.pars[1]=par.Isc;
	P.pars[2]=par.Vmpp;
	P.pars[3]=par.Impp;
	P.pars[4]=par.Vqmp1;
	P.pars[5]=par.Iqmp1;
	P.pars[6]=par.Vqmp2;
	P.pars[7]=par.Iqmp2;
	P.pars[8]=par.Gsc;
	P.pars[9]=par.Roc;
	
	
	Sd=Generic2APWLPars(P,S);
	spar=InitSmooth(Sd);
	(*Efit)=SPWLRMSErr_i(spar, C,S);
	FreeSmooth(&spar);
	return P;
}
