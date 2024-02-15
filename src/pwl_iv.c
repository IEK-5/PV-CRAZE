/*
    Piece-Wise Linear IV model.
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
#include "curve.h"
#include "fitmodels.h"
#include "curvemodel.h"

// Exported Function
// return curvepars according to the model fit.
curvepars PWLPars(generic_pars P, MODELSET *S)
{
	curvepars par;
	
	if (P.TYPE!=M_PWL)
	{
		fprintf(stderr,"Error: wrong parameter type in PWLPars\n");
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

point PWL_coord(double t, int i, generic_pars P)
{
	double V0, I0, V1, I1;
	point p={0,0};
	switch(i)
	{
		case 0:
			if (t>0)
			{
				V0=0;
				I0=P.pars[1];
				V1=P.pars[4];
				I1=P.pars[5];
				p.x=(1-t)*V0+t*V1;
				p.y=(1-t)*I0+t*I1;
			}
			else
			{
				I0=P.pars[1];
				p.x=t;
				p.y=I0-t*P.pars[8];
			}
			return p;
		case 1:
			V0=P.pars[4];
			I0=P.pars[5];
			V1=P.pars[2];
			I1=P.pars[3];
			p.x=(1-t)*V0+t*V1;
			p.y=(1-t)*I0+t*I1;
			return p;
		case 2:
			V0=P.pars[2];
			I0=P.pars[3];
			V1=P.pars[6];
			I1=P.pars[7];
			p.x=(1-t)*V0+t*V1;
			p.y=(1-t)*I0+t*I1;
			return p;
		case 3:
			if (t<=1)
			{
				V0=P.pars[6];
				I0=P.pars[7];
				V1=P.pars[0];
				I1=0;
				p.x=(1-t)*V0+t*V1;
				p.y=(1-t)*I0+t*I1;
			}
			else
			{
				V1=P.pars[0];
				p.x=V1+(t-1)*P.pars[9];
				p.y=-(t-1);
			}
			return p;
		default:
			return p;
	}
} 

point PWL_deriv(double t, int i, generic_pars P)
{
	double V0, I0, V1, I1;
	point p={1,-1};
	switch(i)
	{
		case 0:
			if (t>0)
			{
				V0=0;
				I0=P.pars[1];
				V1=P.pars[4];
				I1=P.pars[5];
				p.x=V1-V0;
				p.y=I1-I0;
			}
			else
			{
				p.y=-P.pars[8];
				p.x=1;
			}
			return p;
		case 1:
			V0=P.pars[4];
			I0=P.pars[5];
			V1=P.pars[2];
			I1=P.pars[3];
			p.x=V1-V0;
			p.y=I1-I0;
			return p;
		case 2:
			V0=P.pars[2];
			I0=P.pars[3];
			V1=P.pars[6];
			I1=P.pars[7];
			p.x=V1-V0;
			p.y=I1-I0;
			return p;
		case 3:
			if (t<1)
			{
				V0=P.pars[6];
				I0=P.pars[7];
				V1=P.pars[0];
				I1=0;
			}
			else
			{
				p.y=-1;
				p.x=P.pars[9];
			}
			return p;
		default:
			return p;
	}
} 

// pars are generic_pars
point PWL_i(double t, point *dp, void * pars)
{
	int i;
	point p;
	generic_pars par;
	par=*((generic_pars *)pars);
	// for t from 0-1 we go through 4 sections, each with a t' from 0 to 1
	t*=4.0;
	i=0;
	i+=(t>1.0);
	i+=(t>2.0);
	i+=(t>3.0);
	// i now points to the knot index before the interval
	t=t-(double)i;
	// compute the time for this interval
	
	p=PWL_coord(t, i, par);
	
	if (dp)
		(*dp)=PWL_deriv(t,i, par);
	return p;
}

// Exported Function
curve PWLToCurve(generic_pars P, int N, MODELSET *S)
{
	curve C;
	C=ModelToCurve((void *) &P, &PWL_i, N,S);
	return C;
}

// Exported Function
double PWLRMSErr(generic_pars P, curve C, MODELSET *S)
{
	return RMSErr(P.pars[0], P.pars[1], (void *) &P, &PWL_i, C,S);
}

// Exported Function
curve PWLDistCurve(generic_pars P, curve C, MODELSET *S)
{
	return DistCurve(P.pars[0], P.pars[1], (void *) &P, &PWL_i, C,S);
}

// Exported Function
generic_pars PWLFit(curvepars par, curve C, double *Efit, MODELSET *S)
{	
	generic_pars P;
	
	P.TYPE=M_PWL;
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
	
	
	(*Efit)=PWLRMSErr(P, C,S);
	return P;
}
