/*
	Routines to aid building new curve models for PV-CRAZE.

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
#include "curve.h"
#include "fitmodels.h"


/* Compute curves and curveparameters for a model
 * a model must provide routines with the following template:
 * 
 * point <model>(double var, point *dp, void *pars)
 * 		return:	point x,y for voltage/current
 * 		var: 	variable on the interval [0,1]. the end points must be 
 * 				open or short circuit
 * 		dp:		derivative point to var.
 * 		pars	anonymous pointer to parameter struct 
 * 
 */ 	

// most routines take Isc and Voc as normalization constants to compute errors and such.


#define POINTDIST2(p1,p2) ((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y))
// model prototype of the form
// point=model(double var, point dp, void *pars)
double CurvePointDist(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), point target, MODELSET *S)
{
	double vmin=0, vmax=1, var, dvar, nvar;
	point p, dp;
	int i=0;
	var=(vmin+vmax)/2;
	// normalize
	target.x/=Voc;
	target.y/=Isc;
	while ((vmax-vmin>S->meps)&&(i<S->miter))
	{
		// compute the tangent line and find the voltage for which the tangent line is closest to the point
		// Use this voltage as our next estimate 
		p=model(var,&dp,pars);
		// normalize
		p.x/=Voc;
		p.y/=Isc;
		dp.x/=Voc;
		dp.y/=Isc;
		// newton
		dvar=(dp.y*target.y + dp.x*target.x - dp.y*p.y - dp.x*p.x)/(dp.y*dp.y+dp.x*dp.x);	
		nvar=var+dvar;
		if (isfinite(nvar))
		{
			if (i%2) // in case we always find the same v we are at the minimum, make sure we bracket the minimum
			{
				if (nvar>var)
					vmin=var;
				else
					vmax=var;
			}
			else
			{
				if (nvar<var)
					vmax=var;
				else
					vmin=var;
			}
			var=nvar;
			if ((var<vmin)||(var>vmax)) // keep v within brackets
				var=(vmin+vmax)/2;
		}	
		else
			var=(vmin+vmax)/2;	
		i++;	
	}
	p=model(var,NULL,pars);
	p.x/=Voc;
	p.y/=Isc;
	return POINTDIST2(target,p);
}

// RMS distance of a curve to a curve model
double RMSErr(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), curve C, MODELSET *S)
{
	int i;
	double d=0;
	if (C.N>0)
	{
		for (i=0;i<C.N;i++)
		{
			d+=CurvePointDist(Voc, Isc, pars, model, C.P[i], S);
		}
		return sqrt(d/(double)C.N);
	}
	return 0.0;
}


// RMS distance of a curve resoved by voltage
curve DistCurve(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), curve C, MODELSET *S)
{
	int i;
	curve D;
	if (C.N>0)
	{
		D.N=C.N;
		if((D.P=Emalloc(D.N*sizeof(point)))==NULL)
		{
			D.N=0;
			return D;
		}
		for (i=0;i<C.N;i++)
		{
			D.P[i].x=C.P[i].x;
			D.P[i].y=CurvePointDist(Voc, Isc, pars, model, C.P[i], S);
		}
	}
	else
	{
		D.N=0;
		D.P=NULL;
	}
	return D;
}
// solve the current for a specified voltage
double Solve_t_IV(double V, double Voc, void * pars, point (*model)(double, point *, void *), MODELSET *S)
{
	double vmin=0, vmax=1,var, dvar, nvar;
	point p, dp;
	int i=0;
	var=V/Voc;
	if (var<0)
		var=0.5;
	if (var>1)
		var=0.5;
	
	p.x=V+Voc;
	while (((vmax-vmin>S->meps)||(fabs(p.x-V)/Voc<S->meps))&&(i<S->miter))
	{
		p=model(var,&dp,pars);
		// newton-like optimization
		dvar=(V - p.x)/dp.x;	
		nvar=var+dvar;
		if (isfinite(nvar))
		{
			if (i%2) // in case we always find the same v we are at the minimum, make sure we bracket the minimum
			{
				if (nvar>var)
					vmin=var;
				else
					vmax=var;
			}
			else
			{
				if (nvar<var)
					vmax=var;
				else
					vmin=var;
			}
			var=nvar;
			if ((var<vmin)||(var>vmax)) // keep v within brackets
				var=(vmin+vmax)/2;
		}	
		else
			var=(vmin+vmax)/2;	
		i++;	
	}
	return var;
}

// solve the voltage for a specified current
double Solve_t_VI(double I, double Isc,  void * pars, point (*model)(double, point *, void *), MODELSET *S)
{
	double vmin=0, vmax=1,var, dvar, nvar;
	point p, dp;
	int i=0;
	var=1-I/Isc;
	if (var<0)
		var=0.5;
	if (var>1)
		var=0.5;
	
	p.y=I+Isc;
	while (((vmax-vmin>S->meps)||(fabs(p.y-I)/Isc<S->meps))&&(i<S->miter))
	{
		p=model(var,&dp,pars);
		// newton
		dvar=(I - p.y)/dp.y;	
		nvar=var+dvar;
		if (isfinite(nvar))
		{
			if (i%2) // in case we always find the same v we are at the minimum, make sure we bracket the minimum
			{
				if (nvar>var)
					vmin=var;
				else
					vmax=var;
			}
			else
			{
				if (nvar<var)
					vmax=var;
				else
					vmin=var;
			}
			var=nvar;
			if ((var<vmin)||(var>vmax)) // keep v within brackets
				var=(vmin+vmax)/2;
		}	
		else
			var=(vmin+vmax)/2;	
		i++;	
	}
	return var;
}


// solve the current for a specified voltage
double Solve_IV(double V, double Voc, point *dp, void * pars, point (*model)(double, point *, void *), MODELSET *S)
{
	point p;
	p=model(Solve_t_IV(V,  Voc, pars, model, S), dp, pars);
	return p.y;
}

// solve the voltage for a specified current
double Solve_VI(double I, double Isc,  point *dp, void * pars, point (*model)(double, point *, void *), MODELSET *S)
{
	point p;
	p=model(Solve_t_VI(I,  Isc, pars, model, S), dp, pars);
	return p.x;
}


point FindMp(point off, void * pars, point (*model)(double, point *, void *), MODELSET *S)
{
	double vmin=0, vmax=1, v, vn=0.5;
	point p, dp;
	int iter=0;
	// offset voltage and current can be used to determine quasi mpp's
	do
	{
		if ((vn>vmax)||(vn<vmin))
			v=(vmin+vmax)/2;
		else
			v=vn;
		p=model(v,&dp,pars);	
		vn=v-(dp.x*(p.y-off.y)+dp.y*(p.x-off.x))/(2*dp.x*dp.y); // go to maximum power point of the tangent line	
		//if (-(p.y-off.y)/(p.x-off.x)<dp.y/dp.x)
		if (vn>v)
			vmin=v;
		else
			vmax=v;
		iter++;
	} while ((fabs(v-vn)>S->meps)&&(iter<S->miter));
	
	if ((vn>vmax)||(vn<vmin))
		v=(vmin+vmax)/2;
	else
		v=vn;
	return model(v,NULL,pars);	
}

curvepars SolPar(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), MODELSET *S)
{
	curvepars par;
	point p, off;
	par.Isc=Solve_IV(0, Voc, &p, pars, model, S);
	par.Gsc=-p.y/p.x;
	par.Voc=Solve_VI(0, Isc, &p, pars, model, S);
	par.Roc=-p.x/p.y;
	off.x=0;
	off.y=0;
	p=FindMp(off, pars, model, S);
	par.Vmpp=p.x;
	par.Impp=p.y;
	
	off.x=0;
	off.y=par.Impp;
	p=FindMp(off, pars, model, S);
	par.Vqmp1=p.x;
	par.Iqmp1=p.y;
	
	off.x=par.Vmpp;
	off.y=0;
	p=FindMp(off, pars, model, S);
	par.Vqmp2=p.x;
	par.Iqmp2=p.y;
	return par;
}


// would be nice to make a good mesh, i.e. more points where the second derivative is high
curve ModelToCurve(void * pars, point (*model)(double, point *, void *), int N, MODELSET *S)
{
	double vmin=0, vmax=1, v, dv, Voc, Isc;
	point p;
	curve C;
	curvepars par;
	double Vp[3], Ip[3];
	int i, j, k;
	if (N<=1)
	{
		C.P=NULL;
		C.N=0;
		return C;
	}
	N=N+3; // make space for the mpp and qmpp's
	if((C.P=Emalloc(N*sizeof(point)))==NULL)
	{
		C.N=0;
		return C;
	}
	C.N=N;
	// determine start and end of the parametric variable
	vmax=Solve_t_VI(0,1,pars,model, S);
	vmin=Solve_t_IV(0,1,pars,model, S);
	
	// I want to have the model solar cell parameters always included
	p=model(vmin,NULL,pars);
	Isc=p.y;
	p=model(vmax,NULL,pars);
	Voc=p.x;
	
	par=SolPar(Voc, Isc, pars, model, S);
	Vp[0]=par.Vqmp1;
	Vp[1]=par.Vmpp;
	Vp[2]=par.Vqmp2;
	Ip[0]=par.Iqmp1;
	Ip[1]=par.Impp;
	Ip[2]=par.Iqmp2;
	j=0;
	v=vmin;
	dv=(vmax-vmin)/(double)(N-1);
	p=model(v,NULL,pars);
	C.P[0].x=p.x;
	C.P[0].y=p.y;
	v+=dv;
	k=1;
	for (i=1;i<N;i++)
	{
		p=model(v,NULL,pars);
		if (j<3)
			if((C.P[i-1].x<Vp[j])&&(p.x>Vp[j]))
			{
				// insert point
				C.P[i].x=Vp[j];
				C.P[i].y=Ip[j];
				j++;
			}
		C.P[k].x=p.x;
		C.P[k].y=p.y;
		k++;
		v+=dv;
	}
	if (k>C.N) // oops, cannot be right, not all mpp/qmpps in there?
		C.N=k;
	return C;
}
