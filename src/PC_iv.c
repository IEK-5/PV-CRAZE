/*
    Pindado-Cubas IV model.
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



typedef struct {
	double Voc, Isc, Vmpp, Impp, eta;
} PCpars;


// convert general to internal parameter struct
PCpars Generic2PCPars(generic_pars P)
{
	PCpars S;
	
	if (P.TYPE!=M_PC)
	{
		fprintf(stderr,"Error: wrong parameter type in Generic2PCPars\n");
		S.Voc=0;
		S.Isc=0;
		S.Vmpp=0;
		S.Impp=0;
		S.eta=0;
		return S;
	}
	S.Voc=P.pars[0];
	S.Isc=P.pars[1];
	S.Vmpp=P.pars[2];
	S.Impp=P.pars[3];
	S.eta=P.pars[4];
	return S;
}



// convert internal to general parameter struct
generic_pars PC2GenericPars(PCpars S)
{
	generic_pars P;
	P.TYPE=M_PC;
	P.pars[0]=S.Voc;
	P.pars[1]=S.Isc;
	P.pars[2]=S.Vmpp;
	P.pars[3]=S.Impp;
	P.pars[4]=S.eta;
	P.pars[5]=0;
	P.pars[6]=0;
	P.pars[7]=0;
	P.pars[8]=0;
	P.pars[9]=0;
	return P;
}	

// Warning: I take the reciprocal value of m as then both gamma and m 
// are between 0 and 1 which makes it easier to do a contraint optimization
// in the generic pars the meaning of m is, however, as it is in literature
#define KHj(v,g,m) (1.0-(1.0-(g))*(v)-(g)*pow((v),1/(m))) 
#define dKHj(v,g,m) ((g)-1.0-(g)*pow((v),(1/(m)-1.0))/(m))
// pars are KHpars
point PCIV_i (double t, point *dp, void *pars)
{
	point p;
	PCpars par;
	
	par=*((PCpars *)pars);
	p.x=t*par.Voc; // for t 0-1 we cover a bit wider voltage range to be sure
	
	if (p.x<=par.Vmpp)
		p.y=par.Isc*(1-(1-par.Impp/par.Isc)*pow(p.x/par.Vmpp,par.Impp/(par.Isc-par.Impp)));
	else
		p.y=(par.Impp*par.Vmpp/p.x)*(1-pow((p.x-par.Vmpp)/(par.Voc-par.Vmpp),par.eta));
	if (dp)
	{
		if (p.x<=par.Vmpp)
			dp->y=-(par.Impp*par.Isc*(1-par.Impp/par.Isc)*pow(p.x/par.Vmpp, par.Impp/(par.Isc-par.Impp)))/((par.Isc-par.Impp)*t);
		else
			dp->y=-(par.Impp*par.Vmpp*par.eta*pow((p.x-par.Vmpp)/(par.Voc-par.Vmpp),par.eta))/(t*(p.x-par.Vmpp))-p.y/t;
		dp->x=par.Voc;
	}
	return p;
}

// Exported Function
curve PCToCurve(generic_pars P, int N, MODELSET *S)
{
	PCpars par;
	par=Generic2PCPars(P);
	return ModelToCurve	((void *) &par, &PCIV_i, N, S);
}

double PCRMSErr_i(PCpars P, curve C, MODELSET *S)
{
	return RMSErr(P.Voc, P.Isc, (void *) &P, &PCIV_i, C, S);
}

// Exported Function
double PCRMSErr(generic_pars P, curve C, MODELSET *S)
{
	return PCRMSErr_i(Generic2PCPars(P), C, S);
}

curve PCDistCurve_i(PCpars P, curve C, MODELSET *S)
{
	return DistCurve(P.Voc, P.Isc, (void *) &P, &PCIV_i, C,S);
}

// Exported Function
curve PCDistCurve(generic_pars P, curve C, MODELSET *S)
{
	return PCDistCurve_i(Generic2PCPars(P), C,S);
}
curvepars PCPars(generic_pars P, MODELSET *S)
{
	PCpars pp;
	pp=Generic2PCPars(P);
	
	return SolPar(pp.Voc, pp.Isc, (void *) &pp, &PCIV_i, S);	
}

// Exported Function
double PCIV(double V, generic_pars P, MODELSET *S)
{	
	PCpars par;
	par=Generic2PCPars(P);
	return Solve_IV(V,par.Voc,NULL,(void *) &par, &PCIV_i, S);	
}

// Exported Function
double PCVI(double I, generic_pars P, MODELSET *S)
{	
	PCpars par;
	par=Generic2PCPars(P);
	return Solve_VI(I,par.Isc,NULL,(void *) &par, &PCIV_i, S);	
}

generic_pars PCFit(curvepars par, curve C, double *Efit, MODELSET *S)
{
	PCpars P;
	
	P.Voc=par.Voc;
	P.Isc=par.Isc;
	P.Vmpp=par.Vmpp;
	P.Impp=par.Impp;
	P.eta=log(par.Vmpp*par.Impp-par.Vqmp2*par.Iqmp2)-log(par.Vmpp*par.Impp);
	P.eta/=(log(par.Vqmp2-par.Vmpp)-log(par.Voc-par.Vmpp));
	(*Efit)=PCRMSErr_i(P, C,S);
	return PC2GenericPars(P);
}

