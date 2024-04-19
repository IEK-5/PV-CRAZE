
/*
    one diode IV model.
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
#include <time.h>
// use gsl nonlinear lsq
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "emalloc.h"
#include "curve.h"
#include "fitmodels.h"
#include "curvemodel.h"
#include "neldermead.h"
#include "diffevol.h"
#include "quickselect.h"

typedef struct {
	double I0,nkT,Iph,Rs, Rsh, Voc, Isc;
} diodepars;

#define Dj(v,io,nkt,rsh) ((io)*(exp((v)/(nkt))-1)+(v)/(rsh))
#define dDj(v,io,nkt,rsh) ((io)*exp((v)/(nkt))/(nkt)+1/rsh)
#define ddDj(v,io,nkt) ((io)*exp((v)/(nkt))/(nkt)/(nkt))
#define NEWTON_EPS 1e-6


/* in principle Halley's method should be faster
 * than newton. However, I find nothing significant 
 * It seems Halleys is slightly less robust so I stick 
 * with newton for now
 */
//#define HALLEY

#ifdef HALLEY
#define Diode_I(V,dIo,P,S) (Diode_I_halley((V),(dIo),(P),(S)))
#define SolveVoc(P,S) (SolveVoc_halley((P),(S)))
#else
#define Diode_I(V,dIo,P,S) (Diode_I_newton((V),(dIo),(P),(S)))
#define SolveVoc(P,S) (SolveVoc_newton((P),(S)))
#endif


// Solve the diode current using newton
double Diode_I_newton(double V, double *dIo, diodepars P, MODELSET *S)
{
	double v1, I, dE;
	double E;
	int i=0;
	static double Vinit=-1; // keep track of previous solutions to aid vonvergence
	E=1;
	if (Vinit<0)
		v1=V;
	else
		v1=Vinit;
	
	while ((fabs(E)>S->meps)&&(i<S->miter))
	{
		
		I=Dj(v1,P.I0,P.nkT,P.Rsh)-P.Iph;
		dE=dDj(v1,P.I0,P.nkT,P.Rsh)+1/P.Rs;
		E=(I-(V-v1)/P.Rs);
		v1-=E/dE; /* newton-rapson */
		i++;
	}
	I=Dj(v1,P.I0,P.nkT,P.Rsh)-P.Iph;
	if (dIo)
		(*dIo)=-1/(1/dDj(v1,P.I0,P.nkT,P.Rsh)+P.Rs);
	if (fabs(E)<S->meps)
		Vinit=v1;
	else
		Vinit=-1;
	return -I;
}
	
// newton solver for the Voc
// enforce boundaries for stability
double SolveVoc_newton(diodepars P, MODELSET *S)
{
	double V, I, dI;
	double Vmin=0, Vmax;
	int i=0;
	Vmax=P.nkT*log(1+P.Iph/P.I0);
	V=(Vmin+Vmax)/2;
	do
	{
		I=Dj(V,P.I0,P.nkT,P.Rsh)-P.Iph;
		dI=dDj(V,P.I0,P.nkT,P.Rsh);
		if (I>0)
			Vmax=V;
		else
			Vmin=V;
		V-=I/dI;
		if ((V<Vmin)||(V>Vmax))
			V=(Vmin+Vmax)/2;
		i++;
	} while ((fabs(I/P.Isc)>S->meps)&&(i<S->miter));
	return V;
}	
// Solve the diode current using Halley's method (should be faster, did not test it yet)
double Diode_I_halley(double V, double *dIo, diodepars P, MODELSET *S)
{
	double v1, I, dE, ddE;
	double E;
	int i=0;
	static double Vinit=-1; // keep track of previous solutions to aid vonvergence
	E=1;
	if (Vinit<0)
		v1=V;
	else
		v1=Vinit;
	
	while ((fabs(E)>S->meps)&&(i<S->miter))
	{
		
		I=Dj(v1,P.I0,P.nkT,P.Rsh)-P.Iph;
		dE=dDj(v1,P.I0,P.nkT,P.Rsh)+1/P.Rs;
		ddE=ddDj(v1,P.I0,P.nkT);
		E=(I-(V-v1)/P.Rs);
		v1-=2*E*dE/(2*dE*dE-E*ddE); 
		i++;
	}
	I=Dj(v1,P.I0,P.nkT,P.Rsh)-P.Iph;
	if (dIo)
		(*dIo)=-1/(1/dDj(v1,P.I0,P.nkT,P.Rsh)+P.Rs);
	if (fabs(E)<S->meps)
		Vinit=v1;
	else
		Vinit=-1;
	return -I;
}
// halley solver for the Voc
// enforce boundaries for stability
double SolveVoc_halley(diodepars P, MODELSET *S)
{
	double V, I, dI, ddI;
	double Vmin=0, Vmax;
	int i=0;
	Vmax=P.nkT*log(1+P.Iph/P.I0);
	V=(Vmin+Vmax)/2;
	do
	{
		I=Dj(V,P.I0,P.nkT,P.Rsh)-P.Iph;
		dI=dDj(V,P.I0,P.nkT,P.Rsh);
		ddI=ddDj(V,P.I0,P.nkT);
		if (I>0)
			Vmax=V;
		else
			Vmin=V;
		V-=2*I*dI/(2*dI*dI-I*ddI);
		if ((V<Vmin)||(V>Vmax))
			V=(Vmin+Vmax)/2;
		i++;
	} while ((fabs(I/P.Isc)>S->meps)&&(i<S->miter));
	return V;
}	

	

// convert general to internal parameter struct
diodepars Generic2DiodePars(generic_pars P, MODELSET *S)
{
	diodepars Sd;
	
	if (P.TYPE!=M_DIODE)
	{
		fprintf(stderr,"Error: wrong parameter type in Generic2DiodePars\n");
		Sd.I0=0;
		Sd.nkT=1;
		Sd.Iph=0;
		Sd.Rs=0;
		Sd.Rsh=1;
		Sd.Voc=0;
		Sd.Isc=0;
		return Sd;
	}
	Sd.I0=P.pars[0];
	Sd.nkT=P.pars[1];
	Sd.Iph=P.pars[2];
	Sd.Rs=P.pars[3];
	Sd.Rsh=P.pars[4];
	// solve the special points
	Sd.Isc=Diode_I(0, NULL, Sd, S);
	Sd.Voc=SolveVoc(Sd, S);
	return Sd;
}


// convert internal to general parameter struct
generic_pars Diode2GenericPars(diodepars S)
{
	generic_pars P;
	P.TYPE=M_DIODE;
	P.pars[0]=S.I0;
	P.pars[1]=S.nkT;
	P.pars[2]=S.Iph;
	P.pars[3]=S.Rs;
	P.pars[4]=S.Rsh;
	P.pars[5]=0;
	P.pars[6]=0;
	P.pars[7]=0;
	P.pars[8]=0;
	P.pars[9]=0;
	return P;
}	


// pars are diodepars
point diodeIV_i (double t, point *dp, void *pars)
{
	point p;
	double *di=NULL;
	diodepars par;
	
	if (dp)
		di=&(dp->y);
	par=*((diodepars *)pars);
	p.x=(2*t-0.5)*par.Voc; // for t 0-1 we cover a bit wider voltage range to be sure
	p.y=Diode_I(p.x, di, par,&DefaultMODSET);
	if (dp)
	{
		// take care of derivatives we have to compute dx/dt and dy/dt
		dp->x=2*par.Voc;
		dp->y*=2*par.Voc;
	}
	return p;
}


// Exported Function
curve DiodeToCurve(generic_pars P, int N, MODELSET *S)
{
	diodepars par;
	par=Generic2DiodePars(P, S);
	return ModelToCurve	((void *) &par, &diodeIV_i, N, S);
}

double DiodeRMSErr_i(diodepars P, curve C, MODELSET *S)
{
	return RMSErr(P.Voc, P.Isc, (void *) &P, &diodeIV_i, C,S);
}

// Exported Function
double DiodeRMSErr(generic_pars P, curve C, MODELSET *S)
{
	return DiodeRMSErr_i(Generic2DiodePars(P, S), C,S);
}


curve DiodeDistCurve_i(diodepars P, curve C, MODELSET *S)
{
	return DistCurve(P.Voc, P.Isc, (void *) &P, &diodeIV_i, C,S);
}

// Exported Function
curve DiodeDistCurve(generic_pars P, curve C, MODELSET *S)
{
	return DiodeDistCurve_i(Generic2DiodePars(P, S), C,S);
}

// Exported Function
curvepars DiodePars(generic_pars P, MODELSET *S)
{
	diodepars pp;
	pp=Generic2DiodePars(P, S);
	
	return SolPar(pp.Voc, pp.Isc, (void *) &pp, &diodeIV_i, S);	
}

// Exported Function
double DiodeIV(double V, generic_pars P, MODELSET *S)
{	
	diodepars par;
	par=Generic2DiodePars(P,S);
	return Solve_IV(V,par.Voc,NULL,(void *) &par, &diodeIV_i, S);	
}

// Exported Function
double DiodeVI(double I, generic_pars P, MODELSET *S)
{	
	diodepars par;
	par=Generic2DiodePars(P,S);
	return Solve_VI(I,par.Isc,NULL,(void *) &par, &diodeIV_i, S);	
}


// function to select support points in the oblicque asymptote parametrization routines
void Select(curvepars P, int i, double *V, double *I)
{
	switch (i)
	{
		case 0:
			(*V)=0;
			(*I)=P.Isc;
			break;
		case 1:
			(*V)=P.Vqmp1;
			(*I)=P.Iqmp1;
			break;
		case 2:
			(*V)=P.Vmpp;
			(*I)=P.Impp;
			break;
		case 3:
			(*V)=P.Vqmp2;
			(*I)=P.Iqmp2;
			break;
		case 4:
			(*V)=P.Voc;
			(*I)=0;
			break;
		default:
			(*V)=0;
			(*I)=0;
	}	
}

/* the oblique asymptote method of parameter estimation comes from:
 * F.J. Toledo, Jose M. Blanes, "Geometric properties of the 
 * single-diode photovoltaic model and a new very simple method 
 * for parameters extraction", Renewable Energy 72 (2014) 125-133
 */
diodepars InitParsOA(curvepars par, int *err, int i1, int i2, int i3)
{
	/* oblique asymptote approx */
	diodepars P;
	double V1, I1, V2, I2, V3, I3;
	double A,B,E, lD, lC, Isc, Gsc, f1, f2, f3, d;
	(*err)=0;
	// some basic checks whether a diode model could fit
	if ((par.Voc<1e-6)||(par.Isc<1e-6))
	{
		P.Rs  = 0;
		P.Rsh = 0;
		P.Iph = 0;	
		P.Voc = 0;
		P.Isc = 0;
		P.nkT = 0;
		P.I0  = 0;	
		(*err)=1;
		return P;
	}
	
	if (par.Voc<par.Vmpp)
	{
		// this is a mess 
		d=par.Vmpp;
		par.Vmpp=par.Voc;
		par.Voc=d;
	}
	if (par.Isc<par.Impp)
	{
		// this is a mess 
		d=par.Impp;
		par.Impp=par.Isc;
		par.Isc=d;
	}	
	
	if ((par.Isc-par.Impp)/par.Gsc<par.Vmpp)
	{
		// current  s-shape 
		par.Gsc=(par.Isc-par.Impp)/2/par.Vmpp;
	}
	
	
	if (par.Impp*par.Vmpp/par.Isc/par.Voc<=0.25)
	{
		// return a photovoltaic resistor
		P.Rs  = 0;
		P.Rsh = par.Voc/par.Isc;
		P.Iph = par.Isc;	
		P.Voc = par.Voc;
		P.Isc = par.Isc;
		P.nkT = 1;
		P.I0  = 0;
		(*err)=1;
		return P;
	}
	Isc=par.Isc;
	Gsc=par.Gsc;
	/* here we choose which points to use. Possible combinations are:
	 * 1,2,3 (qmp1, mpp,  qmp2)
	 * 1,2,4 (qmp1, mpp,  oc)
	 * 1,3,4 (qmp1, qmp2, oc)
	 * 2,3,4 (mpp,  qmp2, oc)
	 */
	
	Select(par, i1, &V1, &I1);
	Select(par, i2, &V2, &I2);
	Select(par, i3, &V3, &I3);
	f1=(Isc-V1*Gsc-I1);
	if (f1<0)
	{
		Isc-=f1;
		Gsc=(Isc-I1)/V1; // correct Gsc such that f1=0
	}
	f2=(Isc-V2*Gsc-I2);
	if (f2<0)
	{
		Isc-=f1;
		Gsc=(Isc-I2)/V2;
	}
	f3=(Isc-V3*Gsc-I3);
	if (f3<0)
	{
		Isc-=f1;
		Gsc=(Isc-I3)/V3;
	}
	
	f1=log(Isc-V1*Gsc-I1);
	f2=log(Isc-V2*Gsc-I2);
	f3=log(Isc-V3*Gsc-I3);
	

	lD=((f1-f2)*(V2-V3)-(f2-f3)*(V1-V2))/((I1-I2)*(V2-V3)-I2*(V1-V2));
	lC=(f2-f3-I2*lD)/(V2-V3);
	B=exp(f1-V1*lC-I1*lD);
	A=par.Isc-B;
	E=-par.Gsc;
		
	P.Voc = par.Voc;
	P.Isc = par.Isc;
	P.Iph=A*lC/(lC-E*lD);
	P.I0=B*lC/(lC-E*lD);
	P.nkT=1/lC;
	P.Rs=lD/lC;
	P.Rsh=1/E-lD/lC;
	
	return P;
}


#define SCATTER 0.02

#define RandFloat(a,b) ((a)+((b)-(a))*((double)rand())/RAND_MAX)
#define RandInt(a,b) ((a)+rand()%((b)-(a)+1))
// oblique asumptote with random deviations to prime the diff-evol population
diodepars InitParsOA4DE(curvepars par, int *err)
{
	/* oblique asymptote approx */
	diodepars P;
	double V1, I1, V2, I2, V3, I3;
	double A,B,E, lD, lC, Isc, f1, f2, f3, d;
	double Gsc, Gmax;
	(*err)=0;
	// some basic checks whether a diode model could fit
	if ((par.Voc<1e-6)||(par.Isc<1e-6))
	{
		P.Rs  = 0;
		P.Rsh = 0;
		P.Iph = 0;	
		P.Voc = 0;
		P.Isc = 0;
		P.nkT = 0;
		P.I0  = 0;	
		(*err)=1;
		return P;
	}
	
	if (par.Voc<par.Vmpp)
	{
		// this is a mess 
		d=par.Vmpp;
		par.Vmpp=par.Voc;
		par.Voc=d;
	}
	if (par.Isc<par.Impp)
	{
		// this is a mess 
		d=par.Impp;
		par.Impp=par.Isc;
		par.Isc=d;
	}	
	
	if ((par.Isc-par.Impp)/par.Gsc<par.Vmpp)
	{
		// current  s-shape 
		par.Gsc=(par.Isc-par.Impp)/2/par.Vmpp;
	}
	
	
	if (par.Impp*par.Vmpp/par.Isc/par.Voc<=0.25)
	{
		// return a photovoltaic resistor
		P.Rs  = 0;
		P.Rsh = par.Voc/par.Isc;
		P.Iph = par.Isc;	
		P.Voc = par.Voc;
		P.Isc = par.Isc;
		P.nkT = 1;
		P.I0  = 0;
		(*err)=1;
		return P;
	}
	
	
	Isc=par.Isc;
	/*switch (RandInt(1,4)) // vary selected points
	{
		case 1:
			Select(par, 1, &V1, &I1);
			Select(par, 2, &V2, &I2);
			Select(par, 3, &V3, &I3);
			break;
		case 2:*/
			Select(par, 2, &V1, &I1);
			Select(par, 3, &V2, &I2);
			Select(par, 4, &V3, &I3);
		/*	break;
		case 3:
			Select(par, 1, &V1, &I1);
			Select(par, 3, &V2, &I2);
			Select(par, 4, &V3, &I3);
			break;
		default:
			Select(par, 2, &V1, &I1);
			Select(par, 3, &V2, &I2);
			Select(par, 4, &V3, &I3);
	}*/
	
	// some random disturbances in voltages and currents
	V1*=RandFloat(1-SCATTER,1+SCATTER);
	V2*=RandFloat(1-SCATTER,1+SCATTER);
	V3*=RandFloat(1-SCATTER,1+SCATTER);
	I1*=RandFloat(1-SCATTER,1+SCATTER);
	I2*=RandFloat(1-SCATTER,1+SCATTER);
	I3*=RandFloat(1-SCATTER,1+SCATTER);
	
	
	
	Gmax=(Isc-I1)/V1;
	if ((Isc-I2)/V2>Gmax)
		Gmax=(Isc-I2)/V2;
	if ((Isc-I3)/V3>Gmax)
		Gmax=(Isc-I3)/V3;
		
	Gsc=RandFloat(0,Gmax);
	
	f1=log(Isc-V1*Gsc-I1);
	f2=log(Isc-V2*Gsc-I2);
	f3=log(Isc-V3*Gsc-I3);
	
	
	lD=((f1-f2)*(V2-V3)-(f2-f3)*(V1-V2))/((I1-I2)*(V2-V3)-I2*(V1-V2));
	lC=(f2-f3-I2*lD)/(V2-V3);
	B=exp(f1-V1*lC-I1*lD);
	A=par.Isc-B;
	E=-par.Gsc;
		
	P.Voc = par.Voc;
	P.Isc = par.Isc;
	P.Iph=A*lC/(lC-E*lD);
	P.I0=B*lC/(lC-E*lD);
	P.nkT=1/lC;
	P.Rs=lD/lC;
	P.Rsh=(1/E-lD/lC);
	
	return P;
}

// oblique Init
diodepars InitObliqueAsymptote(curvepars par, int *err, double dI, double G, int i1, int i2, int i3)
{
	/* oblique asymptote approx */
	diodepars P;
	double V1, I1, V2, I2, V3, I3;
	double A,B,E, lD, lC, Isc, Gsc, f1, f2, f3, d;
	(*err)=0;
	// some basic checks whether a diode model could fit
	if ((par.Voc<1e-6)||(par.Isc<1e-6))
	{
		P.Rs  = 0;
		P.Rsh = 0;
		P.Iph = 0;	
		P.Voc = 0;
		P.Isc = 0;
		P.nkT = 0;
		P.I0  = 0;	
		(*err)=1;
		return P;
	}
	
	if (par.Voc<par.Vmpp)
	{
		// this is a mess 
		d=par.Vmpp;
		par.Vmpp=par.Voc;
		par.Voc=d;
		(*err)=2;
	}
	if (par.Isc<par.Impp)
	{
		// this is a mess 
		d=par.Impp;
		par.Impp=par.Isc;
		par.Isc=d;
		(*err)=3;
	}	
	
	if ((par.Isc-par.Impp)/par.Gsc<par.Vmpp)
	{
		// current  s-shape 
		par.Gsc=(par.Isc-par.Impp)/2/par.Vmpp;
		(*err)=4;
	}
	
	
	if (par.Impp*par.Vmpp/par.Isc/par.Voc<=0.25)
	{
		// return a photovoltaic resistor
		P.Rs  = 0;
		P.Rsh = par.Voc/par.Isc;
		P.Iph = par.Isc;	
		P.Voc = par.Voc;
		P.Isc = par.Isc;
		P.nkT = 1;
		P.I0  = 0;
		(*err)=5;
		return P;
	}
	Isc=par.Isc+fabs(dI);
	Gsc=fmod(fabs(G),par.Gsc);
	/* here we choose which points to use. Possible combinations are:
	 * 1,2,3 (qmp1, mpp,  qmp2)
	 * 1,2,4 (qmp1, mpp,  oc)
	 * 1,3,4 (qmp1, qmp2, oc)
	 * 2,3,4 (mpp,  qmp2, oc)
	 */
	
	Select(par, i1, &V1, &I1);
	Select(par, i2, &V2, &I2);
	Select(par, i3, &V3, &I3);
	f1=(Isc-V1*Gsc-I1);
	if (f1<0)
	{
		Isc-=f1;
		Gsc=(Isc-I1)/V1; // correct Gsc such that f1=0
	}
	f2=(Isc-V2*Gsc-I2);
	if (f2<0)
	{
		Isc-=f1;
		Gsc=(Isc-I2)/V2;
	}
	f3=(Isc-V3*Gsc-I3);
	if (f3<0)
	{
		Isc-=f1;
		Gsc=(Isc-I3)/V3;
	}
	
	f1=log(Isc-V1*Gsc-I1);
	f2=log(Isc-V2*Gsc-I2);
	f3=log(Isc-V3*Gsc-I3);
	

	lD=((f1-f2)*(V2-V3)-(f2-f3)*(V1-V2))/((I1-I2)*(V2-V3)-I2*(V1-V2));
	lC=(f2-f3-I2*lD)/(V2-V3);
	B=exp(f1-V1*lC-I1*lD);
	A=par.Isc-B;
	E=-par.Gsc;
		
	P.Voc = par.Voc;
	P.Isc = par.Isc;
	P.Iph=A*lC/(lC-E*lD);
	P.I0=B*lC/(lC-E*lD);
	P.nkT=1/lC;
	P.Rs=lD/lC;
	P.Rsh=1/E-lD/lC;
	
	return P;
}


double DPNLRMSErr(double x[], void *args);
// simple 2d optimizer for diode parameters
double ObliqueError(double x[], void *args)
{
	diodepars P;
	int err;
	double I,G, xx[5];
	curvepars *par;
	par=((curvepars *)args);
	I=exp(x[0]);
	G=(par->Isc+I-par->Iqmp1)/par->Vqmp1; // max G, min G=0;
	G=G/x[1];
	
	P=InitObliqueAsymptote(*par, &err, I, G, 1, 2, 4);
	
	xx[0]=P.I0*exp(par->Voc/P.nkT);
	xx[1]=P.nkT;
	xx[2]=P.Iph;
	xx[3]=P.Rs;
	xx[4]=P.Rsh;
	return DPNLRMSErr(xx, args);
/*	
	I=Diode_I(0, &G, P);
	a=(I-par->Isc)/par->Isc;
	E=a*a;
	a=par->Voc*(G-par->Gsc)/par->Isc;
	E+=a*a;
	I=Diode_I(par->Vmpp, &G, P);
	a=par->Voc*(G-par->Impp/par->Vmpp)/par->Isc;
	E+=a*a;
	if (err)
		E*=100;
	
	//printf("%e %e %e\n", exp(x[0])/par->Isc, x[1]/par->Gsc,E);
	return E;*/
}
	

diodepars OAInitializer(curvepars *par, double eps, double *Efit, MODELSET *S)
{
	double x[2];
	double step[2];
	double xmin[2];
	double I,G;
	int NCALLS, NRESTART, Err;
	x[0]=par->Isc/100;
	x[1]=2;
	if (par->Isc+x[0]<par->Impp)
		x[0]=par->Impp*1.01;
	x[0]=log(x[0]);	
		
	(*Efit)=ObliqueError(x, (void *)par);
	//printf("%e\t%e %e \t", (*Efit), x[0], x[1]);
		
	step[0]=fabs(x[0])/10;
	step[1]=1;
	nelmin (ObliqueError, (void *) par, 2, x, xmin, Efit, eps, step, S->fititer/20, S->fititer, &NCALLS,&NRESTART, &Err);
	//printf("%e\t%d %e %e\n", (*Efit), NCALLS, xmin[0], xmin[1]);
	
	I=exp(xmin[0]);
	G=(par->Isc+I-par->Iqmp1)/par->Vqmp1; // max G, min G=0;
	G=G/xmin[1];
	return InitObliqueAsymptote(*par, &Err, I, G, 1, 2, 4);
}	


//    0   1   2     3  4 
// x: Ix  nkT Iph, Rs, Rsh
// convert a neldel mead x to a x of the RHS routine below
void NMXtoX(double x[], double xx[], curvepars par)
{
	xx[0]=fabs(x[0]);
	xx[1]=fabs(x[1]);
	xx[2]=fabs(x[2]);
	xx[3]=fabs(x[3]);
	xx[4]=fabs(x[4]);
}

#define VEC_ERROR(b) (b[0]*b[0]+b[1]*b[1]+b[2]*b[2]+b[3]*b[3]+b[4]*b[4]);
// temporary defines to keep code readable
#define Ix   x[0]
#define nkT  x[1]
#define Iph  x[2]
#define Rs   x[3]
#define Rsh  x[4]
#define Isc  par.Isc
#define Impp par.Impp
#define Voc  par.Voc
#define Vmpp par.Vmpp
#define Gsc  par.Gsc

void RHS(double b[5], double x[5], curvepars par)
{
	int i;
	i=0;
	b[i++]= ((-Ix*(exp((Isc*Rs)/nkT)-1)*exp(-Voc/nkT))-(Isc*Rs)/Rsh-Isc+Iph)/Isc;
	b[i++]= ((-Ix*(exp((Vmpp+Impp*Rs)/nkT)-1)*exp(-Voc/nkT))+((-Vmpp)-Impp*Rs)/Rsh+Iph-Impp)/Isc;
	b[i++]= ((-Ix*exp(-Voc/nkT)*(exp(Voc/nkT)-1))-Voc/Rsh+Iph)/Isc;
	b[i++]= (Isc*((Rsh*nkT)/(Ix*Rsh*exp((Vmpp+Impp*Rs)/nkT-Voc/nkT)+nkT)-Vmpp/Impp+Rs))/Voc;
	b[i++]= (Voc*((Ix*Rsh*exp((Isc*Rs)/nkT-Voc/nkT)+nkT)/(Ix*Rs*Rsh*exp((Isc*Rs)/nkT-Voc/nkT)+(Rsh+Rs)*nkT)-Gsc))/Isc;
}

#undef Ix   
#undef nkT  
#undef Iph  
#undef Rs   
#undef Rsh  
#undef Isc
#undef Impp
#undef Voc
#undef Vmpp
#undef Gsc
double DPNLRMSErr(double x[], void *args) // neldel mead error function, converts 10 parameter double vector to BPar
{
	double b[5], xx[5], E;
	curvepars par;
	par=*((curvepars *) args);
	NMXtoX(x, xx, par);
	RHS(b, xx, par);
	E=VEC_ERROR(b);
	if (!isfinite(E))
		E=1e99;
	return E;
}
#define NPOP 100
#define NPRIM 90
// use differntial evolution
diodepars DEParameterize(curvepars par, diodepars P, MODELSET *S)
{
	double **prime;
	int Nprime=0, Err;
	double lb[5], ub[5];
	double xmin[5];
	de_settings Se=de_default;
	double EPS2=S->fiteps;
	time_t t;
	srand((unsigned) time(&t));
	
	/* square EPS to save some sqrts */
	EPS2*=EPS2*5;
	
	Se.max_iter=S->fititer;
	if (Se.max_iter<1000)
		Se.max_iter=1000;
	
	Se.Npop=S->npop;
	Se.CR=0.3;
	Se.F=0.8;
	Se.eps=S->fiteps;
	Se.prime=1;
	
	P=InitParsOA(par, &Err, 2, 3, 4);
	if((prime=Emalloc(NPRIM*sizeof(double *)))==NULL)
		return P;
	if((prime[Nprime]=Emalloc(5*sizeof(double)))==NULL)
	{
		free(prime);
		return P;
	}
	prime[Nprime][0]=fabs(P.I0)*exp(par.Voc/fabs(P.nkT)); // Ix=I0*exp(Voc/nkT)
	prime[Nprime][1]=fabs(P.nkT);
	prime[Nprime][2]=fabs(P.Iph);
	prime[Nprime][3]=fabs(P.Rs);
	prime[Nprime][4]=fabs(P.Rsh);
	Nprime++;
	if (!Err)
	{
		for (;Nprime<NPRIM;Nprime++)
		{
			P=InitParsOA4DE(par, &Err);
			if((prime[Nprime]=Emalloc(5*sizeof(double)))==NULL)
			{
				Nprime--;
				for(;Nprime>=0;Nprime--)
					free(prime[Nprime]);
				free(prime);
				return P;
			}
			prime[Nprime][0]=fabs(P.I0)*exp(par.Voc/fabs(P.nkT)); // Ix=I0*exp(Voc/nkT)
			prime[Nprime][1]=fabs(P.nkT);
			prime[Nprime][2]=fabs(P.Iph);
			prime[Nprime][3]=fabs(P.Rs);
			prime[Nprime][4]=fabs(P.Rsh);
		}
	}
	
	ub[0]=2*par.Isc; // Ix cannot exceed Isc
	ub[1]=(par.Voc-par.Vmpp)/log(par.Isc/(par.Isc-par.Impp)); // nkT cannot exceed the value without shunt or series resistance
	ub[2]=2*par.Isc;
	ub[3]=par.Roc;
	ub[4]=par.Gsc;
	lb[0]=par.Isc-par.Gsc*par.Voc;
	lb[1]=ub[1]/100;
	lb[2]=0.9*par.Isc;
	lb[3]=par.Roc/1000;
	lb[4]=1e-12;
	
	if (isfinite(diffevol_prime(DPNLRMSErr, (void *) (&par), 5, prime, Nprime, xmin, lb, ub, NULL, &Se)))
	{
		P.I0=xmin[0]*exp(-par.Voc/xmin[1]);
		P.nkT=xmin[1];
		P.Iph=xmin[2];
		P.Rs=xmin[3];
		P.Rsh=xmin[4];
	}
	
	for (;Nprime>0;Nprime--)
		free(prime[Nprime-1]);
		
	free(prime);
	
	return P;
}



/********************************* GSL optimizer ************************/

#define Ix   fabs(gsl_vector_get(x,0))
#define nkT  fabs(gsl_vector_get(x,1))
#define Iph  fabs(gsl_vector_get(x,2))
#define Rs   fabs(gsl_vector_get(x,3))
#define Rsh  fabs(gsl_vector_get(x,4))
#define Isc  par.Isc
#define Impp par.Impp
#define Voc  par.Voc
#define Vmpp par.Vmpp
#define Gsc  par.Gsc
void errhandler (const char * reason, const char * file,  int line, int gsl_errno)
{
	printf("%s\n", reason);
}

int RHS_gsl(const gsl_vector * x, void *data, gsl_vector * b)
{
	int i;
	curvepars par;
	par=*(curvepars *)data;
	
	/* RHS */
	i=0;
	gsl_vector_set(b,i++,((-Ix*(exp((Isc*Rs)/nkT)-1)*exp(-Voc/nkT))-(Isc*Rs)/Rsh-Isc+Iph)/Isc);
	gsl_vector_set(b,i++,((-Ix*(exp((Vmpp+Impp*Rs)/nkT)-1)*exp(-Voc/nkT))+((-Vmpp)-Impp*Rs)/Rsh+Iph-Impp)/Isc);
	gsl_vector_set(b,i++,((-Ix*exp(-Voc/nkT)*(exp(Voc/nkT)-1))-Voc/Rsh+Iph)/Isc);
	gsl_vector_set(b,i++,(Isc*((Rsh*nkT)/(Ix*Rsh*exp((Vmpp+Impp*Rs)/nkT-Voc/nkT)+nkT)-Vmpp/Impp+Rs))/Voc);
	gsl_vector_set(b,i++,(Voc*((Ix*Rsh*exp((Isc*Rs)/nkT-Voc/nkT)+nkT)/(Ix*Rs*Rsh*exp((Isc*Rs)/nkT-Voc/nkT)+(Rsh+Rs)*nkT)-Gsc))/Isc);
	return GSL_SUCCESS;
}

int System_gsl(const gsl_vector * x, void *data, gsl_matrix * A)
{
	int i,j;
	curvepars par;
	par=*(curvepars *)data;
	/* A */
	i=0;
	j=0;	
	
	gsl_matrix_set(A,i++,j,-((exp((Isc*Rs)/nkT)-1)*exp(-Voc/nkT))/Isc);
	gsl_matrix_set(A,i++,j,-(exp(-Voc/nkT)*(exp(Vmpp/nkT+(Impp*Rs)/nkT)-1))/Isc);
	gsl_matrix_set(A,i++,j,-(exp(-Voc/nkT)*(exp(Voc/nkT)-1))/Isc);
	gsl_matrix_set(A,i++,j,-(Isc*(Rsh*Rsh)*nkT*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT))/(2*Ix*Rsh*Voc*nkT*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT)+(Ix*Ix)*(Rsh*Rsh)*Voc*exp((2*Vmpp)/nkT+(2*Impp*Rs)/nkT)+Voc*(nkT*nkT)*exp((2*Voc)/nkT)));
	gsl_matrix_set(A,i++,j,((Rsh*Rsh)*Voc*nkT*exp(Voc/nkT+(Isc*Rs)/nkT))/((2*Isc*Ix*Rs*(Rsh*Rsh)+2*Isc*Ix*(Rs*Rs)*Rsh)*nkT*exp(Voc/nkT+(Isc*Rs)/nkT)+(Isc*(Rsh*Rsh)+2*Isc*Rs*Rsh+Isc*(Rs*Rs))*(nkT*nkT)*exp((2*Voc)/nkT)+Isc*(Ix*Ix)*(Rs*Rs)*(Rsh*Rsh)*exp((2*Isc*Rs)/nkT)));


	j++;
	i=0;
	gsl_matrix_set(A,i++,j,-(((Ix*Voc-Isc*Ix*Rs)*exp((Isc*Rs)/nkT)-Ix*Voc)*exp(-Voc/nkT))/(Isc*(nkT*nkT)));
	gsl_matrix_set(A,i++,j,-(exp(-Voc/nkT)*((Ix*Voc-Ix*Vmpp-Impp*Ix*Rs)*exp(Vmpp/nkT+(Impp*Rs)/nkT)-Ix*Voc))/(Isc*(nkT*nkT)));
	gsl_matrix_set(A,i++,j,(Ix*Voc*exp(-Voc/nkT))/(Isc*(nkT*nkT)));
	gsl_matrix_set(A,i++,j,((Isc*Ix*(Rsh*Rsh)*nkT-Isc*Ix*(Rsh*Rsh)*Voc+Isc*Ix*(Rsh*Rsh)*Vmpp+Impp*Isc*Ix*Rs*(Rsh*Rsh))*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT))/(2*Ix*Rsh*Voc*(nkT*nkT)*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT)+(Ix*Ix)*(Rsh*Rsh)*Voc*nkT*exp((2*Vmpp)/nkT+(2*Impp*Rs)/nkT)+Voc*(nkT*nkT*nkT)*exp((2*Voc)/nkT)));
	gsl_matrix_set(A,i++,j,-((Ix*(Rsh*Rsh)*Voc*nkT-Ix*(Rsh*Rsh)*(Voc*Voc)+Isc*Ix*Rs*(Rsh*Rsh)*Voc)*exp(Voc/nkT+(Isc*Rs)/nkT))/((2*Isc*Ix*Rs*(Rsh*Rsh)+2*Isc*Ix*(Rs*Rs)*Rsh)*(nkT*nkT)*exp(Voc/nkT+(Isc*Rs)/nkT)+(Isc*(Rsh*Rsh)+2*Isc*Rs*Rsh+Isc*(Rs*Rs))*(nkT*nkT*nkT)*exp((2*Voc)/nkT)+Isc*(Ix*Ix)*(Rs*Rs)*(Rsh*Rsh)*nkT*exp((2*Isc*Rs)/nkT)));


	j++;
	i=0;
	gsl_matrix_set(A,i++,j,1/Isc);
	gsl_matrix_set(A,i++,j,1/Isc);
	gsl_matrix_set(A,i++,j,1/Isc);
	gsl_matrix_set(A,i++,j,0);
	gsl_matrix_set(A,i++,j,0);


	j++;
	i=0;
	gsl_matrix_set(A,i++,j,-(exp(-Voc/nkT)*(nkT*exp(Voc/nkT)+Ix*Rsh*exp((Isc*Rs)/nkT)))/(Rsh*nkT));
	gsl_matrix_set(A,i++,j,-(exp(-Voc/nkT)*(Impp*Ix*Rsh*exp(Vmpp/nkT+(Impp*Rs)/nkT)+Impp*nkT*exp(Voc/nkT)))/(Isc*Rsh*nkT));
	gsl_matrix_set(A,i++,j,0);
	gsl_matrix_set(A,i++,j,((2*Isc*Ix*Rsh*nkT-Impp*Isc*Ix*(Rsh*Rsh))*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT)+Isc*(Ix*Ix)*(Rsh*Rsh)*exp((2*Vmpp)/nkT+(2*Impp*Rs)/nkT)+Isc*(nkT*nkT)*exp((2*Voc)/nkT))/(2*Ix*Rsh*Voc*nkT*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT)+(Ix*Ix)*(Rsh*Rsh)*Voc*exp((2*Vmpp)/nkT+(2*Impp*Rs)/nkT)+Voc*(nkT*nkT)*exp((2*Voc)/nkT)));
	gsl_matrix_set(A,i++,j,-((2*Ix*Rsh*Voc*nkT-Isc*Ix*(Rsh*Rsh)*Voc)*exp(Voc/nkT+(Isc*Rs)/nkT)+Voc*(nkT*nkT)*exp((2*Voc)/nkT)+(Ix*Ix)*(Rsh*Rsh)*Voc*exp((2*Isc*Rs)/nkT))/((2*Isc*Ix*Rs*(Rsh*Rsh)+2*Isc*Ix*(Rs*Rs)*Rsh)*nkT*exp(Voc/nkT+(Isc*Rs)/nkT)+(Isc*(Rsh*Rsh)+2*Isc*Rs*Rsh+Isc*(Rs*Rs))*(nkT*nkT)*exp((2*Voc)/nkT)+Isc*(Ix*Ix)*(Rs*Rs)*(Rsh*Rsh)*exp((2*Isc*Rs)/nkT)));


	j++;
	i=0;
	gsl_matrix_set(A,i++,j,Rs/(Rsh*Rsh));
	gsl_matrix_set(A,i++,j,(Vmpp+Impp*Rs)/(Isc*(Rsh*Rsh)));
	gsl_matrix_set(A,i++,j,Voc/(Isc*(Rsh*Rsh)));
	gsl_matrix_set(A,i++,j,(Isc*(nkT*nkT)*exp((2*Voc)/nkT))/(2*Ix*Rsh*Voc*nkT*exp(Voc/nkT+Vmpp/nkT+(Impp*Rs)/nkT)+(Ix*Ix)*(Rsh*Rsh)*Voc*exp((2*Vmpp)/nkT+(2*Impp*Rs)/nkT)+Voc*(nkT*nkT)*exp((2*Voc)/nkT)));
	gsl_matrix_set(A,i++,j,-(Voc*(nkT*nkT)*exp((2*Voc)/nkT))/((2*Isc*Ix*Rs*(Rsh*Rsh)+2*Isc*Ix*(Rs*Rs)*Rsh)*nkT*exp(Voc/nkT+(Isc*Rs)/nkT)+(Isc*(Rsh*Rsh)+2*Isc*Rs*Rsh+Isc*(Rs*Rs))*(nkT*nkT)*exp((2*Voc)/nkT)+Isc*(Ix*Ix)*(Rs*Rs)*(Rsh*Rsh)*exp((2*Isc*Rs)/nkT)));

	return GSL_SUCCESS;
}
#undef Ix  
#undef nkT  
#undef Iph  
#undef Rs   
#undef Rsh  
#undef Isc
#undef Impp
#undef Voc
#undef Vmpp
#undef Gsc

void _PrinPars(curvepars par)
{
	fprintf(stderr,"Error:\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n\n", par.Voc, par.Isc, par.Vmpp, par.Impp, par.Vqmp1, par.Iqmp1,par.Vqmp2, par.Iqmp2, par.Gsc, par.Roc);
}

/* GSL nonlinear least squares */
diodepars LevenbergMarquardtOptimizer(curvepars par, diodepars P, MODELSET *S)
{
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *w;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	const size_t n = 5;
	const size_t p = 5;
	int status;

	gsl_vector *b;
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	
	double x_init[5], weights[5];
	gsl_vector_view x;
	gsl_vector_view wts;
	
	double chisq0;
	int info;
	size_t i;

	const double xtol = S->fiteps;
	const double gtol = S->fiteps;
	const double ftol = 0.0;	
	gsl_set_error_handler_off();
	//gsl_set_error_handler (&errhandler);
	
	x_init[0]=P.I0*exp(par.Voc/P.nkT);
	x_init[1]=P.nkT;
	x_init[2]=P.Iph;
	x_init[3]=P.Rs;
	x_init[4]=P.Rsh;
	
	x = gsl_vector_view_array (x_init, p);
	wts = gsl_vector_view_array(weights, n);
  

	/* define the function to be minimized */
	fdf_params.h_fvv=1e-6;
	fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel; // with geodesic accell
	//fdf_params.trs = gsl_multifit_nlinear_trs_lm; // without geodesic accell
	fdf.f = RHS_gsl;
	fdf.df = System_gsl;   
	fdf.fvv = NULL;     
	fdf.n = n;
	fdf.p = p;
	fdf.params = &par;
	for (i = 0; i < n; i++)
		weights[i] = 1.0;
	

	/* allocate workspace with default parameters */
	w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

	/* initialize solver with starting point and weights */
	status=gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);
	if (status)
		return P;

	/* compute initial cost function */
	b = gsl_multifit_nlinear_residual(w);
	status=gsl_blas_ddot(b, b, &chisq0);

	if (status)
		return P;

	/* solve the system with a maximum of 100 iterations */
	status=gsl_multifit_nlinear_driver( S->fititer, xtol, gtol, ftol, NULL, NULL, &info, w);
	
	if (status)
		return P;
	P.I0=fabs(gsl_vector_get(w->x,0));
	P.nkT=fabs(gsl_vector_get(w->x,1));
	P.Iph=fabs(gsl_vector_get(w->x,2));
	P.Rs=fabs(gsl_vector_get(w->x,3));
	P.Rsh=fabs(gsl_vector_get(w->x,4));
	P.I0*=exp(-par.Voc/P.nkT);
	P.Voc=par.Voc;
	P.Isc=par.Isc;
	gsl_multifit_nlinear_free (w);
	gsl_matrix_free (covar);

	gsl_set_error_handler (NULL);
	return P;
}



// fit using differential evolution
generic_pars DiodeFitDE(curvepars par, curve C, double *Efit, MODELSET *S)
{
	diodepars P;
	P=DEParameterize(par,P,S);
	P=LevenbergMarquardtOptimizer(par,P,S);
	(*Efit)=DiodeRMSErr_i(P, C,S);
	return Diode2GenericPars(P);
}

// fit using only non-linear least squares (faster than usinfg differential evolution but less robust)
generic_pars DiodeFit(curvepars par, curve C, double *Efit, MODELSET *S)
{
	diodepars P;
	int err;
	
	//P=InitParsOA(par, &err, 1,2,4);
	P=InitParsOA(par, &err, 2,3,4); 
	if (err)
		return Diode2GenericPars(P);
	P=LevenbergMarquardtOptimizer(par,P,S);
	//P=OAInitializer(&par, 1e-2, Efit);
	(*Efit)=DiodeRMSErr_i(P, C,S);
	if ((!isfinite(*Efit)||(*Efit>S->oeps))&&(S->fititer))
	{
		P=DEParameterize(par,P,S);
		P=LevenbergMarquardtOptimizer(par,P,S);
		(*Efit)=DiodeRMSErr_i(P, C,S);
	}
	return Diode2GenericPars(P);
}

// fit using only non-linear least squares (faster than usinfg differential evolution but less robust)
generic_pars DiodeFitGqmp(curvepars par, curve C, double *Efit, MODELSET *S)
{
	diodepars P;
	int err;
	// change Gsc value to reflect a more coarsely determined 
	// sc conductivity
	par.Gsc=(par.Isc-par.Iqmp1)/par.Vqmp1;
	
	//P=InitParsOA(par, &err, 1,2,4);
	P=InitParsOA(par, &err, 2,3,4); 
	if (err)
		return Diode2GenericPars(P);
	P=LevenbergMarquardtOptimizer(par,P,S);
	//P=OAInitializer(&par, 1e-2, Efit);
	(*Efit)=DiodeRMSErr_i(P, C,S);
	if ((!isfinite(*Efit)||(*Efit>S->oeps))&&(S->fititer))
	{
		P=DEParameterize(par,P,S);
		P=LevenbergMarquardtOptimizer(par,P,S);
		(*Efit)=DiodeRMSErr_i(P, C,S);
	}
	return Diode2GenericPars(P);
}

/* fitting routines to get a full temperature and irradiance dependent 
 * diode model
 */
 
/* from a dataset (esp vs T and G), determine the median diode properties
 * and conditions
 */
void MedianDiode(curvepars *P, double *G, double *T, int N, 
                        MODELSET *S, double *I0, double *NnkT0, 
                        double *Iph0, double *Rs0, double *Rsh0, 
                        double *G0, double *T0)
{
	double err;
	int i;
	double *I0d, *NnkTd, *Iphd, *Rsd, *Rshd;
	generic_pars gp;
	double *GG, *TT;
	curve C;
	
	(*I0)=-1;
	(*NnkT0)=-1;
	(*Iph0)=-1;
	(*Rs0)=-1;
	(*Rsh0)=-1;	
	
	// base parameters are all median values
	// copy G and T arrays to compute median
	if((GG=Emalloc(N*sizeof(double)))==NULL)
		return;
	if((TT=Emalloc(N*sizeof(double)))==NULL)
	{
		free(GG);
		return;
	}
	for (i=0;i<N;i++)
	{
		GG[i]=G[i];
		TT[i]=T[i];
	}	
	(*G0)=quick_median(GG,N);
	(*T0)=quick_median(TT,N);
	free(GG);
	free(TT);
	
	
	if((I0d=Emalloc(N*sizeof(double)))==NULL)
		return;
	if((NnkTd=Emalloc(N*sizeof(double)))==NULL)
	{
		free(I0d);
		return;
	}
	if((Iphd=Emalloc(N*sizeof(double)))==NULL)
	{
		free(I0d);
		free(NnkTd);
		return;
	}
	if((Rsd=Emalloc(N*sizeof(double)))==NULL)
	{
		free(I0d);
		free(NnkTd);
		free(Iphd);
		return;
	}
	if((Rshd=Emalloc(N*sizeof(double)))==NULL)
	{
		free(I0d);
		free(NnkTd);
		free(Iphd);
		free(Rsd);
		return;
	}
	C.N=5;
	if((C.P=Emalloc(C.N*sizeof(point)))==NULL)
	{
		free(I0d);
		free(NnkTd);
		free(Iphd);
		free(Rsd);
		free(Rshd);
		return;
	}
	C.P[0].x=0.0;
	C.P[4].y=0.0;
	
	for (i=0;i<N;i++)
	{
		C.P[0].y=P[i].Isc;
		C.P[1].x=P[i].Vqmp1;
		C.P[1].y=P[i].Iqmp1;
		C.P[2].x=P[i].Vmpp;
		C.P[2].y=P[i].Impp;
		C.P[3].x=P[i].Vqmp2;
		C.P[3].y=P[i].Iqmp2;
		C.P[4].x=P[i].Voc;
			
		gp=DiodeFit(P[i], C, &err, S);
		I0d[i]   = gp.pars[0];
		NnkTd[i] = gp.pars[1];
		Iphd[i]  = gp.pars[2];
		Rsd[i]   = gp.pars[3];
		Rshd[i]  = gp.pars[4];
	}	
	(*I0)=quick_median(I0d,N);
	(*NnkT0)=quick_median(NnkTd,N);
	(*Iph0)=quick_median(Iphd,N);
	(*Rs0)=quick_median(Rsd,N);
	(*Rsh0)=quick_median(Rshd,N);
	free(I0d);
	free(Iphd);
	free(NnkTd);
	free(Rsd);
	free(Rshd);
	free(C.P);
}

/* the optimization parameter set:
 * [Vx,Nnk,alpha,Rs,Rsh,Ea]
 * Voc0  = Voc at G0,T0 (better parameter to optimize than I0)
 * Iph0  = Iph at G0, T0 
 * NnkT0 = thermal voltage prefactor at T0
 * alpha = temp coeff Isc: Isc = Isc0 [1 + α(T2-T1)] G / G0
 * Rs    = series resistance
 * Gsh   = shunt conductance (better parameter to optimize than Rsh)
 * Ea    = activation energy (Eg/n) 
 */
double XtoDiodepars(double x[7],double G0,double T0, diodepars *P, MODELSET *S);
diodepars DiodeParsTranslate(double alpha, double Ea, double I00, double G0, double T0, double G, double T, diodepars *P, MODELSET *S);

void InitParameters(curvepars *P, double *G, double *T, int N, 
                    MODELSET *S, double x[7], double *T0, double *G0)
{
	double I0, NnkT0, Iph0, Rs0, Rsh0, Voc0;
	double Ea=1.12, alpha=5e-3; // rough values should be OK
	MedianDiode(P,G,T,N,S,&I0,&NnkT0,&Iph0,&Rs0,&Rsh0,G0,T0);
	
	if (I0<0)
	{
		x[0]=-1;
		x[1]=-1;
		x[2]=-1;
		x[3]=-1;
		x[4]=-1;
		x[5]=-1;
		x[6]=-1;
		return;
	}
	Voc0=NnkT0*log(Iph0/I0+1.0); // we could directly take the median of the Voc but I suppose it does not matter
	x[0]=Voc0;
	x[1]=Iph0;
	x[2]=NnkT0;
	x[3]=Rs0;
	x[4]=1/Rsh0;
	x[5]=alpha;
	x[6]=Ea;
}
#define K_B 8.6173899875662609e-05
double XtoDiodepars(double x[7],double G0,double T0, diodepars *P, MODELSET *S)
{
	double I00;
	// make sure parameters are positive!
	P->Iph=fabs(x[1]);
	
	I00=fabs(x[1]/(exp(x[0]/x[2])-1));
	P->I0=I00;
	
	P->nkT=fabs(x[2]);
	P->Rs=fabs(x[3]);
	P->Rsh=fabs(1/x[4]); 
	
	// solve the special points
	P->Isc=Diode_I(0, NULL, *P, S);
	P->Voc=SolveVoc(*P, S);
	return I00;
}

diodepars DiodeParsTranslate(double alpha, double Ea, double I00, double G0, double T0, double G, double T, diodepars *P, MODELSET *S)
{
	double kTT0, TT;
	diodepars Po;
	Po=*P;
	Po.Iph=P->Iph*(1+alpha*(T-T0))*G/G0;
	kTT0=K_B*T*T0;
	TT=T/T0;
	// negative Ea is unphysical
	Po.I0=I00*TT*TT*TT*exp(fabs(Ea)*(T-T0)/kTT0);
	Po.nkT=P->nkT*TT;
	// solve the special points
	Po.Isc=Diode_I(0, NULL, Po, S);
	Po.Voc=SolveVoc(Po, S);	
	return Po;
}


typedef struct {
	int N;
	double G0,T0;
	double *G,*T;
	curvepars *P;
	MODELSET *S;
} DiodeGTErrStruct;
double DiodeGTError(double x[7], void *args)
{
	DiodeGTErrStruct *S;
	int i;
	double E=0,e, I00, Voc0, Isc0;
	S=(DiodeGTErrStruct *)args;
	curvepars MP;
	diodepars D0,D;
	
	I00=XtoDiodepars(x,S->G0,S->T0,&D0, S->S);
	Voc0=D0.Voc;
	Isc0=D0.Isc;
	for (i=0;i<S->N;i++)
	{
		D=DiodeParsTranslate(x[5],x[6],I00,S->G0,S->T0,(S->G)[i],(S->T)[i],&D0,S->S);
		MP=SolPar(D.Voc, D.Isc, (void *) &D, &diodeIV_i, S->S);
		e=CurveParDist(MP, (S->P)[i],Voc0, Isc0);
		if (isfinite(e))
			E+=e;
		else
			E+=1.0; // punishment: 100% error (times 5 to match CurveParDist bahavior)
	}
	E/=(double)(S->N);// mean square error
	return E;	
}

#define SIMPLSIZE 0.2 // simplex size
#define MINSIMPLSIZE 1e-4 // simplex size
generic_pars DiodeGTFit(curvepars *C, double *G, double *T, int N, 
                        double *alpha, double *Ea, double *T0, 
                        double *G0, MODELSET *S)
{
	DiodeGTErrStruct ARG={0};
	double EPS2=S->fiteps, E;
	int NCALLS, NRESTART, Err;	
	double x[7], xmin[7], step[7];
	diodepars D={0};
	
	EPS2*=EPS2;
	// initialize S
	ARG.G=G;
	ARG.T=T;
	ARG.P=C;
	ARG.N=N;
	ARG.S=S;
	
	InitParameters(ARG.P,ARG.G,ARG.T,ARG.N,ARG.S, xmin, &ARG.T0, &ARG.G0); 
	if (xmin[0]<0)
	{
		generic_pars O={0};
		return O;
	}
	
	/* fill the x vector */
	x[0]=xmin[0];
	x[1]=xmin[1];
	x[2]=xmin[2];
	x[3]=xmin[3];
	x[4]=xmin[4];
	x[5]=xmin[5];
	x[6]=xmin[6];
	
	step[0]=SIMPLSIZE*x[0]+MINSIMPLSIZE;
	step[1]=SIMPLSIZE*x[1]+MINSIMPLSIZE;
	step[2]=SIMPLSIZE*x[2]+MINSIMPLSIZE;
	step[3]=SIMPLSIZE*x[3]+MINSIMPLSIZE;
	step[4]=SIMPLSIZE*x[4]+MINSIMPLSIZE;
	step[5]=SIMPLSIZE*x[5]+MINSIMPLSIZE;
	step[6]=SIMPLSIZE*x[6]+MINSIMPLSIZE;
	
	E=DiodeGTError(x, (void *) (&ARG));
	
	//printf("Initial Error: %e\n", E/(double)N);
	if (E>EPS2)
		nelmin (DiodeGTError, (void *) (&ARG), 7, x, xmin, &E, EPS2, step, 50, 1000, &NCALLS,&NRESTART, &Err);
	
	XtoDiodepars(xmin,*G0,*T0,&D, S);
	(*alpha) = xmin[5];
	(*T0)    = ARG.T0;
	(*G0)    = ARG.G0;
	(*Ea)    = fabs(xmin[6]);	
	return Diode2GenericPars(D);
}
#undef SIMPLSIZE
#undef MINSIMPLSIZE

generic_pars DiodeGTTranslate(generic_pars P0, double alpha, double Ea, 
                        double G, double T, double G0, double T0, 
                        MODELSET *S)
{
	diodepars D;
	double kTT0,TT;
	
	D=Generic2DiodePars(P0, S);
	
	// translate D	
	D.Iph=D.Iph*(1+alpha*(T-T0))*G/G0;
	kTT0=K_B*T*T0;
	TT=T/T0;
	D.I0=D.I0*TT*TT*TT*exp(Ea*(T-T0)/kTT0);
	D.nkT=D.nkT*TT;
	D.Isc=D.Iph/(1+D.Rs/D.Rsh);
	D.Voc=D.nkT*log(D.Iph/D.I0+1);
	return Diode2GenericPars(D);	
}
 
