/*
    Cubic Spline IV model.
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
#include "splinefit.h"
#include "curve.h"
#include "fitmodels.h"
#include "curvemodel.h"
#include "neldermead.h"
#include "diffevol.h"
#include "spline_iv.h"
#include "bezier3_iv.h"


// internal struct
typedef struct {
	double V[13];
	double I[13];
	double dI[13]; // 13 is the maximum size, 5 key points pluts 4*2 extra dummy knots
	int N;
	int index[5]; // index of the 5 key points (including sc, which is redundant as it is always 0)
} splinepars;

// Exported Function
// return curvepars according to the model fit.
curvepars SplinePars(generic_pars P, MODELSET *S)
{
	curvepars par;
	
	if ((P.TYPE!=M_SPLINE)&&(P.TYPE!=M_PVSPLINE))
	{
		fprintf(stderr,"Error: wrong parameter type in SplinePars\n");
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

// test splinepars whether they are monotonously decreasing
// index refers to one of 4 sections (between the 5 points on the curve)
int MonoDecrPars(splinepars par, int index)
{
	if ((index<0)||(index>3))
		return -1;
	
	// test if any slope goes up
	if (par.dI[par.index[index]]>0)
		return 0;
	if (par.dI[par.index[index+1]]>0)
		return 0;
	// test if value goes up
	if (par.I[par.index[index]]-par.I[par.index[index+1]]<0)
		return 0;
	// mono decr.
	return 1;
}

int MonoDecrInterp(splinepars par, int index)
{
	double m, a, b;
	if ((index<0)||(index>3))
		return -1;
	
	// test monotonous
	m=(par.I[par.index[index+1]]-par.I[par.index[index]])/(par.V[par.index[index+1]]-par.V[par.index[index]]);
	a=par.dI[par.index[index]]/m;
	b=par.dI[par.index[index+1]]/m;
	
	if (a+b-2>0)
	{
		// one of three conditions must be met
		if (a-pow(2*a+b-3,2)/(3*(a+b-2))>=0)
			return 1;
		if (a+2*b-3<=0)
			return 1;
		if (2*a+b-3<=0)
			return 1;
	}
	else
		return 1;
	// curve is non monotonous
	return 0;
}


int MonoDecreasing(splinepars par, int index)
{
	int r;
	if ((r=MonoDecrPars(par, index))==1)
		return MonoDecrInterp(par, index);
	return r;
}


#define MS_EPS 1e-10 

void MakeMonoInterp(splinepars *S, int index) // index 0,1,2,3
{
	int i;
	double dy, y, dy1, dy2, x1, x2, y1, y2, m;
	

	
	if (MonoDecrPars(*S, index)==0) // not possible to make mono
		return;
	if (MonoDecrInterp(*S, index)!=0) // nothing to do
		return;
	x1=S->V[S->index[index]];
	x2=S->V[S->index[index+1]];
	y1=S->I[S->index[index]];
	y2=S->I[S->index[index+1]];
	dy1=S->dI[S->index[index]];
	dy2=S->dI[S->index[index+1]];
	
	m=(y2-y1)/(x2-x1);

	dy=3*m-(dy1+dy2)/2;	
	
	if ((dy>-MS_EPS)&&((dy1-m)*(dy2-m)>0)) // add two knots
	{
		for (i=S->N-1;i>S->index[index];i--)
		{
			S->V[i+2]=S->V[i];
			S->I[i+2]=S->I[i];
			S->dI[i+2]=S->dI[i];
		}
		S->N+=2;
		y=(dy2*y1+dy1*y2)/(dy1+dy2);
		
		S->V[S->index[index]+1]=x1+3*(y-y1)/dy1;
		S->I[S->index[index]+1]=y;
		S->dI[S->index[index]+1]=-MS_EPS;
		
		S->V[S->index[index]+2]=x2+3*(y-y2)/dy2;
		S->I[S->index[index]+2]=y;
		S->dI[S->index[index]+2]=-MS_EPS;
		for (i=4;i>index;i--)
			S->index[i]+=2;
	}
	else // add one knot
	{
		for (i=S->N-1;i>S->index[index];i--)
		{
			S->V[i+1]=S->V[i];
			S->I[i+1]=S->I[i];
			S->dI[i+1]=S->dI[i];
		}
		S->N+=1;
		
		if ((dy1-m)*(dy2-m)<0) // in this case we go for a piece wise quadratic solution
		{
			S->V[S->index[index]+1]=(2*(y1-y2)+(m+dy2)*x2-(m+dy1)*x1)/(dy2-dy1);
			S->I[S->index[index]+1]=(2*((m+dy2)*y1-(m+dy1)*y2)+(m+dy1)*(m+dy2)*(x2-x1))/(2*(dy2-dy1));
			S->dI[S->index[index]+1]=m;
		}
		else
		{
			S->V[S->index[index]+1]=(x1+x2)/2;
			S->I[S->index[index]+1]=y1+((6*m+dy1-dy2)*(x2-x1))/12;
			S->dI[S->index[index]+1]=dy;
		}	
		for (i=4;i>index;i--)
			S->index[i]+=1;
	}	
}


void MakeMono(splinepars *S)
{
	int i;
	for (i=0;i<4;i++)
		MakeMonoInterp(S, i);
	
}


// convert general to internal parameter struct
splinepars Generic2SplinePars(generic_pars P)
{
	splinepars S;
	// default positions for the solar cell parameters
	S.N=5;
	S.index[0]=0;
	S.index[1]=1;
	S.index[2]=2;
	S.index[3]=3;
	S.index[4]=4;

	if ((P.TYPE!=M_SPLINE)&&(P.TYPE!=M_PVSPLINE))
	{
		int i;
		fprintf(stderr,"Error: wrong parameter type in Generic2SplinePars\n");
		for (i=0;i<5;i++)
		{
			S.V[i]=0;
			S.I[i]=0;
			S.dI[i]=0;
		}
		return S;
	}
	
	
	S.V[0]=0;
	S.V[1]=P.pars[4];
	S.V[2]=P.pars[2];
	S.V[3]=P.pars[6];
	S.V[4]=P.pars[0];
	
	S.I[0]=P.pars[1];
	S.I[1]=P.pars[5];
	S.I[2]=P.pars[3];
	S.I[3]=P.pars[7];
	S.I[4]=0;
	
	S.dI[0]=-P.pars[8];
	S.dI[1]=-(P.pars[5]-P.pars[3])/P.pars[4];
	S.dI[2]=-P.pars[3]/P.pars[2];
	S.dI[3]=-P.pars[7]/(P.pars[6]-P.pars[2]);
	S.dI[4]=-1/P.pars[9];
	if (P.TYPE==M_PVSPLINE)
		MakeMono(&S); // insert dummy knots to make it monotone
	return S;
}
	

// convert internal to general parameter struct
generic_pars Spline2GenericPars(splinepars S, MODEL_TYPE T)
{
	generic_pars P;
	P.TYPE=T;
	P.pars[0]=S.V[S.index[4]];
	P.pars[1]=S.I[S.index[0]];
	P.pars[2]=S.V[S.index[2]];
	P.pars[3]=S.I[S.index[2]];
	P.pars[4]=S.V[S.index[1]];
	P.pars[5]=S.I[S.index[1]];
	P.pars[6]=S.V[S.index[3]];
	P.pars[7]=S.I[S.index[3]];
	P.pars[8]=-S.dI[S.index[0]];
	P.pars[9]=-1.0/(S.dI[S.index[4]]);
	return P;
}	


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



/* return the index of the largest V value below V */
static int binsearch(double V, const splinepars *par)
{
	int i, imin, imax;
	
	if (par->V[0]>V)
		return -1;
	if (par->V[par->N-1]<V)
		return par->N-1;
	imin=0;
	imax=par->N-1;
	while (imax-imin>1)
	{
		i=(imin+imax)/2;
		if (par->V[i]<V)
			imin=i;
		else
			imax=i;
	}
	return imin;	
}

static double Spline_I(double V, double *dI, const splinepars *par)
{
	int i;
	double vmin, vplus, h, t;
	double imin, iplus, mmin, mplus;
	
	i=binsearch(V, par);
	if (i<0)
	{
		/* linear extrapolation */
		if (dI)
			(*dI)=par->dI[0];
		return par->I[0]+(V-par->V[0])*par->dI[0];
	}
	if (i==par->N-1)
	{
		/* linear extrapolation */
		if (dI)
			(*dI)=par->dI[par->N-1];
		return par->I[par->N-1]+(V-par->V[par->N-1])*par->dI[par->N-1];
	}
	vmin=par->V[i];
	vplus=par->V[i+1];
	imin=par->I[i];
	iplus=par->I[i+1];
	mmin=par->dI[i];
	mplus=par->dI[i+1];
	h=vplus-vmin;
	t=(V-vmin)/h;
	if (dI)
		(*dI)=(imin*H00dot(t)+h*mmin*H10dot(t)+iplus*H01dot(t)+h*mplus*H11dot(t))/h;
	return imin*H00(t)+h*mmin*H10(t)+iplus*H01(t)+h*mplus*H11(t);
}


// pars are splinepars
point splineIV_i (double t, point *dp, void *pars)
{
	point p;
	double *di=NULL;
	const splinepars *par;
	
	if (dp)
		di=&(dp->y);
	par=(const splinepars *)pars;
	p.x=(2*t-0.5)*par->V[par->index[4]]; // for t 0-1 we cover a bit wider voltage range to be sure
	p.y=Spline_I(p.x, di, par);
	if (dp)
	{
		// take care of derivatives we have to compute dx/dt and dy/dt
		dp->x=2*par->V[par->index[4]];
		dp->y*=2*par->V[par->index[4]];
	}
	return p;
}

// Exported Function
curve SplineToCurve(generic_pars P, int N, MODELSET *S)
{
	splinepars par;
	curve C;
	par=Generic2SplinePars(P);
	C=ModelToCurve((void *) &par, &splineIV_i, N, S);
	return C;
}

double SplineRMSErr_i(splinepars P, curve C, MODELSET *S)
{
	return RMSErr(P.V[P.index[4]], P.I[P.index[0]], (void *) &P, &splineIV_i, C,S);
}

// Exported Function
double SplineRMSErr(generic_pars P, curve C, MODELSET *S)
{
	splinepars par;
	double E;
	par=Generic2SplinePars(P);
	E=SplineRMSErr_i(par, C,S);
	return E;
}


curve SplineDistCurve_i(splinepars P, curve C, MODELSET *S)
{
	return DistCurve(P.V[P.index[4]], P.I[P.index[0]], (void *) &P, &splineIV_i, C,S);
}

// Exported Function
curve SplineDistCurve(generic_pars P, curve C, MODELSET *S)
{
	splinepars par;
	curve d;
	par=Generic2SplinePars(P);
	d=SplineDistCurve_i(par, C,S);
	return d;
}


// Exported Function
double SplineIV(double V, generic_pars P, MODELSET *S)
{	
	splinepars par;
	par=Generic2SplinePars(P);
	return Solve_IV(V,P.pars[0],NULL,(void *) &par, &splineIV_i, S);	
}

// Exported Function
double SplineVI(double I, generic_pars P, MODELSET *S)
{	
	splinepars par;
	par=Generic2SplinePars(P);
	return Solve_VI(I,P.pars[0],NULL,(void *) &par, &splineIV_i, S);	
}

/* the following routines are to fit a spline to data
 * It should only be used for plain splines and not PV-Splines which 
 * may or may not add dummy knots!
 */
// sort current/voltage arrays
int dcomp(const void *a, const void *b)
{ 
	if (*((double *)a)<*((double *)b))
		return -1;
	if (*((double *)a)>*((double *)b))
		return 1;
	return 0;
} 
void dsort(double *x, int Nt)
{ 
	qsort(x, Nt, sizeof(double), dcomp);
}


splinepars splinesolpar_i(splinepars P, double *E, MODELSET *S)
{
	curvepars par;
	splinepars P2;
	int i;
	
	par=SolPar(1.0, 1.0, (void *) &P, &splineIV_i, S);	
	
	P2.V[0]=0;
	P2.V[1]=par.Vqmp1;
	P2.V[2]=par.Vmpp;
	P2.V[3]=par.Vqmp2;
	P2.V[4]=par.Voc;
	
	P2.I[0]=par.Isc;
	P2.I[1]=par.Iqmp1;
	P2.I[2]=par.Impp;
	P2.I[3]=par.Iqmp2;
	P2.I[4]=0;
	
	P2.dI[0]=-par.Gsc;
	P2.dI[1]=-(par.Iqmp1-par.Impp)/par.Vqmp1;
	P2.dI[2]=-par.Impp/par.Vmpp;
	P2.dI[3]=-par.Iqmp2/(par.Vqmp2-par.Vmpp);
	P2.dI[4]=-1/par.Roc;
	
	
	
	if (E)
	{
		(*E)=0;
		for (i=1;i<5;i++)
			(*E)+=(P2.V[i]-P.V[i])*(P2.V[i]-P.V[i]);
		(*E)=sqrt(*E);
	}
	
	return P2;
}


// mathc90 sc2fit wrapper
// TODO: use GSL instead
#define NB 6 // number of knots (sc, mmp1, mpp, mmp2, oc, max(vcurve))
#define NW 11 // NW must be greater or equal to the number of knots+4, i.e. >= 6+4
#define KNEPS 1e-12
// following routine assumes a normalized curve (currents divided by approx Isc). The scale of the voltages does not matter
splinepars SplineLSQ_step(splinepars S, double maxx, double *x, double *y, int NXY, double *err, int *cerr)
{
	splinepars Sout;
	long int i, ierr;
	double sigfac, w[NB*NW], yknot[NB], ypknot[NB];
	double sd[1];;
	double b[NB];
	int nb;
	(*cerr)=0;
	for (i=0;i<5;i++)
		b[i]=S.V[i];
	if (S.V[4]-maxx>-KNEPS)
	{
		S.V[4]=maxx+KNEPS;
		nb=NB-1;
	}
	else
	{
		b[NB-1]=maxx;
		nb=NB;
	}
	dsort(b, nb);
	for (i=0;i<nb-1;i++)
		if (!(b[i]<b[i+1]))
			b[i+1]+=KNEPS;
	sd[0]=-1.0;
	sc2fit( x, y, sd, NXY, b, nb, w, NW, yknot, ypknot, &sigfac, &ierr );
	if (err)
		(*err)=sigfac;
	if (ierr)
	{
		(*cerr)=1;
		return S;
	}
	for (i=0;i<NB-1;i++)
	{
		Sout.V[i]=S.V[i];
		Sout.I[i]=yknot[i];
		Sout.dI[i]=ypknot[i];
	}
	return Sout;
}

#define f0 0.5

splinepars SplineLSQ_driver(splinepars par, curve C, MODELSET *S)
{
	splinepars Sp, Sn={0}, Sl;
	double *x, *y, maxx, err;
	double E=1, f=f0;
	int i, NXY=0, iter=0, cerr;
	if((x=Emalloc(C.N*sizeof(double)))==NULL)
	{
		return Sn;
	}
	if((y=Emalloc(C.N*sizeof(double)))==NULL)
	{
		free(x);
		return Sn;
	}
	
	for (i=0;i<C.N;i++)
	{
		if ((C.P[i].x>=0)&&(C.P[i].y>=0))
		{
			x[NXY]=(C.P[i].x/par.V[4]);
			y[NXY]=(C.P[i].y/par.I[0]);
			NXY++;
		}
	}
	maxx=x[NXY-1];
	
	Sp=par;
	// normalize
	for (i=0;i<5;i++)
	{
		Sp.V[i]/=par.V[4];
		Sp.I[i]/=par.I[0];
		Sp.dI[i]/=(par.I[0]/par.V[4]);
	}
	Sn=Sp;
	Sl=Sp;

	do
	{
		Sp.V[1]=f*Sl.V[1]+(1-f)*Sn.V[1];
		Sp.V[2]=f*Sl.V[2]+(1-f)*Sn.V[2];
		Sp.V[3]=f*Sl.V[3]+(1-f)*Sn.V[3];
		Sp.V[4]=f*Sl.V[4]+(1-f)*Sn.V[4];
		Sn=SplineLSQ_step(Sp, maxx, x, y, NXY, &err, &cerr);
		if (!cerr)
			Sl=Sp;
		Sn=splinesolpar_i(Sn, &E, S);
		iter++;
	} while ((E>S->fiteps)&&(iter<S->fititer));
	free(x);
	free(y);
	for (i=0;i<5;i++)
	{
		Sn.V[i]*=par.V[4];
		Sn.I[i]*=par.I[0];
		Sn.dI[i]*=(par.I[0]/par.V[4]);
	}
	
	return Sn;	
}

/* end of spline fitting routines */

generic_pars PlainSplineFit(curvepars par, curve C, double *Efit, MODELSET *S)
{	
	splinepars Sp;
	generic_pars P, P0;
	double e0;
	
	P.TYPE=M_SPLINE;
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
	
	Sp=Generic2SplinePars(P);
	P0=P;
	e0=SplineRMSErr_i(Sp, C,S);
	(*Efit)=e0;
	if (S->RefinePars&&S->fititer) // plain splines may be fitted to data
	{
		Sp=SplineLSQ_driver(Sp, C, S);
		(*Efit)=SplineRMSErr_i(Sp, C, S);
		if (((*Efit)>e0)||(!isfinite(*Efit)))
		{
			(*Efit)=e0;
			P=P0;
		}
		else
			P=Spline2GenericPars(Sp, M_SPLINE);
	}
	return P;
}


generic_pars PVSplineFit(curvepars par, curve C, double *Efit, MODELSET *S)
{	
	splinepars Sp;
	generic_pars P;
	
	P.TYPE=M_PVSPLINE;
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
	
	Sp=Generic2SplinePars(P);
	(*Efit)=SplineRMSErr_i(Sp, C, S);
	// we do not fit PV-Splines as they must match 1 to 1 with the given curvepars
	return P;
}

/* IEC 60891 v2 (2009)
* 
* 
*                        G2
* I2 = I1 [1 + α(T2-T1)] ──
*                        G1
* 
*                ⎡               ⎛G2⎞⎤
* V2 = V1 + Voc1 ⎢β(T2-T1) + a ln⎜──⎟⎥ - Rs (I2-I1) - κ I2 (T2-T1)
*                ⎣               ⎝G1⎠⎦
* 
* 
* (V1,I1)@(T1,G1) -> (V2,I2)@(T2,G2) -> (V3,I3)@(T1,G1)
* 
* makes V1 ≠ V3 and I1 ≠ I3 (if we are lucky V1 ≈ V3 and I1 ≈ I3)
* 
*/ 
/* when provided an arbitrary set iof IV's at various T and G values 
 * it is best to find an appropriate values for Tstd and Gstd such 
 * that these values are well defined within the dataset.
 * 
 * Here we use the mean T and G values provided.
 */
 

// defines for the old norm
#define IEC_v2_I2(I1,alpha,DT,GR) ((I1)*(1+(alpha)*(DT))*(GR))
#define IEC_v2_V2(V1,I1,Voc1,I2,alpha,beta,a,Rs,kappa,DT,GR) ((V1)+(Voc1)*((beta)*(DT)+(a)*log((GR)))-(Rs)*((I2)-(I1))-(kappa)*((I2)*(DT)))
#define IEC_v2_dIdV2(dIdV1,alpha,beta,a,Rs,kappa,DT,GR) (-(dIdV1)*(GR)*((alpha)*(DT)+1)/((dIdV1)*(GR)*(DT)*((alpha)*(DT)+1)*(kappa)-1+(GR)*(Rs)*(alpha)*(dIdV1)*(DT)+((GR)-1)*(Rs)*(dIdV1)))

void IEC_v2(double V1, double I1, double Voc1, double alpha, double beta, double a, double Rs, double kappa,
            double T1, double T2, double G1, double G2, double *V2o, double *I2o)
{
	double DT=(T2-T1), GR=(G2/G1);
	(*I2o)=IEC_v2_I2(I1,alpha,DT,GR);
	(*V2o)=IEC_v2_V2(V1,I1,Voc1,(*I2o),alpha,beta,a,Rs,kappa,DT,GR);
}

generic_pars IEC60891_v2(generic_pars P, double alpha, double beta, double a, double Rs, double kappa, double T1, double T2, double G1, double G2, MODELSET *S)
{
	splinepars Sp, Sn;
	curvepars par;
	int i;
	double DT, GR;
	if (P.TYPE==M_PVBEZIER)
		return B3IEC60891_v2(P, alpha, beta, a, Rs, kappa, T1, T2, G1, G2, S);
	Sp=Generic2SplinePars(P);
	Sn=Sp;
	
	DT=(T2-T1);
	GR=(G2/G1);
	
	for (i=0;i<Sp.N;i++)
	{
		Sn.I[i]=IEC_v2_I2(Sp.I[i],alpha,DT,GR);
		Sn.V[i]=IEC_v2_V2(Sp.V[i],Sp.I[i],Sp.V[4],Sn.I[i],alpha,beta,a,Rs,kappa,DT,GR);
		Sn.dI[i]=IEC_v2_dIdV2(Sp.dI[i],alpha,beta,a,Rs,kappa,DT,GR);
	}
	
	par=SolPar(Sn.V[Sn.index[4]], Sn.I[Sn.index[0]], (void *) &Sn, &splineIV_i, S);
	
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
	return P;
}
/*  IEC 60891 v3 (2021):
* 
* 
*         G2 ⎛1 + α(T2-Tstd) ⎞
* I2 = I1 ── ⎜────────────── ⎟
*         G1 ⎝1 + α(T1-Tstd) ⎠
*         
* 
* V2 = V1 - Rs1 (I2 - I1) - κ I2 (T2 - T1)
* 
*               ⎛                                             1       1  ⎞
*     + Voc,std ⎜β[f(G2) (T2 - Tstd) - f(G1) (T1 - Tstd)] + ───── - ─────⎟
*               ⎝                                           f(G2)   f(G1)⎠
* 
* Rs1 = Rs + κ(T1 - Tstd)
* 
*        Voc,std        2⎛Gstd⎞        ⎛Gstd⎞
* f(G) = ─────── = B2 ln ⎜────⎟ + B1 ln⎜────⎟ + 1
*         Voc(G)         ⎝ G  ⎠        ⎝ G  ⎠
* 
* 
* A lot more complicated than before but at least the round trip consistency is now given
* 
* 
* 
* 
* 
I2(I1):=(I1*(G2/G1)*((1+alpha*(T2-Ts))/(1+alpha*(T1-Ts))));
V2(V1,I1):=V1-(Rs+kappa*(T1-Ts))*((I1*(G2/G1)*((1+alpha*(T2-Ts))/(1+alpha*(T1-Ts))))-I1)-kappa*((I1*(G2/G1)*((1+alpha*(T2-Ts))/(1+alpha*(T1-Ts)))))*(T2-T1)+Vocs*(beta*((B1*log(Gs/G2)+B2*(log(Gs/G2))^2+1)*(T2-Ts)-(B1*log(Gs/G1)+B2*(log(Gs/G1))^2+1)*(T1-Ts))+1/(B1*log(Gs/G2)+B2*(log(Gs/G2))^2+1)-1/(B1*log(Gs/G1)+B2*(log(Gs/G1))^2+1));
dI2(h):=(I2(I1+h*dIdV1)-I2(I1-h*dIdV1))/2;
dV2(h):=(V2(V1+h,I1+h*dIdV1)-V2(V1-h,I1-h*dIdV1))/2;
dIdV2(h):=dI2(h)/dV2(h);
limit(dIdV2(h),h,0);

dIdV2=(((G2*Ts-G2*T2)*alpha-G2)*dIdV1)/((((G2-G1)*Ts^2+(2*G1*T1-2*G2*T2)*Ts+G2*T2^2-G1*T1^2)*alpha+(G1-G2)*Ts+G2*T2-G1*T1)*dIdV1*kappa+(((G1-G2)*Rs*Ts+G2*Rs*T2-G1*Rs*T1)*alpha+(G2-G1)*Rs)*dIdV1+(G1*Ts-G1*T1)*alpha-G1)

*/

void IEC_v3(double V1, double Vs, double I1, double dIdV1,
            double G1, double G2, double Gs, 
            double T1, double T2, double Ts,
            double alpha, double beta, double B1, double B2, double Rs, double kappa,
            double *V2o, double *I2o, double *dIdV2o)
{
	double I2,V2,dIdV2;
	double Rs1, l, FG1, FG2;
	
	if ((fabs((G1-G2)/Gs)<1e-9)&&(fabs((T1-T2)/Ts)<1e-9))
	{
		(*V2o)=V1;
		(*I2o)=I1;
		(*dIdV2o)=dIdV1;
		return;
	}
	
	
	I2=I1*(G2/G1)*(1+alpha*(T2-Ts))/(1+alpha*(T1-Ts));	
	Rs1=Rs+kappa*(T1-Ts);
	
	l=log(Gs/G1);
	FG1=(B2*l+B1)*l+1.0;
	
	l=log(Gs/G2);
	FG2=(B2*l+B1)*l+1.0;
	
	V2=V1-Rs1*(I2-I1)-kappa*I2*(T2-T1);
	V2+=Vs*(beta*(FG2*(T2-Ts)-FG1*(T1-Ts))+1/FG2-1/FG1);
	
	dIdV2=(((G2*Ts-G2*T2)*alpha-G2)*dIdV1)/((((G2-G1)*Ts*Ts+(2*G1*T1-2*G2*T2)*Ts+G2*T2*T2-G1*T1*T1)*alpha+(G1-G2)*Ts+G2*T2-G1*T1)*dIdV1*kappa+(((G1-G2)*Rs*Ts+G2*Rs*T2-G1*Rs*T1)*alpha+(G2-G1)*Rs)*dIdV1+(G1*Ts-G1*T1)*alpha-G1);
	(*I2o)=I2;
	(*V2o)=V2;
	(*dIdV2o)=dIdV2;
}

generic_pars IEC60891_v3(generic_pars P, double alpha, double beta, double B1, double B2, double Rs, double kappa, double Vocstd, double T1, double T2, double Tstd, double G1, double G2, double Gstd, MODELSET *S)
{
	splinepars Sp, Sn;
	curvepars par;
	int i;
	if (P.TYPE==M_PVBEZIER)
		return B3IEC60891_v3(P, alpha, beta, B1, B2, Rs, kappa, Vocstd, T1, T2, Tstd, G1, G2, Gstd, S);
	Sp=Generic2SplinePars(P);
	Sn=Sp;
	
	for (i=0;i<Sp.N;i++)
		IEC_v3(Sp.V[i], Vocstd, Sp.I[i], Sn.dI[i],G1, G2, Gstd,T1, T2, Tstd,alpha,beta,B1,B2,Rs,kappa,Sn.V+i,Sn.I+i,Sn.dI+i);
	
	par=SolPar(Sn.V[Sn.index[4]], Sn.I[Sn.index[0]], (void *) &Sn, &splineIV_i, S);
	
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
	return P;
}

/* ***********************************************************************
 * Fitting the IEC 60891 (v2 and v3) procedure 2 coefficients
 *
 * Given ESP measured at N conditions (T, G), find the coefficients that
 * best predict the measured ESP from one reference ESP at (Tstd, Gstd).
 * For a set of coefficients the reference ESP is the mean of the measured
 * ESP translated to Tstd, Gstd. It is then translated back to every
 * measured condition, and the fit minimizes the sum over the curves of the
 * RMS distance of the predicted curve to the 5 key points (sc, qmpp-, mpp,
 * qmpp+, oc) of the measured ESP.
 * Judging the coefficients in this forward direction matters for v2, which
 * is not invertible: translating back with the true coefficients does not
 * give identical ESP. The coefficients are defined at Tstd, Gstd, as the
 * translation routines use them (v3 is not invariant to the reference).
 *
 * Start values come from linear least squares fits of the v2 relations
 * for Isc, Voc and the mpp. A downhill simplex refines them. For large
 * data sets the simplex uses an evenly spaced subset of IEC_NFIT curves
 * (the cost is linear in the number of curves while 6 coefficients are
 * well determined by a few hundred); the reference ESP is computed from
 * all curves.
 *
 * The coefficients are kept in the v3 layout
 *   c[0] alpha, c[1] beta, c[2] a (v2) or B1 (v3), c[3] B2 (v3), c[4] Rs, c[5] kappa
 * ***********************************************************************/

extern void dgelsd_( int* m, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, double* s, double* rcond, int* rank,
                double* work, int* lwork, int* iwork, int* info );

/* least squares solution x (m) of A x = b with A (n x m, column major).
 * A and b are overwritten. A rank deficient system (e.g. all curves at the
 * same temperature) gets the minimum norm solution. Returns 1 on failure. */
static int LSQ(double *A, double *b, int n, int m, double *x)
{
	int nrhs=1, lda=n, ldb=n, lwork=-1, rank, info, i;
	double rcond=-1.0, wkopt, *s, *work=NULL;
	int *iwork;
	if (n<m)
		return 1;
	s=Emalloc(m*sizeof(double));
	iwork=Emalloc((11*m+32)*sizeof(int));
	if (s && iwork)
	{
		dgelsd_(&n, &m, &nrhs, A, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork, iwork, &info);
		lwork=(int)wkopt;
		if ((info==0)&&((work=Emalloc(lwork*sizeof(double)))!=NULL))
			dgelsd_(&n, &m, &nrhs, A, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
	}
	if (!work)
		info=-1;
	for (i=0;(info==0)&&(i<m);i++)
		x[i]=b[i];
	free(s);
	free(iwork);
	free(work);
	return (info!=0);
}

/* start values c (v3 layout, B2=0) from linear least squares fits of the v2 relations
 *   Isc_i Gstd/G_i (1 + alpha (Tstd-T_i)) = Iscstd
 *   Voc_i (1 + beta (Tstd-T_i) + a ln(Gstd/G_i)) = Vocstd
 *   Vmpp_i + Voc_i (beta (Tstd-T_i) + a ln(Gstd/G_i)) - Rs (Imppstd-Impp_i)
 *          - kappa Imppstd (Tstd-T_i) = Vmppstd
 * with Imppstd the mean Impp translated with alpha. Returns 1 on failure. */
static int IECStart(const curvepars *C, const double *G, const double *T, int N, double Tstd, double Gstd, double c[6])
{
	double *A, *b, x[3], Imppstd=0;
	int i, err=1;
	A=Emalloc(3*N*sizeof(double));
	b=Emalloc(N*sizeof(double));
	if (!A || !b)
		goto IECSTARTFINISH;
	for (i=0;i<N;i++)
	{
		A[i]=Gstd*C[i].Isc*(Tstd-T[i])/G[i];
		A[N+i]=-1.0;
		b[i]=-Gstd*C[i].Isc/G[i];
	}
	if (LSQ(A, b, N, 2, x))
		goto IECSTARTFINISH;
	c[0]=x[0];
	for (i=0;i<N;i++)
	{
		A[i]=-1.0;
		A[N+i]=C[i].Voc*(Tstd-T[i]);
		A[2*N+i]=C[i].Voc*log(Gstd/G[i]);
		b[i]=-C[i].Voc;
	}
	if (LSQ(A, b, N, 3, x))
		goto IECSTARTFINISH;
	c[1]=x[1];
	c[2]=x[2];
	for (i=0;i<N;i++)
		Imppstd+=IEC_v2_I2(C[i].Impp, c[0], Tstd-T[i], Gstd/G[i]);
	Imppstd/=(double)N;
	for (i=0;i<N;i++)
	{
		A[i]=C[i].Impp-Imppstd;
		A[N+i]=Imppstd*(T[i]-Tstd);
		A[2*N+i]=-1.0;
		b[i]=-C[i].Vmpp-C[i].Voc*(c[1]*(Tstd-T[i])+c[2]*log(Gstd/G[i]));
	}
	if (LSQ(A, b, N, 3, x))
		goto IECSTARTFINISH;
	c[3]=0;
	c[4]=fabs(x[0]);
	c[5]=x[1];
	err=0;
	for (i=0;i<6;i++)
		err|=!isfinite(c[i]);
IECSTARTFINISH:
	free(A);
	free(b);
	return err;
}

typedef struct {
	int version;           // IEC 60891 version, 2 or 3
	int nc;                // number of fitted coefficients
	const int *idx;        // their positions in the v3 layout
	const curvepars *C;    // the ESP to fit
	const double *G, *T;   // and their conditions
	int N;
	double Tstd, Gstd;     // reference conditions
	MODEL_TYPE M;          // interpolation model
	MODELSET *S;
	double c0[6], scale[6];// coefficients c = c0 + scale y
	generic_pars *tr;      // work space: the translated ESP
	generic_pars ref;      // mean translated ESP
	curve *K;              // the 5 key points of each measured ESP
} IECFit;

static void IECCoef(const IECFit *F, const double *y, double c[6])
{
	int k;
	for (k=0;k<6;k++)
		c[k]=F->c0[k];
	for (k=0;k<F->nc;k++)
		c[F->idx[k]]+=F->scale[F->idx[k]]*y[k];
	c[4]=fabs(c[4]); // Rs >= 0
}

/* v3: the Voc at the reference conditions equals the mean translated Voc.
 * At I=0 the translation gives Voc_std,i = Voc_i + Vocstd X_i, so
 * Vocstd = <Voc_i>/(1-<X_i>) */
static double IECVocstd(const IECFit *F, const double c[6])
{
	double l, f, X=0, V=0;
	int i;
	for (i=0;i<F->N;i++)
	{
		l=log(F->Gstd/F->G[i]);
		f=(c[3]*l+c[2])*l+1.0;
		X+=-c[1]*f*(F->T[i]-F->Tstd)+1.0-1.0/f;
		V+=F->C[i].Voc;
	}
	X/=(double)F->N;
	V/=(double)F->N;
	if (fabs(1-X)<1e-3) // not sensible, fall back to the mean Voc
		return V;
	return V/(1-X);
}

static int IECFinite(const generic_pars *P)
{
	int k;
	for (k=0;k<10;k++)
		if (!isfinite(P->pars[k]))
			return 0;
	return 1;
}

/* translate all ESP to the reference conditions and set F->ref to their
 * mean. Returns the number of finite translated ESP. */
static int IECTranslateAll(IECFit *F, const double c[6])
{
	generic_pars P;
	double Vocstd=0;
	int i, k, n=0;
	if (F->version==3)
		Vocstd=IECVocstd(F, c);
	for (k=0;k<10;k++)
		F->ref.pars[k]=0;
	F->ref.TYPE=F->M;
	P.TYPE=F->M;
	for (i=0;i<F->N;i++)
	{
		P.pars[0]=F->C[i].Voc;
		P.pars[1]=F->C[i].Isc;
		P.pars[2]=F->C[i].Vmpp;
		P.pars[3]=F->C[i].Impp;
		P.pars[4]=F->C[i].Vqmp1;
		P.pars[5]=F->C[i].Iqmp1;
		P.pars[6]=F->C[i].Vqmp2;
		P.pars[7]=F->C[i].Iqmp2;
		P.pars[8]=F->C[i].Gsc;
		P.pars[9]=F->C[i].Roc;
		if (F->version==2)
			F->tr[i]=IEC60891_v2(P, c[0], c[1], c[2], c[4], c[5], F->T[i], F->Tstd, F->G[i], F->Gstd, F->S);
		else
			F->tr[i]=IEC60891_v3(P, c[0], c[1], c[2], c[3], c[4], c[5], Vocstd, F->T[i], F->Tstd, F->Tstd, F->G[i], F->Gstd, F->Gstd, F->S);
		if (IECFinite(F->tr+i))
		{
			for (k=0;k<10;k++)
				F->ref.pars[k]+=F->tr[i].pars[k];
			n++;
		}
	}
	for (k=0;(n>0)&&(k<10);k++)
		F->ref.pars[k]/=(double)n;
	return n;
}

// RMS distance of the (translated) curve to the key points
static double IEC_RMSErr(generic_pars P, curve C, MODELSET *S)
{
	if (P.TYPE==M_PVBEZIER)
		return B3RMSErr(P, C, S);
	return SplineRMSErr(P, C, S);
}

/* objective: sum of the RMS distances of the reference ESP translated to
 * the measured conditions to the measured key points. A curve that fails to
 * translate counts as 1 (100%) */
static double IECErr(double y[], void *args)
{
	IECFit *F=(IECFit *)args;
	generic_pars Q;
	double c[6], e, E=0;
	int i;
	IECCoef(F, y, c);
	if (!IECTranslateAll(F, c))
		return (double)F->N;
	for (i=0;i<F->N;i++)
	{
		if (F->version==2)
			Q=IEC60891_v2(F->ref, c[0], c[1], c[2], c[4], c[5], F->Tstd, F->T[i], F->Gstd, F->G[i], F->S);
		else
			Q=IEC60891_v3(F->ref, c[0], c[1], c[2], c[3], c[4], c[5], F->ref.pars[0], F->Tstd, F->T[i], F->Tstd, F->Gstd, F->G[i], F->Gstd, F->S);
		e=IECFinite(&Q) ? IEC_RMSErr(Q, F->K[i], F->S) : 1.0;
		E+=isfinite(e) ? e : 1.0;
	}
	return E;
}

/* a curve the fit can use: finite ESP with a positive mpp at a positive irradiance */
static int IECValid(const curvepars *C, double G, double T)
{
	return isfinite(C->Voc) && isfinite(C->Isc) && isfinite(C->Vmpp) && isfinite(C->Impp)
	    && isfinite(C->Vqmp1) && isfinite(C->Iqmp1) && isfinite(C->Vqmp2) && isfinite(C->Iqmp2)
	    && isfinite(C->Gsc) && isfinite(C->Roc) && isfinite(G) && isfinite(T)
	    && (G>0) && (C->Voc>0) && (C->Isc>0) && (C->Vmpp>0) && (C->Impp>0);
}

#define IEC_NFIT 128   // maximum number of curves used by the simplex
#define IEC_STEP 0.2   // initial simplex size relative to the coefficient scale
#define IEC_POLISH 0.05 // simplex size of the restart from the optimum
/* fit IEC 60891 version 2 or 3. Returns the ESP at Tstd, Gstd and sets the
 * coefficients c (v3 layout). On failure all values are NaN. */
static generic_pars IEC60891Fit(int version, curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double c[6], MODEL_TYPE M, MODELSET *S)
{
	static const int idx2[5]={0, 1, 2, 4, 5}, idx3[6]={0, 1, 2, 3, 4, 5};
	generic_pars P;
	IECFit F={0};
	double y[6]={0}, ymin[6]={0}, step[6], floor[6], E0, E, Emin, R0=0, reqmin;
	point *kp=NULL;
	curvepars *Cs=NULL, *Cv=NULL;
	double *Gs=NULL, *Ts=NULL, *Gv=NULL, *Tv=NULL;
	int i, k, n, icount, numres, ifault, run;

	P.TYPE=M;
	for (k=0;k<10;k++)
		P.pars[k]=NAN;
	for (k=0;k<6;k++)
		c[k]=NAN;
	/* skip the curves the fit cannot use */
	Cv=Emalloc(N*sizeof(curvepars));
	Gv=Emalloc(N*sizeof(double));
	Tv=Emalloc(N*sizeof(double));
	if (!Cv || !Gv || !Tv)
		goto IECFITFINISH;
	for (i=0, k=0;i<N;i++)
		if (IECValid(C+i, G[i], T[i]))
		{
			Cv[k]=C[i];
			Gv[k]=G[i];
			Tv[k]=T[i];
			k++;
		}
	C=Cv;
	G=Gv;
	T=Tv;
	N=k;
	if ((N<3) || IECStart(C, G, T, N, Tstd, Gstd, F.c0))
		goto IECFITFINISH;

	F.version=version;
	F.nc=(version==2) ? 5 : 6;
	F.idx=(version==2) ? idx2 : idx3;
	/* the curves for the simplex: all, or an evenly spaced subset */
	n=(N>IEC_NFIT) ? IEC_NFIT : N;
	if (n<N)
	{
		Cs=Emalloc(n*sizeof(curvepars));
		Gs=Emalloc(n*sizeof(double));
		Ts=Emalloc(n*sizeof(double));
		if (!Cs || !Gs || !Ts)
			goto IECFITFINISH;
		for (i=0;i<n;i++)
		{
			k=(int)(((long)i*N)/n);
			Cs[i]=C[k];
			Gs[i]=G[k];
			Ts[i]=T[k];
		}
	}
	F.C=(n<N) ? Cs : C;
	F.G=(n<N) ? Gs : G;
	F.T=(n<N) ? Ts : T;
	F.N=n;
	F.Tstd=Tstd;
	F.Gstd=Gstd;
	F.M=M;
	F.S=S;
	F.tr=Emalloc(N*sizeof(generic_pars));
	F.K=Emalloc(n*sizeof(curve));
	kp=Emalloc(5*n*sizeof(point));
	if (!F.tr || !F.K || !kp)
		goto IECFITFINISH;
	for (i=0;i<n;i++)
	{
		F.K[i].N=5;
		F.K[i].P=kp+5*i;
		F.K[i].P[0].x=0;
		F.K[i].P[0].y=F.C[i].Isc;
		F.K[i].P[1].x=F.C[i].Vqmp1;
		F.K[i].P[1].y=F.C[i].Iqmp1;
		F.K[i].P[2].x=F.C[i].Vmpp;
		F.K[i].P[2].y=F.C[i].Impp;
		F.K[i].P[3].x=F.C[i].Vqmp2;
		F.K[i].P[3].y=F.C[i].Iqmp2;
		F.K[i].P[4].x=F.C[i].Voc;
		F.K[i].P[4].y=0;
	}

	/* scale of each coefficient: its start value, at least a typical size.
	 * R0 is a resistance scale for Rs and kappa */
	for (k=0;k<N;k++)
		R0+=C[k].Voc/C[k].Isc;
	R0/=(double)N;
	floor[0]=1e-4;    // alpha [1/K]
	floor[1]=1e-4;    // beta [1/K]
	floor[2]=1e-2;    // a, B1
	floor[3]=1e-2;    // B2
	floor[4]=1e-2*R0; // Rs
	floor[5]=1e-4*R0; // kappa [Ohm/K]
	for (k=0;k<6;k++)
		F.scale[k]=fabs(F.c0[k])+floor[k];

	Emin=E0=IECErr(y, &F);
	/* first a simplex from the start values, then a smaller one restarted
	 * from the optimum as the simplex may collapse prematurely */
	for (run=0;run<2;run++)
	{
		for (k=0;k<F.nc;k++)
			step[k]=run ? IEC_POLISH : IEC_STEP;
		reqmin=S->fiteps*Emin;
		reqmin*=reqmin;
		nelmin(IECErr, &F, F.nc, y, ymin, &E, reqmin, step, 10, S->fititer, &icount, &numres, &ifault);
		if (!isfinite(E) || (E>=Emin))
			break;
		Emin=E;
		for (k=0;k<F.nc;k++)
			y[k]=ymin[k];
	}
	// y holds the best point (the start values if the simplex did not improve on them)
	IECCoef(&F, y, c);
	// the reference ESP from all curves
	F.C=C;
	F.G=G;
	F.T=T;
	F.N=N;
	if (IECTranslateAll(&F, c)==0)
	{
		for (k=0;k<6;k++)
			c[k]=NAN;
		goto IECFITFINISH;
	}
	if (version==2)
		c[3]=0;
	P=F.ref;
IECFITFINISH:
	free(F.tr);
	free(F.K);
	free(kp);
	free(Cs);
	free(Gs);
	free(Ts);
	free(Cv);
	free(Gv);
	free(Tv);
	return P;
}
#undef IEC_NFIT
#undef IEC_STEP
#undef IEC_POLISH

/* fit the IEC 60891 v2 coefficients, returns the ESP at Tstd, Gstd (NaN on failure) */
generic_pars IEC60891Fit_v2(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double *alpha, double *beta, double *a, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S)
{
	double c[6];
	generic_pars P=IEC60891Fit(2, C, G, T, N, Tstd, Gstd, c, M, S);
	(*alpha)=c[0];
	(*beta)=c[1];
	(*a)=c[2];
	(*Rs)=c[4];
	(*kappa)=c[5];
	return P;
}

/* fit the IEC 60891 v3 coefficients, returns the ESP at Tstd, Gstd (NaN on failure) */
generic_pars IEC60891Fit_v3(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double *alpha, double *beta, double *B1, double *B2, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S)
{
	double c[6];
	generic_pars P=IEC60891Fit(3, C, G, T, N, Tstd, Gstd, c, M, S);
	(*alpha)=c[0];
	(*beta)=c[1];
	(*B1)=c[2];
	(*B2)=c[3];
	(*Rs)=c[4];
	(*kappa)=c[5];
	return P;
}

// some utility functions to test plausibility of a set of curvepars
/* hermite cubic spline base equations for the second derivative */
#define H00dotdot(t) (12*(t)-6)
#define H10dotdot(t) (6*(t)-4)
#define H01dotdot(t) (-12*(t)+6)
#define H11dotdot(t) (6*(t)-2)

// measure the convexity of a spline
double SplineConvexity(splinepars par, int index)
{
	double t, vmin, vplus, imin, iplus, mmin, mplus, h;
	double dd1, dd2, C=0;
	int istart, iend, i;
	if ((index<0)||(index>3))
		return 0;
	switch(index)
	{
		case 0:
			istart=0;
			iend=par.index[index+1];
			break;
		default:
			istart=par.index[index];
			iend=par.index[index+1];
			break;
	}	
	for (i=istart;i<iend;i++)
	{		
		vmin=par.V[i];
		vplus=par.V[i+1];
		imin=par.I[i];
		iplus=par.I[i+1];
		mmin=par.dI[i];
		mplus=par.dI[i+1];
		h=vplus-vmin;
	
		// we only need to check the second derivative at the interval endpoints as the second derivative is linear
		t=0;
		dd1=(imin*H00dotdot(t)+h*mmin*H10dotdot(t)+iplus*H01dotdot(t)+h*mplus*H11dotdot(t))/h/h;
		t=1;
		dd2=(imin*H00dotdot(t)+h*mmin*H10dotdot(t)+iplus*H01dotdot(t)+h*mplus*H11dotdot(t))/h/h;
	
		//integrate the convex parts only (i.e. concave parts are 0) 
		if ((dd1<0)&&(dd2>0))
			C+=(-dd2/dd1)*h*dd2/2;
		else if ((dd1>0)&&(dd2<0))
			C+=(-dd1/dd2)*h*dd1/2;
		else if ((dd1<0)&&(dd2<0))
			C+=h*(dd1+dd2)/2;
	}
	return C;
}


void Convexity(curvepars par, double *convex1, double *convex2, double *convex3, double *convex4)
{
	// measure the convexity of the curve in two parts: sc-mpp and mpp-voc
	splinepars S;
	generic_pars P;
	
	P.TYPE=M_PVSPLINE; // do not use SPLINES, splines may introduce convexity
					   // as splines may become non-monotonous under circumstances
					   // even if the actual curvepars are strictly monotonous
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
	
	S=Generic2SplinePars(P);
	
	(*convex1)=SplineConvexity(S,0)*par.Voc/par.Isc;
	(*convex2)=SplineConvexity(S,1)*par.Voc/par.Isc;
	(*convex3)=SplineConvexity(S,2)*par.Voc/par.Isc;
	(*convex4)=SplineConvexity(S,3)*par.Voc/par.Isc;
}
