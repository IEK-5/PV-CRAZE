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
// version 2 fitting utilities
extern void dgelsd_( int* m, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, double* s, double* rcond, int* rank,
                double* work, int* lwork, int* iwork, int* info );


int FitIsc(double *Isc, double *G, double *T, int N, double Gstd, double Tstd, double *Iscstd, double *alpha)
{
	double *M, *b, *s;
	int i;
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	double wkopt;
	double* work;
    int *iwork;
		
	if((M=Emalloc(N*2*sizeof(double)))==NULL)
	{
		return 1;
	}
	if((b=Emalloc(N*sizeof(double)))==NULL)
	{
		free(M);
		return 1;
	}
	for (i=0;i<N;i++)
	{
		M[i]=Gstd*Isc[i]*(Tstd-T[i])/G[i];
		M[N+i]=-1.0;
		b[i]=-Gstd*Isc[i]/G[i];
	}
		
	n=N;
	m=2;
	if((s=Emalloc(m*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		return 1;
	}
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(11*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		return 1;
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
		return 1;
	}
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*alpha)=0;
		(*Iscstd)=0;
	}
	else
	{
		(*alpha)=b[0];
		(*Iscstd)=b[1];
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
	if (info>0)
		return 1;
	else
		return 0;
}

int FitVoc(double *Voc, double *G, double *T, int N, double Gstd, double Tstd, double *Vocstd, double *beta, double *a)
{
	double *M, *b, *s;
	int i;
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	double wkopt;
	double* work;
    int *iwork;
			
	if((M=Emalloc(N*3*sizeof(double)))==NULL)
	{
		return 1;
	}
	if((b=Emalloc(N*sizeof(double)))==NULL)
	{
		free(M);
		return 1;
	}
	for (i=0;i<N;i++)
	{
		M[i]=-1.0;
		M[N+i]=Voc[i]*(Tstd-T[i]);
		M[2*N+i]=Voc[i]*log(Gstd/G[i]);
		b[i]=-Voc[i];
	}
		
	n=N;
	m=3;
	if((s=Emalloc(m*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		return 1;
	}
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(11*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		return 1;
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
		return 1;
	}
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*a)=0;
		(*Vocstd)=0;
		(*beta)=0;
	}
	else
	{
		(*Vocstd)=b[0];
		(*beta)=b[1];
		(*a)=b[2];
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
	if (info>0)
		return 1;
	else
		return 0;
}

int FitMpp(double *Vmpp, double *Impp, double *Voc, double *G, double *T, int N, double Gstd, double Tstd, double alpha, double beta, double a, double *Imppstd, double *Vmppstd, double *Rs, double *kappa)
{
	double *M, *b, *s, Istd=0;
	int i;
	int m, n, nrhs = 1, lda, ldb, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	double wkopt;
	double* work;
    int *iwork;
	for (i=0;i<N;i++)
		Istd+=IEC_v2_I2(Impp[i],alpha,Tstd-T[i],Gstd/G[i]);
	(*Imppstd)=Istd/(double)N;
	
	if((M=Emalloc(N*3*sizeof(double)))==NULL)
	{
		return 1;
	}
	if((b=Emalloc(N*sizeof(double)))==NULL)
	{
		free(M);
		return 1;
	}
	
	for (i=0;i<N;i++)
	{
		M[i]=(Impp[i]-(*Imppstd));
		M[N+i]=(*Imppstd)*(T[i]-Tstd);
		M[2*N+i]=-1.0;
		b[i]=-Vmpp[i]-Voc[i]*(beta*(Tstd-T[i])+a*log(Gstd/G[i]));
	}
		
	n=N;
	m=3;
	if((s=Emalloc(m*sizeof(double)))==NULL)
	{
		free(M);
		free(b);
		return 1;
	}
	lda=n;
	ldb=n;
	//print_matrix(n, m, M);
	
	if((iwork=Emalloc(11*m*sizeof(int)))==NULL)
	{
		free(M);
		free(b);
		free(s);
		return 1;
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
		return 1;
	}
	/* Solve the equations A*X = B */
	dgelsd_( &n, &m, &nrhs, M, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
					iwork, &info );
	/* Check for convergence */
	if( info > 0 )
	{
		(*Rs)=0;
		(*kappa)=0;
		(*Vmppstd)=0;
	}
	else
	{	
		(*Rs)=b[0];
		(*kappa)=b[1];
		(*Vmppstd)=b[2];
	}
	free(b);
	free(M);
	free(s);
	free(work);
	free(iwork);
	if (info>0)
		return 1;
	else
		return 0;
}
typedef struct {
	curvepars *C;
	curve P;// a curve containing only the 5 key points
	double *G, *T;
	int N;
	double Tstd, Gstd;
	MODEL_TYPE M; // interpolation model used for the translation
	MODELSET *S;
} errpars;

// RMS distance of the (translated) curve to the key points
static double IEC_RMSErr(generic_pars P, curve C, MODELSET *S)
{
	if (P.TYPE==M_PVBEZIER)
		return B3RMSErr(P, C, S);
	return SplineRMSErr(P, C, S);
}

/* optimization vector x preprocessing macros
 * x[0]: alpha
 * x[1]: beta
 * x[2]: a
 * x[3]: Rs
 * x[4]: kappa
 */
#define PRE_NONE(a)  (a) 		// preprocessor: do nothing
#define IPRE_NONE(a) (a) 		// inverse thereof
#define PRE_ABS(a)   (fabs(a))  // preprocessor: absolute value
#define IPRE_ABS(a)  (a) 		// inverse thereof
#define PRE_EXP(a)   (exp(a))   // preprocessor: absolute value
#define IPRE_EXP(a)  (log(a))	// inverse thereof

#define FX0(a)  (PRE_NONE(a))
#define IFX0(a) (IPRE_NONE(a))

#define FX1(a)  (PRE_NONE(a))
#define IFX1(a) (IPRE_NONE(a))

#define FX2(a)  (PRE_ABS(a))
#define IFX2(a) (IPRE_ABS(a))

#define FX3(a)  (PRE_ABS(a))
#define IFX3(a) (IPRE_ABS(a))

#define FX4(a)  (PRE_NONE(a))
#define IFX4(a) (IPRE_NONE(a))


double FitErr(double x[5], void *args)
{
	double E=0;
	errpars *pp;
	int i;
	generic_pars P2;
	pp=(errpars *) args;
	P2.TYPE=pp->M;
	for (i=0;i<pp->N;i++)
	{
		double e;
		P2.pars[0]=pp->C[i].Voc;
		P2.pars[1]=pp->C[i].Isc;
		P2.pars[2]=pp->C[i].Vmpp;
		P2.pars[3]=pp->C[i].Impp;
		P2.pars[4]=pp->C[i].Vqmp1;
		P2.pars[5]=pp->C[i].Iqmp1;
		P2.pars[6]=pp->C[i].Vqmp2;
		P2.pars[7]=pp->C[i].Iqmp2;
		P2.pars[8]=pp->C[i].Gsc;
		P2.pars[9]=pp->C[i].Roc;
		P2=IEC60891_v2(P2, FX0(x[0]), FX1(x[1]), FX2(x[2]), FX3(x[3]), FX4(x[4]), pp->T[i], pp->Tstd, pp->G[i], pp->Gstd, pp->S); // from stc to actual
		e=IEC_RMSErr(P2, pp->P, pp->S);
		if (!isfinite(e))
			e=1; // 100% error
		E+=e;
	}
	return E;
}

#define SIMPLSIZE 0.2 // simplex size
#define MINSIMPLSIZE 1e-4 // simplex size
void OptimizeIEC_v2(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, curvepars *Cstd, double *alpha, double *beta, double *a, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S)
{
	double x[5];
	double step[5];
	double xmin[5];
	int NCALLS, NRESTART, Err;	
	double EPS2=S->fiteps, E;
	errpars pp;
	
	
	/* square EPS to save some sqrts */
	EPS2*=EPS2*N*10;	
	
	/* fill the x vector */
	xmin[0]=x[0]=IFX0(*alpha);
	xmin[1]=x[1]=IFX1(*beta);
	xmin[2]=x[2]=IFX2(*a);
	xmin[3]=x[3]=IFX3(*Rs);
	xmin[4]=x[4]=IFX4(*kappa);
	
	
	step[0]=SIMPLSIZE*x[0]+MINSIMPLSIZE;
	step[1]=SIMPLSIZE*x[1]+MINSIMPLSIZE;
	step[2]=SIMPLSIZE*x[2]+MINSIMPLSIZE;
	step[3]=SIMPLSIZE*x[3]+MINSIMPLSIZE;
	step[4]=SIMPLSIZE*x[4]+MINSIMPLSIZE;
	
	pp.P.N=5;
	if((pp.P.P=Emalloc(5*sizeof(point)))==NULL)	
		return;
	pp.P.P[0].x=0;
	pp.P.P[0].y=Cstd->Isc;
	pp.P.P[1].x=Cstd->Vqmp1;
	pp.P.P[1].y=Cstd->Iqmp1;
	pp.P.P[2].x=Cstd->Vmpp;
	pp.P.P[2].y=Cstd->Impp;
	pp.P.P[3].x=Cstd->Vqmp2;
	pp.P.P[3].y=Cstd->Iqmp2;
	pp.P.P[4].x=Cstd->Voc;
	pp.P.P[4].y=0;
	pp.C=C;
	pp.G=G;
	pp.T=T;
	pp.N=N;
	pp.Tstd=Tstd;
	pp.Gstd=Gstd;
	pp.M=M;
	pp.S=S;
	
	
	E=FitErr(x, (void *) (&pp));
	//printf("Initial Error: %e\n", E/(double)N);
	if (E>EPS2)
		nelmin (FitErr, (void *) (&pp), 5, x, xmin, &E, EPS2, step, 50, 1000, &NCALLS,&NRESTART, &Err);
	//printf("Final Error:   %e\n", E/(double)N);
	free(pp.P.P);
	(*alpha)=FX0(xmin[0]);
	(*beta)=FX0(xmin[1]);
	(*a)=FX0(xmin[2]);
	(*Rs)=FX0(xmin[3]);
	(*kappa)=FX0(xmin[4]);
}

/* fit the IEC norm, return spline pars at Tstd,Gstd and set alpha, beta, ... etc. */
generic_pars IEC60891Fit_v2(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double *alpha, double *beta, double *a, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S)
{
	int i;
	double *arr, Ts, Gs;
	generic_pars P={0}, P2;
	curvepars Cstd;

	if((arr=Emalloc(4*N*sizeof(double)))==NULL)
		return P;
	
	P.TYPE=M;
	// compute temporary standard temperatures and irradiances
	// as the mean temperature and irradiance
	Ts=0;
	Gs=0;
	for (i=0;i<N;i++)
	{
		Ts+=T[i];
		Gs+=G[i];
		
		arr[i]=C[i].Voc;
		arr[N+i]=C[i].Isc;
		arr[2*N+i]=C[i].Vmpp;
		arr[3*N+i]=C[i].Impp;
	}
	Ts/=(double)N;
	Gs/=(double)N;	
	//printf("Tstd,Gstd: %f, %f\n", Ts, Gs);
	
	
	// initial fit using a sequence of OLS-like solutions
	(*alpha)=0;
	(*beta)=0;
	(*a)=0;
	(*Rs)=0;
	(*kappa)=0;
	
	if (FitIsc(arr+N, G, T, N, Gs, Ts, &(P.pars[1]), alpha))
		goto IECFITFINISH; // give up
	if (FitVoc(arr, G, T, N, Gs, Ts, &(P.pars[0]), beta, a))
		goto IECFITFINISH; // give up
	if(FitMpp(arr+2*N, arr+3*N, arr, G, T, N, Gs, Ts, *alpha, *beta, *a, &(P.pars[3]), &(P.pars[2]), Rs, kappa))
		goto IECFITFINISH; // give up
		
	// compute the STC curve
	for (i=0;i<N;i++)
	{
		/* determine std conditions by translating every curve to std and determine mean parameter values */
		P2.TYPE=M;
		P2.pars[0]=C[i].Voc;
		P2.pars[1]=C[i].Isc;
		P2.pars[2]=C[i].Vmpp;
		P2.pars[3]=C[i].Impp;
		P2.pars[4]=C[i].Vqmp1;
		P2.pars[5]=C[i].Iqmp1;
		P2.pars[6]=C[i].Vqmp2;
		P2.pars[7]=C[i].Iqmp2;
		P2.pars[8]=C[i].Gsc;
		P2.pars[9]=C[i].Roc;
		P2=IEC60891_v2(P2, *alpha, *beta, *a, *Rs, *kappa, T[i], Ts, G[i], Gs, S);
		// all parameters beyond the standard solpars are not initialized yet
		P.pars[4]+=P2.pars[4];
		P.pars[5]+=P2.pars[5];
		P.pars[6]+=P2.pars[6];
		P.pars[7]+=P2.pars[7];
		P.pars[8]+=P2.pars[8];
		P.pars[9]+=P2.pars[9];
	}
	P.pars[4]/=(double)N;
	P.pars[5]/=(double)N;
	P.pars[6]/=(double)N;
	P.pars[7]/=(double)N;
	P.pars[8]/=(double)N;
	P.pars[9]/=(double)N;
	
	Cstd.Voc=P.pars[0];
	Cstd.Isc=P.pars[1];
	Cstd.Vmpp=P.pars[2];
	Cstd.Impp=P.pars[3];
	Cstd.Vqmp1=P.pars[4];
	Cstd.Iqmp1=P.pars[5];
	Cstd.Vqmp2=P.pars[6];
	Cstd.Iqmp2=P.pars[7];
	Cstd.Gsc=P.pars[8];
	Cstd.Roc=P.pars[9];
	
	
	// optimize solution using a downhill simplex
	OptimizeIEC_v2(C, G, T, N, Ts, Gs, &Cstd, alpha, beta, a, Rs, kappa, M, S);
	
	// rebase the std condition from Ts,Gs to Tstd,Gstd
	P.pars[0]=0;
	P.pars[1]=0;
	P.pars[2]=0;
	P.pars[3]=0;
	P.pars[4]=0;
	P.pars[5]=0;
	P.pars[6]=0;
	P.pars[7]=0;
	P.pars[8]=0;
	P.pars[9]=0;
	for (i=0;i<N;i++)
	{
		/* determine std conditions by translating every curve to std and determine mean parameter values */
		P2.TYPE=M;
		P2.pars[0]=C[i].Voc;
		P2.pars[1]=C[i].Isc;
		P2.pars[2]=C[i].Vmpp;
		P2.pars[3]=C[i].Impp;
		P2.pars[4]=C[i].Vqmp1;
		P2.pars[5]=C[i].Iqmp1;
		P2.pars[6]=C[i].Vqmp2;
		P2.pars[7]=C[i].Iqmp2;
		P2.pars[8]=C[i].Gsc;
		P2.pars[9]=C[i].Roc;
		P2=IEC60891_v2(P2, *alpha, *beta, *a, *Rs, *kappa, T[i], Tstd, G[i], Gstd, S);
		// all parameters beyond the standard solpars are not initialized yet
		P.pars[0]+=P2.pars[0];
		P.pars[1]+=P2.pars[1];
		P.pars[2]+=P2.pars[2];
		P.pars[3]+=P2.pars[3];
		P.pars[4]+=P2.pars[4];
		P.pars[5]+=P2.pars[5];
		P.pars[6]+=P2.pars[6];
		P.pars[7]+=P2.pars[7];
		P.pars[8]+=P2.pars[8];
		P.pars[9]+=P2.pars[9];
	}
	P.pars[0]/=(double)N;
	P.pars[1]/=(double)N;
	P.pars[2]/=(double)N;
	P.pars[3]/=(double)N;
	P.pars[4]/=(double)N;
	P.pars[5]/=(double)N;
	P.pars[6]/=(double)N;
	P.pars[7]/=(double)N;
	P.pars[8]/=(double)N;
	P.pars[9]/=(double)N;
	
	
IECFITFINISH:
	free(arr);
	return P;
}


#undef FX0
#undef IFX0
#undef FX1
#undef IFX1
#undef FX2
#undef IFX2
#undef FX3
#undef IFX3
#undef FX4
#undef IFX4

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

/* optimization vector x preprocessing macros
 * x[0]: alpha
 * x[1]: beta
 * x[2]: B1
 * x[3]: B2
 * x[4]: Rs
 * x[5]: kappa
 */

#define FX0(a)  (PRE_NONE(a))
#define IFX0(a) (IPRE_NONE(a))

#define FX1(a)  (PRE_NONE(a))
#define IFX1(a) (IPRE_NONE(a))

#define FX2(a)  (PRE_NONE(a))
#define IFX2(a) (IPRE_NONE(a))

#define FX3(a)  (PRE_NONE(a))
#define IFX3(a) (IPRE_NONE(a))

#define FX4(a)  (PRE_ABS(a))
#define IFX4(a) (IPRE_ABS(a))

#define FX5(a)  (PRE_NONE(a))
#define IFX5(a) (IPRE_NONE(a))

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

double FitErr2(double x[6], void *args)
{
	double E=0, e;
	errpars *pp;
	int i;
	generic_pars P2;
	pp=(errpars *) args;
	for (i=0;i<pp->N;i++)
	{
		P2.TYPE=pp->M;
		P2.pars[0]=pp->C[i].Voc;
		P2.pars[1]=pp->C[i].Isc;
		P2.pars[2]=pp->C[i].Vmpp;
		P2.pars[3]=pp->C[i].Impp;
		P2.pars[4]=pp->C[i].Vqmp1;
		P2.pars[5]=pp->C[i].Iqmp1;
		P2.pars[6]=pp->C[i].Vqmp2;
		P2.pars[7]=pp->C[i].Iqmp2;
		P2.pars[8]=pp->C[i].Gsc;
		P2.pars[9]=pp->C[i].Roc;
		P2=IEC60891_v3(P2, FX0(x[0]), FX1(x[1]), FX2(x[2]), FX3(x[3]), FX4(x[4]), FX5(x[5]), pp->P.P[4].x, pp->T[i], pp->Tstd, pp->Tstd, pp->G[i], pp->Gstd, pp->Gstd, pp->S); // from actual to std
		e=IEC_RMSErr(P2, pp->P, pp->S);
		if (!isfinite(e))
			e=1;
		E+=e;
	}

	//printf("%e %e %e %e %e %e %e\n", x[0], x[1], x[2], x[3], x[4], x[5], E);
	return E;
}

// part errors
// only Voc
double FitErrV(double x[6], void *args)
{
	double E=0, f, e;
	errpars *pp;
	int i;
	pp=(errpars *) args;
	for (i=0;i<pp->N;i++)
	{
		f=log(pp->Gstd/pp->G[i]);
		f=1+FX2(x[1])*f+FX3(x[2])*f*f;
		e=pp->P.P[4].x-pp->C[i].Voc*(1/f+FX1(x[0])*(pp->Tstd - pp->T[i])*f);
		E+=e*e;
		//printf("%e \n", e);
	}
	//printf("%e \n\n", E);
	return E;
}

// only Mpp
double FitErrM(double x[6], void *args)
{
	double E=0;
	errpars *pp;
	int i;
	generic_pars P2;
	pp=(errpars *) args;
	for (i=0;i<pp->N;i++)
	{
		P2.TYPE=pp->M;
		P2.pars[0]=pp->C[i].Voc;
		P2.pars[1]=pp->C[i].Isc;
		P2.pars[2]=pp->C[i].Vmpp;
		P2.pars[3]=pp->C[i].Impp;
		P2.pars[4]=pp->C[i].Vqmp1;
		P2.pars[5]=pp->C[i].Iqmp1;
		P2.pars[6]=pp->C[i].Vqmp2;
		P2.pars[7]=pp->C[i].Iqmp2;
		P2.pars[8]=pp->C[i].Gsc;
		P2.pars[9]=pp->C[i].Roc;
		P2=IEC60891_v3(P2, FX0(x[-4]), FX1(x[-3]), FX2(x[-2]), FX3(x[-1]), FX4(x[0]), FX5(x[1]), pp->P.P[4].x, pp->T[i], pp->Tstd, pp->Tstd, pp->G[i], pp->Gstd, pp->Gstd, pp->S); // from actual to std
	
		E+=IEC_RMSErr(P2, pp->P, pp->S);	
	}
	//printf("%e %e %e %e %e %e %e\n", x[0], x[1], x[2], x[3], x[4], x[5], E);
	return E;
}
// only Isc
double FitErrI(double x[6], void *args)
{
	double E=0, e;
	errpars *pp;
	int i;
	pp=(errpars *) args;
	for (i=0;i<pp->N;i++)
	{
		e=pp->P.P[0].y-pp->C[i].Isc*pp->Gstd/(pp->G[i]*(1+FX1(x[1])*(pp->T[i]-pp->Tstd)));
		E+=e*e;
		//printf("%e \n", e);
	}
	//printf("%e \n\n", E);
	return E;
}

void OptimizeIEC_v3(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, curvepars *Cstd, double *alpha, double *beta, double *B1, double *B2, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S)
{
	double x[6], x0[6];
	double step[6];
	double xmin[6];
	int NCALLS, NRESTART, Err;	
	double EPS2=S->fiteps, E,e, E0;
	errpars pp;
	
	
	//printf("Settings: eps=%e N=%d\n", S->fiteps, S->fititer);
	/* square EPS to save some sqrts */
	EPS2*=EPS2;	
	
	/* fill the x vector */
	x[0]=x0[0]=xmin[0]=IFX0(*alpha);
	x[1]=x0[1]=xmin[1]=IFX1(*beta);
	x[2]=x0[2]=xmin[2]=IFX2(*B1);// initially should be equal to a in the version 2 norm!
	x[3]=x0[3]=xmin[3]=IFX3(*B1);
	x[4]=x0[4]=xmin[4]=IFX4(*Rs);
	x[5]=x0[5]=xmin[5]=IFX5(*kappa);
	/*
	printf("alpha:   %e\n", FX0(xmin[0]));
	printf("beta :   %e\n", FX1(xmin[1]));
	printf("B1   :   %e\n", FX2(xmin[2]));
	printf("B2   :   %e\n", FX3(xmin[3]));
	printf("Rs   :   %e\n", FX4(xmin[4]));
	printf("kappa:   %e\n", FX5(xmin[5]));*/
	step[0]=SIMPLSIZE*x[0]+MINSIMPLSIZE;
	step[1]=SIMPLSIZE*x[1]+MINSIMPLSIZE;
	step[2]=SIMPLSIZE*x[2]+MINSIMPLSIZE;
	step[3]=SIMPLSIZE*x[2]+MINSIMPLSIZE; // note size related to x2
	step[4]=SIMPLSIZE*x[4]+MINSIMPLSIZE;
	step[5]=SIMPLSIZE*x[5]+MINSIMPLSIZE;
	
	pp.P.N=5;
	if((pp.P.P=Emalloc(5*sizeof(point)))==NULL)
	{
		return;
	}	
	pp.P.P[0].x=0;
	pp.P.P[0].y=Cstd->Isc;
	pp.P.P[1].x=Cstd->Vqmp1;
	pp.P.P[1].y=Cstd->Iqmp1;
	pp.P.P[2].x=Cstd->Vmpp;
	pp.P.P[2].y=Cstd->Impp;
	pp.P.P[3].x=Cstd->Vqmp2;
	pp.P.P[3].y=Cstd->Iqmp2;
	pp.P.P[4].x=Cstd->Voc;
	pp.P.P[4].y=0;
	pp.C=C;
	pp.G=G;
	pp.T=T;
	pp.N=N;
	pp.Tstd=Tstd;
	pp.Gstd=Gstd;
	pp.M=M;
	pp.S=S;
	
	
	E=FitErr2(x, (void *) (&pp));
	E0=E;
	//printf("Initial Error: %e\n", E/(double)N);	
	/* only fit Voc
	nelmin (FitErrV, (void *) (&pp), 3, x+1, xmin+1, &E, EPS2, step+1, 10, S->fititer, &NCALLS,&NRESTART, &Err);
	e=FitErr2(xmin, (void *) (&pp));
	if (e<E0)
	{
		E0=e;
		x0[3]=xmin[3];
		x0[2]=xmin[2];
		x0[1]=xmin[1];	
	}
	printf("OC Error: %e\n", E0/(double)N);
	e=FitErr2(x0, (void *) (&pp));
	printf("Error0: %e\n", e/(double)N);
	
	printf("beta :   %e\n", FX1(x0[1]));
	printf("B1   :   %e\n", FX2(x0[2]));
	printf("B2   :   %e\n", FX3(x0[3]));
	x[0]=xmin[0]=x0[0];
	x[1]=xmin[1]=x0[1];
	x[2]=xmin[2]=x0[2];
	x[3]=xmin[3]=x0[3];
	x[4]=xmin[4]=x0[4];
	x[5]=xmin[5]=x0[5]; */
	/* fit Rs and Kappa 
	nelmin (FitErrM, (void *) (&pp), 2, x+4, xmin+4, &E, EPS2, step+4, 10, S->fititer, &NCALLS,&NRESTART, &Err);
	e=FitErr2(xmin, (void *) (&pp));
	if (e<E0)
	{
		E0=e;
		x0[4]=xmin[4];
		x0[5]=xmin[5];
	}
	e=FitErr2(x0, (void *) (&pp));
	printf("MP Error: %e\n", E0/(double)N);
	printf("Error0: %e\n", e/(double)N);
	printf("Rs   :   %e\n", FX1(x0[4]));
	printf("kappa:   %e\n", FX2(x0[5]));
	
	x[0]=xmin[0]=x0[0];
	x[1]=xmin[1]=x0[1];
	x[2]=xmin[2]=x0[2];
	x[3]=xmin[3]=x0[3];
	x[4]=xmin[4]=x0[4];
	x[5]=xmin[5]=x0[5]; */
	
	// fit all
	nelmin (FitErr2, (void *) (&pp), 6, x, xmin, &E, EPS2, step, 10, S->fititer, &NCALLS,&NRESTART, &Err);
	e=FitErr2(xmin, (void *) (&pp));
	if (e<E0)
	{
		x0[0]=xmin[0];
		x0[1]=xmin[1];
		x0[2]=xmin[2];
		x0[3]=xmin[3];
		x0[4]=xmin[4];
		x0[5]=xmin[5];
	}
	free(pp.P.P);
	/*printf("alpha:   %e\n", FX0(x0[0]));
	printf("beta :   %e\n", FX1(x0[1]));
	printf("B1   :   %e\n", FX2(x0[2]));
	printf("B2   :   %e\n", FX3(x0[3]));
	printf("Rs   :   %e\n", FX4(x0[4]));
	printf("kappa:   %e\n", FX5(x0[5])); */
	(*alpha)= FX0(x0[0]);
	(*beta) = FX1(x0[1]);
	(*B1)   = FX2(x0[2]);
	(*B2)   = FX3(x0[3]);
	(*Rs)   = FX4(x0[4]);
	(*kappa)= FX5(x0[5]);
}


generic_pars IEC60891Fit_v3(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double *alpha, double *beta, double *B1, double *B2, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S)
{
	int i;
	double *arr, Ts=0, Gs=0;
	generic_pars P={0}, P2;
	curvepars Cstd;
	 
	/* we start with the version 2 of the IEC60891, it is fast and easily
	 * fitted with 3 OLS-like fits
	 */
	 
	if((arr=Emalloc(4*N*sizeof(double)))==NULL)
		return P;
	P.TYPE=M;
	
	// compute temporary std G and T
	for (i=0;i<N;i++)
	{
		Ts+=T[i];
		Gs+=G[i];
		arr[i]=C[i].Voc;
		arr[N+i]=C[i].Isc;
		arr[2*N+i]=C[i].Vmpp;
		arr[3*N+i]=C[i].Impp;
	}
	Ts/=(double)N;
	Gs/=(double)N;
	
	// initial fit
	(*alpha)=0;
	(*beta)=0;
	(*B1)=0;
	(*B2)=0;
	(*Rs)=0;
	(*kappa)=0;
	
	if (FitIsc(arr+N, G, T, N, Gs, Ts, &(P.pars[1]), alpha))
		goto IECFITFINISH; // give up
	if (FitVoc(arr, G, T, N, Gs, Ts, &(P.pars[0]), beta, B1)) // initially we put the value of a in B1
		goto IECFITFINISH; // give up
	if(FitMpp(arr+2*N, arr+3*N, arr, G, T, N, Gs, Ts, *alpha, *beta, *B1, &(P.pars[3]), &(P.pars[2]), Rs, kappa))
		goto IECFITFINISH; // give up
		
	// compute STC curve
	for (i=0;i<N;i++)
	{
		/* determine std conditions by translating every curve to std and determine mean parameter values */
		P2.TYPE=M;
		P2.pars[0]=C[i].Voc;
		P2.pars[1]=C[i].Isc;
		P2.pars[2]=C[i].Vmpp;
		P2.pars[3]=C[i].Impp;
		P2.pars[4]=C[i].Vqmp1;
		P2.pars[5]=C[i].Iqmp1;
		P2.pars[6]=C[i].Vqmp2;
		P2.pars[7]=C[i].Iqmp2;
		P2.pars[8]=C[i].Gsc;
		P2.pars[9]=C[i].Roc;
		P2=IEC60891_v2(P2, *alpha, *beta, *B1, *Rs, *kappa, T[i], Ts, G[i], Gs, S);
		// all parameters beyond the standard solpars are not initialized yet
		P.pars[4]+=P2.pars[4];
		P.pars[5]+=P2.pars[5];
		P.pars[6]+=P2.pars[6];
		P.pars[7]+=P2.pars[7];
		P.pars[8]+=P2.pars[8];
		P.pars[9]+=P2.pars[9];
	}
	P.pars[4]/=(double)N;
	P.pars[5]/=(double)N;
	P.pars[6]/=(double)N;
	P.pars[7]/=(double)N;
	P.pars[8]/=(double)N;
	P.pars[9]/=(double)N;
	
	Cstd.Voc=P.pars[0];
	Cstd.Isc=P.pars[1];
	Cstd.Vmpp=P.pars[2];
	Cstd.Impp=P.pars[3];
	Cstd.Vqmp1=P.pars[4];
	Cstd.Iqmp1=P.pars[5];
	Cstd.Vqmp2=P.pars[6];
	Cstd.Iqmp2=P.pars[7];
	Cstd.Gsc=P.pars[8];
	Cstd.Roc=P.pars[9];
	
	// optimize using a downhill simplex
	OptimizeIEC_v3(C, G, T, N, Ts, Gs, &Cstd, alpha, beta, B1, B2, Rs, kappa, M, S);
	
	// rebase the std condition from Ts,Gs to Tstd,Gstd
	P=IEC60891_v3(P, *alpha, *beta, *B1, *B2, *Rs, *kappa, Cstd.Voc, Ts, Tstd, Ts, Gs, Gstd, Gs, S);
	
IECFITFINISH:
	free(arr);
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
