/*
    Cubic Bézier IV model.
    Copyright (C) 2021  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2021                                   *
 *                                                               *             
 *****************************************************************/ 
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "splinefit.h"
#include "curve.h"
#include "fitmodels.h"
#include "curvemodel.h"
#include "spline_iv.h"
// internal struct
typedef struct {
	double V[13];
	double I[13];
} B3pars;

// Exported Function
// return curvepars according to the model fit.
curvepars B3Pars(generic_pars P, MODELSET *S)
{
	curvepars par;
	
	if (P.TYPE!=M_PVBEZIER)
	{
		fprintf(stderr,"Error: wrong parameter type in B3Pars\n");
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

// helper function to construct the Bézier
B3pars InitFitBezier(B3pars B, double dIsc, double dIqmp1, double dImpp, double dIqmp2, double dIoc, MODELSET *S);
// convert general to internal parameter struct
B3pars Generic2B3Pars(generic_pars P, MODELSET *S)
{
	B3pars Sd;
	double dIsc, dIoc, dIqmp1, dImpp, dIqmp2;

	if (P.TYPE!=M_PVBEZIER)
	{
		int i;
		fprintf(stderr,"Error: wrong parameter type in Generic2B3Pars\n");
		for (i=0;i<13;i++)
		{
			Sd.V[i]=0;
			Sd.I[i]=0;
		}
		return Sd;
	}
	/* procedure EPAR 2 Bézier 
	 * 1 populate the Bézier arrays sparcely, using only the knots at sc mpp and oc
	 * 2 Fit controll points for the sparse Bézier curve
	 * 3 Split the sparse Bézier to obtain all the knots
	 * 4 Correct the quasi mpp points
	 */
	// initialize sparcely populated Bézier arrays (step 1)
	Sd.V[0] =0;
	Sd.V[3] =P.pars[4];
	Sd.V[6] =P.pars[2];
	Sd.V[9] =P.pars[6];
	Sd.V[12]=P.pars[0];
	
	Sd.I[0] =P.pars[1];
	Sd.I[3] =P.pars[5];
	Sd.I[6] =P.pars[3];
	Sd.I[9] =P.pars[7];
	Sd.I[12]=0;

	dIsc  =-P.pars[8];
	dIqmp1=-(P.pars[5]-P.pars[3])/P.pars[4];
	dImpp =-P.pars[3]/P.pars[2];
	dIqmp2=-P.pars[7]/(P.pars[6]-P.pars[2]);
	dIoc  =-1/P.pars[9];
	// Make the full parametrization (steps 2-4)
	Sd=InitFitBezier(Sd, dIsc, dIqmp1, dImpp, dIqmp2, dIoc,S);
	/*
	printf("--\n");
	for(i=0;i<13;i++)
		printf("%e %e\n", S.V[i],S.I[i]);
		
	printf("--\n");*/
	return Sd;
}
	

// convert internal to general parameter struct
generic_pars B32GenericPars(B3pars S)
{
	generic_pars P;
	P.TYPE=M_PVBEZIER;
	P.pars[0]=S.V[12];
	P.pars[1]=S.I[0];
	P.pars[2]=S.V[6];
	P.pars[3]=S.I[6];
	P.pars[4]=S.V[3];
	P.pars[5]=S.I[3];
	P.pars[6]=S.V[9];
	P.pars[7]=S.I[9];
	P.pars[8]=(S.I[0]-S.I[1])/S.V[1];
	P.pars[9]=(S.V[12]-S.V[11])/S.I[11];
	return P;
}	

/* Bézier Base Functions */
#define B30(t) ((1-(t))*(1-(t))*(1-(t)))
#define B31(t) (3*(t)*(1-(t))*(1-(t)))
#define B32(t) (3*(t)*(t)*(1-(t)))
#define B33(t) ((t)*(t)*(t))

#define DB30(t) (3*(1-(t))*(1-(t)))
#define DB31(t) (6*(t)*(1-(t)))
#define DB32(t) (3*(t)*(t))

double B3_coord(double t, double points[])
{
	return B30(t)*points[0]+B31(t)*points[1]+B32(t)*points[2]+B33(t)*points[3];
} 

double B3_velocity(double t, double points[])
{
	return DB30(t)*(points[1]-points[0])+DB31(t)*(points[2]-points[1])+DB32(t)*(points[3]-points[2]);
} 

// pars are B3pars
point B3_i(double t, point *dp, void * pars)
{
	int i;
	point p;
	B3pars par;
	par=*((B3pars *)pars);
	
	
	
	// for t from 0-1 we go through 4 sections, each with a t' from 0 to 1
	t*=4.0;
	i=0;
	i+=(t>1.0);
	i+=(t>2.0);
	i+=(t>3.0);
	// i now points to the knot index before the interval
	t=t-(double)i;
	// compute the time for this interval
	
	if ((t>=0)&&(t<=1))
	{
		p.x=B3_coord(t, par.V+3*i);
		p.y=B3_coord(t, par.I+3*i);
		
		if (dp)
		{
			dp->x=B3_velocity(t, par.V+3*i);
			dp->y=B3_velocity(t, par.I+3*i);
		}
		return p;
	}
	
	// linear extrapolation
	if ((t<0)&&(i==0))
	{
		p.x=par.V[0]+t*(par.V[1]-par.V[0]);
		p.y=par.I[0]+t*(par.I[1]-par.I[0]);;
	
		if (dp)
		{
			dp->x=(par.V[1]-par.V[0]);
			dp->y=(par.I[1]-par.I[0]);
		}
	}
	else if((t>1)&&(i==3))
	{
		p.x=par.V[12]+t*(par.V[12]-par.V[11]);
		p.y=par.I[12]+t*(par.I[12]-par.I[11]);;
	
		if (dp)
		{
			dp->x=(par.V[12]-par.V[11]);
			dp->y=(par.I[12]-par.I[11]);
		}
	}
	else
	{
		// cannot happen
		p.x=0;
		p.y=0;		
		if (dp)
		{
			dp->x=1;
			dp->y=1;
		}
		
	}
				
	return p;
}

// Exported Function
curve B3ToCurve(generic_pars P, int N, MODELSET *S)
{
	B3pars par;
	curve C;
	par=Generic2B3Pars(P,S);
	C=ModelToCurve((void *) &par, &B3_i, N,S);
	return C;
}

double B3RMSErr_i(B3pars P, curve C, MODELSET *S)
{
	return RMSErr(P.V[12], P.I[0], (void *) &P, &B3_i, C,S);
}

// Exported Function
double B3RMSErr(generic_pars P, curve C, MODELSET *S)
{
	B3pars par;
	double E;
	par=Generic2B3Pars(P,S);
	E=B3RMSErr_i(par, C,S);
	return E;
}

curve B3DistCurve_i(B3pars P, curve C, MODELSET *S)
{
	return DistCurve(P.V[12], P.I[0], (void *) &P, &B3_i, C,S);
}

// Exported Function
curve B3DistCurve(generic_pars P, curve C, MODELSET *S)
{
	B3pars par;
	curve d;
	par=Generic2B3Pars(P,S);
	d=B3DistCurve_i(par, C,S);
	return d;
}


// Exported Function
double B3IV(double V, generic_pars P, MODELSET *S)
{	
	B3pars par;
	par=Generic2B3Pars(P,S);
	return Solve_IV(V,P.pars[0],NULL,(void *) &par, &B3_i, S);	
}

// Exported Function
double B3VI(double I, generic_pars P, MODELSET *S)
{	
	B3pars par;
	par=Generic2B3Pars(P,S);
	return Solve_VI(I,P.pars[0],NULL,(void *) &par, &B3_i, S);	
}


/* 
 * We solve controll points such that we cross the coordinate (X,Y) at t=t
 * Startpoint is (x0,y0) endpoint is (x1,y2)
 * Start derivatives (dx0,dy0) End derivatives (dy1,dy1)
 * Controll coordinates are:
 * (x0+a*dx0,y0+a*dy0) and (x1-b*dy1, y1-b*dy1)
 * 
 * We can solve in maxima with:
 * B(t, p0, p1, d0, d1):=(1-t)^3*p0+3*(1-t)^2*t*(p0+a*d0)+3*(1-t)*t^2*(p1-b*d1)+t^3*p1;
 * eq1: B(t,x0,x1, dx0, dx1)=X;
 * eq2: B(t,y0,y1, dy0, dy1)=Y;
 * solve([eq1,eq2],[a,b]);
 * 
 * I define D as the cross product of the two derivative vectors:
 * D=(dx1*dy0-dx0*dy1)
 * 
 * Note the solution requires that the cross product of the derivatives is not 0, i.e. the derivatives must span the 2D space
 * 
 * We then find: 
 * a = (dx1*(2*t^3*(y1-y0)-3*t^2*(y1-y0)+(Y-y0))-dy1*(2*t^3*(x1-x0)-3*t^2*(x1-x0)+(X-x0))/(3*D*(t^3-2*t^2+t))
 * b = (dx0*(2*t^3*(y0-y1)-3*t^2*(y0-y1)+(y0-Y))+dy0*(2*t^3*(x1-x0)-3*t^2*(x1-x0)+(X-x0)))/(3*D*(t^3-t^2))
 */
void ControllAB(double t, double x0, double x1, double y0, double y1, double dx0, double dx1, double dy0, double dy1, double X, double Y, double *a, double *b)
{
	double D=(dx1*dy0-dx0*dy1);
	double t3=t*t*t, t2=t*t;
	(*a)=(dx1*(2*t3*(y1-y0)-3*t2*(y1-y0)+(Y-y0))-dy1*(2*t3*(x1-x0)-3*t2*(x1-x0)+(X-x0)))/(3*D*(t3-2*t2+t));
	(*b)=(dx0*(2*t3*(y0-y1)-3*t2*(y0-y1)+(y0-Y))+dy0*(2*t3*(x1-x0)-3*t2*(x1-x0)+(X-x0)))/(3*D*(t3-t2));
}

void SetControll(double t, double *x, double *y, double dy0, double dy1)
{
	// points is bezier section with length 7 as not all elements are used (yet)
	double a,b;
	ControllAB(t, x[0], x[6], y[0], y[6], 1.0, 1.0, dy0, dy1, x[3], y[3], &a, &b);
	a=fabs(a);
	x[1]=x[0]+a;
	y[1]=y[0]+a*dy0;
	b=fabs(b);
	x[5]=x[6]-b;
	y[5]=y[6]-b*dy1;
}

// before ther bezier is fully parameterized we use a sparse bezier array
// meaning we do not use all points yet, for now we do not touch point indexes 2-4 of 6
// i.2. an arry with 6 points is used for a single cubic Bézier section
double B3_init_velocity(double t, double points[])
{
	return DB30(t)*(points[1]-points[0])+DB31(t)*(points[5]-points[1])+DB32(t)*(points[6]-points[5]);
} 
double B3_init_coord(double t, double points[])
{
	return B30(t)*points[0]+B31(t)*points[1]+B32(t)*points[5]+B33(t)*points[6];
} 



/* check whether a bezier section is monotonic
 * The function takes the 4 coordinates of a bezier curve. It takes only one dimension at a time
 * i.e. the 2D curve is monotonic if both x(t) and y(t) are and requires two calls to this function
 * 
 * The optional Error argument E computes a measure of how non-monotonic the curve is (0 if the curve is monotonic)
 * Note, for a meaningful error value one should normalize the x and y coordinates when you check a 2D curve */
int Monotone(double x0, double x1, double x2, double x3, double *E)
{
	double t, t2, t3, it, it2, it3, num, denom, d0, d1, d;
	double xinfl;
	int nr, NR=0;
	if ((x1>x0)&&(x1<x2)&&(x2<x3))// check coordinates are sorted
		return 1;
	d0=x1-x0;
	d1=x3-x2;
	denom=x3-3*x2+3*x1-x0;
	num= -(x2-2*x1+x0);
	if (E)
		(*E)=0;		
	if ((d0*d1>=0)&&((fabs(denom)<fabs(num))||(denom*num<0))) // slopes at t=0 and t=1 of same sign and no inflection point on 0<=t<=1
		return 1;
		
	if (fabs(denom)<fabs(num)||(denom*num<0))
	{
		// no inflection point but slopes at start and end of section of unequal sign, there must be exactly one extremum on interval 0<=t<=1
		if (!E)
			return 0;		
		xinfl=(x0+x3)/2;
		NR=1; // we have exactly one root
	}
	else
	{
		// there is an inflection point and slopes at t=0 and t=1 are of same sign or 0
		t=num/denom;
		t2=t*t;
		t3=t2*t;
		it=(1-t);
		it2=it*it;
		it3=it*it2;		
		xinfl=it3*x0+3*t*it2*x1+3*t2*it*x2+t3*x3;
		
		d=3*t2*x3-3*t2*x2+6*it*t*x2-6*it*t*x1+3*it2*x1-3*it2*x0;		
		if ((d*d0>=0)&&(d*d1>=0)) // monotonic
			return 1;
		if (!E)			// not monotonic
			return 0;	// fixme, If we do not compute the error we do not consider cases where numerical noise falsly makes us think we are non-monotonic
		NR=2; // we have 1 or two roots
	}
	nr=0;	
	// find maxima on interval and determine error
	//t=((+/-)sqrt((y0-y1)*y3+y2^2+((-y1)-y0)*y2+y1^2)-y2+2*y1-y0)/(y3-3*y2+3*y1-y0)
	d=(x0-x1)*x3+x2*x2-(x1+x0)*x2+x1*x1;
	if (d<0)
	{
		// impossible error, there is no real root of the derivative and thus no extremum
		// we are probably facing numerical noise, let us give the curve the benefit of the doubt
		return 1;
	}
	num=sqrt(d)-x2+2*x1-x0; // try first root
	if (fabs(denom)>=fabs(num)&&(denom*num>=0)) // we have found a root on the interval
	{ 
		t=num/denom;
		t2=t*t;
		t3=t2*t;
		it=(1-t);
		it2=it*it;
		it3=it*it2;
		(*E)+=fabs(it3*x0+3*t*it2*x1+3*t2*it*x2+t3*x3-xinfl);
		nr=1;
	}
	num=-sqrt(d)-x2+2*x1-x0;
	if (fabs(denom)>=fabs(num)&&(denom*num>=0)) // we have found a root on the interval
	{ 
		t=num/denom;
		t2=t*t;
		t3=t2*t;
		it=(1-t);
		it2=it*it;
		it3=it*it2;
		(*E)+=fabs(it3*x0+3*t*it2*x1+3*t2*it*x2+t3*x3-xinfl);
		nr++;
	}
	/*if ((NR==1)&&(nr!=1))
	{ 
		// impossible error
		//fprintf(stderr,"Warning: found %d extrema where there should be one\n", nr);
		(*E)+=0;
	}*/
	if ((NR==2)&&(nr==0))
	{ 
		// impossible error
		//fprintf(stderr,"Warning: found no extremum where there should be at least one\n");
		return 1; // assume numerical errors lead to this, pretend all is fine
	}
	return (nr==0); /* monotonic if no maxima are present */	
}


#define CEPS 1e-6 // this is a bit ugly, This value puts more weight on hitting the correct coordinate than to get the derivative right
				  // I hope the exact value is not all that important as we do a correction anyway afterwards but it is a problem when the coordinate shifts alot
				  // it makes for funky interpollations
double FitControll(double *x, double *y, double dy0, double dy1, double dyc, MODELSET *S)
{

	double tmin=0, topt, tmax=1, t, Eopt, E, Em1=0, Em2=0;
	double dy,dx, xx, yy;
	int MONO, MONO_OPT=0;
	
	
	t=(tmin+tmax)/2;
	SetControll(t, x, y, dy0, dy1);
	dx=B3_init_velocity(t, x);
	dy=B3_init_velocity(t, y);
	if ((Monotone(x[0], x[1], x[5], x[6], &Em1))&&(Monotone(y[0], y[1], y[5], y[6], &Em2)))
	{
		MONO_OPT=1;
		xx=B3_init_coord(t, x);
		yy=B3_init_coord(t, y);
		
		E=(dy/dx-dyc)*(dy/dx-dyc)+(xx-x[3])*(xx-x[3])+(yy-y[3])*(yy-y[3])/CEPS;
	}
	else
	{
		xx=B3_init_coord(t, x);
		yy=B3_init_coord(t, y);
		E=Em1+Em2;
		E=E*E+(xx-x[3])*(xx-x[3])+(yy-y[3])*(yy-y[3])/CEPS;
	}
	topt=t;
	Eopt=E;
	while ((tmax-tmin>S->fiteps))
	{
		// alternatingly we bi-sect tmin-topt and topt-tmax
		if (topt-tmin>tmax-topt)
			t=(topt+tmin)/2;
		else
			t=(topt+tmax)/2;
		
		SetControll(t, x, y, dy0, dy1);
		dx=B3_init_velocity(t, x);
		dy=B3_init_velocity(t, y);
		MONO=((Monotone(x[0], x[1], x[5], x[6], &Em1))&&(Monotone(y[0], y[1], y[5], y[6], &Em2)));
		
		if (MONO)
		{
			xx=B3_init_coord(t, x);
			yy=B3_init_coord(t, y);
			E=(dy/dx-dyc)*(dy/dx-dyc)+(xx-x[3])*(xx-x[3])+(yy-y[3])*(yy-y[3])/CEPS;
		}
		else
		{
			xx=B3_init_coord(t, x);
			yy=B3_init_coord(t, y);
			E=Em1+Em2;
			E=E*E+(xx-x[3])*(xx-x[3])+(yy-y[3])*(yy-y[3])/CEPS;
		}
		if (MONO && (!MONO_OPT))
		{
			MONO_OPT=MONO;
			if (t<topt)
			{
				tmax=topt;
				topt=t;
				Eopt=E;
			}
			else
			{
				tmin=topt;
				topt=t;
				Eopt=E;
			}
		}
		else if (!MONO && MONO_OPT)
		{
			if (t<topt)
				tmin=t;
			else
				tmax=t;
		}
		else
		{
			if (E<Eopt)
			{
				if (t<topt)
				{
					tmax=topt;
					topt=t;
					Eopt=E;
				}
				else
				{
					tmin=topt;
					topt=t;
					Eopt=E;
				}
			}
			else
			{
				if (t<topt)
					tmin=t;
				else
					tmax=t;
			}
		}
	}

	SetControll(t, x, y, dy0, dy1);
	
	if (!((Monotone(x[0], x[1], x[5], x[6], NULL))&&(Monotone(y[0], y[1], y[5], y[6], NULL))))
	{
		t=0.5;
		SetControll(t, x, y, dy0, dy1);
	}
		
	return t;
}

int ParMono(double x0, double y0, double dy0, double x1, double y1, double dy1)
{
	if (y1>y0)
		return 0;
	if (x0>x1)
		return 0;
	if (dy0*dy1<0)
		return 0;
	return 1;
}

void MonotonizeSection(double *x, double *y)
{
	if (!((Monotone(x[0], x[1], x[2], x[3], NULL))&&(Monotone(y[0], y[1], y[2], y[3], NULL))))
	{
		// only do things when there is a problem
		double tx, ty, t;
		/* fix wrong way lines (derivative I/V is right but V(t)/t is not) */
		if (x[1]<x[0])
		{
			x[1]=2*x[0]-x[1];
			y[1]=2*y[0]-y[1];
		}
		if (x[2]>x[3])
		{
			x[2]=2*x[3]-x[2];
			y[2]=2*y[3]-y[2];
		}
		/* to enforce a monotone curve we make the control points monotone in both x and y coordinate
		 * i.e. x0<x1<x2<x3 and y0>y1>y2>y3		 */
		// parametric curves along start/end-control-point-vector
		tx=(x[3]-x[0])/(x[3]-x[0]+x[1]-x[2]);
		ty=(y[3]-y[0])/(y[3]-y[0]+y[1]-y[2]);
		// limit time such that one controll points does not lie beyond the other in either x or y 
		t=1;
		if (tx<ty)
		{
			if ((tx<1)&&(tx>=0))
				t=tx;
		}
		else
		{
			if ((ty<1)&&(ty>=0))
				t=ty;
		}
		x[1]=x[0]+t*(x[1]-x[0]);
		y[1]=y[0]+t*(y[1]-y[0]);
		x[2]=x[3]+t*(x[2]-x[3]);
		y[2]=y[3]+t*(y[2]-y[3]);
	}
}
/* Keep the controll points of a section inside the triangle spanned by the two
 * end tangents (only if the tangents intersect within the section).
 * This is applied to the sections adjacent to the mpp only: limiting the
 * sections towards sc and oc as well made the fit to measured IV curves
 * noticeably worse while those sections do not cause the mpp kink.
 * The distance of each controll point to its knot is limited to lie between
 * TRI_FMIN and 1 times the distance of the knot to the tangent intersection.
 * The upper limit ensures the section is monotone and does not cross its end
 * tangents, i.e. around the mpp the power does not exceed Pmpp. The lower limit
 * prevents degenerate (zero length) handles which lead to kinks.
 * Only the handle lengths change, the knots and slopes remain as they are. */
#define TRI_FMIN 0.3
void TriangleSection(double *x, double *y, double s0, double s3)
{
	double xi, L0, L3, h0, h3;
	h0=x[1]-x[0];
	h3=x[3]-x[2];
	if (s0==s3)
		return;
	// voltage of the tangent intersection
	xi=((y[3]-s3*x[3])-(y[0]-s0*x[0]))/(s0-s3);
	if ((xi<=x[0])||(xi>=x[3]))
		return;
	L0=xi-x[0];
	L3=x[3]-xi;
	if (h0>L0) h0=L0;
	if (h0<TRI_FMIN*L0) h0=TRI_FMIN*L0;
	if (h3>L3) h3=L3;
	if (h3<TRI_FMIN*L3) h3=TRI_FMIN*L3;
	x[1]=x[0]+h0;
	y[1]=y[0]+h0*s0;
	x[2]=x[3]-h3;
	y[2]=y[3]-h3*s3;
}

/* mpp_at_end: 1 if the section ends at the mpp (the sc->mpp half),
 *             0 if it starts at the mpp (the mpp->oc half) */
void SplitSection(double *x, double *y, double dy0, double dyc, double dy1, double t, int mpp_at_end)
{
	// split bezier section of 4 points in 2 connected sections of 7 points total (one point is shared by both)
	double m1x, m2x, m3x, m1y, m2y, m3y;
	double mm1x, mm2x, mm1y, mm2y, xc, yc;
	double l;
	m1x=(1-t)*x[0]+t*x[1];
	m1y=(1-t)*y[0]+t*y[1];
	
	m2x=(1-t)*x[1]+t*x[5];
	m2y=(1-t)*y[1]+t*y[5];
	
	m3x=(1-t)*x[5]+t*x[6];
	m3y=(1-t)*y[5]+t*y[6];
	
	mm1x=(1-t)*m1x+t*m2x;
	mm1y=(1-t)*m1y+t*m2y;
	mm2x=(1-t)*m2x+t*m3x;
	mm2y=(1-t)*m2y+t*m3y;
	
	/* compute central point */
	xc=B30(t)*x[0] + B31(t)*x[1] + B32(t)*x[5] + B33(t)*x[6];
	yc=B30(t)*y[0] + B31(t)*y[1] + B32(t)*y[5] + B33(t)*y[6];
	
	x[1]=m1x;
	y[1]=m1y;
	x[5]=m3x;
	y[5]=m3y;	

	// compute the controll points with the same length and the gicen dyc
	l=(mm1x-xc)*(mm1x-xc)+(mm1y-yc)*(mm1y-yc);
	x[2]=x[3]-sqrt(l)/sqrt(dyc*dyc+1);
	y[2]=y[3]-dyc*sqrt(l)/sqrt(dyc*dyc+1);
	
	l=(mm2x-xc)*(mm2x-xc)+(mm2y-yc)*(mm2y-yc);
	x[4]=x[3]+sqrt(l)/sqrt(dyc*dyc+1);
	y[4]=y[3]+dyc*sqrt(l)/sqrt(dyc*dyc+1);
	
	// keep the controll points of the sub-section at the mpp within its tangent triangle
	if (mpp_at_end)
		TriangleSection(x+3, y+3, dyc, dy1);
	else
		TriangleSection(x, y, dy0, dyc);
	
	// ensure monotone interpolation
	if (ParMono(x[0], y[0], (y[1]-y[0])/(x[1]-x[0]), x[3], y[3], dyc))
		MonotonizeSection(x, y);
	
	if (ParMono(x[3], y[3], dyc, x[6], y[6],(y[6]-y[5])/(x[6]-x[5])))
		MonotonizeSection(x+3, y+3);
}

// do initial parametrization
// the parameters are at the start a sparcely occupied bezier array where we 
// have 6 points per section of which we first only use the outer 4
B3pars InitFitBezier(B3pars B, double dIsc, double dIqmp1, double dImpp, double dIqmp2, double dIoc, MODELSET *S)
{
	double t1, t2;
	// first we fit the controll points
	t1=FitControll(B.V, B.I, dIsc, dImpp, dIqmp1,S);
	t2=FitControll(B.V+6, B.I+6, dImpp, dIoc, dIqmp2,S);
	SplitSection(B.V, B.I, dIsc, dIqmp1, dImpp, t1, 1);
	SplitSection(B.V+6, B.I+6, dImpp, dIqmp2, dIoc, t2, 0);
	return B;
}


generic_pars B3Fit(curvepars par, curve C, double *Efit, MODELSET *S)
{	
	B3pars Sd;
	generic_pars P;
	
	P.TYPE=M_PVBEZIER;
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
	
	Sd=Generic2B3Pars(P,S);
	(*Efit)=B3RMSErr_i(Sd, C,S);
	return P;
}


/* IEC60891 procedure 2 translation of a PV-Bézier
 * For a given curve the IEC60891 translation maps every point (V1,I1) with an
 * affine map (I2 is proportional to I1 and V2 is linear in V1 and I1).
 * Bézier curves are invariant under affine maps, i.e. mapping the controll
 * points translates the curve exactly. The PV-spline in contrast is a
 * function I(V), which is not closed under the shear of the V-axis by the
 * series resistance term, so there only the knots and slopes are translated. */
static generic_pars B3SolPars(B3pars *B, MODELSET *S)
{
	generic_pars P;
	curvepars par;
	par=SolPar(B->V[12], B->I[0], (void *) B, &B3_i, S);
	if ((B->V[0]>0)&&(B->V[1]>B->V[0]))
	{
		/* The translated curve starts at V>0. The solver only searches
		 * 0<=t<=1, so extrapolate linearly to short circuit here, as B3_i
		 * does for t<0 */
		double s=(B->I[1]-B->I[0])/(B->V[1]-B->V[0]);
		par.Isc=B->I[0]-s*B->V[0];
		par.Gsc=-s;
	}
	P.TYPE=M_PVBEZIER;
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

generic_pars B3IEC60891_v2(generic_pars P, double alpha, double beta, double a, double Rs, double kappa, double T1, double T2, double G1, double G2, MODELSET *S)
{
	B3pars B;
	double Voc1;
	int i;
	B=Generic2B3Pars(P,S);
	Voc1=B.V[12];
	for (i=0;i<13;i++)
		IEC_v2(B.V[i], B.I[i], Voc1, alpha, beta, a, Rs, kappa, T1, T2, G1, G2, B.V+i, B.I+i);
	return B3SolPars(&B, S);
}

generic_pars B3IEC60891_v3(generic_pars P, double alpha, double beta, double B1, double B2, double Rs, double kappa, double rho, double Vocstd, double T1, double T2, double Tstd, double G1, double G2, double Gstd, MODELSET *S)
{
	B3pars B;
	double dIdV;
	int i;
	B=Generic2B3Pars(P,S);
	for (i=0;i<13;i++)
		IEC_v3(B.V[i], Vocstd, B.I[i], 0, G1, G2, Gstd, T1, T2, Tstd, alpha, beta, B1, B2, Rs, kappa, rho, B.V+i, B.I+i, &dIdV);
	return B3SolPars(&B, S);
}
