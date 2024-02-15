/*
    Quadratic Bézier IV model.
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
	double Vk[5], Ik[5]; // 5 key points
	double Vi[4], Ii[4]; // intersection points
} B2pars;


/* quadratic Bézier base equations */
#define BB0(t) ((1.0-(t))*(1.0-(t)))
#define BB1(t) (2.0*(1.0-(t))*(t))
#define BB2(t) ((t)*(t))

/* derivatives */
#define DBB0(t) (2.0*(1.0-(t)))
#define DBB1(t) (2.0*(t))

B2pars Generic2B2Pars(generic_pars P)
{
	B2pars S;
	double I0, I1, V0, V1, dI0, dI1;
	int i;
	
	if (P.TYPE!=M_B2)
	{
		fprintf(stderr,"Error: wrong parameter type in Generic2B2Pars\n");
		S.Vk[0]=0;
		S.Vk[1]=0;
		S.Vk[2]=0;
		S.Vk[3]=0;
		S.Vk[4]=0;
		
		S.Ik[0]=0;
		S.Ik[1]=0;
		S.Ik[2]=0;
		S.Ik[3]=0;
		S.Ik[4]=0;
		
		S.Vi[0]=0;
		S.Vi[1]=0;
		S.Vi[2]=0;
		S.Vi[3]=0;
		
		S.Ii[0]=0;
		S.Ii[1]=0;
		S.Ii[2]=0;
		S.Ii[3]=0;
		
		return S;
	}
	// the knots are directly available
	S.Vk[0]=0;
	S.Vk[1]=P.pars[4];
	S.Vk[2]=P.pars[2];
	S.Vk[3]=P.pars[6];
	S.Vk[4]=P.pars[0];
	
	S.Ik[0]=P.pars[1];
	S.Ik[1]=P.pars[5];
	S.Ik[2]=P.pars[3];
	S.Ik[3]=P.pars[7];
	S.Ik[4]=0;
	
	
	// the intermediate points need to be computed from the knot coordinates and derivatives
	for (i=0;i<4;i++)
	{
		I0=S.Ik[i];
		V0=S.Vk[i];
		I1=S.Ik[i+1];
		V1=S.Vk[i+1];
		switch (i)
		{
			case 0:
				dI0=-P.pars[8];
				dI1=(S.Ik[2]-S.Ik[1])/S.Vk[1];
				break;
			case 1:
				dI0=(S.Ik[2]-S.Ik[1])/S.Vk[1];
				dI1=-S.Ik[2]/S.Vk[2];
				break;
			case 2:
				dI0=-S.Ik[2]/S.Vk[2];
				dI1=S.Ik[3]/(S.Vk[3]-S.Vk[4]);
				break;
			case 3:
			default:
				dI0=S.Ik[3]/(S.Vk[3]-S.Vk[4]);
				dI1=-1/P.pars[9];
				break;
		}
		if (fabs(dI0-dI1)>1e-10*fabs(I1-I0))
		{
			S.Vi[i]=(V1*dI1-V0*dI0-I1+I0)/(dI1-dI0);
			S.Ii[i]=I0+dI0*(S.Vi[i]-V0);
		}
		else
		{
			// either the derivatives are 0, or of opposite sign and equal magnitude
			if (fabs(dI0)<1e-10*fabs(I1-I0))
			{
				// panic! lines do not cross
				S.Vi[i]=(V0+V1)/2;
				S.Ii[i]=(I0+I1)/2;
			}
			else
			{
				S.Vi[i]=(V0+V1)/2;
				S.Ii[i]=I0+dI0*(S.Vi[i]-V0);
			}
		}
	}
	
	return S;
}



// convert internal to general parameter struct
generic_pars B22GenericPars(B2pars S)
{
	generic_pars P;
	P.TYPE=M_B2;
	P.pars[0]=S.Vk[4];
	P.pars[1]=S.Ik[0];
	P.pars[2]=S.Vk[2];
	P.pars[3]=S.Ik[2];
	P.pars[4]=S.Vk[1];
	P.pars[5]=S.Ik[1];
	P.pars[6]=S.Vk[3];
	P.pars[7]=S.Ik[3];
	P.pars[8]=(S.Ik[0]-S.Ii[0])/S.Vk[1];
	P.pars[9]=(S.Vk[4]-S.Vi[3])/S.Ii[3];
	return P;
}	

point B2IV_i (double t, point *dp, void *pars)
{
	point p;
	B2pars S;
	int i;
	
	S=*((B2pars *)pars);
	
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
		// quadratic Bézier interpolation
		p.x=BB0(t)*S.Vk[i]+BB1(t)*S.Vi[i]+BB2(t)*S.Vk[i+1];
		p.y=BB0(t)*S.Ik[i]+BB1(t)*S.Ii[i]+BB2(t)*S.Ik[i+1];
	
		if (dp)
		{
			dp->x=DBB0(t)*(S.Vi[i]-S.Vk[i])+DBB1(t)*(S.Vk[i+1]-S.Vi[i]);
			dp->y=DBB0(t)*(S.Ii[i]-S.Ik[i])+DBB1(t)*(S.Ik[i+1]-S.Ii[i]);
		}
	}
	else
	{
		// linear extrapolation
		if ((t<0)&&(i==0))
		{
			p.x=S.Vk[0]+t*(S.Vi[0]-S.Vk[0]);
			p.y=S.Ik[0]+t*(S.Ii[0]-S.Ik[0]);;
		
			if (dp)
			{
				dp->x=(S.Vi[0]-S.Vk[0]);
				dp->y=(S.Ii[0]-S.Ik[0]);
			}
		}
		else if((t>1)&&(i==3))
		{
			p.x=S.Vk[4]+t*(S.Vk[4]-S.Vi[3]);
			p.y=S.Ik[4]+t*(S.Ik[4]-S.Ii[3]);;
		
			if (dp)
			{
				dp->x=(S.Vk[4]-S.Vi[3]);
				dp->y=(S.Ik[4]-S.Ii[3]);
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
	}			
	return p;
}

// Exported Function
curve B2ToCurve(generic_pars P, int N, MODELSET *S)
{
	B2pars par;
	par=Generic2B2Pars(P);
	return ModelToCurve	((void *) &par, &B2IV_i, N, S);
}

double B2RMSErr_i(B2pars P, curve C, MODELSET *S)
{
	return RMSErr(P.Vk[4], P.Ik[0], (void *) &P, &B2IV_i, C,S);
}

// Exported Function
double B2RMSErr(generic_pars P, curve C, MODELSET *S)
{
	return B2RMSErr_i(Generic2B2Pars(P), C,S);
}

curve B2DistCurve_i(B2pars P, curve C, MODELSET *S)
{
	return DistCurve(P.Vk[4], P.Ik[0], (void *) &P, &B2IV_i, C,S);
}

// Exported Function
curve B2DistCurve(generic_pars P, curve C, MODELSET *S)
{
	return B2DistCurve_i(Generic2B2Pars(P), C,S);
}
curvepars B2Pars(generic_pars P, MODELSET *S)
{
	curvepars par;
	
	if (P.TYPE!=M_B2)
	{
		fprintf(stderr,"Error: wrong parameter type in B2Pars\n");
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


generic_pars B2Fit(curvepars par, curve C, double *Efit, MODELSET *S)
{
	generic_pars P;
	B2pars Sd;
	P.TYPE=M_B2;
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
	
	Sd=Generic2B2Pars(P);
	(*Efit)=B2RMSErr_i(Sd, C, S);
	return P;
}


