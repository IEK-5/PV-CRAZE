/*
    Karmalkar-Haneefa IV model.
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
#include "diffevol.h"

// Warning: I take the reciprocal value of m as then both gamma and m 
// are between 0 and 1 which makes it easier to do a contraint optimization
// in the generic pars the meaning of m is, however, as it is in literature
#define KHj(v,g,m) (1.0-(1.0-(g))*(v)-(g)*pow((v),1/(m))) 
#define dKHj(v,g,m) ((g)-1.0-(g)*pow((v),(1/(m)-1.0))/(m))


typedef struct {
	double gamma,m,Voc,Isc;
} KHpars;


// convert general to internal parameter struct
KHpars Generic2KHPars(generic_pars P)
{
	KHpars S;
	
	if (P.TYPE!=M_KH)
	{
		fprintf(stderr,"Error: wrong parameter type in Generic2KHPars\n");
		S.gamma=0;
		S.m=1;
		S.Voc=0;
		S.Isc=0;
		
		return S;
	}
	S.Voc=P.pars[0];
	S.Isc=P.pars[1];
	S.gamma=P.pars[2];
	S.m=1/P.pars[3];
	return S;
}



// convert internal to general parameter struct
generic_pars KH2GenericPars(KHpars S)
{
	generic_pars P;
	P.TYPE=M_KH;
	P.pars[0]=S.Voc;
	P.pars[1]=S.Isc;
	P.pars[2]=S.gamma;
	P.pars[3]=1/S.m;
	P.pars[4]=0;
	P.pars[5]=0;
	P.pars[6]=0;
	P.pars[7]=0;
	P.pars[8]=0;
	P.pars[9]=0;
	return P;
}	


// pars are KHpars
point KHIV_i (double t, point *dp, void *pars)
{
	point p;
	KHpars par;
	
	par=*((KHpars *)pars);
	p.x=t*par.Voc; // for t 0-1 we cover a bit wider voltage range to be sure
	p.y=KHj(t,par.gamma,par.m)*par.Isc;
	if (dp)
	{
		dp->y=dKHj(t,par.gamma,par.m)*par.Isc;
		dp->x=par.Voc;
	}
	return p;
}

// Exported Function
curve KHToCurve(generic_pars P, int N, MODELSET *S)
{
	KHpars par;
	par=Generic2KHPars(P);
	return ModelToCurve	((void *) &par, &KHIV_i, N, S);
}

double KHRMSErr_i(KHpars P, curve C, MODELSET *S)
{
	return RMSErr(P.Voc, P.Isc, (void *) &P, &KHIV_i, C,S);
}

// Exported Function
double KHRMSErr(generic_pars P, curve C, MODELSET *S)
{
	return KHRMSErr_i(Generic2KHPars(P), C,S);
}

curve KHDistCurve_i(KHpars P, curve C, MODELSET *S)
{
	return DistCurve(P.Voc, P.Isc, (void *) &P, &KHIV_i, C, S);
}

// Exported Function
curve KHDistCurve(generic_pars P, curve C, MODELSET *S)
{
	return KHDistCurve_i(Generic2KHPars(P), C, S);
}
curvepars KHPars(generic_pars P, MODELSET *S)
{
	KHpars pp;
	pp=Generic2KHPars(P);
	
	return SolPar(pp.Voc, pp.Isc, (void *) &pp, &KHIV_i, S);	
}


/* I built 3 fitting routines:
 * 1 newton solver
 * 1 Levenberg Marquard
 * 1 Differential evolution
 * As far as I know all work fine. The differential evolution one, however
 * is a constrained solver which only produces monotonic solutions.
 */

// Newton solver
KHpars KHparsameterize(curvepars par, double *Epar, MODELSET *S)
{
	KHpars P;
	int iter=0;
	double vp=par.Vmpp/par.Voc;
	double ip=par.Impp/par.Isc;
	double dg, dm, e1, e2;
	double D=1;
	if (ip<1-vp)
	{
		// TODO: add warning system, FF under 25%
		P.gamma=0;
		P.m=1;
		return P;
	}
	P.Voc=par.Voc;
	P.Isc=par.Isc;	
	
	
	P.gamma=-(par.Isc-par.Impp)*P.Voc/(P.Isc*par.Vmpp)+1;
	P.m=P.gamma/(P.Voc*par.Impp/((par.Voc-par.Vmpp)*P.Isc)+1-P.gamma);
	while ((D>S->fiteps)&&(iter<S->fititer))
	{
		e1=KHj(vp,P.gamma,P.m)-ip;
		e2=ip/vp+dKHj(vp,P.gamma,P.m);
		
		
		dg = (P.m*vp*log(vp)*e2+((-log(vp))-P.m)*e1)/((P.m-1)*vp*log(vp)+P.m*pow(vp,1/P.m)-P.m*vp);
		dm = ((pow(P.m,3)*pow(vp,1/P.m+1)-pow(P.m,3)*vp*vp)*e2+(pow(P.m,3)*vp-pow(P.m,2)*pow(vp,1/P.m))*e1)/((P.gamma*P.m-P.gamma)*pow(vp,1/P.m+1)*log(vp)-P.gamma*P.m*pow(vp,1/P.m+1)+P.gamma*P.m*pow(vp,2/P.m));
		
		D=sqrt((dg*dg+dm*dm)/2);
		if (isfinite(dm))
			P.m-=dm;
		else
			break;
		if (isfinite(dg))
			P.gamma-=dg;
		else
			break;	
		iter++;	
	}
	(*Epar)=D;
	return P;
}



/********************************* GSL optimizer ************************/

#define Impp (par.Impp/par.Isc)
#define Vmpp  (par.Vmpp/par.Voc)
#define Pg (gsl_vector_get(x,0))
#define Pm (gsl_vector_get(x,1))
int KH_RHS_gsl(const gsl_vector * x, void *data, gsl_vector * b)
{
	int i;
	curvepars par;
	par=*(curvepars *)data;
	
	/* RHS */
	i=0;
	gsl_vector_set(b,i++,(-Pg*pow(Vmpp,(1/Pm)))-(1-Pg)*Vmpp-Impp+1);
	gsl_vector_set(b,i++,(-(Pg*pow(Vmpp,(1/Pm-1)))/Pm)+Impp/Vmpp+Pg-1);
	return GSL_SUCCESS;
}

int KH_System_gsl(const gsl_vector * x, void *data, gsl_matrix * A)
{
	int i,j;
	curvepars par;
	par=*(curvepars *)data;
	/* A */
	i=0;
	j=0;	
	
	gsl_matrix_set(A,i++,j,Vmpp-pow(Vmpp,(1/Pm)));
	gsl_matrix_set(A,i++,j,-(pow(Vmpp,(1/Pm))-Pm*Vmpp)/(Pm*Vmpp));
	
	j++;
	i=0;
	gsl_matrix_set(A,i++,j,(Pg*pow(Vmpp,(1/Pm))*log(Vmpp))/pow(Pm,2));
	gsl_matrix_set(A,i++,j,(Pg*pow(Vmpp,(1/Pm))*log(Vmpp)+Pg*Pm*pow(Vmpp,(1/Pm)))/(pow(Pm,3)*Vmpp));

	return GSL_SUCCESS;
}
#undef Impp
#undef Vmpp
#undef Pg
#undef Pm


/* GSL nonlinear least squares */
KHpars KH_LevenbergMarquardtOptimizer(curvepars par, KHpars P, MODELSET *S)
{
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *w;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	const size_t n = 2;
	const size_t p = 2;

	gsl_vector *b;
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	
	double x_init[2], weights[2];
	gsl_vector_view x;
	gsl_vector_view wts;
	
	double chisq0;
	int info;
	size_t i;

	const double xtol = S->fiteps;
	const double gtol = S->fiteps;
	const double ftol = 0.0;	
	
	x_init[0]=P.gamma;
	x_init[1]=P.m;
	
	x = gsl_vector_view_array (x_init, p);
	wts = gsl_vector_view_array(weights, n);
  

	/* define the function to be minimized */
	fdf_params.h_fvv=1e-6;
	fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
	fdf.f = KH_RHS_gsl;
	fdf.df = KH_System_gsl;   /* set to NULL for finite-difference Jacobian */
	fdf.fvv = NULL;     /* not using geodesic acceleration */
	fdf.n = n;
	fdf.p = p;
	fdf.params = &par;
	for (i = 0; i < n; i++)
		weights[i] = 1.0;
	

	/* allocate workspace with default parameters */
	w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

	/* initialize solver with starting point and weights */
	gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

	/* compute initial cost function */
	b = gsl_multifit_nlinear_residual(w);
	gsl_blas_ddot(b, b, &chisq0);

	/* solve the system with a maximum of 100 iterations */
	gsl_multifit_nlinear_driver( S->fititer, xtol, gtol, ftol, NULL, NULL, &info, w);
	
	P.gamma=fabs(gsl_vector_get(w->x,0));
	P.m=fabs(gsl_vector_get(w->x,1));
	gsl_multifit_nlinear_free (w);
	gsl_matrix_free (covar);

	return P;
}

double KH_DE_Err(double x[], void *args)
{
	curvepars pars;
	double E, vp, ip, e;
	pars=*((curvepars *)args);
	// constraints implemented as returning inf for invalid parameters
	// g>=-1/(m-1)=-m'/(1-m')
	if (x[0]<-x[1]/(1-x[1]))
	{
		return INFINITY;
	}
	vp=pars.Vmpp/pars.Voc;
	ip=pars.Impp/pars.Isc;
	
	e=KHj(vp,x[0],x[1])-ip;
	E=e*e;
	e=ip/vp+dKHj(vp,x[0],x[1]);
	E+=e*e;
	return E;
}
#define RandFloat(a,b) ((a)+((b)-(a))*((double)rand())/RAND_MAX)
KHpars KH_DE(curvepars par, MODELSET *S)
{
	double **prime;
	int Nprime=0;
	de_settings Se=de_default;
	double m, g, gmin=1.0;
	KHpars P={0};
	double x[2], lb[2], ub[2];
	time_t t;
	srand((unsigned) time(&t));

	Se.eps=S->fiteps;
	Se.max_iter=5000;
	Se.Npop=90;
	// constrained parameter initialization for monotonous solutions
	// gamma<1
	// m >= 1 (m' = 1/m < 1)
	// g >= -1/(m-1) = -m'/(1-m')
	if((prime=Emalloc(Se.Npop*sizeof(double *)))==NULL)
		return P;
		
	for (Nprime=0;Nprime<Se.Npop;Nprime++)
	{
		if((prime[Nprime]=Emalloc(2*sizeof(double)))==NULL)
		{
			Nprime--;
			for(;Nprime>=0;Nprime--)
				free(prime[Nprime]);
			free(prime);
			return P;
		}
		m=RandFloat(0,1); // first choose m'
		g=RandFloat(-m/(1.0-m),1);// constrain gamma
		if (g<gmin)
			gmin=g;
		prime[Nprime][0]=g;
		prime[Nprime][1]=m;
	}
	lb[0]=gmin;
	lb[1]=0;
	
	ub[0]=1;
	ub[1]=1;
	
	if (isfinite(diffevol_prime(&KH_DE_Err, (void *) (&par), 2, prime, Nprime, x, lb, ub, NULL, &Se)))
	{
		P.Voc=par.Voc;
		P.Isc=par.Isc;
		P.gamma=x[0];
		P.m=x[1];
	}
	else
	{
		P.Voc=par.Voc;
		P.Isc=par.Isc;
		P.gamma=0.5;
		P.m=0.5;
	}
	for (;Nprime>0;Nprime--)
		free(prime[Nprime-1]);
	free(prime);
		
	return P;
}

int MonoPars(curvepars par)
{
	if (par.Isc<0)
		return 1;
	if (par.Voc<0)
		return 1;
				
	if (par.Isc<par.Impp)
		return 1;
	if (par.Voc<par.Vmpp)
		return 1;
		
	if (par.Vmpp<0)
		return 1;
	if (par.Impp<0)
		return 1;
	return 0;
}



generic_pars KHFit(curvepars par, curve C, double *Efit, MODELSET *S)
{
	KHpars P;
	
	P.Voc=par.Voc;
	P.Isc=par.Isc;
	P.gamma=-(par.Isc-par.Impp)*P.Voc/(P.Isc*par.Vmpp)+1;
	P.m=P.gamma/(P.Voc*par.Impp/((par.Voc-par.Vmpp)*P.Isc)+1-P.gamma);
	if (P.gamma<1e-10)
		P.gamma=0.5;
	if (P.gamma>1-1e-10)
		P.gamma=0.5;
	if (P.m<1e-10)
		P.m=0.5;
	if (P.m>1-1e-10)
		P.m=0.5;		
	P=KH_LevenbergMarquardtOptimizer(par, P, S);
	// several conditions are no good
	// if P.m>1 (m<1) the derivative at v=0 is not finite
	// if gamma<0 && P.m<P.gamma/(P.gamma-1)
	if ((P.gamma>1)||(P.m>1)||((P.gamma<0)&&(P.m<P.gamma/(P.gamma-1))))
	{
		// optimum is not monotonous!
		// in case the curve parameters *are* monotonous we
		// do a new constrained optimization
		if (!MonoPars(par)) 
			P=KH_DE(par, S);
	}
	//P=KH_LevenbergMarquardtOptimizer(par, P);
	//P=KHparsameterize(par, &e);
	(*Efit)=KHRMSErr_i(P, C,S);
	return KH2GenericPars(P);
}

