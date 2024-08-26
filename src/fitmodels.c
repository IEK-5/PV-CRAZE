/*
    fitmodels
    Framework interface to the various models 
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
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "emalloc.h"
#include "tictoc.h"
#include "curiosity.h"
#include "curve.h"
#include "fitmodels.h"
#include "model_table.h" // only for fitmodels.h to figure out how to treat a model

MODELSET DefaultMODSET = {
	M_NONE,
	50,					// max iteration for basic fitting algorithm (if applicale)
	1e-5,				// eps for basic fitting optimizer (if applicale)
	50,					// iterations for model evaluation (if applicale)
	100,				// population size (if applicale)
	1e-5,				// model evaluation eps (if applicale)
	1e-2,				// oeps eps threshold for expensive optim,ization (if applicable)
	false,				// refine solar cell parameter during fit (if applicable)
};


typedef enum {
	KW_MODEL,			// model
	KW_FITIT,			// fititer
	KW_FITEPS,			// fiteps
	KW_MIT,				// miter
	KW_NPOP,			// npop
	KW_MEPS,			// meps
	KW_OEPS,			// oeps
	KW_REFINE,			// RefinePars
	KW_NONE
} KEYWORDS;

typedef struct {
	char *keyword;
	KEYWORDS K;
} key_entry;


key_entry parsetab[] = {
	{"model",   		KW_MODEL},
	{"fitit",   		KW_FITIT},
	{"npop",   			KW_NPOP},
	{"fiteps",  		KW_FITEPS},
	{"mit",     		KW_MIT},
	{"meps",    		KW_MEPS},
	{"oeps",    		KW_OEPS},
	{"refine",  		KW_REFINE},
	{NULL,      		KW_NONE}
};

KEYWORDS LookupKey(char *word)
{
	int i=0;
	while (parsetab[i].keyword)
	{
		if (strncmp(parsetab[i].keyword,word,50)==0)
			return parsetab[i].K;
		i++;
	}
	return KW_NONE;
}

MODEL_TYPE LookupModel(char *word)
{
	int i=0;
	while (MODEL_TABLE[i].TYPE!=M_NONE)
	{
		if (strncmp(MODEL_TABLE[i].name,word,10)==0)
			return MODEL_TABLE[i].TYPE;
		i++;
	}
	if (strncmp(word,"kitten",10)==0)
	{
		/* 25.03.24, it is almost easter. So,
		 * if you are looking for a kitten,
		 * you may find one on Mars...
		 */
		Curiosity(5, 5);
	}
	return M_NONE;
}

char *pv_craze_reverselookup(MODEL_TYPE T)
{
	int i=0;
	while (MODEL_TABLE[i].TYPE!=M_NONE)
	{
		if (MODEL_TABLE[i].TYPE==T)
			return MODEL_TABLE[i].name;
		i++;
	}
	return NULL;
}
int IsInStr(char c, char *string)
{
	char *p;
	if (c=='\0')
		return 1; // all strings should be null terminated so null is always in every string
	p=string;
	while (*p)
	{
		if (*p==c)
			return 1;
		p++;
	}
	return 0;
}
#define MAXSTRING 1000
char *NextWord(const char *start, char *word)
{
	char *p, *w;
	w=word;
	p=(char *)start;
	while (*p)
	{
		if (IsInStr(*p, ", \t()"))
			break;
		*w=*p;
		w++;
		if (w-word==MAXSTRING-2)
			break;
		p++;
	}
	*w='\0';
	while ((*p)&&(IsInStr(*p, ", \t()")))
		p++;
	if (*p)
		return p;
	else
		return NULL;
}
/* the spec string:
 * comma separated, white space is ignored
 * fields have the form
 * <keyword>(<arg>)
 */


#define MINITER 0
#define MAXITER 10000
#define MINEPS 1e-12
#define MAXEPS 1e-1
#define MINPOP 10
#define MAXPOP 1000


int ParseModel(const char *spec, char *word, MODELSET *S)
{
	char *next;
	int inum;
	double dnum;
	next=(char *)spec;
	if (S==NULL)
		return 1;
	// parse model specs
	while (next)
	{
		next=NextWord(next, word);
		switch (LookupKey(word))
		{
			case KW_MODEL:
				if (S->M!=M_NONE)// model
				{
					fprintf(stderr,"Warning: cannot specify two models in one go. Ignoring %s.\n", word);
					next=NextWord(next, word); // skip the model specifier
				}
				else
				{
					next=NextWord(next, word);
					S->M=LookupModel(word);
					if (S->M==M_NONE)
						fprintf(stderr,"Warning: no model %s available.\n", word);
				}
				break;
			case KW_FITIT:			// fititer
				next=NextWord(next, word);
				inum=atoi(word);
				if (inum<MINITER)
					fprintf(stderr,"Warning: cannot specify a negative number of fit iterations\n");
				else if (inum>MAXITER)
					fprintf(stderr,"Warning: you are going overboard with the number of fit iterations\n");
				else
					S->fititer=inum; 
				break;		
			case KW_NPOP:			// npop
				next=NextWord(next, word);
				inum=atoi(word);
				if (inum<MINPOP)
					fprintf(stderr,"Warning: cannot specify a population size less than %d\n", MINPOP);
				else if (inum>MAXPOP)
					fprintf(stderr,"Warning: you are going overboard with the population size\n");
				else
					S->npop=inum; 
				break;		
			case KW_FITEPS:			// fiteps
				next=NextWord(next, word);
				dnum=atof(word);
				if (dnum<MINEPS)
					fprintf(stderr,"Warning: cannot specify an epsilon smaller than %e\n", MINEPS);
				else if (dnum>MAXEPS)
					fprintf(stderr,"Warning: an epsilon value should be small (not larger than %e)\n", MAXEPS);
				else 
					S->fiteps=dnum;
				break;	
			case KW_MIT:				// miter
				next=NextWord(next, word);
				inum=atoi(word);
				if (inum<MINITER)
					fprintf(stderr,"Warning: cannot specify a negative number of model iterations\n");
				else if (inum>MAXITER)
					fprintf(stderr,"Warning: you are going overboard with the number of model iterations\n");
				else 
					S->miter=inum;
				break;		
			case KW_MEPS:			// meps
				next=NextWord(next, word);
				dnum=atof(word);
				if (dnum<MINEPS)
					fprintf(stderr,"Warning: cannot specify an epsilon smaller than %e\n", MINEPS);
				else if (dnum>MAXEPS)
					fprintf(stderr,"Warning: an epsilon value should be small (not larger than %e)\n", MAXEPS);
				else 
					S->meps=dnum;
				break;					
			case KW_OEPS:			// meps
				next=NextWord(next, word);
				dnum=atof(word);
				if (dnum<MINEPS)
					fprintf(stderr,"Warning: cannot specify an epsilon smaller than %e\n", MINEPS);
				else if (dnum>MAXEPS)
					fprintf(stderr,"Warning: an epsilon value should be small (not larger than %e)\n", MAXEPS);
				else 
					S->oeps=dnum;
				break;	
			case KW_REFINE:
				printf("       refine solar cell parameters\n");
				S->RefinePars=true;
				break;		
			case KW_NONE:
			default:
				fprintf(stderr,"Error: unknown keyword %s\n", word);
				break;
		}
	}
	if (S->M==M_NONE)
		return 1;
	else
		return 0;
}


void FreeList(MODEL_LIST *L)
{
	if (L->N)
	{
		free(L->List);
		free(L->pars);
		free(L->Efit);
		free(L->tm);
		L->List=NULL;
		L->Efit=NULL;
		L->pars=NULL;
		L->tm=NULL;
		L->N=0;
	}
}


generic_pars init_gen_pars={{0,0,0,0,0,0,0,0,0,0},M_NONE};
MODEL_LIST InitList()
{
	MODEL_LIST L={0};
	if((L.List=Emalloc(sizeof(MODELSET)))==NULL)
		return L;
	if((L.pars=Emalloc(sizeof(generic_pars)))==NULL)
	{
		free(L.List);
		L.List=NULL;
		return L;
	}
	if((L.Efit=Emalloc(sizeof(double)))==NULL)
	{
		free(L.List);
		L.List=NULL;
		free(L.pars);
		L.pars=NULL;
		return L;
	}
	if((L.tm=Emalloc(sizeof(double)))==NULL)
	{
		free(L.List);
		L.List=NULL;
		free(L.pars);
		L.pars=NULL;
		free(L.Efit);
		L.Efit=NULL;
		return L;
	}
	L.List[0]=DefaultMODSET;
	L.pars[0]=init_gen_pars;
	L.Efit[0]=0;
	L.tm[0]=0;
	L.N=1;
	return L;
}

void AddModel(MODEL_LIST *L, const char *spec)
{
	char *word;
	void *tmp;
	
	if (L->N<=0)
		(*L)=InitList();
	else
	{
		if((tmp=Erealloc(L->List, (L->N+1)*sizeof(MODELSET)))==NULL)
			return;
		else
			L->List=(MODELSET *)tmp;
			
		if((tmp=Erealloc(L->pars, (L->N+1)*sizeof(generic_pars)))==NULL)
			return;
		else
			L->pars=(generic_pars *)tmp;
			
		if((tmp=Erealloc(L->Efit, (L->N+1)*sizeof(double)))==NULL)
			return;
		else
			L->Efit=(double *)tmp;
			
		if((tmp=Erealloc(L->tm, (L->N+1)*sizeof(double)))==NULL)
			return;
		else
			L->tm=(double *)tmp;
		
		L->List[L->N]=DefaultMODSET;
		L->pars[L->N]=init_gen_pars;
		L->Efit[L->N]=0;
		L->tm[L->N]=0;
		L->N++;
	}
	if((word=Emalloc(MAXSTRING*sizeof(char)))==NULL)
	{
		L->N--;
		return;
	}
	
	if (ParseModel(spec,word, L->List+L->N-1))
	{
		fprintf(stderr,"Error: could not add model %s\n", spec);
		L->N--;		
	}
	free(word);
}


int ModelIndex(MODEL_TYPE T)
{
	int i=0;
	while (MODEL_TABLE[i].TYPE!=M_NONE)
	{
		if (MODEL_TABLE[i].TYPE==T)
			return i;
		i++;
	}
	return i;
}


curve _ToCurve(generic_pars P, int N, MODELSET *S)
{
	int i;
	if (P.TYPE==M_NONE)
	{
		curve C;
		C.N=0;
		C.P=NULL;
		return C;
	}
	i=ModelIndex(P.TYPE);
	return MODEL_TABLE[i].CURVE(P,N,S);
}
double 	_RMSErr(generic_pars P, curve C, MODELSET *S)
{
	int i;
	if (P.TYPE==M_NONE)
		return NAN;
	i=ModelIndex(P.TYPE);
	return MODEL_TABLE[i].RMS(P,C,S);
}

curve 	_DistCurve(generic_pars P, curve C, MODELSET *S)
{
	int i;
	if (P.TYPE==M_NONE)
	{
		curve C;
		C.N=0;
		C.P=NULL;
		return C;
	}
	i=ModelIndex(P.TYPE);
	return MODEL_TABLE[i].DIST(P,C,S);
}
curvepars _Pars(generic_pars P, MODELSET *S)
{
	int i;
	if (P.TYPE==M_NONE)
		return nullpars;
	i=ModelIndex(P.TYPE);
	return MODEL_TABLE[i].PARS(P,S);
}

double _IV(double V, generic_pars P, MODELSET *S)
{
	int i;
	if (P.TYPE==M_NONE)
		return NAN;
	i=ModelIndex(P.TYPE);
	return MODEL_TABLE[i].IV(V,P,S);
}
double _VI(double I, generic_pars P, MODELSET *S)
{
	int i;
	if (P.TYPE==M_NONE)
		return NAN;
	i=ModelIndex(P.TYPE);
	return MODEL_TABLE[i].VI(I,P,S);
}

generic_pars _Fit(MODEL_TYPE T, curvepars par, curve C, double *Efit, MODELSET *S)
{
	int i;
	if (T==M_NONE)
	{
		generic_pars p;
		for (i=0;i<10;i++)
			p.pars[i]=0;
		p.TYPE=M_NONE;
		return p;
	}
	i=ModelIndex(T);
	
	return MODEL_TABLE[i].FIT(par,C,Efit,S);
}

curve pv_craze_g_tocurve(generic_pars P, int N)
{
	return _ToCurve(P, N, &DefaultMODSET);
}
double 	pv_craze_g_rmse(generic_pars P, curve C)
{
	return _RMSErr(P, C, &DefaultMODSET);
}
curve pv_craze_g_distcurve(generic_pars P, curve C)
{
	return _DistCurve(P, C, &DefaultMODSET);
}
curvepars pv_craze_g_pars(generic_pars P)
{
	return _Pars(P, &DefaultMODSET);
}
generic_pars pv_craze_g_fit(MODEL_TYPE T, curvepars par, curve C, double *Efit)		
{
	return _Fit(T, par, C, Efit, &DefaultMODSET);
}

double pv_craze_IV(double V, generic_pars P)
{
	return _IV(V,P, &DefaultMODSET);
}
double pv_craze_VI(double I, generic_pars P)
{
	return _VI(I,P, &DefaultMODSET);
}

void FitModels(MODEL_LIST *L, curvepars par, curve C, int R)
{
	int i;
	curve Cs;
	curve Ca;
	double t;
	if (C.N==0)
	{
		// if no curve is given but we do have ESP, we make a minimal curve out of those
		Ca.N=5;
		Ca.P=malloc(5*sizeof(point));
		Ca.P[0].x=0;
		Ca.P[0].y=par.Isc;
		Ca.P[1].x=par.Vqmp1;
		Ca.P[1].y=par.Iqmp1;
		Ca.P[2].x=par.Vmpp;
		Ca.P[2].y=par.Impp;
		Ca.P[3].x=par.Vqmp2;
		Ca.P[3].y=par.Iqmp2;
		Ca.P[4].x=par.Voc;
		Ca.P[4].y=0;
		R=0; // resample makes no sense here
	}
	else
		Ca=PreprocessCurve(C,0);
	if (Ca.N==0)
		Cs=Ca;
	else 
	{
		if (R)
			Cs=Resample(Ca, R, par.Voc, par.Isc);
		else
			Cs=Ca;
	}
	for (i=0;i<L->N;i++)
	{
		TIC(&t);
		L->pars[i]=_Fit(L->List[i].M, par, Cs, L->Efit+i,L->List+i);
		L->tm[i]+=TOC(&t);
	}
	if ((Ca.N)&&(R))
		free(Cs.P);
	free(Ca.P);
	return;
}

int PropIndex(MODEL_TYPE T)
{
	int i=0;
	while (MODEL_PROP[i].TYPE!=M_NONE)
	{
		if (MODEL_PROP[i].TYPE==T)
			return i;
		i++;
	}
	return i;
}

int pv_craze_numfields(MODEL_TYPE T)
{
	int i;
	i=PropIndex(T);
	return MODEL_PROP[i].N;
}
