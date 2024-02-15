/*
 * Interface to pv-craze:
 * A top level interface to pv-craze routines
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
#include "Error.h"
#include "curve.h"
#include "tictoc.h"
#include "fitmodels.h"
#include "spline_iv.h" // convexity
#include "lib-pv-craze-structs.h"

typedef struct {
	int flag;
	int count;
} FLAGCOUNT;

typedef struct {
	FLAGCOUNT *F;
	int N;
} FLAGLIST;


#define VERB_TIMING 1		// print timing
#define VERB_REPORT 2		// print report
#define VERB_PROGRESS 4		// print progress bar
#define VERB_HEAD 8			// print function header
 
int pvcraze_verbose=VERB_TIMING|VERB_REPORT|VERB_PROGRESS|VERB_HEAD;

void pv_craze_quiet()
{
	pvcraze_verbose=0;
}
void pv_craze_verbose()
{
	pvcraze_verbose=VERB_TIMING|VERB_REPORT|VERB_PROGRESS|VERB_HEAD;
}
void VerboseControll(int set, int verbose)
{
	if (set)
		pvcraze_verbose|=verbose;
	else
		pvcraze_verbose&=(~verbose);
}

void pv_craze_verbosetiming(int set)
{
	VerboseControll(set, VERB_TIMING);
}
void pv_craze_verbosereport(int set)
{
	VerboseControll(set, VERB_REPORT);
}
void pv_craze_verboseprogress(int set)
{
	VerboseControll(set, VERB_PROGRESS);
}
void pv_craze_verbosehead(int set)
{
	VerboseControll(set, VERB_HEAD);
}
// default settings for curve analysis
analyze_settings def_settings = {
	0,								// SpikeFilter
	0,								// Ident Filter
	0,								// resample
	1e-12,							// Ident-Eps
	0.5,							// Min. Rel. Range
	0.12, 0.12,						// MPP_f_V, MPP_f_I
	0.1, 0.1,						// OC_f_I, SC_f_V
	's',							// am
};

MODEL_TYPE pv_craze_lookupmodel(char *word)
{
	return LookupModel(word);
}
#define LISTBLOCK 16
int IsInList(FLAGLIST *list, int id)
{
	int i;
	if (list->N==0)
		return -1;
	i=0;
	while(i<list->N)
	{
		if (list->F[i].flag==id)
			return i;
		i++;
	}
	return -1;
}

FLAGLIST InitFlagList()
{
	FLAGLIST F;
	F.N=0;
	F.F=(FLAGCOUNT *)Emalloc((3+LISTBLOCK)*sizeof(FLAGCOUNT));
	return F;
}

int AddToList(FLAGLIST *list, int id)
{
	int index;
	if (!list->F)
		return 1;
	if ((index=IsInList(list, id))<0)
	{
		if ((list->N+3)%LISTBLOCK==0)
		{
			FLAGCOUNT *tmp;
			tmp=(FLAGCOUNT *)Erealloc(list->F,(list->N+3+LISTBLOCK)*sizeof(FLAGCOUNT));
			if (tmp)
				list->F=tmp;
			else
				return 1;
		}	
		list->F[list->N].flag=id;
		list->F[list->N].count=1;
		list->N++;
		return 0;
	}
	list->F[index].count++;
	return 0;
}

int CompareErr(const void *a, const void * b)
{
	FLAGCOUNT A, B;
	A=*((FLAGCOUNT *)a);
	B=*((FLAGCOUNT *)b);
	if (A.count<B.count)
		return 1;
	if (A.count>B.count)
		return -1;
	return 0;
} 

void SortErrWarn(FLAGLIST *errwarn)
{
	qsort((void *)errwarn->F, errwarn->N, sizeof(FLAGCOUNT), &CompareErr);
}

// ProgressBar
void BarScale(int len, int tics)
{	
	int l, tic;
	if (!(pvcraze_verbose&VERB_PROGRESS))
		return;
		
	l=0;
	tic=len/tics;
	while (l<len)
	{							
		l++;
		if (l==1)
		{
			printf("0 %%-");
			l+=2;
		}
		else if (l==len)
		{
			printf("100 %%");
			l+=4;
		}
		else if (l%tic==0)
		{
			printf("+");
		}
		else
			printf("-");	
	}
	printf("\n|");	
	fflush(stdout);
}

	
int ProgressBar(int pcn, int pco, int len, int tics)
/* pcn: new percentage complete */
/* pco: old percentage complete */
/* len: length of the progress bar */
/* tics: number of tics */
{
	int pc_n=(pcn*len)/100;
	static int del=0;
	int pc=pco, tic;
	int i,j;
	pco=(pco*len)/100;
	if (!(pvcraze_verbose&VERB_PROGRESS))
		return pc;
	if (pco==len)
		return pc;
	
	tic=len/tics;
	
	for (i=pco+1;i<=pc_n;i++)
	{
		if(i%tic == 0)
		{
			if (i==0)
				printf("|");
			else
			{
				for (j=0;j<del;j++)
					printf("\b");
				printf("|");
			}
		}
		else
		{
			for (j=0;j<del;j++)
				printf("\b");
			printf("=");
		}
		pc=pcn;
		printf("%3i%%",pc);
		del=4;
	}
	if (pcn>pc)
	{
		for (j=0;j<del;j++)
			printf("\b");
		printf("%3i%%",pcn);
		del=4;
		pc=pcn;	
	}
		
	if (pc_n==len)
	{
		printf("\n");	
		del=0;
	}
	fflush(stdout);
	return pc;
}
#define LLEN 65
#define TICS 5


// initialize settings
ProcessSettings pv_craze_initsettings()
{
	ProcessSettings S;
	
	S.analyze=def_settings; // curve analysis settings (see above and curve.c/h)
	S.L.N=0;			// fit model list length
	return S;
}

// add a fit model
void pv_craze_addmodeltosettings(ProcessSettings *S, const char *str)
{
	int n;
	n=S->L.N;
	AddModel(&(S->L),str);	// see fitmodels.c/h
	if (S->L.N>n)
	{
		if (!(pvcraze_verbose&VERB_HEAD))
			return;
		printf("Add Model:\n");
		printf("Model Nr. %d = %s\n", n, str);
	}
	else
		fprintf(stderr, "**Failed to add Model: %s**\n", str);
}

// free settings
void pv_craze_freesettings(ProcessSettings *S)
{
	FreeList(&(S->L));
}

// top level analysis routine
// takes a single curve and returns a CurveProp given some settings
CurveProp pv_craze_analyze(curve C, ProcessSettings *S, double *tana)
{
	CurveProp Prop={0};
	int j, k;
	curve Cf;
	double t;

	if (S->L.N>0)
	{
		Prop.E=Emalloc(S->L.N*sizeof(double));
		if (Prop.E)
		{
			if((Prop.P=Emalloc(S->L.N*sizeof(generic_pars)))==NULL)
			{				
				free(Prop.E);
				Prop.E=NULL;
				Prop.N=0;
			}
			else
				Prop.N=S->L.N;
		}
		else
		{
			Prop.E=NULL;
			Prop.P=NULL;
			Prop.N=0;
		}
	}
	else
	{
		Prop.E=NULL;
		Prop.P=NULL;
		Prop.N=0;
	}
	Prop.A_FLAGS=FLAG_CODE(FLAG_NOFLAG);
	TIC(&t);
	Prop.par=AnalyzeLSQ(C, &Cf, &Prop.A_FLAGS, S->analyze);
	
	(*tana)+=TOC(&t);
	// fit the various models
	FitModels(&(S->L), Prop.par,Cf, 0); 
	free(Cf.P);
	Cf.P=NULL;
	Cf.N=0;
	for (j=0;j<Prop.N;j++)
	{
		for (k=0;k<pv_craze_numfields(S->L.List[j].M);k++)
			Prop.P[j].pars[k]=S->L.pars[j].pars[k];
		Prop.P[j].TYPE=S->L.pars[j].TYPE;
		Prop.E[j]=S->L.Efit[j];
	}
	return Prop;
}

// free curve property struct
void pv_craze_freeprop(CurveProp *P)
{
	free(P->E);
	free(P->P);
}

#define NROWLINE 4
#define NCOLLINE 4

// summarize analysis results for a list of curve properties
void Report(CurveProp *Prop, int N)
{	
	int i, j, k, l, e;
	FLAGLIST FlagList;
	char *WLAB="|Warnings";
	char *ELAB="Errors";
	char *s=WLAB;
	if (!(pvcraze_verbose&VERB_REPORT))
		return;
	
	FlagList=InitFlagList();
	if (FlagList.F==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in Report, no report can be produced\n");
		return;
	}
	for (i=0;i<N;i++)
	{
		if (AddToList(&FlagList, Prop[i].A_FLAGS))	
		{
			free(FlagList.F);
			fprintf(stderr,"INFO: Memory allocation failure in Report, no report can be produced\n");
			return;
		}
	}	
	SortErrWarn(&FlagList);
	/* create result table */
	printf("\n************************************************************************\n");
	printf(  "*                            Analysis Report:                          *\n");
	printf(  "************************************************************************\n\n");
	
	printf(  "%d IV characteristics were analyzed\n", N);
	
	printf(  "The Flag Table (all available flags):\n");
	printf(  "------------------------------------------------------------------------\n");
	printf(  "Flag | Code       | Err |                 Error Message                 \n");
	printf(  "-----+------------+-----+-----------------------------------------------\n");
	for (i=0;i<NIVF;i++)
	{
		char c;
		int m;
		m=i%52;
		if (m<26)
			c='a'+m;
		else
			c='A'+m-26;
		if (i>=NIVWARN)
			printf("%c    | 0x%08x | yes | %s\n", c, FLAG_CODE(i), Fmessage(i));
		else
			printf("%c    | 0x%08x | no  | %s\n", c, FLAG_CODE(i), Fmessage(i));
	}
	for (;i<NIVF;i++)
	{
		char c;
		int m;
		m=i%52;
		if (m<26)
			c='a'+m;
		else
			c='A'+m-26;
		printf("%c    | 0x%08x | no  | %s\n", c, FLAG_CODE(i), Fmessage(i));
	}
		
	printf(  "-----+------------+---+-------------------------------------------------\n\n");
	printf(  "Result Table:          ");
	i=0;
	for (j=1;j<NIVF;j++)
	{
		if (j==NIVWARN)
		{
			printf("||");
			s=ELAB;
			i=0;
		}
		else if (j%NCOLLINE==1)
		{
			if (s[i])
			{
				putchar(s[i]);
				i++;
			}
			else
				putchar(' ');
				
		}	
		if (s[i])
		{
			putchar(s[i]);
			i++;
		}
		else
			putchar(' ');
	}
	putchar('\n');
	
	printf(  "---------------------+-");
	for (j=1;j<NIVF;j++)
	{
		if (j==NIVWARN)
			printf("++");
		else if (j%NCOLLINE==1)
			printf("+");	
		printf("-");		
	}
	printf(  "\n");
	printf(  " Count  | Code       ");
	for (i=0;i<NIVF;i++)
	{
		char c;
		int m;
		m=i%52;
		if (m<26)
			c='a'+m;
		else
			c='A'+m-26;	
		if (i==NIVWARN)
			printf("||");
		else if ((i%NCOLLINE==1)||(i==0))
			printf("|");
		printf("%c", c);
	}
	printf("\n");
	k=0; // flaged curve counter
	l=0; // unique flag counter
	e=0; // error counter
	for (i=0;i<FlagList.N;i++)
	{
		if (l%NROWLINE==0)
		{
			printf("--------+------------");
			printf("+-");
			for (j=1;j<NIVF;j++)
			{
				if (j==NIVWARN)
					printf("++");
				else if (j%NCOLLINE==1)
					printf("+");	
				printf("-");		
			}
			printf("\n");
		}		
		printf("%8d| 0x%08x ",FlagList.F[i].count, FlagList.F[i].flag);
		if (FlagList.F[i].flag)
		{
			k+=FlagList.F[i].count;
			if (IsError(FlagList.F[i].flag))
				e+=FlagList.F[i].count;
		}
		if (FlagList.F[i].flag==0)
			printf("|X");
		else
			printf("| ");
		for (j=1;j<NIVF;j++)
		{
			if (j==NIVWARN)
				printf("||");
			else if (j%NCOLLINE==1)
				printf("|");	
			if (ISFLAG(FlagList.F[i].flag,j))
			{
				if (j>=NIVWARN)
					printf("!");
				else
					printf("X");
			}
			else
				printf(" ");
		}
		printf("\n");	
		l++;	
	}
	printf("--------+------------+-");
	for (i=1;i<NIVF;i++)
	{
		if (i==NIVWARN)
			printf("++");
		else if (i%NCOLLINE==1)
			printf("+");	
		printf("-");		
	}		
	printf("\n");
	printf("Total %d unique combinations of flags\n", l);			
	printf("Total %d IV curves with flags (%.1f %%)\n", k, 100.0*((double)k)/((double)N));
	printf("Total %d IV curves with errors (%.1f %%)\n", e, 100.0*((double)e)/((double)N));
	printf("------------------------------------------------------------------------\n\n");
	free(FlagList.F);
}

#define HLLEN 70
void pvcraze_PrintCentered(int LL, char *s)
{
	int l, left, right;
	l=strlen(s);
	if (l==0)
	{
		putchar('\n');
		return;
	}	
	if (l>=LL)
	{
		printf("%s\n",s);
		return;
	}
	left=(LL-l)/2;
	right=LL-left-l;
	putchar('>');
	for (l=1;l<left;l++)
		putchar(' ');
	printf("%s",s);
	for (l=0;l<right-1;l++)
		putchar(' ');
	putchar('<');
	putchar('\n');
	
}
void pvcraze_PrintHead(int N, char *title)
{
	char line[HLLEN+1];
	char task[HLLEN+1];
	int i;
	for (i=0;i<HLLEN;i++)
		line[i]='*';
	line[i]='\0';
	snprintf(task,HLLEN,"Processing %d IV Characteristics", N);
	pvcraze_PrintCentered(HLLEN, line);
	pvcraze_PrintCentered(HLLEN, title);
	pvcraze_PrintCentered(HLLEN, task);
	pvcraze_PrintCentered(HLLEN, line);
	putchar('\n');
}

// do a batch analysis for a list of curves
CurveProp * pv_craze_analyze_list(curve *C, int N, ProcessSettings *S)
{
	int i, pc=0;
	CurveProp *Prop;
	double ttotal=0, trest;
	double tana=0;
	
	for (i=0;i<S->L.N;i++)
		S->L.tm[i]=0; // reset timers
	
	if((Prop=Emalloc(N*sizeof(CurveProp)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in analyze list, no curves can be analyzed\n");
		return Prop;
	}
		
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE)");
		
	BarScale(LLEN, TICS);
	
	TIC(&ttotal);
	for (i=0;i<N;i++)
	{
		Prop[i]=pv_craze_analyze(C[i], S, &tana);
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}	
	pc=ProgressBar(100, pc, LLEN, TICS);
	ttotal=TOC(&ttotal);
	
	if (pvcraze_verbose&VERB_TIMING)
	{
		printf(  "------------------------------------------------------------------------\n");
		printf(  "Curves Analyzed: %d\n", N);
		printf(  "Total Time (s):  %.2f\n", ttotal);
		printf(  "Time/Datum (us): %.2f\n", 1000000.0*ttotal/(double)N);
		printf(  "Breakdown:\n");
		printf(  "------------------------------------------------------------------------\n");
		printf(  "Curve Analysis:                                  %12.3e s (%4.1f %%)\n", tana, 100*tana/ttotal);
		trest=ttotal-tana;
		for (i=0;i<S->L.N;i++)
		{
			printf(  "Model %2d:                                        %12.3e s (%4.1f %%)\n", i, S->L.tm[i], 100*S->L.tm[i]/ttotal);
			trest-=S->L.tm[i];
		}
		printf(  "Rest:                                            %12.3e s (%4.1f %%)\n", trest, 100*trest/ttotal);
		printf(  "------------------------------------------------------------------------\n");
	}
	Report(Prop, N);
	return Prop;
}

CurveProp pv_craze_modelfit(curvepars par, ProcessSettings *S)
{
	CurveProp Prop;
	int j, k;
	curve Cf={0};
	
	if (S->L.N>0)
	{
		Prop.E=Emalloc(S->L.N*sizeof(double));
		if (Prop.E)
		{
			if((Prop.P=Emalloc(S->L.N*sizeof(generic_pars)))==NULL)
			{				
				free(Prop.E);
				Prop.E=NULL;
				Prop.N=0;
			}
			else
				Prop.N=S->L.N;
		}
		else
		{
			Prop.E=NULL;
			Prop.P=NULL;// is it useful to continue, now that memory allocation fails?
			Prop.N=0;
		}
	}
	else
	{
		Prop.E=NULL;
		Prop.P=NULL;
		Prop.N=0; // is it useful to continue?
	}
	Prop.A_FLAGS=FLAG_CODE(FLAG_NOFLAG);
	Prop.par=par;
	// fit the various models
	FitModels(&(S->L), Prop.par,Cf, 0); 
	for (j=0;j<Prop.N;j++)
	{
		for (k=0;k<pv_craze_numfields(S->L.List[j].M);k++)
			Prop.P[j].pars[k]=S->L.pars[j].pars[k];
			
		Prop.P[j].TYPE=S->L.pars[j].TYPE;
		Prop.E[j]=S->L.Efit[j];
	}
	return Prop;
}


// do a batch analysis for a list of epars
CurveProp * pv_craze_modelfit_list(curvepars *par, int N, ProcessSettings *S)
{
	int i, pc=0;
	CurveProp *Prop;
	double ttotal=0, trest;
	
	
	for (i=0;i<S->L.N;i++)
		S->L.tm[i]=0; // reset timers
	
	if((Prop=Emalloc(N*sizeof(CurveProp)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in model fit list, no curves can be analyzed\n");
		return Prop;		
	}
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - Model Fit");
		
	BarScale(LLEN, TICS);
	
	TIC(&ttotal);
	for (i=0;i<N;i++)
	{
		Prop[i]=pv_craze_modelfit(par[i], S);
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}	
	pc=ProgressBar(100, pc, LLEN, TICS);
	ttotal=TOC(&ttotal);
	
	if (pvcraze_verbose&VERB_TIMING)
	{
		printf(  "------------------------------------------------------------------------\n");
		printf(  "Curves Analyzed: %d\n", N);
		printf(  "Total Time (s):  %.2f\n", ttotal);
		printf(  "Time/Datum (us): %.2f\n", 1000000.0*ttotal/(double)N);
		printf(  "Breakdown:\n");
		printf(  "------------------------------------------------------------------------\n");
		trest=ttotal;
		for (i=0;i<S->L.N;i++)
		{
			printf(  "Model %2d:                                        %12.3e s (%4.1f %%)\n", i, S->L.tm[i], 100*S->L.tm[i]/ttotal);
			trest-=S->L.tm[i];
		}
		printf(  "Rest:                                            %12.3e s (%4.1f %%)\n", trest, 100*trest/ttotal);
		printf(  "------------------------------------------------------------------------\n");
	}
	return Prop;
}

/* Compute an error 
 * computes the RMS euclidean distance of the curve points to the model
 */
double *pv_craze_error(generic_pars *pars, curve *C, int N)
{
	int i, pc=0;
	double *d;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - Error");
	BarScale(LLEN, TICS);
	if((d=Emalloc(N*sizeof(double)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in error, no errors can be computed\n");
		return d;
	}
	for (i=0;i<N;i++)
	{
		d[i]=pv_craze_g_rmse(pars[i], C[i]);
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}
	pc=ProgressBar(100, pc, LLEN, TICS);
	return d;
}

/* Compute an error curve 
 * computed the euclidean distance of the curve points to the model
 * curve in the normalized current-voltage space
 * Produces a voltage dependent distance
 */
curve *pv_craze_errorcurve(generic_pars *pars, curve *C, int N)
{
	int i, pc=0;
	curve *Co;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - Error Curve");
	BarScale(LLEN, TICS);
	if((Co=Emalloc(N*sizeof(curve)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in error curve, no errors can be computed\n");
		return Co;
	}
	for (i=0;i<N;i++)
	{
		Co[i]=pv_craze_g_distcurve(pars[i], C[i]);
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}
	pc=ProgressBar(100, pc, LLEN, TICS);
	return Co;
}
// reconstruct IV curves from curve model parameters
curve *pv_craze_curvemodel(generic_pars *pars, int N, int M)
{
	int i, pc=0;
	curve *C;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - Curve Model");
		
	BarScale(LLEN, TICS);
	if((C=Emalloc(N*sizeof(curve)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in curve model, no curves can be computed\n");
		return C;
	}
	for (i=0;i<N;i++)
	{
		C[i]=pv_craze_g_tocurve(pars[i], M);
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}
	pc=ProgressBar(100, pc, LLEN, TICS);
	return C;
}
// translate curves using the IEC60891 procedure 2 
/* takes N curves at N conditions and translates them to M conditions
 * This can be used to
 * 	1: set N=n and M=1
 * 		translate n IVs from n constions to 1 standard condition
 *  2: set N=1 and M=m
 * 		translate one IV to m conditions
 */
curvepars *pv_craze_translate_IEC60891_v2(curvepars *pars, double *G1, double *T1, int N, double *G2, double *T2, int M, 
                   double alpha, double beta, double a, double Rs, double kappa)
{
	int i, pc=0, nn;
	curvepars *Cout;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - IEC60891 v2 (2009) Translate");
		
	BarScale(LLEN, TICS);
	nn=N;
	if (M>nn)
		nn=M;
	if((Cout=Emalloc(nn*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in EC60891 v2 translation, no curves can be translated\n");
		return Cout;
	}
	for (i=0;i<nn;i++)
	{
		generic_pars P, O;
		P.TYPE=M_PVSPLINE;
		P.pars[0]=pars[i%N].Voc;
		P.pars[1]=pars[i%N].Isc;
		P.pars[2]=pars[i%N].Vmpp;
		P.pars[3]=pars[i%N].Impp;
		P.pars[4]=pars[i%N].Vqmp1;
		P.pars[5]=pars[i%N].Iqmp1;
		P.pars[6]=pars[i%N].Vqmp2;
		P.pars[7]=pars[i%N].Iqmp2;
		P.pars[8]=pars[i%N].Gsc;
		P.pars[9]=pars[i%N].Roc;
		O=IEC60891_v2(P, alpha, beta, a, Rs, kappa, T1[i%N], T2[i%M], G1[i%N], G2[i%M], &DefaultMODSET);
		Cout[i].Voc=O.pars[0];
		Cout[i].Isc=O.pars[1];
		Cout[i].Vmpp=O.pars[2];
		Cout[i].Impp=O.pars[3];
		Cout[i].Vqmp1=O.pars[4];
		Cout[i].Iqmp1=O.pars[5];
		Cout[i].Vqmp2=O.pars[6];
		Cout[i].Iqmp2=O.pars[7];
		Cout[i].Gsc=O.pars[8];
		Cout[i].Roc=O.pars[9];	
		/*printf("---------------------------------\n");	
		printf("%e -> %e %e -> %e\n", T1[i%N], T2[i%M], G1[i%N], G2[i%M]);	
		printf(" %e -> %e\n", pars[i%N].Voc, Cout[i].Voc);
		printf(" %e -> %e\n", pars[i%N].Isc, Cout[i].Isc);
		printf(" %e -> %e\n", pars[i%N].Vmpp, Cout[i].Vmpp);
		printf(" %e -> %e\n", pars[i%N].Impp, Cout[i].Impp);
		printf(" %e -> %e\n", pars[i%N].Vqmp1, Cout[i].Vqmp1);
		printf(" %e -> %e\n", pars[i%N].Iqmp1, Cout[i].Iqmp1);
		printf(" %e -> %e\n", pars[i%N].Vqmp2, Cout[i].Vqmp2);
		printf(" %e -> %e\n", pars[i%N].Iqmp2, Cout[i].Iqmp2);
		printf(" %e -> %e\n", pars[i%N].Gsc, Cout[i].Gsc);
		printf(" %e -> %e\n", pars[i%N].Roc, Cout[i].Roc);
		printf("---------------------------------\n");*/
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}
	pc=ProgressBar(100, pc, LLEN, TICS);
	return Cout;
}



curvepars *pv_craze_translate_IEC60891_v3(curvepars *pars, double *G1, double *T1, int N, double *G2, double *T2, int M, 
                   double Vocstd, double Tstd, double Gstd, double alpha, double beta, double B1, double B2, double Rs, double kappa)
{
	int i, pc=0, nn;
	curvepars *Cout;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - IEC60891 v3 (2021) Translate");
		
	BarScale(LLEN, TICS);
	nn=N;
	if (M>nn)
		nn=M;
	if((Cout=Emalloc(nn*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr,"INFO: Memory allocation failure in EC60891 v3 translation, no curves can be translated\n");
		return Cout;
	}
	for (i=0;i<nn;i++)
	{
		generic_pars P, O;
		P.TYPE=M_PVSPLINE;
		P.pars[0]=pars[i%N].Voc;
		P.pars[1]=pars[i%N].Isc;
		P.pars[2]=pars[i%N].Vmpp;
		P.pars[3]=pars[i%N].Impp;
		P.pars[4]=pars[i%N].Vqmp1;
		P.pars[5]=pars[i%N].Iqmp1;
		P.pars[6]=pars[i%N].Vqmp2;
		P.pars[7]=pars[i%N].Iqmp2;
		P.pars[8]=pars[i%N].Gsc;
		P.pars[9]=pars[i%N].Roc;
		O=IEC60891_v3(P, alpha, beta, B1, B2, Rs, kappa, Vocstd, T1[i%N], T2[i%M], Tstd, G1[i%N], G2[i%M], Gstd, &DefaultMODSET);
		Cout[i].Voc=O.pars[0];
		Cout[i].Isc=O.pars[1];
		Cout[i].Vmpp=O.pars[2];
		Cout[i].Impp=O.pars[3];
		Cout[i].Vqmp1=O.pars[4];
		Cout[i].Iqmp1=O.pars[5];
		Cout[i].Vqmp2=O.pars[6];
		Cout[i].Iqmp2=O.pars[7];
		Cout[i].Gsc=O.pars[8];
		Cout[i].Roc=O.pars[9];	
		/*printf("---------------------------------\n");	
		printf("%e -> %e %e -> %e\n", T1[i%N], T2[i%M], G1[i%N], G2[i%M]);	
		printf(" %e -> %e\n", pars[i%N].Voc, Cout[i].Voc);
		printf(" %e -> %e\n", pars[i%N].Isc, Cout[i].Isc);
		printf(" %e -> %e\n", pars[i%N].Vmpp, Cout[i].Vmpp);
		printf(" %e -> %e\n", pars[i%N].Impp, Cout[i].Impp);
		printf(" %e -> %e\n", pars[i%N].Vqmp1, Cout[i].Vqmp1);
		printf(" %e -> %e\n", pars[i%N].Iqmp1, Cout[i].Iqmp1);
		printf(" %e -> %e\n", pars[i%N].Vqmp2, Cout[i].Vqmp2);
		printf(" %e -> %e\n", pars[i%N].Iqmp2, Cout[i].Iqmp2);
		printf(" %e -> %e\n", pars[i%N].Gsc, Cout[i].Gsc);
		printf(" %e -> %e\n", pars[i%N].Roc, Cout[i].Roc);
		printf("---------------------------------\n");*/
		pc=ProgressBar((100*i)/N, pc, LLEN, TICS);
	}
	pc=ProgressBar(100, pc, LLEN, TICS);
	return Cout;
}

MODELSET FitIECMODSET = {
	M_NONE,
	1000,				// max iteration for basic fitting algorithm (if applicale)
	1e-5,				// eps for basic fitting optimizer (if applicale)
	50,					// iterations for model evaluation (if applicale)
	100,				// population size (if applicale)
	1e-5,				// model evaluation eps (if applicale)
	1e-2,				// oeps eps threshold for expensive optim,ization (if applicable)
	false,				// refine solar cell parameter during fit (if applicable)
};

curvepars pv_craze_fit_IEC60891_v2(curvepars *pars, double *G, double *T, int N, double Tstd, double Gstd,
                   double *alpha, double *beta, double *a, double *Rs, double *kappa)
{
	curvepars Out;
	generic_pars P;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - IEC60891 v2 (2009) Fit");
	P=IEC60891Fit_v2(pars, G, T, N, Tstd, Gstd, alpha, beta, a, Rs, kappa, &FitIECMODSET);
	Out.Voc  =P.pars[0];
	Out.Isc  =P.pars[1];
	Out.Vmpp =P.pars[2];
	Out.Impp =P.pars[3];
	Out.Vqmp1=P.pars[4];
	Out.Iqmp1=P.pars[5];
	Out.Vqmp2=P.pars[6];
	Out.Iqmp2=P.pars[7];
	Out.Gsc  =P.pars[8];
	Out.Roc  =P.pars[9];	

	return Out;
}


curvepars pv_craze_fit_IEC60891_v3(curvepars *pars, double *G, double *T, int N, double Tstd, double Gstd,
                   double *alpha, double *beta, double *B1, double *B2, double *Rs, double *kappa)
{
	curvepars Out;
	generic_pars P;
	if (pvcraze_verbose&VERB_HEAD)
		pvcraze_PrintHead(N, "Photo-Voltaic CuRve AnalyZEr (PV-CRAZE) - IEC60891 v3 (2021) Fit");
	P=IEC60891Fit_v3(pars, G, T, N, Tstd, Gstd, alpha, beta, B1, B2, Rs, kappa, &FitIECMODSET);
	Out.Voc  =P.pars[0];
	Out.Isc  =P.pars[1];
	Out.Vmpp =P.pars[2];
	Out.Impp =P.pars[3];
	Out.Vqmp1=P.pars[4];
	Out.Iqmp1=P.pars[5];
	Out.Vqmp2=P.pars[6];
	Out.Iqmp2=P.pars[7];
	Out.Gsc  =P.pars[8];
	Out.Roc  =P.pars[9];	

	return Out;
}
