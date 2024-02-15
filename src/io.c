/*
    io module for the commandline PV-CRAZE interface
    Copyright (C) 2020  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2020                                   *
 *                                                               *             
 *****************************************************************/                                                                             

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "emalloc.h"
#include "curve.h"
#include "fitmodels.h"
#include "io.h"

typedef enum {ILIST, IRANGE, IORANGE, INONE} ITYPE;
typedef struct {
	ITYPE T;
	int *list;
	int N;
	int Start, Step, End;
} ISpec;

char *StripSep(char *str, char sep, int *e)
// separate string at seperation characters, set e to 1 if a separation character is found
{
	char *p;
	p=str;
	(*e)=0;
	while (*p)
	{
		if (*p==sep)
		{
			*p='\0';
			if (sep=='\0')
				return p;
			else
			{
				(*e)=1;
				return p+1;
			}
		}
		p++;
	}
	return p;
}
void FreeMat(MatrixData *M)
{
	int i;
	for (i=0;i<M->n;i++)
		free(M->data[i]);
	free(M->data);
	M->n=0;
	M->m=0;
	M->data=NULL;
}
		

MatrixData ReadMatrix(char *fn, char sep)
{
	int n=0,m=0;
	FILE *f;
	char *word, *cur, *next, ch;
	int Nac=50, Nar=50, i, e;
	char ssep[6];
	MatrixData M={0};
	
	M.m=0;
	M.n=0;
	if((M.data=Emalloc(Nar*sizeof(double *)))==NULL)
		return M;
	M.data[0]=Emalloc(Nac*sizeof(double));
	{
		free(M.data);
		M.data=NULL;
		return M;
	}	
	if ((f=fopen(fn, "r"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for reading\n",fn);
		free(M.data[0]);
		free(M.data);
		M.data=NULL;
		return M;
	}
	
	snprintf(ssep,6,"%%1[%c]", sep);
	
    while ((feof(f)==0))
    {
		if (fscanf(f,"%[#]",&ch)==1)
			while ((feof(f)==0)&&(ch!='\n'))
				fscanf(f,"%c",&ch);
		else
		{
			if (fscanf(f,"%ms", &word)==1)
			{
				cur=word;
				printf("word %s\n", word);
				while ((*cur)!='\0')
				{
					next=StripSep(cur, sep, &e);
					if (*cur=='\0')
						M.data[n][m++]=NAN;
					else
						M.data[n][m++]=atof(cur);
					cur=next;
					if (m==Nac)
					{
						Nac+=50;
						/* reallocate all lines (including this one)*/
						for (i=0;i<=n;i++)
						{
							double *tmp;
							if((tmp=Erealloc(M.data[i],Nac*sizeof(double)))==NULL)
							{									
								for (i=0;i<=n;i++) // clean all up
									free(M.data[i]);
								free(M.data);
								M.data=NULL;
								M.n=0;
								M.m=0;
								fclose(f);
								return M;
							}
							M.data[i]=tmp;
						}
					}
					if (m>M.m)
					{
						M.m=m;
						/* add NAN's to previous lines */
						for (i=0;i<n;i++)
							M.data[i][m-1]=NAN;
					}
				}
				free(word);
			}
			if (!e&&sep) /* we need to read white space till the seperation sign */
			{
				while (fscanf(f,"%1[ \t]",&ch)==1); /* only read separation characters */	
				fscanf(f,ssep,&ch) ;// if nothing is read then the white space will substitute the separation
			} 
			
			printf("scanning whistepace\n");
			while (fscanf(f,"%1[ \t\r]",&ch)==1); /* only read separationb characters */				
			fscanf(f,"%1[\n]",&ch); /* only read endlines */
			if (ch=='\n')
			{
				printf("found an endline\n");
				if (m>0)
				{
					for (;m<M.m;m++) /* pad line with NAN's */
						M.data[n][m]=NAN;
					n++;
					if((M.data[n]=Emalloc(Nac*sizeof(double)))==NULL)
						return M;	
					
					if (n>M.n)
					{
						M.n=n;
						if (n==Nar-1)
						{
							double **tmp;
							Nar+=50;
							if((tmp=Erealloc(M.data, Nar*sizeof(double *)))==NULL)
							{
								M.n--;
								return M;
							}
							M.data=tmp;
						}
					}
					m=0;
				}
			}
		}
	}
	if (m==0)
		free(M.data[n]);
	else
	{
		for (;m<M.m;m++) /* pad line with NAN's */
			M.data[n][m]=NAN;
		n++;
	}
	fclose(f);
	return M;
}

int IsInSpec(int i, ISpec S)
{
	int s,e, st;
	if (S.End>S.Start)
	{
		s=S.Start;
		e=S.End;
	}
	else
	{
		s=S.End;
		e=S.Start;
	}
	st=abs(S.Step);
	
	switch(S.T)
	{
		case ILIST:
		{
			int j;
			for (j=0;j<S.N;j++)
				if (S.list[j]==i)
					return 1;
			return 0;
		}
		case IRANGE:
			if (i>e)
				return 0; // else go to IORANGE
		case IORANGE:
			i-=s;
			if (i>=0)
				return (i%st==0);
			return 0;
		case INONE:
		default:
			return 0;
	}
}

void PrintSpec(ISpec S)
{
	printf("Spec\n");
	switch(S.T)
	{
		case ILIST:
		{
			int j;
			printf("list: ");
			for (j=0;j<S.N;j++)
				printf("%d, ", S.list[j]);
			printf("\n");
			return;
		}
		case IRANGE:
			printf("range: %d:%d:%d\n", S.Start,S.Step,S.End);
			return;
		case IORANGE:
			printf("open range: %d:%d:end\n", S.Start,S.Step);
			return;
		case INONE:
		default:
			printf("none\n");
			return;
	}
}

ISpec MakeSpec(const char *spec)
{
	ISpec S={0};
	char *p;
	int i, j;
	char *s, *p2;
	
	if((s=Emalloc(1024*sizeof(char)))==NULL)
		return S;
	i=strlen(spec);
	if (i>1023)
	{
		fprintf(stderr, "Error: colspec string too long\n");
		exit(1);
	}
	strncpy(s,spec,1023);
	s[1023]='\0';
	
	S.list=NULL;
	S.N=0;
	S.Start=0;
	S.Step=1;
	S.End=0;
	
	// Test pattern Start:Step:End
	p=s;
	i=0;
	j=0;
	while (*p)
	{
		if (*p==':')
		{
			i++;
			p++;
		}
		else if (*p==',')
		{
			j++;
			p++;
		}
		else
		{
			if (strncmp(p,"end",3)==0)
				p+=3;
			else
			{
				if (((*p>='0')&&(*p<='9'))||(*p=='-'))
					p++;
				else
				{
					fprintf(stderr, "Error: invalid ISpec string 0 %s\n", s);
					free(s);
					return S;
				}
			}
		}
	}
	if (i==0)
	{
		S.T=ILIST;
		S.N=j+1;
		if((S.list=Emalloc(S.N*sizeof(int)))==NULL)
		{
			free(s);
			return S;
		}
		p=s;
		p2=p;
		i=0;
		while (*p2)
		{
			if (*p2==',')
			{
				*p2='\0';
				S.list[i++]=atoi(p);
				p=p2+1;
			}
			p2++;
		}
		S.list[i++]=atoi(p);
		if (S.list[i-1]<0)
		{
			free(S.list);
			S.list=NULL;
			S.N=0;
			S.Start=0;
			S.Step=1;
			S.End=0;
			fprintf(stderr, "Error: a row/column index is negative\n");
		}
		free(s);
		return S;
	}
	p=s;
	i=0;
	while (*p)
	{
		if (*p==':')
		{
			*p='\0';
			S.Start=atoi(s);
			if (S.Start<0)
			{
				S.list=NULL;
				S.N=0;
				S.Start=0;
				S.Step=1;
				S.End=0;
				fprintf(stderr, "Error: Start row/column index is negative\n");
				free(s);
				return S;
			}
			
			p++;
			break;
		}
		if ((*p>='0')&&(*p<='9')) // no negative start index allowed
			p++;
		else
		{
			fprintf(stderr, "Error: invalid ISpec string 1 %s\n",s);
			free(s);
			return S;
		}
	}
	if (!*p)
	{
		fprintf(stderr, "Error: invalid ISpec string 2 %s\n",s);
		free(s);
		return S;
	}
	p2=p;
	while (*p)
	{
		if (*p==':')
		{
			*p='\0';
			S.Step=atoi(p2);
			if (S.Step==0)
			{
				S.list=NULL;
				S.N=0;
				S.Start=0;
				S.Step=1;
				S.End=0;
				fprintf(stderr, "Error: indes step is 0\n");
				return S;
			}
			p++;
			break;
		}
		if (((*p>='0')&&(*p<='9'))||(*p=='-'))
			p++;
		else
		{
			fprintf(stderr, "Error: invalid ISpec string 3 %s\n",s);
			free(s);
			return S;
		}
	}
	if (!*p)
	{
		fprintf(stderr, "Error: invalid ISpec string 4 %s\n",s);
		free(s);
		return S;
	}
	
	if (strncmp(p,"end",3)==0)
	{
		S.T=IORANGE;
		S.End=S.Start+1;
		free(s);
		return S;
	}
	S.T=IRANGE;
	S.End=atoi(p);
	if (S.End<0)
	{
		S.list=NULL;
		S.N=0;
		S.Start=0;
		S.Step=1;
		S.End=0;
		fprintf(stderr, "Error: End row/column index is negative\n");
		free(s);
		return S;
	}
	if (S.Start>S.End) // direction only controlled by the step
	{
		i=S.Start;
		S.Start=S.End;
		S.End=i;
	}
	free(s);
	return S;
}


typedef struct {
		int list_i;
		int list_v;
} IndexList;

int IndexCompare (const void * a, const void * b)
{
    return ( (*(IndexList*)a).list_v - (*(IndexList*)b).list_v );
}

void ReorderMatrix(MatrixData *M, ISpec Col, ISpec Row)
{
	int i, j;
	// first do columns then rows;
	double dc; // dummy column value
	double *dr; // dummy row value
	
	switch (Col.T)
	{
		case ILIST:
		{
			// swap according to list order
			IndexList *L;
			if((L=Emalloc(Col.N*sizeof(IndexList)))==NULL)
				return;
			for (i=0;i<Col.N;i++)
			{
				L[i].list_i=i;
				L[i].list_v=Col.list[i];
			}
			qsort (L, Col.N, sizeof(IndexList), IndexCompare);
			for (i=0;i<M->n;i++)
				for (j=0;j<Col.N/2;j++)
				{
					dc=M->data[i][j];
					M->data[i][L[j].list_i]=M->data[i][L[Col.N-j-1].list_i];
					M->data[i][L[Col.N-j-1].list_i]=dc;
				}
			free(L);
			
			break;
		}
		case IRANGE:
		case IORANGE:
			if (Col.Step<0)
			{
				for (i=0;i<M->n;i++)
				{
					for (j=0;j<M->m/2;j++)
					{
						dc=M->data[i][j];
						M->data[i][j]=M->data[i][M->m-j-1];
						M->data[i][M->m-j-1]=dc;
					}
				}
			}
			break;
		default:
			fprintf(stderr, "Error: unknown error\n");
			exit(1);
	}
	switch (Row.T)
	{
		case ILIST:
		{
			// swap according to list order
			IndexList *L;
			if((L=Emalloc(Row.N*sizeof(IndexList)))==NULL)
				return;
			for (i=0;i<Col.N;i++)
			{
				L[i].list_i=i;
				L[i].list_v=Row.list[i];
			}
			qsort (L, Row.N, sizeof(IndexList), IndexCompare);
			for (i=0;i<M->n;i++)
				for (j=0;j<Row.N/2;j++)
				{
					dc=M->data[i][j];
					M->data[i][L[j].list_i]=M->data[i][L[Row.N-j-1].list_i];
					M->data[i][L[Row.N-j-1].list_i]=dc;
				}
			free(L);
			
			break;
		}
		case IRANGE:
		case IORANGE:
			if (Row.Step<0)
			{
				for (i=0;i<M->n/2;i++)
				{
					dr=M->data[i];
					M->data[i]=M->data[M->n-i-1];
					M->data[M->n-i-1]=dr;
				}
			}
			break;
		default:
			fprintf(stderr, "Error: unknown error\n");
			exit(1);
	}
	
}

//MatrixData ReadSpecMatrix(char *fn, char sep, ISpec Col, ISpec Row)
MatrixData ReadSpecMatrix(char *fn, char sep, const char *cspec, const char *rspec)
{
	int n=0,m=0, row, col;
	FILE *f;
	char *word, *cur, *next, ch;
	int Nac=50, Nar=50, i, e;
	char ssep[6];
	MatrixData M={0};
	ISpec Col, Row;
	
	Col=MakeSpec(cspec);
	Row=MakeSpec(rspec);
	// PrintSpec(Col);
	row=0;
	col=0;
	M.m=0;
	M.n=0;
	if((M.data=Emalloc(Nar*sizeof(double *)))==NULL)
		return M;
	if((M.data[0]=Emalloc(Nac*sizeof(double)))==NULL)
	{
		free(M.data);
		M.data=NULL;
		return M;
	}
	
	if ((f=fopen(fn, "r"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for reading\n",fn);
		return M;
	}
	
	snprintf(ssep,6,"%%1[%c]", sep);
	
    while ((feof(f)==0))
    {
		if (fscanf(f,"%[#]",&ch)==1)
			while ((feof(f)==0)&&(ch!='\n'))
				fscanf(f,"%c",&ch);
		else
		{
			if (fscanf(f,"%ms", &word)==1)
			{
				cur=word;
				while ((*cur)!='\0')
				{
					next=StripSep(cur, sep, &e);
					if (IsInSpec(col, Col))
					{
						if (*cur=='\0')
							M.data[n][m]=NAN;
						else
							M.data[n][m]=atof(cur);
						m++;						
						if (m==Nac)
						{
							Nac+=50;
							/* reallocate all lines (including this one)*/
							for (i=0;i<=n;i++)
							{
								double *tmp;
								if((tmp=Erealloc(M.data[i],Nac*sizeof(double)))==NULL)
								{									
									for (i=0;i<=n;i++) // clean all up
										free(M.data[i]);
									free(M.data);
									M.data=NULL;
									M.n=0;
									M.m=0;
									fclose(f);
									return M;
								}
								M.data[i]=tmp;
							}
						}
						if (m>M.m)
						{
							M.m=m;
							/* add NAN's to previous lines */
							for (i=0;i<n;i++)
								M.data[i][m-1]=NAN;
						}
					}
					cur=next;
					col++;
				}
				free(word);
			}
			if (!e&&sep) /* we need to read white space till the seperation sign */
			{
				while (fscanf(f,"%1[ \t]",&ch)==1); /* only read separation characters */	
				fscanf(f,ssep,&ch) ;// if nothing is read then the white space will substitute the separation
			} 
			
			while (fscanf(f,"%1[ \r\t]",&ch)==1); /* only read separationb characters */				
			fscanf(f,"%1[\n]",&ch); /* only read endlines */
			if (ch=='\n')
			{
				if (m>0)
				{
					for (;m<M.m;m++) /* pad line with NAN's */
						M.data[n][m]=NAN;
					if (IsInSpec(row, Row))
					{
						n++;
						if((M.data[n]=Emalloc(Nac*sizeof(double)))==NULL)
							return M;
							
						if (n>M.n)
						{
							M.n=n;
							if (n==Nar-1)
							{
								double **tmp;
								Nar+=50;
								if((tmp=Erealloc(M.data, Nar*sizeof(double *)))==NULL)
								{
									M.n--;
									return M;
								}
								M.data=tmp;
							}
						}
					}
					row++;
					m=0;
					col=0;
				}
			}
		}
	}
	if (m==0)
		free(M.data[n]);
	else
	{
		for (;m<M.m;m++) /* pad line with NAN's */
			M.data[n][m]=NAN;
		n++;
	}
	fclose(f);
	ReorderMatrix(&M, Col, Row);
	
	if (Col.list)
		free(Col.list);
	if (Row.list)
		free(Row.list);
	return M;
}
curve MatrixRows2Curve(MatrixData *V, MatrixData *I, int rv, int ri)
{
	curve C={0};
	int i, k=0;
	C.N=0;
	if (V->m!=I->m)
	{
		C.P=Emalloc(sizeof(point));// we do not want to free a non allocated array		
		fprintf(stderr, "Error: lengths of voltage and current matrices does not match\n");
		return C;
	}
	if ((rv<0)||(rv>=V->n))
	{
		C.P=Emalloc(sizeof(point));// we do not want to free a non allocated array		
		fprintf(stderr, "Error: voltage matrix row index out of bounds\n");
		return C;
	}
	if ((ri<0)||(ri>=I->n))
	{
		C.P=Emalloc(sizeof(point));// we do not want to free a non allocated array		
		fprintf(stderr, "Error: current matrix row index out of bounds\n");
		return C;
	}
	
	if((C.P=Emalloc(V->m*sizeof(point)))==NULL)
		return C;
		
	for (i=0;i<V->m;i++)
	{
		double vv,ii;
		vv=V->data[rv][i];
		ii=I->data[ri][i];
		if (isfinite(vv)&&isfinite(ii)) // skip all NaN's and such
		{
			C.P[k].x=vv;
			C.P[k].y=ii;
			k++;
		}
	}
	C.N=k;
	return C;
}

curve MatrixCols2Curve(MatrixData *V, MatrixData *I, int cv, int ci)
{
	curve C;
	int i, k=0;
	C.N=0;
	if (V->n!=I->n)
	{
		C.P=Emalloc(sizeof(point));// we do not want to free a non allocated array		
		fprintf(stderr, "Error: lengths of voltage and current matrices does not match\n");
		return C;
	}
	if ((cv<0)||(cv>=V->m))
	{
		C.P=Emalloc(sizeof(point));// we do not want to free a non allocated array		
		fprintf(stderr, "Error: voltage matrix column index out of bounds\n");
		return C;
	}
	if ((ci<0)||(ci>=I->m))
	{
		C.P=Emalloc(sizeof(point));// we do not want to free a non allocated array		
		fprintf(stderr, "Error: current matrix column index out of bounds\n");
		return C;
	}
	
	if((C.P=Emalloc(V->n*sizeof(point)))==NULL)
		return C;
	for (i=0;i<V->n;i++)
	{
		double vv,ii;
		vv=V->data[i][cv];
		ii=I->data[i][ci];
		if (isfinite(vv)&&isfinite(ii)) // skip all NaN's and such
		{
			C.P[k].x=vv;
			C.P[k].y=ii;
			k++;
		}
	}
	C.N=k;
	return C;
}	
#define MAXLEN 2048
void PadString(FILE *f, char *str, int N)
{
	int l;
	l=strlen(str);
	if (l>MAXLEN)
		l=MAXLEN;
	fprintf(f,"%s", str);
	// always separate
	fprintf(f," ");
	l++;
	while(l<N)
	{
		fprintf(f," ");
		l++;
	}
}



void PrintLabels(FILE *f, MODEL_TYPE T, int N, int start)
{
	int i,j;
	i=PropIndex(T);
	
	for (j=0;j<NumFields(T);j++)
	{
		if ((j==0)&&(start)&&(N>2))
			PadString(f, MODEL_PROP[i].labels[j], N-2);
		else
			PadString(f, MODEL_PROP[i].labels[j], N);
	}
}


void WriteExtrParMatrix(char *fn, MatrixData par)
{
	FILE *f;
	int i, j;
	
	if (par.m!=16)
	{
		printf("Provided data is not extracted solar cell parameters\n");
		return;
	} 
	if ((f=fopen(fn, "w"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for writing\n",fn);
		return;
	}
	fprintf(f, "# ");
	PrintLabels(f, M_NONE, 14, 1);
	PadString(f, "E-Flags", 14);
	PadString(f, "M-Flags", 14);
	PadString(f, "Convexity 1", 14);
	PadString(f, "Convexity 2", 14);
	PadString(f, "Convexity 3", 14);
	PadString(f, "Convexity 4", 14);
	fprintf(f, "\n");
	
	for (i=0;i<par.n;i++)
	{
		fprintf(f, "%12e", par.data[i][0]);	
		for (j=1;j<10;j++)
			fprintf(f, "  %12e", par.data[i][j]);	
		fprintf(f, "  %12d", (int)par.data[i][10]);
		fprintf(f, "  %12d", (int)par.data[i][11]);
		fprintf(f, "  %12e", par.data[i][12]);
		fprintf(f, "  %12e", par.data[i][13]);
		fprintf(f, "  %12e", par.data[i][14]);
		fprintf(f, "  %12e", par.data[i][15]);
		fprintf(f,"\n");
	}
			
	fclose(f);
}

void WriteModelParsMatrix(char *fn, MatrixData par, MODEL_LIST L)
{
	FILE *f;
	int i, j, k, l;
	
	if ((f=fopen(fn, "w"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for writing\n",fn);
		return;
	}
	
	fprintf(f, "# col: %-5d  ",1);
	for(i=1;i<par.m+L.N;i++)
		fprintf(f,"col: %-7d  ",i+1);
	fprintf(f, "\n# ");
	
	if (L.N)
	{
		PadString(f,"TYPE", 12);
		PrintLabels(f, L.List[0].M, 14, 0);
		PadString(f, "Err", 14);
	}
	for (i=1;i<L.N;i++)
	{
		PadString(f,"TYPE", 14);
		PrintLabels(f, L.List[i].M, 14, 0);
		PadString(f, "Err", 14);
	}
	fprintf(f, "\n");	
	
	for (i=0;i<par.n;i++)
	{
		l=0;
		for (j=0;j<L.N;j++)
		{
			fprintf(f, "%-13d ", L.List[j].M);
			for (k=0;k<NumFields(L.List[j].M);k++)
				fprintf(f, "%-13e ", par.data[i][l++]);
			fprintf(f, "%-13e ", par.data[i][l++]);//error	
		}		
		fprintf(f,"\n");
	}
			
	fclose(f);
}

void WriteParMatrix(char *fn, MatrixData par, MODEL_LIST L)
{
	FILE *f;
	int i, j;
	
	if (par.m%11!=0)
	{
		printf("Provided data is not solar cell parameters\n");
		return;
	} 
	 
	if ((f=fopen(fn, "w"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for writing\n",fn);
		return;
	}
	fprintf(f, "# col: %-5d  ",1);
	for(i=1;i<par.m;i++)
		fprintf(f,"col: %-7d  ",i+1);
	fprintf(f, "\n# ");
	
	if (L.N)
	{
		PrintLabels(f, M_NONE, 14, 1);
		PadString(f, "Err", 14);
	}
	for (i=1;i<L.N;i++)
	{
		PrintLabels(f, M_NONE, 14, 0);
		PadString(f, "Err", 14);
	}
	fprintf(f, "\n");
	
	for (i=0;i<par.n;i++)
	{
		fprintf(f, "%-13e ", par.data[i][0]);	
		for (j=1;j<par.m;j++)
			fprintf(f, "%-13e ", par.data[i][j]);
		fprintf(f,"\n");
	}
			
	fclose(f);
}

MatrixData ReadModelParMatrix(char *fn, char sep, int *M)
{
	MatrixData par;
	int k, m;
	par=ReadSpecMatrix(fn, sep, "0:1:end", "0:1:end");
	// test format
	k=0;
	m=0;
	while (k<par.m)
	{
		k+=NumFields((MODEL_TYPE)par.data[0][k])+2;
		m++;
	}
	if (k!=par.m)
		fprintf(stderr, "Warning: data-file seems no PV-CRAZE Model Parameter Matrix\n");
	(*M)=m;	
	return par;
}

generic_pars GetNextPars(MatrixData par, int r /* row */, int *c /* column */)
{
	int n, j;
	generic_pars P={{0,0,0,0,0,0,0,0,0,0}, M_NONE};
	MODEL_TYPE T;
	if (*c>=par.m)
	{
		(*c)=0;
		return P;
	}
	// c points to the TYPE column
	T=(MODEL_TYPE)par.data[r][*c];
	if ((T<=0)||(T>=M_NONE))
	{
		(*c)=0;
		fprintf(stderr, "Warning: data-matrix seems no PV-CRAZE Model Parameter Matrix\n");
		return P;
	}
	j=(*c)+1; // j points to firt parameter
	if (j+NumFields(T)>=par.m)
	{
		(*c)=0;
		fprintf(stderr, "Warning: data-matrix seems no PV-CRAZE Model Parameter Matrix\n");
		return P;
	}
	
	P.TYPE=T;
	
	for (n=0;n<NumFields(T);n++)
		P.pars[n]=par.data[r][j++];
	// j points to err field
	j++;// j points to next type field
	if (j>=par.m)
		(*c)=0;
	else
		(*c)=j;
	return P;
}

void WriteMatrix(char *fn, MatrixData M)
{
	FILE *f;
	int i, j;
	
	 
	if ((f=fopen(fn, "w"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for writing\n",fn);
		return;
	}
	
	for (i=0;i<M.n;i++)
	{
		fprintf(f, "%12e", M.data[i][0]);	
		for (j=1;j<M.m;j++)
			fprintf(f, "  %12e", M.data[i][j]);	
		fprintf(f, "\n");
	}
			
	fclose(f);
}

curve FileToCurve(char *fn, int colx, int coly)
{
	FILE *f;
	char *word, ch;
	int col=0;
	int Na=50, Nb=50, np=0;
	curve C={0};
	C.N=0;
	if((C.P=Emalloc(Na*sizeof(point)))==NULL)
		return C;
	if ((f=fopen(fn, "r"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for reading\n",fn);
		return C;
	}
    while ((feof(f)==0))
    {
		if (fscanf(f,"%[#]",&ch)==1)
			while ((feof(f)==0)&&(ch!='\n'))
				fscanf(f,"%c",&ch);
		else
		{
			if (fscanf(f,"%ms", &word)==1)
			{
				col++;
				if (col==colx)
				{
					C.P[C.N].x=atof(word);
					np++;
				}
				if (col==coly)
				{
					C.P[C.N].y=atof(word);
					np++;
				}
				free(word);
			}
			if (np<2)
			{
				while (fscanf(f,"%1[ \t]",&ch)==1); /* only read characters that are whitespace */
				fscanf(f,"%1[\n]",&ch); /* only read endlines */
				if ((ch=='\n')&&(np==1))
				{
					fprintf(stderr, "Warning: incomplete data point\n");
					col=0;
					np=0;
				}
			}
			else
			{
				C.N+=1;
				col=0;
				np=0;
				if (C.N==Na)
				{
					point *tmp;
					Na+=Nb;
					if((tmp=Erealloc(C.P, Na*sizeof(point)))==NULL)
						return C; // return what we collcted so far	
					C.P=tmp;
				}
				while ((feof(f)==0)&&(ch!='\n'))
					fscanf(f,"%c",&ch);
			}
		}
	}
	fclose(f);
	return C;
}

void WriteCurve(char *fn, curve C)
{
	FILE *f;
	int i;
	 
	if ((f=fopen(fn, "w"))==NULL)
	{
		fprintf(stderr, "Error: cannot open file %s for writing\n",fn);
		return;
	}
	
	for (i=0;i<C.N;i++)
		fprintf(f, "%e %e\n", C.P[i].x, C.P[i].y);
	fclose(f);
}
