/*
 * pvcraze._core: low level Python bindings to the PV-CRAZE library
 * Copyright (C) 2024-2026  B. E. Pieters,
 * IEK-5 Photovoltaik, Forschunszentrum Juelich
 *
 * This module is private. It deliberately does as little as possible: it
 * converts numpy arrays to the library's C structs, calls the library and
 * converts the results back. Validation of user input, broadcasting and all
 * user friendly wrapping lives in the pure Python layer (python/pvcraze/).
 *
 * Conventions used throughout:
 *   esp      (N, 10) float64 array, columns in curvepars order
 *   types    (N,)    int64 array with model IDs (MODEL_TYPE)
 *   pars     (N, 10) float64 array with model parameters
 *   curves   sequence of (n, 2) float64 arrays with (V, I) columns
 * Temperatures are always in degrees Celsius.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <limits.h>
#include <string.h>

/* library headers first: numpy pulls in complex.h, which defines the macro I
 * that clashes with parameter names in the library prototypes */
#include "lib-pv-craze.h"
#include "Error.h"

#include <numpy/arrayobject.h>

extern int pvcraze_verbose; /* verbosity bit mask, see lib-pv-craze.c */

/* The library prints with printf. Flush after each call so its output is not
 * held back in the C stdio buffer and interleaves properly with Python's. */
#define FLUSH_LIBRARY_OUTPUT() do { if (pvcraze_verbose) fflush(stdout); } while (0)

/* **************************************************
 * conversion helpers
 * **************************************************/

/* Convert obj to a C-contiguous float64 array with ndim dimensions.
 * If ncols>=0 the last dimension must have that length.
 * Returns a new reference or NULL with an exception set. */
static PyArrayObject *as_double(PyObject *obj, int ndim, npy_intp ncols, const char *name)
{
	PyArrayObject *a;
	a=(PyArrayObject *)PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (a==NULL)
		return NULL;
	if (PyArray_NDIM(a)!=ndim)
	{
		PyErr_Format(PyExc_ValueError, "%s must be %d-dimensional, got %d dimensions", name, ndim, PyArray_NDIM(a));
		Py_DECREF(a);
		return NULL;
	}
	if ((ncols>=0)&&(PyArray_DIM(a, ndim-1)!=ncols))
	{
		PyErr_Format(PyExc_ValueError, "%s must have %zd columns, got %zd", name, (Py_ssize_t)ncols, (Py_ssize_t)PyArray_DIM(a, ndim-1));
		Py_DECREF(a);
		return NULL;
	}
	if (PyArray_SIZE(a)>INT_MAX)
	{
		PyErr_Format(PyExc_ValueError, "%s is too large", name);
		Py_DECREF(a);
		return NULL;
	}
	return a;
}

/* 1D float64 array of length n */
static PyArrayObject *as_vector(PyObject *obj, npy_intp n, const char *name)
{
	PyArrayObject *a=as_double(obj, 1, -1, name);
	if (a&&(PyArray_DIM(a, 0)!=n))
	{
		PyErr_Format(PyExc_ValueError, "%s must have length %zd, got %zd", name, (Py_ssize_t)n, (Py_ssize_t)PyArray_DIM(a, 0));
		Py_DECREF(a);
		return NULL;
	}
	return a;
}

static void pars_to_row(const curvepars *p, double *r)
{
	r[0]=p->Voc;
	r[1]=p->Isc;
	r[2]=p->Vmpp;
	r[3]=p->Impp;
	r[4]=p->Vqmp1;
	r[5]=p->Iqmp1;
	r[6]=p->Vqmp2;
	r[7]=p->Iqmp2;
	r[8]=p->Gsc;
	r[9]=p->Roc;
}

static curvepars row_to_pars(const double *r)
{
	curvepars p;
	p.Voc=r[0];
	p.Isc=r[1];
	p.Vmpp=r[2];
	p.Impp=r[3];
	p.Vqmp1=r[4];
	p.Iqmp1=r[5];
	p.Vqmp2=r[6];
	p.Iqmp2=r[7];
	p.Gsc=r[8];
	p.Roc=r[9];
	return p;
}

/* esp (N,10) -> malloc'ed curvepars array. Sets *n. */
static curvepars *esp_from_obj(PyObject *obj, int *n)
{
	PyArrayObject *a;
	curvepars *p;
	const double *d;
	int i;
	if ((a=as_double(obj, 2, 10, "esp"))==NULL)
		return NULL;
	*n=(int)PyArray_DIM(a, 0);
	if ((p=malloc(((*n)>0 ? *n : 1)*sizeof(curvepars)))==NULL)
	{
		Py_DECREF(a);
		PyErr_NoMemory();
		return NULL;
	}
	d=(const double *)PyArray_DATA(a);
	for (i=0; i<*n; i++)
		p[i]=row_to_pars(d+10*i);
	Py_DECREF(a);
	return p;
}

static PyArrayObject *esp_to_array(const curvepars *p, int n)
{
	npy_intp dims[2]={n, 10};
	PyArrayObject *a;
	double *d;
	int i;
	if ((a=(PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE))==NULL)
		return NULL;
	d=(double *)PyArray_DATA(a);
	for (i=0; i<n; i++)
		pars_to_row(p+i, d+10*i);
	return a;
}

static int valid_type(long long t)
{
	return (t>=0)&&(t<(long long)M_NONE);
}

/* types (N,) + pars (N,10) -> malloc'ed generic_pars array. Sets *n. */
static generic_pars *generic_from_obj(PyObject *types_obj, PyObject *pars_obj, int *n)
{
	PyArrayObject *t, *a;
	generic_pars *g;
	const double *d;
	const npy_int64 *ti;
	int i, j;

	if ((a=as_double(pars_obj, 2, 10, "pars"))==NULL)
		return NULL;
	t=(PyArrayObject *)PyArray_FROM_OTF(types_obj, NPY_INT64, NPY_ARRAY_IN_ARRAY);
	if (t==NULL)
	{
		Py_DECREF(a);
		return NULL;
	}
	if ((PyArray_NDIM(t)!=1)||(PyArray_DIM(t, 0)!=PyArray_DIM(a, 0)))
	{
		PyErr_SetString(PyExc_ValueError, "types must be 1-dimensional with one entry per row in pars");
		Py_DECREF(a);
		Py_DECREF(t);
		return NULL;
	}
	*n=(int)PyArray_DIM(a, 0);
	if ((g=malloc(((*n)>0 ? *n : 1)*sizeof(generic_pars)))==NULL)
	{
		Py_DECREF(a);
		Py_DECREF(t);
		PyErr_NoMemory();
		return NULL;
	}
	d=(const double *)PyArray_DATA(a);
	ti=(const npy_int64 *)PyArray_DATA(t);
	for (i=0; i<*n; i++)
	{
		if (!valid_type(ti[i]))
		{
			PyErr_Format(PyExc_ValueError, "unknown model id %lld", (long long)ti[i]);
			free(g);
			Py_DECREF(a);
			Py_DECREF(t);
			return NULL;
		}
		g[i].TYPE=(MODEL_TYPE)ti[i];
		for (j=0; j<10; j++)
			g[i].pars[j]=d[10*i+j];
	}
	Py_DECREF(a);
	Py_DECREF(t);
	return g;
}

static generic_pars single_generic(int type, PyObject *pars_obj, int *ok)
{
	generic_pars g={0};
	PyArrayObject *a;
	int j;
	*ok=0;
	if (!valid_type(type))
	{
		PyErr_Format(PyExc_ValueError, "unknown model id %d", type);
		return g;
	}
	if ((a=as_vector(pars_obj, 10, "pars"))==NULL)
		return g;
	g.TYPE=(MODEL_TYPE)type;
	for (j=0; j<10; j++)
		g.pars[j]=((const double *)PyArray_DATA(a))[j];
	Py_DECREF(a);
	*ok=1;
	return g;
}

/* (types (N,), pars (N,10)) arrays from generic_pars, zeroing unused slots */
static PyObject *generic_to_tuple(const generic_pars *g, int n)
{
	npy_intp dims[2]={n, 10};
	PyArrayObject *t, *a;
	double *d;
	npy_int64 *ti;
	int i, j, nf;
	t=(PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT64);
	a=(PyArrayObject *)PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
	if ((t==NULL)||(a==NULL))
	{
		Py_XDECREF(t);
		Py_XDECREF(a);
		return NULL;
	}
	d=(double *)PyArray_DATA(a);
	ti=(npy_int64 *)PyArray_DATA(t);
	for (i=0; i<n; i++)
	{
		ti[i]=(npy_int64)g[i].TYPE;
		nf=pv_craze_numfields(g[i].TYPE);
		for (j=0; j<nf && j<10; j++)
			d[10*i+j]=g[i].pars[j];
	}
	return Py_BuildValue("(NN)", t, a);
}

static void free_curves(curve *C, int n)
{
	int i;
	if (C==NULL)
		return;
	for (i=0; i<n; i++)
		free(C[i].P);
	free(C);
}

/* sequence of (n,2) arrays -> malloc'ed curve array. Sets *n. */
static curve *curves_from_obj(PyObject *obj, int *n)
{
	PyObject *seq;
	curve *C;
	Py_ssize_t N, i;
	int j;

	if ((seq=PySequence_Fast(obj, "curves must be a sequence of (n, 2) arrays"))==NULL)
		return NULL;
	N=PySequence_Fast_GET_SIZE(seq);
	if (N>INT_MAX)
	{
		Py_DECREF(seq);
		PyErr_SetString(PyExc_ValueError, "too many curves");
		return NULL;
	}
	if ((C=calloc(N>0 ? N : 1, sizeof(curve)))==NULL)
	{
		Py_DECREF(seq);
		PyErr_NoMemory();
		return NULL;
	}
	for (i=0; i<N; i++)
	{
		PyArrayObject *a;
		const double *d;
		a=as_double(PySequence_Fast_GET_ITEM(seq, i), 2, 2, "each curve");
		if (a==NULL)
		{
			free_curves(C, (int)i);
			Py_DECREF(seq);
			return NULL;
		}
		C[i].N=(int)PyArray_DIM(a, 0);
		if ((C[i].P=malloc((C[i].N>0 ? C[i].N : 1)*sizeof(point)))==NULL)
		{
			Py_DECREF(a);
			free_curves(C, (int)i);
			Py_DECREF(seq);
			PyErr_NoMemory();
			return NULL;
		}
		d=(const double *)PyArray_DATA(a);
		for (j=0; j<C[i].N; j++)
		{
			C[i].P[j].x=d[2*j];
			C[i].P[j].y=d[2*j+1];
		}
		Py_DECREF(a);
	}
	Py_DECREF(seq);
	*n=(int)N;
	return C;
}

static PyObject *curve_to_array(curve C)
{
	npy_intp dims[2]={C.N, 2};
	PyArrayObject *a;
	double *d;
	int j;
	if ((a=(PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE))==NULL)
		return NULL;
	d=(double *)PyArray_DATA(a);
	for (j=0; j<C.N; j++)
	{
		d[2*j]=C.P[j].x;
		d[2*j+1]=C.P[j].y;
	}
	return (PyObject *)a;
}

/* list of arrays from a curve array; frees the curves */
static PyObject *curves_to_list(curve *C, int n)
{
	PyObject *list, *item;
	int i;
	if ((list=PyList_New(n))==NULL)
	{
		free_curves(C, n);
		return NULL;
	}
	for (i=0; i<n; i++)
	{
		if ((item=curve_to_array(C[i]))==NULL)
		{
			Py_DECREF(list);
			free_curves(C, n);
			return NULL;
		}
		PyList_SET_ITEM(list, i, item);
	}
	free_curves(C, n);
	return list;
}

/* add model specification strings to the settings */
static int add_models(ProcessSettings *S, PyObject *models)
{
	PyObject *seq;
	Py_ssize_t i, n;
	if ((seq=PySequence_Fast(models, "models must be a sequence of strings"))==NULL)
		return 0;
	n=PySequence_Fast_GET_SIZE(seq);
	for (i=0; i<n; i++)
	{
		PyObject *s=PySequence_Fast_GET_ITEM(seq, i);
		const char *spec;
		int before=S->L.N;
		if (!PyUnicode_Check(s))
		{
			PyErr_SetString(PyExc_TypeError, "model specifications must be strings");
			Py_DECREF(seq);
			return 0;
		}
		if ((spec=PyUnicode_AsUTF8(s))==NULL)
		{
			Py_DECREF(seq);
			return 0;
		}
		pv_craze_addmodeltosettings(S, spec);
		if (S->L.N==before)
		{
			PyErr_Format(PyExc_ValueError, "invalid model specification %R", s);
			Py_DECREF(seq);
			return 0;
		}
	}
	Py_DECREF(seq);
	return 1;
}

/* Collect per model results from CurveProp into (types (M,), pars (M,N,10), err (N,M)) */
static PyObject *props_to_fit_tuple(const ProcessSettings *S, const CurveProp *P, int n)
{
	int M=S->L.N, i, k, j, nf;
	npy_intp tdims[1]={M};
	npy_intp pdims[3]={M, n, 10};
	npy_intp edims[2]={n, M};
	PyArrayObject *types, *pars, *err;
	npy_int64 *ti;
	double *pd, *ed;

	types=(PyArrayObject *)PyArray_SimpleNew(1, tdims, NPY_INT64);
	pars=(PyArrayObject *)PyArray_ZEROS(3, pdims, NPY_DOUBLE, 0);
	err=(PyArrayObject *)PyArray_SimpleNew(2, edims, NPY_DOUBLE);
	if ((types==NULL)||(pars==NULL)||(err==NULL))
	{
		Py_XDECREF(types);
		Py_XDECREF(pars);
		Py_XDECREF(err);
		return NULL;
	}
	ti=(npy_int64 *)PyArray_DATA(types);
	pd=(double *)PyArray_DATA(pars);
	ed=(double *)PyArray_DATA(err);
	for (k=0; k<M; k++)
		ti[k]=(npy_int64)S->L.List[k].M;
	for (i=0; i<n; i++)
		for (k=0; k<M; k++)
		{
			if (k<P[i].N)
			{
				nf=pv_craze_numfields(P[i].P[k].TYPE);
				for (j=0; j<nf && j<10; j++)
					pd[(k*(npy_intp)n+i)*10+j]=P[i].P[k].pars[j];
				ed[i*M+k]=P[i].E[k];
			}
			else
				ed[i*M+k]=NAN; /* allocation failure inside the library */
		}
	return Py_BuildValue("(NNN)", types, pars, err);
}

/* **************************************************
 * library information
 * **************************************************/

static PyObject *core_model_table(PyObject *self, PyObject *noargs)
{
	PyObject *list, *labels, *item;
	int t, i, idx, nf;
	char *name;
	if ((list=PyList_New(0))==NULL)
		return NULL;
	for (t=0; t<(int)M_NONE; t++)
	{
		if ((name=pv_craze_reverselookup((MODEL_TYPE)t))==NULL)
			continue;
		idx=PropIndex((MODEL_TYPE)t);
		nf=MODEL_PROP[idx].N;
		if ((labels=PyTuple_New(nf))==NULL)
		{
			Py_DECREF(list);
			return NULL;
		}
		for (i=0; i<nf; i++)
			PyTuple_SET_ITEM(labels, i, PyUnicode_FromString(MODEL_PROP[idx].labels[i]));
		item=Py_BuildValue("(isN)", t, name, labels);
		if ((item==NULL)||(PyList_Append(list, item)<0))
		{
			Py_XDECREF(item);
			Py_DECREF(list);
			return NULL;
		}
		Py_DECREF(item);
	}
	return list;
}

static PyObject *core_flag_table(PyObject *self, PyObject *noargs)
{
	PyObject *list, *item;
	int i;
	if ((list=PyList_New(0))==NULL)
		return NULL;
	for (i=1; i<NIVF; i++)
	{
		item=Py_BuildValue("(kssO)", (unsigned long)FLAG_CODE(i), Fname(i), Fmessage(i), (i>=NIVWARN) ? Py_True : Py_False);
		if ((item==NULL)||(PyList_Append(list, item)<0))
		{
			Py_XDECREF(item);
			Py_DECREF(list);
			return NULL;
		}
		Py_DECREF(item);
	}
	return list;
}

static PyObject *core_default_analysis_settings(PyObject *self, PyObject *noargs)
{
	ProcessSettings S=pv_craze_initsettings();
	analyze_settings *a=&S.analyze;
	char am[2]={a->am, '\0'};
	return Py_BuildValue("{s:O,s:d,s:i,s:s,s:d,s:d,s:d,s:d,s:d}",
		"spike_filter", a->SpikeFilter ? Py_True : Py_False,
		"identical_eps", a->Identicals ? a->IdentEps : 0.0,
		"resample", a->Resample,
		"mpp_method", am,
		"mpp_dv", a->MPP_f_V,
		"mpp_di", a->MPP_f_I,
		"sc_dv", a->SC_f_V,
		"oc_di", a->OC_f_I,
		"min_rel_range", a->MinRelRange);
}

static PyObject *core_get_verbosity(PyObject *self, PyObject *noargs)
{
	return PyLong_FromLong(pvcraze_verbose);
}

static PyObject *core_set_verbosity(PyObject *self, PyObject *args)
{
	int mask;
	if (!PyArg_ParseTuple(args, "i", &mask))
		return NULL;
	pvcraze_verbose=mask;
	Py_RETURN_NONE;
}

/* **************************************************
 * analysis and model fitting
 * **************************************************/

static PyObject *core_analyze(PyObject *self, PyObject *args, PyObject *kwargs)
{
	static char *kwlist[]={"curves", "models", "spike_filter", "identical_eps", "resample",
		"mpp_method", "mpp_dv", "mpp_di", "sc_dv", "oc_di", "min_rel_range", NULL};
	PyObject *curves_obj, *models_obj, *fits, *ret;
	PyArrayObject *esp, *flags;
	int spike, resample, N=0, i;
	double ident_eps, mpp_dv, mpp_di, sc_dv, oc_di, min_rel;
	const char *method;
	ProcessSettings S;
	CurveProp *P;
	curve *C;
	npy_intp dims[1];
	npy_int64 *fd;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO$pdisddddd", kwlist, &curves_obj, &models_obj,
			&spike, &ident_eps, &resample, &method, &mpp_dv, &mpp_di, &sc_dv, &oc_di, &min_rel))
		return NULL;
	if (strlen(method)!=1)
	{
		PyErr_SetString(PyExc_ValueError, "mpp_method must be a single character");
		return NULL;
	}

	S=pv_craze_initsettings();
	S.analyze.SpikeFilter=spike;
	S.analyze.Identicals=(ident_eps>0);
	if (ident_eps>0)
		S.analyze.IdentEps=ident_eps;
	S.analyze.Resample=resample;
	S.analyze.am=method[0];
	S.analyze.MPP_f_V=mpp_dv;
	S.analyze.MPP_f_I=mpp_di;
	S.analyze.SC_f_V=sc_dv;
	S.analyze.OC_f_I=oc_di;
	S.analyze.MinRelRange=min_rel;
	if (!add_models(&S, models_obj))
	{
		pv_craze_freesettings(&S);
		return NULL;
	}
	if ((C=curves_from_obj(curves_obj, &N))==NULL)
	{
		pv_craze_freesettings(&S);
		return NULL;
	}

	P=NULL;
	if (N>0)
	{
		P=pv_craze_analyze(C, N, &S);
		FLUSH_LIBRARY_OUTPUT();
		if (P==NULL)
		{
			free_curves(C, N);
			pv_craze_freesettings(&S);
			return PyErr_NoMemory();
		}
	}
	free_curves(C, N);

	dims[0]=N;
	{
		npy_intp edims[2]={N, 10};
		esp=(PyArrayObject *)PyArray_SimpleNew(2, edims, NPY_DOUBLE);
	}
	flags=(PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT64);
	fits=props_to_fit_tuple(&S, P, N);
	if ((esp!=NULL)&&(flags!=NULL))
	{
		fd=(npy_int64 *)PyArray_DATA(flags);
		for (i=0; i<N; i++)
		{
			pars_to_row(&P[i].par, (double *)PyArray_DATA(esp)+10*i);
			fd[i]=(npy_int64)(unsigned int)P[i].A_FLAGS;
		}
	}
	for (i=0; i<N; i++)
		pv_craze_freeprop(&P[i]);
	free(P);
	pv_craze_freesettings(&S);
	if ((esp==NULL)||(flags==NULL)||(fits==NULL))
	{
		Py_XDECREF(esp);
		Py_XDECREF(flags);
		Py_XDECREF(fits);
		return NULL;
	}
	ret=Py_BuildValue("(NNO)", esp, flags, fits);
	Py_DECREF(fits);
	return ret;
}

static PyObject *core_fit(PyObject *self, PyObject *args)
{
	PyObject *esp_obj, *models_obj, *ret;
	ProcessSettings S;
	CurveProp *P=NULL;
	curvepars *pars;
	int N=0, i;

	if (!PyArg_ParseTuple(args, "OO", &esp_obj, &models_obj))
		return NULL;
	S=pv_craze_initsettings();
	if (!add_models(&S, models_obj))
	{
		pv_craze_freesettings(&S);
		return NULL;
	}
	if ((pars=esp_from_obj(esp_obj, &N))==NULL)
	{
		pv_craze_freesettings(&S);
		return NULL;
	}
	if (N>0)
	{
		P=pv_craze_modelfit(pars, N, &S);
		FLUSH_LIBRARY_OUTPUT();
		if (P==NULL)
		{
			free(pars);
			pv_craze_freesettings(&S);
			return PyErr_NoMemory();
		}
	}
	free(pars);
	ret=props_to_fit_tuple(&S, P, N);
	for (i=0; i<N; i++)
		pv_craze_freeprop(&P[i]);
	free(P);
	pv_craze_freesettings(&S);
	return ret;
}

/* **************************************************
 * model evaluation
 * **************************************************/

static PyObject *core_curves(PyObject *self, PyObject *args)
{
	PyObject *types_obj, *pars_obj;
	generic_pars *g;
	curve *C=NULL;
	int n, npoints;

	if (!PyArg_ParseTuple(args, "OOi", &types_obj, &pars_obj, &npoints))
		return NULL;
	if (npoints<2)
	{
		PyErr_SetString(PyExc_ValueError, "npoints must be at least 2");
		return NULL;
	}
	if ((g=generic_from_obj(types_obj, pars_obj, &n))==NULL)
		return NULL;
	if (n>0)
	{
		C=pv_craze_curvemodel(g, n, npoints);
		FLUSH_LIBRARY_OUTPUT();
		if (C==NULL)
		{
			free(g);
			return PyErr_NoMemory();
		}
	}
	free(g);
	return curves_to_list(C, n);
}

static PyObject *core_model_esp(PyObject *self, PyObject *args)
{
	PyObject *types_obj, *pars_obj, *ret;
	generic_pars *g;
	curvepars *cp;
	int n, i;

	if (!PyArg_ParseTuple(args, "OO", &types_obj, &pars_obj))
		return NULL;
	if ((g=generic_from_obj(types_obj, pars_obj, &n))==NULL)
		return NULL;
	if ((cp=malloc((n>0 ? n : 1)*sizeof(curvepars)))==NULL)
	{
		free(g);
		return PyErr_NoMemory();
	}
	for (i=0; i<n; i++)
		cp[i]=pv_craze_g_pars(g[i]);
	ret=(PyObject *)esp_to_array(cp, n);
	free(cp);
	free(g);
	return ret;
}

static PyObject *core_rmse(PyObject *self, PyObject *args)
{
	PyObject *types_obj, *pars_obj, *curves_obj;
	PyArrayObject *out;
	generic_pars *g;
	curve *C;
	double *E=NULL;
	int n, nc=0;
	npy_intp dims[1];

	if (!PyArg_ParseTuple(args, "OOO", &types_obj, &pars_obj, &curves_obj))
		return NULL;
	if ((g=generic_from_obj(types_obj, pars_obj, &n))==NULL)
		return NULL;
	if ((C=curves_from_obj(curves_obj, &nc))==NULL)
	{
		free(g);
		return NULL;
	}
	if (nc!=n)
	{
		PyErr_Format(PyExc_ValueError, "got %d parameter sets but %d curves", n, nc);
		free(g);
		free_curves(C, nc);
		return NULL;
	}
	dims[0]=n;
	if ((out=(PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE))!=NULL && n>0)
	{
		E=pv_craze_error(g, C, n);
		FLUSH_LIBRARY_OUTPUT();
		if (E==NULL)
		{
			Py_CLEAR(out);
			PyErr_NoMemory();
		}
		else
			memcpy(PyArray_DATA(out), E, n*sizeof(double));
	}
	free(E);
	free(g);
	free_curves(C, nc);
	return (PyObject *)out;
}

static PyObject *core_distance_curves(PyObject *self, PyObject *args)
{
	PyObject *types_obj, *pars_obj, *curves_obj;
	generic_pars *g;
	curve *C, *D=NULL;
	int n, nc=0;

	if (!PyArg_ParseTuple(args, "OOO", &types_obj, &pars_obj, &curves_obj))
		return NULL;
	if ((g=generic_from_obj(types_obj, pars_obj, &n))==NULL)
		return NULL;
	if ((C=curves_from_obj(curves_obj, &nc))==NULL)
	{
		free(g);
		return NULL;
	}
	if (nc!=n)
	{
		PyErr_Format(PyExc_ValueError, "got %d parameter sets but %d curves", n, nc);
		free(g);
		free_curves(C, nc);
		return NULL;
	}
	if (n>0)
	{
		D=pv_craze_errorcurve(g, C, n);
		FLUSH_LIBRARY_OUTPUT();
		if (D==NULL)
		{
			free(g);
			free_curves(C, nc);
			return PyErr_NoMemory();
		}
	}
	free(g);
	free_curves(C, nc);
	return curves_to_list(D, n);
}

/* evaluate I(V) (direction 0) or V(I) (direction 1)
 * types (N,), pars (N,10), x (N,m) -> (N,m) */
static PyObject *evaluate(PyObject *args, int direction)
{
	PyObject *types_obj, *pars_obj, *x_obj;
	PyArrayObject *x, *out;
	generic_pars *g;
	const double *xd;
	double *od;
	int n, i;
	npy_intp j, m;

	if (!PyArg_ParseTuple(args, "OOO", &types_obj, &pars_obj, &x_obj))
		return NULL;
	if ((g=generic_from_obj(types_obj, pars_obj, &n))==NULL)
		return NULL;
	if ((x=as_double(x_obj, 2, -1, direction ? "I" : "V"))==NULL)
	{
		free(g);
		return NULL;
	}
	if (PyArray_DIM(x, 0)!=n)
	{
		PyErr_SetString(PyExc_ValueError, "evaluation points must have one row per parameter set");
		free(g);
		Py_DECREF(x);
		return NULL;
	}
	m=PyArray_DIM(x, 1);
	if ((out=(PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(x), NPY_DOUBLE))==NULL)
	{
		free(g);
		Py_DECREF(x);
		return NULL;
	}
	xd=(const double *)PyArray_DATA(x);
	od=(double *)PyArray_DATA(out);
	for (i=0; i<n; i++)
		for (j=0; j<m; j++)
			od[i*m+j]=direction ? pv_craze_VI(xd[i*m+j], g[i]) : pv_craze_IV(xd[i*m+j], g[i]);
	free(g);
	Py_DECREF(x);
	return (PyObject *)out;
}

static PyObject *core_iv(PyObject *self, PyObject *args)
{
	return evaluate(args, 0);
}

static PyObject *core_vi(PyObject *self, PyObject *args)
{
	return evaluate(args, 1);
}

/* **************************************************
 * irradiance / temperature translation
 * **************************************************/

/* parse esp + four condition vectors of equal length */
typedef struct {
	curvepars *esp;
	PyArrayObject *t1, *g1, *t2, *g2;
	int n;
} translate_args;

static void free_translate_args(translate_args *a)
{
	free(a->esp);
	Py_XDECREF(a->t1);
	Py_XDECREF(a->g1);
	Py_XDECREF(a->t2);
	Py_XDECREF(a->g2);
}

static int get_translate_args(translate_args *a, PyObject *esp, PyObject *t1, PyObject *g1, PyObject *t2, PyObject *g2)
{
	memset(a, 0, sizeof(*a));
	if ((a->esp=esp_from_obj(esp, &a->n))==NULL)
		return 0;
	if (((a->t1=as_vector(t1, a->n, "T_from"))==NULL)||((a->g1=as_vector(g1, a->n, "G_from"))==NULL)||
		((a->t2=as_vector(t2, a->n, "T"))==NULL)||((a->g2=as_vector(g2, a->n, "G"))==NULL))
	{
		free_translate_args(a);
		return 0;
	}
	if (a->n==0)
	{
		PyErr_SetString(PyExc_ValueError, "nothing to translate");
		free_translate_args(a);
		return 0;
	}
	return 1;
}

#define VEC(A) ((double *)PyArray_DATA(A))

static PyObject *core_iec60891_v2_translate(PyObject *self, PyObject *args)
{
	PyObject *esp, *t1, *g1, *t2, *g2, *ret;
	double alpha, beta, a, rs, kappa;
	translate_args A;
	curvepars *out;
	if (!PyArg_ParseTuple(args, "OOOOOddddd", &esp, &t1, &g1, &t2, &g2, &alpha, &beta, &a, &rs, &kappa))
		return NULL;
	if (!get_translate_args(&A, esp, t1, g1, t2, g2))
		return NULL;
	out=pv_craze_translate_IEC60891_v2(A.esp, VEC(A.g1), VEC(A.t1), A.n, VEC(A.g2), VEC(A.t2), A.n, alpha, beta, a, rs, kappa);
	FLUSH_LIBRARY_OUTPUT();
	if (out==NULL)
	{
		free_translate_args(&A);
		return PyErr_NoMemory();
	}
	ret=(PyObject *)esp_to_array(out, A.n);
	free(out);
	free_translate_args(&A);
	return ret;
}

static PyObject *core_iec60891_v3_translate(PyObject *self, PyObject *args)
{
	PyObject *esp, *t1, *g1, *t2, *g2, *ret;
	double voc_ref, t_ref, g_ref, alpha, beta, b1, b2, rs, kappa;
	translate_args A;
	curvepars *out;
	if (!PyArg_ParseTuple(args, "OOOOOddddddddd", &esp, &t1, &g1, &t2, &g2, &voc_ref, &t_ref, &g_ref,
			&alpha, &beta, &b1, &b2, &rs, &kappa))
		return NULL;
	if (!get_translate_args(&A, esp, t1, g1, t2, g2))
		return NULL;
	out=pv_craze_translate_IEC60891_v3(A.esp, VEC(A.g1), VEC(A.t1), A.n, VEC(A.g2), VEC(A.t2), A.n,
		voc_ref, t_ref, g_ref, alpha, beta, b1, b2, rs, kappa);
	FLUSH_LIBRARY_OUTPUT();
	if (out==NULL)
	{
		free_translate_args(&A);
		return PyErr_NoMemory();
	}
	ret=(PyObject *)esp_to_array(out, A.n);
	free(out);
	free_translate_args(&A);
	return ret;
}

/* parse esp + T, G vectors of equal length for the fitting routines */
typedef struct {
	curvepars *esp;
	PyArrayObject *t, *g;
	int n;
} fit_args;

static void free_fit_args(fit_args *a)
{
	free(a->esp);
	Py_XDECREF(a->t);
	Py_XDECREF(a->g);
}

static int get_fit_args(fit_args *a, PyObject *esp, PyObject *t, PyObject *g)
{
	memset(a, 0, sizeof(*a));
	if ((a->esp=esp_from_obj(esp, &a->n))==NULL)
		return 0;
	if (((a->t=as_vector(t, a->n, "T"))==NULL)||((a->g=as_vector(g, a->n, "G"))==NULL))
	{
		free_fit_args(a);
		return 0;
	}
	if (a->n==0)
	{
		PyErr_SetString(PyExc_ValueError, "no data to fit");
		free_fit_args(a);
		return 0;
	}
	return 1;
}

static PyObject *core_iec60891_v2_fit(PyObject *self, PyObject *args)
{
	PyObject *esp, *t, *g;
	double t_ref, g_ref, alpha=0, beta=0, a=0, rs=0, kappa=0;
	int max_points, seed;
	fit_args A;
	curvepars out;
	double row[10];
	npy_intp dims[1]={10};
	PyArrayObject *o;
	if (!PyArg_ParseTuple(args, "OOOddii", &esp, &t, &g, &t_ref, &g_ref, &max_points, &seed))
		return NULL;
	if (!get_fit_args(&A, esp, t, g))
		return NULL;
	out=pv_craze_fit_IEC60891_v2(A.esp, VEC(A.g), VEC(A.t), A.n, t_ref, g_ref, &alpha, &beta, &a, &rs, &kappa, max_points, seed);
	FLUSH_LIBRARY_OUTPUT();
	free_fit_args(&A);
	if ((o=(PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE))==NULL)
		return NULL;
	pars_to_row(&out, row);
	memcpy(PyArray_DATA(o), row, sizeof(row));
	return Py_BuildValue("(Nddddd)", o, alpha, beta, a, rs, kappa);
}

static PyObject *core_iec60891_v3_fit(PyObject *self, PyObject *args)
{
	PyObject *esp, *t, *g;
	double t_ref, g_ref, alpha=0, beta=0, b1=0, b2=0, rs=0, kappa=0;
	int max_points, seed;
	fit_args A;
	curvepars out;
	double row[10];
	npy_intp dims[1]={10};
	PyArrayObject *o;
	if (!PyArg_ParseTuple(args, "OOOddii", &esp, &t, &g, &t_ref, &g_ref, &max_points, &seed))
		return NULL;
	if (!get_fit_args(&A, esp, t, g))
		return NULL;
	out=pv_craze_fit_IEC60891_v3(A.esp, VEC(A.g), VEC(A.t), A.n, t_ref, g_ref, &alpha, &beta, &b1, &b2, &rs, &kappa, max_points, seed);
	FLUSH_LIBRARY_OUTPUT();
	free_fit_args(&A);
	if ((o=(PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE))==NULL)
		return NULL;
	pars_to_row(&out, row);
	memcpy(PyArray_DATA(o), row, sizeof(row));
	return Py_BuildValue("(Ndddddd)", o, alpha, beta, b1, b2, rs, kappa);
}

static PyObject *core_diode_gt_fit(PyObject *self, PyObject *args)
{
	PyObject *esp, *t, *g;
	double t_ref, g_ref, alpha=0, ea=0;
	int max_points, seed;
	fit_args A;
	generic_pars out;
	npy_intp dims[1]={10};
	PyArrayObject *o;
	if (!PyArg_ParseTuple(args, "OOOddii", &esp, &t, &g, &t_ref, &g_ref, &max_points, &seed))
		return NULL;
	if (!get_fit_args(&A, esp, t, g))
		return NULL;
	out=pv_craze_fit_diode_GT(A.esp, VEC(A.g), VEC(A.t), A.n, &alpha, &ea, t_ref, g_ref, 'C', max_points, seed);
	FLUSH_LIBRARY_OUTPUT();
	free_fit_args(&A);
	if ((o=(PyArrayObject *)PyArray_ZEROS(1, dims, NPY_DOUBLE, 0))==NULL)
		return NULL;
	memcpy(PyArray_DATA(o), out.pars, pv_craze_numfields(M_DIODE)*sizeof(double));
	return Py_BuildValue("(Ndd)", o, alpha, ea);
}

static PyObject *core_diode_gt_translate(PyObject *self, PyObject *args)
{
	PyObject *pars_obj, *t_obj, *g_obj, *ret;
	PyArrayObject *t, *g;
	double t_ref, g_ref, alpha, ea;
	generic_pars P, *out;
	npy_intp n;
	int ok;
	if (!PyArg_ParseTuple(args, "OddOOdd", &pars_obj, &t_ref, &g_ref, &t_obj, &g_obj, &alpha, &ea))
		return NULL;
	P=single_generic((int)M_DIODE, pars_obj, &ok);
	if (!ok)
		return NULL;
	if ((t=as_double(t_obj, 1, -1, "T"))==NULL)
		return NULL;
	n=PyArray_DIM(t, 0);
	if ((g=as_vector(g_obj, n, "G"))==NULL)
	{
		Py_DECREF(t);
		return NULL;
	}
	if (n==0)
	{
		Py_DECREF(t);
		Py_DECREF(g);
		PyErr_SetString(PyExc_ValueError, "nothing to translate");
		return NULL;
	}
	out=pv_craze_translate_diode_GT(P, g_ref, t_ref, VEC(g), VEC(t), (int)n, alpha, ea, 'C');
	FLUSH_LIBRARY_OUTPUT();
	Py_DECREF(t);
	Py_DECREF(g);
	if (out==NULL)
		return PyErr_NoMemory();
	ret=generic_to_tuple(out, (int)n);
	free(out);
	return ret;
}

/* **************************************************
 * module definition
 * **************************************************/

static PyMethodDef core_methods[]={
	{"model_table", core_model_table, METH_NOARGS, "model_table() -> list of (id, name, labels)"},
	{"flag_table", core_flag_table, METH_NOARGS, "flag_table() -> list of (code, name, message, is_error)"},
	{"default_analysis_settings", core_default_analysis_settings, METH_NOARGS, "default_analysis_settings() -> dict"},
	{"get_verbosity", core_get_verbosity, METH_NOARGS, "get_verbosity() -> int bit mask"},
	{"set_verbosity", core_set_verbosity, METH_VARARGS, "set_verbosity(mask)"},
	{"analyze", (PyCFunction)(void (*)(void))core_analyze, METH_VARARGS|METH_KEYWORDS,
		"analyze(curves, models, *, spike_filter, identical_eps, resample, mpp_method, mpp_dv, mpp_di, sc_dv, oc_di, min_rel_range)"
		" -> (esp, flags, (types, pars, err))"},
	{"fit", core_fit, METH_VARARGS, "fit(esp, models) -> (types, pars, err)"},
	{"curves", core_curves, METH_VARARGS, "curves(types, pars, npoints) -> list of (n, 2) arrays"},
	{"model_esp", core_model_esp, METH_VARARGS, "model_esp(types, pars) -> esp"},
	{"rmse", core_rmse, METH_VARARGS, "rmse(types, pars, curves) -> (N,) array"},
	{"distance_curves", core_distance_curves, METH_VARARGS, "distance_curves(types, pars, curves) -> list of (n, 2) arrays"},
	{"iv", core_iv, METH_VARARGS, "iv(types, pars, V) -> I, V of shape (N, m)"},
	{"vi", core_vi, METH_VARARGS, "vi(types, pars, I) -> V, I of shape (N, m)"},
	{"iec60891_v2_translate", core_iec60891_v2_translate, METH_VARARGS,
		"iec60891_v2_translate(esp, T_from, G_from, T, G, alpha, beta, a, rs, kappa) -> esp"},
	{"iec60891_v3_translate", core_iec60891_v3_translate, METH_VARARGS,
		"iec60891_v3_translate(esp, T_from, G_from, T, G, voc_ref, t_ref, g_ref, alpha, beta, b1, b2, rs, kappa) -> esp"},
	{"iec60891_v2_fit", core_iec60891_v2_fit, METH_VARARGS,
		"iec60891_v2_fit(esp, T, G, t_ref, g_ref, max_points, seed) -> (esp_ref, alpha, beta, a, rs, kappa)"},
	{"iec60891_v3_fit", core_iec60891_v3_fit, METH_VARARGS,
		"iec60891_v3_fit(esp, T, G, t_ref, g_ref, max_points, seed) -> (esp_ref, alpha, beta, b1, b2, rs, kappa)"},
	{"diode_gt_fit", core_diode_gt_fit, METH_VARARGS,
		"diode_gt_fit(esp, T, G, t_ref, g_ref, max_points, seed) -> (pars, alpha, ea)"},
	{"diode_gt_translate", core_diode_gt_translate, METH_VARARGS,
		"diode_gt_translate(pars, t_ref, g_ref, T, G, alpha, ea) -> (types, pars)"},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef core_module={
	PyModuleDef_HEAD_INIT,
	"pvcraze._core",
	"Private low level bindings to the PV-CRAZE C library. Use the pvcraze package instead.",
	-1,
	core_methods
};

PyMODINIT_FUNC PyInit__core(void)
{
	PyObject *m;
	import_array();
	if ((m=PyModule_Create(&core_module))==NULL)
		return NULL;
	if (PyModule_AddIntConstant(m, "MODEL_NONE", (long)M_NONE)<0)
	{
		Py_DECREF(m);
		return NULL;
	}
	return m;
}
