/*
 * Python interface to PV-CRAZE
 * Copyright (C) 2024  B. E. Pieters, 
 * IEK-5 Photovoltaik, Forschunszentrum Juelich
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
#include <pvcraze.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/* **************************************************
 * utilities 
 * **************************************************/
static PyObject *py_pv_craze_LookupModel(PyObject *self, PyObject *args);
static PyObject *py_pv_craze_ReverseLookup(PyObject *self, PyObject *args);
static PyObject *py_pv_craze_FlagTable(PyObject *self, PyObject *noargs);
static PyObject *py_pv_craze_quiet(PyObject *self, PyObject *noargs);
static PyObject *py_pv_craze_verbose(PyObject *self, PyObject *noargs);
static PyObject *py_pv_craze_verbose_set_report(PyObject *self, PyObject *args);
static PyObject *py_pv_craze_verbose_set_progress(PyObject *self, PyObject *args);
static PyObject *py_pv_craze_verbose_set_timing(PyObject *self, PyObject *args);
static PyObject *py_pv_craze_verbose_set_head(PyObject *self, PyObject *args);

/* **************************************************
 * curve modeling and analysis routines 
 * **************************************************/
static PyObject *py_pv_craze_analyze(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_modelfit(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_curvemodel(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_modelpars(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_error(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_iv(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_vi(PyObject *self, PyObject *args, PyObject *kwargs);
//static PyObject *py_pv_craze_errorcurve(PyObject *self, PyObject *args);

/* **************************************************
 * IEC60891 (version 2 and 3) procedure 2 routines
 * **************************************************/
// IEC60891 version 2 (2009)
static PyObject *py_pv_craze_translate_IEC60891_v2(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_fit_IEC60891_v2(PyObject *self, PyObject *args, PyObject *kwargs);
// IEC60891 version 3 (2021)
static PyObject *py_pv_craze_translate_IEC60891_v3(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_fit_IEC60891_v3(PyObject *self, PyObject *args, PyObject *kwargs);

/* **************************************************
 * Diode Temperature and Irradiance routines
 * **************************************************/
static PyObject *py_pv_craze_translate_diode_GT(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *py_pv_craze_fit_diode_GT(PyObject *self, PyObject *args, PyObject *kwargs);


/* ==== Set up the methods table ====================== */
static PyMethodDef pv_craze_methods[] = {
	{"quiet",  					py_pv_craze_quiet								, METH_NOARGS					,"Silence the pv-craze library"},
	{"verbose",  				py_pv_craze_verbose								, METH_NOARGS					,"Make the pv-craze library verbose"},
	{"verbose_set_report",		py_pv_craze_verbose_set_report					, METH_VARARGS 					,"Enable/Disable printing analysis report"},
	{"verbose_set_progress",	py_pv_craze_verbose_set_progress				, METH_VARARGS 					,"Enable/Disable printing a progress bar"},
	{"verbose_set_timing",		py_pv_craze_verbose_set_timing					, METH_VARARGS 					,"Enable/Disable printing timing information"},
	{"verbose_set_head",		py_pv_craze_verbose_set_head					, METH_VARARGS 					,"Enable/Disable printing a header"},
	{"analyze", 				(PyCFunction)py_pv_craze_analyze				, METH_VARARGS | METH_KEYWORDS	,"PV-CRAZE curve analysis. Extracts the extended solar cell parameters and optionally fits curve models"},
	{"curvemodel", 				(PyCFunction)py_pv_craze_curvemodel				, METH_VARARGS | METH_KEYWORDS	,"Reconstruct IV from curve model parameters"},
	{"modelpars", 				(PyCFunction)py_pv_craze_modelpars				, METH_VARARGS | METH_KEYWORDS	,"Compute extended solar cell parameters from generic parameters"},
	{"lookupmodel", 			py_pv_craze_LookupModel							, METH_VARARGS 					,"Returns a model ID based on a model name"},
	{"reverselookup", 			py_pv_craze_ReverseLookup						, METH_VARARGS 					,"Returns a model name based on a model ID"},
	{"flagtable",  				py_pv_craze_FlagTable							, METH_NOARGS					,"generate a flag table for the used pv-craze library"},
	{"modelfit", 				(PyCFunction)py_pv_craze_modelfit			    , METH_VARARGS | METH_KEYWORDS	,"Takes extended solar cell parameters and fits curve models to them"},
	{"modelerror", 				(PyCFunction)py_pv_craze_error					, METH_VARARGS | METH_KEYWORDS	,"Compute RMS distance between IV-curve and model"},
	{"iv", 						(PyCFunction)py_pv_craze_iv						, METH_VARARGS | METH_KEYWORDS	,"Compute current from voltage with model parameters"},
	{"vi", 						(PyCFunction)py_pv_craze_vi						, METH_VARARGS | METH_KEYWORDS	,"Compute voltage from current with model parameters"},
//    {"pv_craze_errorcurve", py_pv_craze_errorcurve, METH_VARARGS,""},
	{"translate_IEC60891_v2", 	(PyCFunction)py_pv_craze_translate_IEC60891_v2	, METH_VARARGS | METH_KEYWORDS	,"IEC60891 v2 translation routine"},
	{"fit_IEC60891_v2", 		(PyCFunction)py_pv_craze_fit_IEC60891_v2			, METH_VARARGS | METH_KEYWORDS	,"IEC60891 v2 fitting routine"},
	{"translate_IEC60891_v3", 	(PyCFunction)py_pv_craze_translate_IEC60891_v3	, METH_VARARGS | METH_KEYWORDS	,"IEC60891 v3 translation routine"},
	{"fit_IEC60891_v3", 		(PyCFunction)py_pv_craze_fit_IEC60891_v3			, METH_VARARGS | METH_KEYWORDS	,"IEC60891 v3 fitting routine"},
	{"translate_diode_GT",		(PyCFunction)py_pv_craze_translate_diode_GT		, METH_VARARGS | METH_KEYWORDS	,"Diode model irradiance/temperature translation routine"},
	{"fit_diode_GT",	 		(PyCFunction)py_pv_craze_fit_diode_GT			, METH_VARARGS | METH_KEYWORDS	,"Diode model irradiance/temperature fitting routine"},
	{ NULL, NULL, 0, NULL }     /* Sentinel - marks the end of this structure */
};

static struct PyModuleDef python_pvcraze = {
	PyModuleDef_HEAD_INIT,
	"pvcraze",
	"Photo-Voltaic CuRve AnalyZEr (PV-CRAZE)",
	-1,
	pv_craze_methods
};

int IsDoubleArray(PyArrayObject *Array, int n, int m)
{
	int typ;
	int SIZE=0;
	
	if (NULL==Array)
	{
		PyErr_SetString(PyExc_ValueError, "Error: expect a numpy array with double type data");
		return 0;	
	}
	
	if (!PyObject_IsInstance((PyObject *) Array, (PyObject *)&PyArray_Type))
	{
		PyErr_SetString(PyExc_ValueError, "Error: expect a numpy array with double type data");
		return 0;
	}
	if (n>=0)
		SIZE|=1;
	if (m>=0)
		SIZE|=2;
	typ=PyArray_TYPE(Array);
	switch(typ)
	{
		case NPY_DOUBLE:
			switch(SIZE)
			{
				case 0:
					return 1;
				case 1:
					if (Array->dimensions[0]==n)
						return 1;
					else
						PyErr_Format(PyExc_ValueError, "Error: expect a numpy array with dimension 0 length %d",n);
					
					break;
				case 2:
					if (Array->dimensions[1]==m)
						return 1;
					else
						PyErr_Format(PyExc_ValueError, "Error: expect a numpy array with dimension 1 length %d",m);
					break;
				case 3:
					if ((Array->dimensions[0]==n) &&(Array->dimensions[1]==m))
						return 1;
					else
						PyErr_Format(PyExc_ValueError, "Error: expect a numpy array with size %d x %d",n,m);
					break;
				default:
					break;
			}
			return 0;
		default:
			PyErr_SetString(PyExc_ValueError, "Error: expect a numpy array with double type data");
			return 0;
	}
}

void double_into_2darray(PyArrayObject *np_array, npy_intp i, npy_intp j, double value)
{
    if (!PyArray_Check(np_array))
    {
        PyErr_SetString(PyExc_TypeError, "Input must be a numpy.ndarray.");
        return;
    }
    if (PyArray_TYPE(np_array) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "Array must be of type np.float64 (double).");
        return;
    }
    if (PyArray_NDIM(np_array) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "Array must be 2-dimensional.");
        return;
    }
    // Ensure we have a writable array and not a read-only view or such.
    if (!PyArray_ISWRITEABLE(np_array)) {
        PyErr_SetString(PyExc_RuntimeError, "Array is not writable.");
        return;
    }

    // Get a pointer to the element at position (i, j)
    double *elem = (double *)PyArray_GETPTR2(np_array, i, j);

    // Insert the value
    *elem = value;
}

void int_into_2darray(PyArrayObject *np_array, npy_intp i, npy_intp j, int value)
{
    if (!PyArray_Check(np_array))
    {
        PyErr_SetString(PyExc_TypeError, "Input must be a numpy.ndarray.");
        return;
    }
    if (PyArray_TYPE(np_array) != NPY_INT)
    {
        PyErr_SetString(PyExc_TypeError, "Array must be of type np.int32 (int).");
        return;
    }
    if (PyArray_NDIM(np_array) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "Array must be 2-dimensional.");
        return;
    }
    // Ensure we have a writable array and not a read-only view or such.
    if (!PyArray_ISWRITEABLE(np_array)) {
        PyErr_SetString(PyExc_RuntimeError, "Array is not writable.");
        return;
    }

    // Get a pointer to the element at position (i, j)
    int *elem = (int *)PyArray_GETPTR2(np_array, i, j);

    // Insert the value
    *elem = value;
}

void double_into_1darray(PyArrayObject *np_array, npy_intp i, double value)
{
    if (!PyArray_Check(np_array))
    {
        PyErr_SetString(PyExc_TypeError, "Input must be a numpy.ndarray.");
        return;
    }
    if (PyArray_TYPE(np_array) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "Array must be of type np.float64 (double).");
        return;
    }
    if (PyArray_NDIM(np_array) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Array must be 2-dimensional.");
        return;
    }
    // Ensure we have a writable array and not a read-only view or such.
    if (!PyArray_ISWRITEABLE(np_array)) {
        PyErr_SetString(PyExc_RuntimeError, "Array is not writable.");
        return;
    }

    // Get a pointer to the element at position (i, j)
    double *elem = (double *)PyArray_GETPTR1(np_array, i);

    // Insert the value
    *elem = value;
}



PyMODINIT_FUNC PyInit_pvcraze(void)
{
	PyObject *module = PyModule_Create(&python_pvcraze);

	if (module == NULL)
		return NULL;

	PyModule_AddStringConstant(module, "__author__", "Bart E. Pieters and Jan Körtgen");
	PyModule_AddStringConstant(module, "__version__", "0.1.0");
	import_array();  // Must be present for NumPy.  Called first after above line.
	return module;
}


static PyObject *py_pv_craze_quiet(PyObject *self, PyObject *noargs)
{
	pv_craze_quiet();
	return Py_None;
}

static PyObject *py_pv_craze_verbose(PyObject *self, PyObject *noargs)
{
	pv_craze_verbose();
	return Py_None;
}

static PyObject *py_pv_craze_verbose_set_report(PyObject *self, PyObject *args)
{
	int b=-1; 
	
	if(!PyArg_ParseTuple(args, "i", &b))
		return NULL;
	if(b<0)
		return NULL;
	pv_craze_verbosereport(b);
	return Py_None;
}
static PyObject *py_pv_craze_verbose_set_progress(PyObject *self, PyObject *args)
{
	int b=-1; 
	
	if(!PyArg_ParseTuple(args, "i", &b))
		return NULL;
	if(b<0)
		return NULL;
	pv_craze_verboseprogress(b);
	return Py_None;
}
static PyObject *py_pv_craze_verbose_set_timing(PyObject *self, PyObject *args)
{
	int b=-1; 
	
	if(!PyArg_ParseTuple(args, "i", &b))
		return NULL;
	if(b<0)
		return NULL;
	pv_craze_verbosetiming(b);
	return Py_None;
}
static PyObject *py_pv_craze_verbose_set_head(PyObject *self, PyObject *args)
{
	int b=-1; 
	
	if(!PyArg_ParseTuple(args, "i", &b))
		return NULL;
	if(b<0)
		return NULL;
	pv_craze_verbosehead(b);
	return Py_None;
}


static PyObject *py_pv_craze_LookupModel(PyObject *self, PyObject *args)
{
	char * str; 
	MODEL_TYPE T;
	if(!PyArg_ParseTuple(args, "s", &str))
		return NULL;
	if(NULL==str)
		return NULL;
	
	T=pv_craze_lookupmodel(str);
	if (T==M_NONE)
	{
		PyErr_Format(PyExc_ValueError, "Error: model \"%s\"does not exist",str);
		return NULL;
	}
	return PyLong_FromLong((long)T);
}

static PyObject *py_pv_craze_ReverseLookup(PyObject *self, PyObject *args)
{
	char * str; 
	int a=-1;
	if(!PyArg_ParseTuple(args, "i", &a))
		return NULL;
	if(a<0)
	{
		PyErr_SetString(PyExc_ValueError, "Error: model does not exist or no valid argument");
		return NULL;
	}
	str=pv_craze_reverselookup((MODEL_TYPE)a);
	if (str==NULL)
	{
		PyErr_Format(PyExc_ValueError, "Error: model with ID %d does not exist",a);
		return NULL;
	}
	return Py_BuildValue("z", str);
}

static PyObject *py_pv_craze_FlagTable(PyObject *self, PyObject *noargs)
{
	AnalyzeFlag *T;
	int i,N;
	PyObject *messages;
	PyObject *ret;
	PyArrayObject *flags;
	PyArrayObject *err;
	PyObject *X;
	npy_intp ArrayDim[1]={0};
	
	T=pv_craze_flagtable(&N);
	if (!T)
	{
		PyErr_SetString(PyExc_ValueError, "Error: could not generate a flag table, pv-craze appears to be broken!");
		return NULL;
	}
	ArrayDim[0]=N;
	flags=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_INT);
	err=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_BOOL);
	messages=PyList_New(N);
	
	for (i=0;i<N;i++)
	{
		X=PyLong_FromLong((long)T[i].FLAG);
		PyArray_SETITEM(flags, PyArray_GETPTR1(flags, i), X);
		Py_DECREF(X);		
		
		X=PyBool_FromLong((long)T[i].ERR);
		PyArray_SETITEM(err, PyArray_GETPTR1(err, i), X);
		Py_DECREF(X);	
		
		X=PyUnicode_FromString(T[i].message);
		PyList_SetItem(messages, i, X);
	}
	pv_craze_free(T);
	ret=Py_BuildValue("{s:O,s:O,s:O}", "flags", flags,"messages", messages,"err",err);
	Py_DECREF(flags);
	Py_DECREF(messages);
	Py_DECREF(err);	
	return ret;
}

npy_intp PARLEN=11;
/* curve analysis tool:
 * py_pv_craze_analyze takes as input a list of numpy arrays and some optional arguments
 * It returns a dictionary with in it:
 * epars: an Nx10 numpy array with for all N input IV's the extended solar cell parameters
 * flags: flags for each curve
 * genpars: list of N lists with one numpy array per model with parameters
 * E: Fit error values in a numpy array with N rows and one column per fitted model 
 */
static PyObject *py_pv_craze_analyze(PyObject *self, PyObject *args, PyObject *kwargs)
{
	// output objects
	PyArrayObject *flags, *E, *epars; // output arrays
	PyArrayObject *genpars; // generic parameter values for all IV's
	PyObject *genlist;		// list of generic parameter arrays, one for every model
	PyObject *ret;
	
	// input objects
	PyObject *py_curveList;
	PyObject * StringList=NULL; 
	
	// internal
	PyArrayObject *py_curve;
	PyObject * str;  
	ProcessSettings S;
	CurveProp *P;
	curve *C;
	npy_intp ArrayDim[2]={0};
	double eps=-1, mppdv=-1,mppdi=-1,scdv=-1,ocdi=-1;
	int n, Ns=0, i, j, k, spike=0, N, am=-1;
	int resample=-1;

	static char *kwlist[] = {
		"IV",
		"models",
		"mpp-method",
		"I_filter",
		"S_filter",
		"mpp_dv", 
		"mpp_di", 
		"sc_dv", 
		"oc_di", 
		"resample", 
		NULL
	};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|O!Cdiddddi",kwlist,&PyList_Type, &py_curveList, &PyList_Type, &StringList, &am, &eps, &spike, &mppdv, &mppdi, &scdv, &ocdi, &resample)) 
		return NULL;
	if(NULL==py_curveList) 
		return NULL;
		
	S=pv_craze_initsettings();
	if (eps>0)
	{
		S.analyze.Identicals=1;
		S.analyze.IdentEps=eps;
	}
	if (spike)
		S.analyze.SpikeFilter=1;	
	if (mppdv>0)
		S.analyze.MPP_f_V=mppdv;
	if (mppdi>0)
		S.analyze.MPP_f_I=mppdi;
	if (scdv>0)
		S.analyze.SC_f_V=scdv;
	if (ocdi>0)
		S.analyze.OC_f_I=ocdi;
	if (am>0)
		S.analyze.am=(char)am;
	if (resample>0)
		S.analyze.Resample=resample;		
	// parse optional models
	if (StringList)
	{	
		Ns=PyList_Size(StringList);
		for(npy_intp i=0; i<Ns; i++)
		{
			str = PyList_GetItem(StringList, i);
			pv_craze_addmodeltosettings(&S, PyUnicode_AsUTF8(str));
		}
	}
	
	// prepare the curves
	N=PyList_Size(py_curveList);
	if ((C=malloc(N*sizeof(curve)))==NULL)
	{
		pv_craze_freesettings(&S);
		PyErr_SetString(PyExc_ValueError, "Error: pv_craze_analyze could not allocate memory");
		return NULL;
	}
	for(j=0; j<N; j++)
	{
		py_curve=(PyArrayObject *) PyList_GetItem(py_curveList, (Py_ssize_t)j);
		if (!IsDoubleArray(py_curve,-1,2))
		{
			pv_craze_freesettings(&S);
			return NULL;
		}
		n=py_curve->dimensions[0];
		if ((C[j].P=malloc(n*sizeof(point)))==NULL)
		{
			for(i=j-1; i>=0; i--)
				free(C[i].P);	
			pv_craze_freesettings(&S);
			PyErr_SetString(PyExc_ValueError, "Error: pv_craze_analyze could not allocate memory");
			return NULL;
		}
		for(i=0; i<n; i++)
		{
			C[j].P[i].x=*(double *)(py_curve->data+i*py_curve->strides[0]);
			C[j].P[i].y=*(double *)(py_curve->data+i*py_curve->strides[0]+py_curve->strides[1]);
		}
		C[j].N=n;
	}

	P=pv_craze_analyze(C, N, &S);

	ArrayDim[0]=N;
	ArrayDim[1]=10;
	epars=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
	
	
	ArrayDim[1]=S.L.N; // number of models
	E=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
	
	ArrayDim[1]=1; // number of flags
	flags=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_INT);
	
	genlist=PyList_New(S.L.N);
	// allocate space for arrays
	ArrayDim[0]=N;
	ArrayDim[1]=11;
	for (i=0; i<S.L.N; i++)
	{
		genpars=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
		PyList_SetItem(genlist, i, (PyObject *) genpars);
	}
	
	for (i=0; i<N; i++)
	{
		double_into_2darray(epars, i, 0, P[i].par.Voc);
		double_into_2darray(epars, i, 1, P[i].par.Isc);
		double_into_2darray(epars, i, 2, P[i].par.Vmpp);
		double_into_2darray(epars, i, 3, P[i].par.Impp);
		double_into_2darray(epars, i, 4, P[i].par.Vqmp1);
		double_into_2darray(epars, i, 5, P[i].par.Iqmp1);
		double_into_2darray(epars, i, 6, P[i].par.Vqmp2);
		double_into_2darray(epars, i, 7, P[i].par.Iqmp2);
		double_into_2darray(epars, i, 8, P[i].par.Gsc);
		double_into_2darray(epars, i, 9, P[i].par.Roc);
		
		
		int_into_2darray(flags, i, 0, P[i].A_FLAGS);
		for (k=0;k<P[i].N;k++)
		{
			// fetch next model parameters
			genpars=(PyArrayObject *) PyList_GetItem(genlist, (Py_ssize_t)k);	
			double_into_2darray(genpars, i, 0, (double)P[i].P[k].TYPE);
			for(j=0; j<10; j++)
				double_into_2darray(genpars, i, j+1, P[i].P[k].pars[j]);
			double_into_2darray(E, i, k,P[i].E[k]);
		}  
	}
	for(i=0; i<N; i++)
	{
		free(C[i].P);
		pv_craze_freeprop(&P[i]);
	}
	Ns=S.L.N;
	pv_craze_freesettings(&S);
	pv_craze_free(P);
	free(C);
	if (Ns>0)
	{
		ret=Py_BuildValue("{s:O,s:O,s:O,s:O}", "esp", epars,"flags", flags, "modelpar", genlist, "err", E);
		Py_DECREF(genlist);
		Py_DECREF(E);
	}
	else
		ret=Py_BuildValue("{s:O,s:O}", "esp", epars,"flags", flags);
	Py_DECREF(epars);
	Py_DECREF(flags);
	return ret;
	
}

/* epar fitting tool:
 * py_pv_craze_modelfit takes as input a list of epar numpy arrays (nx10)
 * and a list of strings to specify the models to fit.
 * It returns a list with generic parameter numpy arrays (nx11)
 */
static PyObject *py_pv_craze_modelfit(PyObject *self, PyObject *args, PyObject *kwargs)
{
	// output objects
	PyArrayObject *genpars; // generic parameter values for all IV's
	PyObject *genlist;		// list of generic parameter arrays, one for every model
	
	PyArrayObject *fiterr; // model fit error
	PyObject *fiterrlist;  // list of fit errors
	
	// input objects
	PyArrayObject *epars; // epars
	PyObject * StringList=NULL; 
	
	// internal
	PyObject * str;  
	PyObject * ret; 
	ProcessSettings S;
	CurveProp *P;
	curvepars *pars;
	npy_intp ArrayDim[2]={0};
	int n, Ns=0, i, j, k;

	static char *kwlist[] = {
		"esp",
		"models",
		NULL
	};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!",kwlist,&PyArray_Type, &epars, &PyList_Type, &StringList)) 
		return NULL;
	if(NULL==epars) 
		return NULL;	

		
	if (StringList)
	{	
		S=pv_craze_initsettings();
		Ns=PyList_Size(StringList);
		for(npy_intp i=0; i<Ns; i++)
		{
			str = PyList_GetItem(StringList, i);
			pv_craze_addmodeltosettings(&S, PyUnicode_AsUTF8(str));
		}
	}
	else
		return NULL;	
		
	if (S.L.N==0)
	{
		pv_craze_freesettings(&S);
		PyErr_SetString(PyExc_ValueError, "Error: no valid model strings found");
		return NULL;
	}
		
	// prepare the curves
	n=epars->dimensions[0];
	if (!IsDoubleArray(epars,-1,10))
	{
		pv_craze_freesettings(&S);
		return NULL;
	}
	if ((pars=malloc(n*sizeof(curvepars)))==NULL)
	{
		pv_craze_freesettings(&S);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	
	for(j=0; j<n; j++)
	{
		pars[j].Voc=*(double *)(epars->data+j*epars->strides[0]+0*epars->strides[1]);
		pars[j].Isc=*(double *)(epars->data+j*epars->strides[0]+1*epars->strides[1]);
		pars[j].Vmpp=*(double *)(epars->data+j*epars->strides[0]+2*epars->strides[1]);
		pars[j].Impp=*(double *)(epars->data+j*epars->strides[0]+3*epars->strides[1]);
		pars[j].Vqmp1=*(double *)(epars->data+j*epars->strides[0]+4*epars->strides[1]);
		pars[j].Iqmp1=*(double *)(epars->data+j*epars->strides[0]+5*epars->strides[1]);
		pars[j].Vqmp2=*(double *)(epars->data+j*epars->strides[0]+6*epars->strides[1]);
		pars[j].Iqmp2=*(double *)(epars->data+j*epars->strides[0]+7*epars->strides[1]);
		pars[j].Gsc=*(double *)(epars->data+j*epars->strides[0]+8*epars->strides[1]);
		pars[j].Roc=*(double *)(epars->data+j*epars->strides[0]+9*epars->strides[1]);
	}

	P=pv_craze_modelfit(pars, n, &S);
	
	genlist=PyList_New(S.L.N);
	fiterrlist=PyList_New(S.L.N);
	// allocate space for arrays
	ArrayDim[0]=n;
	ArrayDim[1]=11;
	for (i=0; i<S.L.N; i++)
	{
		genpars=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
		fiterr=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_DOUBLE);
			
		for (k=0; k<n; k++)
		{	
			// model parameters
			double_into_2darray(genpars, k, 0, (double)P[k].P[i].TYPE);
			
			for(j=0; j<10; j++)
				double_into_2darray(genpars, k, j+1, P[k].P[i].pars[j]);
				
			// fit error
			double_into_1darray(fiterr, k, P[k].E[i]);
		}  
		
		PyList_SetItem(genlist, i, (PyObject *)genpars);
				
		PyList_SetItem(fiterrlist, i, (PyObject *)fiterr);
	}
	

	for(i=0; i<n; i++)
		pv_craze_freeprop(&P[i]);
	pv_craze_freesettings(&S);
	pv_craze_free(P);
	free(pars);

	ret=Py_BuildValue("{s:O,s:O}", "modelpar", genlist,"err", fiterrlist);
	Py_DECREF(genlist);
	Py_DECREF(fiterrlist);
	return ret;
}


MODEL_TYPE PVC_def_type=M_PVBEZIER;
#define MAX_NIV 1000000
/* curve model:
 * py_pv_craze_curvemodel takes as input a list of list of generic parameters numpy arrays (nx11), 
 * the number of points per IV, and some optional arguments.
 * It returns a list per generic parameter input array, with a list of IV data (one array per row in generic parameter array) 
 */
static PyObject *py_pv_craze_curvemodel(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *py_genparList, *py_curveList_out, *py_curveListList_out;
	PyArrayObject *py_genpars, *py_curve_out;
	curve *curve_out;
	generic_pars *genpars_in;
	int n,m, Na=-1, M, Niv, i, j, k, l;
	int *len;
	npy_intp ArrayDim[2]={0};

	static char *kwlist[] = {
		"modelpar",
		"N",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!i",kwlist,&PyList_Type, &py_genparList,&Niv)) 
		return NULL;
	if(NULL==py_genparList) 
		return NULL;
	if ((Niv<=0) ||(Niv>MAX_NIV))
	{
		PyErr_SetString(PyExc_ValueError, "Error: invalid number of IV points");
		return NULL;	
	}


	M=PyList_Size(py_genparList); // number of models

	if ((genpars_in=malloc(M*sizeof(generic_pars)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((len=malloc(M*sizeof(int)))==NULL)
	{
		free(genpars_in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	Na=M;
	k=0; // counter for the number of genpars
	for(j=0; j<M; j++)
	{
		py_genpars=(PyArrayObject *) PyList_GetItem(py_genparList, (Py_ssize_t)j);
		if (!IsDoubleArray(py_genpars,-1,-1))
		{
			free(genpars_in);
			free(len);
			return NULL;
		}
		n=py_genpars->dimensions[0];
		m=py_genpars->dimensions[1];
		len[j]=n;
		if ((m!=11)&&(m!=10))
		{
			free(genpars_in);
			free(len);
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an n x 11 array");
			return NULL;
		}
		if (k+n>Na)
		{
			generic_pars *dummy;			
			if ((dummy=realloc(genpars_in, (k+n)*sizeof(generic_pars)))==NULL)
			{
				free(genpars_in);
				free(len);
				PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
				return NULL;	
			}
			genpars_in=dummy;
		}
		for (i=0;i<n;i++)
		{
			double a;
			if (m==11)
			{
				a=*(double *)(py_genpars->data+i*py_genpars->strides[0]);
				if (fabs(a-round(a))>1e-15)
				{
					free(genpars_in);
					free(len);
					PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
					return NULL;	
				}
				if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
				{
					free(genpars_in);
					free(len);
					PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
					return NULL;	
				}
				genpars_in[k].TYPE=(MODEL_TYPE)a;
				for(l=0; l<10; l++)
					genpars_in[k].pars[l]=*(double *)(py_genpars->data+i*py_genpars->strides[0]+(l+1)*py_genpars->strides[1]);
			}
			else
			{
				genpars_in[k].TYPE=PVC_def_type;
				for(l=0; l<10; l++)
					genpars_in[k].pars[l]=*(double *)(py_genpars->data+i*py_genpars->strides[0]+l*py_genpars->strides[1]);
			}
			k++;
		}       
	}

	curve_out=pv_craze_curvemodel(genpars_in, k, Niv);

	py_curveListList_out=PyList_New(M);
	ArrayDim[0]=Niv;
	ArrayDim[1]=2;
	k=0;
	for(j=0; j<M; j++)
	{
		py_curveList_out=PyList_New(len[j]);
		for(i=0; i<len[j]; i++)
		{
			ArrayDim[0]=curve_out[k].N;
			py_curve_out=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
			for(l=0; l<curve_out[k].N; l++)
			{
				double_into_2darray(py_curve_out, l, 0, curve_out[k].P[l].x);
				double_into_2darray(py_curve_out, l, 1, curve_out[k].P[l].y);
			}
			pv_craze_free(curve_out[k].P); // use pv-craze free as pv-craze allocated this memory
			k++;
			PyList_SetItem(py_curveList_out, i, (PyObject *)py_curve_out); 
		}
		PyList_SetItem(py_curveListList_out, j, py_curveList_out); 
	}
	free(len);
	free(genpars_in);
	pv_craze_free(curve_out);

	return py_curveListList_out;
}

/* model pars:
 * py_pv_craze_modelpars takes as input a list of list of generic parameters numpy arrays (nx11), 
 * It returns a list per generic parameter input array, with an array of esp
 */
static PyObject *py_pv_craze_modelpars(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *py_genparList, *py_solparsList_out;
	PyArrayObject *py_genpars, *py_solpars_out;
	curvepars solpars;
	generic_pars genpars_in;
	
	
	int n,m, M, i, j, l;
	npy_intp ArrayDim[2]={0};

	static char *kwlist[] = {
		"modelpar",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!",kwlist,&PyList_Type, &py_genparList)) 
		return NULL;
	if(NULL==py_genparList) 
		return NULL;
	
	M=PyList_Size(py_genparList); // number of models	
	py_solparsList_out=PyList_New(M);
	
	for(j=0; j<M; j++)
	{
		py_genpars=(PyArrayObject *) PyList_GetItem(py_genparList, (Py_ssize_t)j);
		if (!IsDoubleArray(py_genpars,-1,-1))
		{
			Py_XDECREF(py_solparsList_out);						
			return NULL;
		}
		n=py_genpars->dimensions[0];
		m=py_genpars->dimensions[1];
		
		
		if ((m!=11)&&(m!=10))
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an n x 11 array");
			Py_XDECREF(py_solparsList_out);						
			return NULL;
		}
		ArrayDim[0]=n;
		ArrayDim[1]=10;
		py_solpars_out=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
		for (i=0;i<n;i++)
		{
			double a;
			if (m==11)
			{
				a=*(double *)(py_genpars->data+i*py_genpars->strides[0]);
				if (fabs(a-round(a))>1e-15)
				{
					PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
					Py_XDECREF(py_solparsList_out);						
					Py_XDECREF(py_solpars_out);						
					return NULL;	
				}
				if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
				{
					PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
					Py_XDECREF(py_solparsList_out);						
					Py_XDECREF(py_solpars_out);						
					return NULL;	
				}
				genpars_in.TYPE=(MODEL_TYPE)a;
				for(l=0; l<10; l++)
					genpars_in.pars[l]=*(double *)(py_genpars->data+i*py_genpars->strides[0]+(l+1)*py_genpars->strides[1]);
			}
			else
			{
				genpars_in.TYPE=PVC_def_type;
				for(l=0; l<10; l++)
					genpars_in.pars[l]=*(double *)(py_genpars->data+i*py_genpars->strides[0]+l*py_genpars->strides[1]);
			}		
			solpars=pv_craze_g_pars(genpars_in);
			double_into_2darray(py_solpars_out, i, 0, solpars.Voc);
			double_into_2darray(py_solpars_out, i, 1, solpars.Isc);
			double_into_2darray(py_solpars_out, i, 2, solpars.Vmpp);
			double_into_2darray(py_solpars_out, i, 3, solpars.Impp);
			double_into_2darray(py_solpars_out, i, 4, solpars.Vqmp1);
			double_into_2darray(py_solpars_out, i, 5, solpars.Iqmp1);
			double_into_2darray(py_solpars_out, i, 6, solpars.Vqmp2);
			double_into_2darray(py_solpars_out, i, 7, solpars.Iqmp2);
			double_into_2darray(py_solpars_out, i, 8, solpars.Gsc);
			double_into_2darray(py_solpars_out, i, 9, solpars.Roc);
		}
		PyList_SetItem(py_solparsList_out, j, (PyObject *)py_solpars_out); 
		       
	}
	return py_solparsList_out;
}
/* py_pv_craze_error
 * Takes a list of generic parameter arrays and a list of IV curves
 * For each array in the generic parameter list, the errors w.r.t. each IV is computed
 * This means all arrays in the generic parameter list must be of size nx11, where n is 
 * also the length of the IV list. It outputs a numpy array with one column per generic 
 * parameter array and n rows
 */
static PyObject *py_pv_craze_error(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *py_genparList;    
	PyObject *py_curveList;  
	PyArrayObject *py_genpars, *py_error, *py_curve_in;
	
	generic_pars *genpars_in;
	curve *curve_in;
	double *E;
	int Nm, Niv, i, j, l, n=0, m;
	npy_intp ArrayDim[2]={0};

	static char *kwlist[] = {
		"modelpar",
		"IV",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!",kwlist,&PyList_Type, &py_genparList,&PyList_Type, &py_curveList)) 
		return NULL;
	if(NULL==py_genparList) 
		return NULL;


	Nm=PyList_Size(py_genparList);
	Niv=PyList_Size(py_curveList);
	
	if ((curve_in=malloc(Niv*sizeof(curve)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}

	for(i=0; i<Niv; i++)
	{
		py_curve_in=(PyArrayObject *) PyList_GetItem(py_curveList, (Py_ssize_t)i);
		if (!IsDoubleArray(py_curve_in,-1,2))
		{
			for(i=i-1; i>=0; i--)
				free(curve_in[i].P);
			free(curve_in);
			return NULL;	
		}
		
		n=py_curve_in->dimensions[0];
		m=py_curve_in->dimensions[1];
		if ((m!=2)||(n<1))
		{
			for(i=i-1; i>=0; i--)
				free(curve_in[i].P);
			free(curve_in);
			PyErr_SetString(PyExc_ValueError, "Error: expect a numpy arrays of dimensions Nx2");
			return NULL;	
		}
		
		if ((curve_in[i].P=malloc(n*sizeof(point)))==NULL)
		{
			for(i=i-1; i>=0; i--)
				free(curve_in[i].P);
			free(curve_in);
			PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
			return NULL;	
		}
		for(j=0;j<n;j++)
		{
			curve_in[i].P[j].x=*(double *)(py_curve_in->data+j*py_curve_in->strides[0]);
			curve_in[i].P[j].y=*(double *)(py_curve_in->data+j*py_curve_in->strides[0]+py_curve_in->strides[1]);
		}
		curve_in[i].N=n;
	} 
	
	if ((genpars_in=malloc(Niv*sizeof(generic_pars)))==NULL)
	{
		for(i=0; i<Niv; i++)
			free(curve_in[i].P);
		free(curve_in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	ArrayDim[0]=Niv;
	ArrayDim[1]=Nm;
	py_error=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
	for(i=0; i<Nm; i++)
	{
		py_genpars=(PyArrayObject *) PyList_GetItem(py_genparList, (Py_ssize_t)i);
		if (!IsDoubleArray(py_genpars,-1,-1))
		{
			for(i=0; i<Niv; i++)
				free(curve_in[i].P);
			free(curve_in);
			free(genpars_in);
			return NULL;
		}
		n=py_genpars->dimensions[0];
		m=py_genpars->dimensions[1];
		if ((n!=Niv)||((m!=11)&&(m!=10)))
		{
			for(i=0; i<Niv; i++)
				free(curve_in[i].P);
			free(curve_in);
			free(genpars_in);
			PyErr_SetString(PyExc_ValueError, "Error: generic parameter array must have dimension nx10 or nx11 with n equal to the IV curve list length");
			return NULL;	
		}
		for (j=0;j<n;j++)
		{
			double a;
			if (m==11)
			{
				a=*(double *)(py_genpars->data+j*py_genpars->strides[0]);
				if (fabs(a-round(a))>1e-15)
				{
					for(i=0; i<Niv; i++)
						free(curve_in[i].P);
					free(curve_in);
					free(genpars_in);
					PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
					return NULL;	
				}
				if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
				{
					for(i=0; i<Niv; i++)
						free(curve_in[i].P);
					free(curve_in);
					free(genpars_in);
					PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
					return NULL;	
				}
				genpars_in[j].TYPE=(MODEL_TYPE)a;
				for(l=0; l<10; l++)
					genpars_in[j].pars[l]=*(double *)(py_genpars->data+j*py_genpars->strides[0]+(l+1)*py_genpars->strides[1]);
			}
			else
			{
				genpars_in[j].TYPE=PVC_def_type;
				for(l=0; l<10; l++)
					genpars_in[j].pars[l]=*(double *)(py_genpars->data+j*py_genpars->strides[0]+l*py_genpars->strides[1]);
			}
		}   
		E=pv_craze_error(genpars_in, curve_in, n);
		for (j=0;j<n;j++)
			double_into_2darray(py_error, j, i, E[j]);
		pv_craze_free(E);// pvcraze allocated E
	}
	for (j=0;j<n;j++)
		free(curve_in[j].P);
	
	free(genpars_in);
	free(curve_in);

	return (PyObject *)py_error;

}

static PyObject *py_pv_craze_iv(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *py_genpars, *py_V, *py_I;
	generic_pars P;
	double a,V;
	int i, j, n, m;
	npy_intp ArrayDim[2]={0};

	static char *kwlist[] = {
		"V",
		"modelpar",
		NULL
	};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!",kwlist, &PyArray_Type, &py_V,  &PyArray_Type, &py_genpars)) 
		return NULL;
	if(NULL==py_genpars) 
		return NULL;
	if(NULL==py_V) 
		return NULL;
	
	if (PyArray_NDIM(py_genpars)==1)
	{
		n=py_genpars->dimensions[0];
		if ((n!=11)&&(n!=10))
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 11");
			return NULL;
		}
		if (n==11)
		{
			a=*(double *)(py_genpars->data);
			if (fabs(a-round(a))>1e-15)
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
				return NULL;	
			}
			if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
				return NULL;	
			}
			P.TYPE=(MODEL_TYPE)a;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i+1)*py_genpars->strides[0]);
		}
		else
		{
			P.TYPE=PVC_def_type;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i)*py_genpars->strides[0]);
		}
	}
	else if (PyArray_NDIM(py_genpars)==2)
	{
		n=py_genpars->dimensions[0];
		m=py_genpars->dimensions[1];
		if (n!=1)
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 1x10 or 1x11");
			return NULL;
		}
		if ((m!=11)&&(m!=10))
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 1x10 or 1x11");
			return NULL;
		}
		if (m==11)
		{
			a=*(double *)(py_genpars->data);
			if (fabs(a-round(a))>1e-15)
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
				return NULL;	
			}
			if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
				return NULL;	
			}
			P.TYPE=(MODEL_TYPE)a;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i+1)*py_genpars->strides[1]);
		}
		else
		{
			P.TYPE=PVC_def_type;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i)*py_genpars->strides[1]);
		}
	}
	else
	{
		PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 1x10 or 1x11");
		return NULL;		
	}
	if (PyArray_NDIM(py_V)==1)
	{
		n=py_V->dimensions[0];
		ArrayDim[0]=n;
		py_I=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_DOUBLE);
		for (i=0;i<n;i++)
		{
			V=*(double *)(py_V->data+(i)*py_V->strides[0]);
			double_into_1darray(py_I, i, pv_craze_IV(V, P));
		}
	}
	else if (PyArray_NDIM(py_V)==2)
	{
		n=py_V->dimensions[0];
		m=py_V->dimensions[1];
		ArrayDim[0]=n;
		ArrayDim[1]=m;
		py_I=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
		for (i=0;i<n;i++)
			for (j=0;j<m;j++)
			{
				V=*(double *)(py_V->data+i*py_V->strides[0]+j*py_V->strides[1]);
				double_into_2darray(py_I, i, j, pv_craze_IV(V, P));
			}
	}
	else
		return NULL;
	
	return (PyObject *)py_I;
}

static PyObject *py_pv_craze_vi(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *py_genpars, *py_V, *py_I;
	
	generic_pars P;
	double a,Ii;
	int i, j, n, m;
	npy_intp ArrayDim[2]={0};

	static char *kwlist[] = {
		"I",
		"modelpar",
		NULL
	};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!",kwlist, &PyArray_Type, &py_I,  &PyArray_Type, &py_genpars)) 
		return NULL;
	if(NULL==py_genpars) 
		return NULL;
	if(NULL==py_I) 
		return NULL;
	
	if (PyArray_NDIM(py_genpars)==1)
	{
		n=py_genpars->dimensions[0];
		if ((n!=11)&&(n!=10))
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 11");
			return NULL;
		}
		if (n==11)
		{
			a=*(double *)(py_genpars->data);
			if (fabs(a-round(a))>1e-15)
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
				return NULL;	
			}
			if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
				return NULL;	
			}
			P.TYPE=(MODEL_TYPE)a;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i+1)*py_genpars->strides[0]);
		}
		else
		{
			P.TYPE=PVC_def_type;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i)*py_genpars->strides[0]);
		}
	}
	else if (PyArray_NDIM(py_genpars)==2)
	{
		n=py_genpars->dimensions[0];
		m=py_genpars->dimensions[1];
		if (n!=1)
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 1x10 or 1x11");
			return NULL;
		}
		if ((m!=11)&&(m!=10))
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 1x10 or 1x11");
			return NULL;
		}
		if (m==11)
		{
			a=*(double *)(py_genpars->data);
			if (fabs(a-round(a))>1e-15)
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
				return NULL;	
			}
			if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
			{
				PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
				return NULL;	
			}
			P.TYPE=(MODEL_TYPE)a;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i+1)*py_genpars->strides[1]);
		}
		else
		{
			P.TYPE=PVC_def_type;
			for(i=0; i<10; i++)
				P.pars[i]=*(double *)(py_genpars->data+(i)*py_genpars->strides[1]);
		}
	}
	else
	{
		PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, expect an array of length 1x10 or 1x11");
		return NULL;		
	}
	if (PyArray_NDIM(py_I)==1)
	{
		n=py_I->dimensions[0];
		ArrayDim[0]=n;
		py_V=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_DOUBLE);
		for (i=0;i<n;i++)
		{
			Ii=*(double *)(py_I->data+i*py_I->strides[0]);
			double_into_1darray(py_V, i, pv_craze_VI(Ii, P));
		}
	}
	else if (PyArray_NDIM(py_I)==2)
	{
		n=py_I->dimensions[0];
		m=py_I->dimensions[1];
		ArrayDim[0]=n;
		ArrayDim[1]=m;
		py_V=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
		for (i=0;i<n;i++)
			for (j=0;j<m;j++)
			{
				Ii=*(double *)(py_I->data+i*py_I->strides[0]+j*py_I->strides[1]);
				double_into_2darray(py_V, i, j, pv_craze_VI(Ii, P));
			}
	}
	else
		return NULL;
	return (PyObject *)py_V;
}




static PyObject *py_pv_craze_translate_IEC60891_v2(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *epars_in, *epars_out;
	PyArrayObject *G1, *G2, *T1, *T2;
	double alpha=0, beta=0, a=0, Rs=0, kappa=0;
	double *_G1, *_G2, *_T1, *_T2;
	int N, M, n, i;
	curvepars *in, *out;
	npy_intp ArrayDim[2]={0};
	
	static char *kwlist[] = {
		"esp",
		"T1",
		"G1",
		"T2",
		"G2",
		"alpha",
		"beta",
		"a",
		"Rs",
		"kappa",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!ddddd",kwlist,&PyArray_Type, &epars_in,&PyArray_Type, &T1,&PyArray_Type, &G1,&PyArray_Type, &T2,&PyArray_Type, &G2, &alpha, &beta, &a, &Rs, &kappa)) 
		return NULL;
	if(NULL==epars_in) 
		return NULL;
	if(NULL==G1) 
		return NULL;
	if(NULL==G2) 
		return NULL;
	if(NULL==T1) 
		return NULL;
	if(NULL==T2) 
		return NULL;
		
	if (!IsDoubleArray(epars_in,-1,10))
		return NULL;
	if (!IsDoubleArray(G1,-1,-1))
		return NULL;
	if (!IsDoubleArray(T1,-1,-1))
		return NULL;
	if (!IsDoubleArray(G2,-1,-1))
		return NULL;
	if (!IsDoubleArray(T2,-1,-1))
		return NULL;
	n=epars_in->dimensions[0];
	N=n;
	n=G1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	n=T1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: T1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	M=G2->dimensions[0];
	n=T2->dimensions[0];
	if (n!=M)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G2 and T2 input arrays must be of equal length");
		return NULL;	
	}
	if ((in=malloc(N*sizeof(curvepars)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
	{
		in[i].Voc=*(double *)(epars_in->data+i*epars_in->strides[0]+0*epars_in->strides[1]);
		in[i].Isc=*(double *)(epars_in->data+i*epars_in->strides[0]+1*epars_in->strides[1]);
		in[i].Vmpp=*(double *)(epars_in->data+i*epars_in->strides[0]+2*epars_in->strides[1]);
		in[i].Impp=*(double *)(epars_in->data+i*epars_in->strides[0]+3*epars_in->strides[1]);
		in[i].Vqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+4*epars_in->strides[1]);
		in[i].Iqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+5*epars_in->strides[1]);
		in[i].Vqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+6*epars_in->strides[1]);
		in[i].Iqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+7*epars_in->strides[1]);
		in[i].Gsc=*(double *)(epars_in->data+i*epars_in->strides[0]+8*epars_in->strides[1]);
		in[i].Roc=*(double *)(epars_in->data+i*epars_in->strides[0]+9*epars_in->strides[1]);
	}
	if ((_G1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
		_G1[i]=*(double *)(G1->data+i*G1->strides[0]);
	if ((_T1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
		_T1[i]=*(double *)(T1->data+i*T1->strides[0]);
		
	if ((_G2=malloc(M*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		free(_T1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<M; i++)
		_G2[i]=*(double *)(G2->data+i*G2->strides[0]);
		
	if ((_T2=malloc(M*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		free(_T1);
		free(_G2);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<M; i++)
		_T2[i]=*(double *)(PyArray_GETPTR1(T2, i));
	
	out=pv_craze_translate_IEC60891_v2(in, _G1, _T1, N, _G2, _T2, M, alpha, beta, a, Rs, kappa);
	/* clean up*/
	free(in);
	free(_G1);	
	free(_T1);
	free(_G2);
	free(_T2);
	n=N;
	if (M>n)
		n=M;
	
	ArrayDim[0]=n;
	ArrayDim[1]=10;
	epars_out=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
	for(i=0; i<n; i++)
	{
		double_into_2darray(epars_out, i, 0, out[i].Voc);
		double_into_2darray(epars_out, i, 1, out[i].Isc);
		double_into_2darray(epars_out, i, 2, out[i].Vmpp);
		double_into_2darray(epars_out, i, 3, out[i].Impp);
		double_into_2darray(epars_out, i, 4, out[i].Vqmp1);
		double_into_2darray(epars_out, i, 5, out[i].Iqmp1);
		double_into_2darray(epars_out, i, 6, out[i].Vqmp2);
		double_into_2darray(epars_out, i, 7, out[i].Iqmp2);
		double_into_2darray(epars_out, i, 8, out[i].Gsc);
		double_into_2darray(epars_out, i, 9, out[i].Roc);
	}
	pv_craze_free(out); // pvcraze allocated out
	return (PyObject *)epars_out;
} 

static PyObject *py_pv_craze_fit_IEC60891_v2(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *epars_in, *epars_out;
	PyArrayObject *G1, *T1;
	PyObject *ret;
	double alpha=1, beta=1, a=1, Rs=1, kappa=1, Tstd=25, Gstd=1000;
	double *_G1, *_T1;
	int N, n, i, fN=0, off=0;
	curvepars *in, out;
	npy_intp ArrayDim[1]={10};
	
	static char *kwlist[] = {
		"esp",
		"T",
		"G",
		"Tstd",
		"Gstd",
		"furthest_N",
		"offset",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!|ddii",kwlist,&PyArray_Type, &epars_in,&PyArray_Type, &T1,&PyArray_Type, &G1, &Tstd, &Gstd, &fN, &off)) 
		return NULL;
	if(NULL==epars_in) 
		return NULL;
	if(NULL==G1) 
		return NULL;
	if(NULL==T1) 
		return NULL;
		
	if (fN<0)
		fN=0;
		
	if (!IsDoubleArray(epars_in,-1,10))
		return NULL;
	n=epars_in->dimensions[0];
	N=n;
	if (!IsDoubleArray(G1,-1,-1))
		return NULL;
	n=G1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	if (!IsDoubleArray(T1,-1,-1))
		return NULL;
	n=T1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: T1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	
	
	if ((in=malloc(N*sizeof(curvepars)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_G1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_T1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
	{
		in[i].Voc=*(double *)(epars_in->data+i*epars_in->strides[0]+0*epars_in->strides[1]);
		in[i].Isc=*(double *)(epars_in->data+i*epars_in->strides[0]+1*epars_in->strides[1]);
		in[i].Vmpp=*(double *)(epars_in->data+i*epars_in->strides[0]+2*epars_in->strides[1]);
		in[i].Impp=*(double *)(epars_in->data+i*epars_in->strides[0]+3*epars_in->strides[1]);
		in[i].Vqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+4*epars_in->strides[1]);
		in[i].Iqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+5*epars_in->strides[1]);
		in[i].Vqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+6*epars_in->strides[1]);
		in[i].Iqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+7*epars_in->strides[1]);
		in[i].Gsc=*(double *)(epars_in->data+i*epars_in->strides[0]+8*epars_in->strides[1]);
		in[i].Roc=*(double *)(epars_in->data+i*epars_in->strides[0]+9*epars_in->strides[1]);
		_G1[i]=*(double *)(G1->data+i*G1->strides[0]);
		_T1[i]=*(double *)(T1->data+i*T1->strides[0]);
	}
	out=pv_craze_fit_IEC60891_v2(in, _G1, _T1, N, Tstd, Gstd,&alpha, &beta, &a, &Rs, &kappa, fN,off);
	/* clean up*/
	free(in);
	free(_G1);	
	free(_T1);
	
	
	epars_out=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_DOUBLE);
	//build dictionary with std pars alpha, beta
	double_into_1darray(epars_out, 0, out.Voc);
	double_into_1darray(epars_out, 1, out.Isc);
	double_into_1darray(epars_out, 2, out.Vmpp);
	double_into_1darray(epars_out, 3, out.Impp);
	double_into_1darray(epars_out, 4, out.Vqmp1);
	double_into_1darray(epars_out, 5, out.Iqmp1);
	double_into_1darray(epars_out, 6, out.Vqmp2);
	double_into_1darray(epars_out, 7, out.Iqmp2);
	double_into_1darray(epars_out, 8, out.Gsc);
	double_into_1darray(epars_out, 9, out.Roc);
	
	
	ret=Py_BuildValue("{s:O,s:d,s:d,s:d,s:d,s:d}","esp_std", epars_out, "alpha", alpha,
		"beta", beta, "a", a, "Rs", Rs, "kappa", kappa);
	Py_DECREF(epars_out);
	return ret;
} 
	
static PyObject *py_pv_craze_translate_IEC60891_v3(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *epars_in, *epars_out;
	PyArrayObject *G1, *G2, *T1, *T2;
	double Vocstd=0, Tstd=25, Gstd=1000, alpha=0, beta=0, B1=0, B2=0, Rs=0, kappa=0;
	double *_G1, *_G2, *_T1, *_T2;
	int N, M, n, i;
	curvepars *in, *out;
	npy_intp ArrayDim[2]={0};
	
	static char *kwlist[] = {
		"esp",
		"T1",
		"G1",
		"T2",
		"G2",
		"Vocstd",
		"alpha",
		"beta",
		"B1",
		"B2",
		"Rs",
		"kappa",
		"Gstd",
		"Tstd",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!dddddd|dd",kwlist,&PyArray_Type, &epars_in,&PyArray_Type, &T1,&PyArray_Type, &G1,&PyArray_Type, &T2,&PyArray_Type, &G2, &Vocstd, &alpha, &beta, &B1, &B2, &Rs, &kappa, &Gstd, &Tstd)) 
		return NULL;
	if(NULL==epars_in) 
		return NULL;
	if(NULL==G1) 
		return NULL;
	if(NULL==G2) 
		return NULL;
	if(NULL==T1) 
		return NULL;
	if(NULL==T2) 
		return NULL;
	
	if (!IsDoubleArray(epars_in, -1,10))
		return NULL;
	if (!IsDoubleArray(G1,-1,-1))
		return NULL;
	if (!IsDoubleArray(T1,-1,-1))
		return NULL;
	if (!IsDoubleArray(G2,-1,-1))
		return NULL;
	if (!IsDoubleArray(T2,-1,-1))
		return NULL;
	n=epars_in->dimensions[0];
	N=n;
	n=G1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	n=T1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: T1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	M=G2->dimensions[0];
	n=T2->dimensions[0];
	if (n!=M)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G2 and T2 input arrays must be of equal length");
		return NULL;	
	}
	if ((in=malloc(N*sizeof(curvepars)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_G1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_T1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
	{
		in[i].Voc=*(double *)(epars_in->data+i*epars_in->strides[0]+0*epars_in->strides[1]);
		in[i].Isc=*(double *)(epars_in->data+i*epars_in->strides[0]+1*epars_in->strides[1]);
		in[i].Vmpp=*(double *)(epars_in->data+i*epars_in->strides[0]+2*epars_in->strides[1]);
		in[i].Impp=*(double *)(epars_in->data+i*epars_in->strides[0]+3*epars_in->strides[1]);
		in[i].Vqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+4*epars_in->strides[1]);
		in[i].Iqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+5*epars_in->strides[1]);
		in[i].Vqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+6*epars_in->strides[1]);
		in[i].Iqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+7*epars_in->strides[1]);
		in[i].Gsc=*(double *)(epars_in->data+i*epars_in->strides[0]+8*epars_in->strides[1]);
		in[i].Roc=*(double *)(epars_in->data+i*epars_in->strides[0]+9*epars_in->strides[1]);
		_G1[i]=*(double *)(G1->data+i*G1->strides[0]);
		_T1[i]=*(double *)(T1->data+i*T1->strides[0]);
	}
		
	if ((_G2=malloc(M*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		free(_T1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_T2=malloc(M*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		free(_T1);
		free(_G2);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<M; i++)
	{
		_G2[i]=*(double *)(G2->data+i*G2->strides[0]);
		_T2[i]=*(double *)(T2->data+i*T2->strides[0]);
	}
	out=pv_craze_translate_IEC60891_v3(in, _G1, _T1, N, _G2, _T2, M, Vocstd, Tstd, Gstd, alpha, beta, B1, B2, Rs, kappa);
	/* clean up*/
	free(in);
	free(_G1);	
	free(_T1);
	free(_G2);
	free(_T2);
	n=N;
	if (M>n)
		n=M;
	
	ArrayDim[0]=n;
	ArrayDim[1]=10;
	epars_out=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
	for(i=0; i<n; i++)
	{
		double_into_2darray(epars_out, i, 0, out[i].Voc);
		double_into_2darray(epars_out, i, 1, out[i].Isc);
		double_into_2darray(epars_out, i, 2, out[i].Vmpp);
		double_into_2darray(epars_out, i, 3, out[i].Impp);
		double_into_2darray(epars_out, i, 4, out[i].Vqmp1);
		double_into_2darray(epars_out, i, 5, out[i].Iqmp1);
		double_into_2darray(epars_out, i, 6, out[i].Vqmp2);
		double_into_2darray(epars_out, i, 7, out[i].Iqmp2);
		double_into_2darray(epars_out, i, 8, out[i].Gsc);
		double_into_2darray(epars_out, i, 9, out[i].Roc);
	}
	pv_craze_free(out);
	return (PyObject *)epars_out;
} 


static PyObject *py_pv_craze_fit_IEC60891_v3(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *epars_in, *epars_out;
	PyArrayObject *G1, *T1;
	PyObject *ret;
	double alpha=0, beta=0, B1=0, B2=0, Rs=0, kappa=0, Tstd=25, Gstd=1000;
	double *_G1, *_T1;
	int N, n, i, fN=0, off=0;
	curvepars *in, out;
	npy_intp ArrayDim[1]={10};
	
	static char *kwlist[] = {
		"esp",
		"T",
		"G",
		"Tstd",
		"Gstd",
		"furthest_N",
		"offset",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!|ddii",kwlist,&PyArray_Type, &epars_in, &PyArray_Type, &T1, &PyArray_Type, &G1, &Tstd, &Gstd,&fN,&off)) 
		return NULL;
	if(NULL==epars_in) 
		return NULL;
	if(NULL==G1) 
		return NULL;
	if(NULL==T1) 
		return NULL;
	if (fN<0)
		fN=0;	
		
	if (!IsDoubleArray(epars_in,-1,10))
		return NULL;
	if (!IsDoubleArray(G1, -1,-1))
		return NULL;
	if (!IsDoubleArray(T1, -1,-1))
		return NULL;
	n=epars_in->dimensions[0];
	N=n;
	n=G1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	n=T1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: T1 input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	if ((in=malloc(N*sizeof(curvepars)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_G1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_T1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
	{
		in[i].Voc=*(double *)(epars_in->data+i*epars_in->strides[0]+0*epars_in->strides[1]);
		in[i].Isc=*(double *)(epars_in->data+i*epars_in->strides[0]+1*epars_in->strides[1]);
		in[i].Vmpp=*(double *)(epars_in->data+i*epars_in->strides[0]+2*epars_in->strides[1]);
		in[i].Impp=*(double *)(epars_in->data+i*epars_in->strides[0]+3*epars_in->strides[1]);
		in[i].Vqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+4*epars_in->strides[1]);
		in[i].Iqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+5*epars_in->strides[1]);
		in[i].Vqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+6*epars_in->strides[1]);
		in[i].Iqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+7*epars_in->strides[1]);
		in[i].Gsc=*(double *)(epars_in->data+i*epars_in->strides[0]+8*epars_in->strides[1]);
		in[i].Roc=*(double *)(epars_in->data+i*epars_in->strides[0]+9*epars_in->strides[1]);
		_G1[i]=*(double *)(G1->data+i*G1->strides[0]);
		_T1[i]=*(double *)(T1->data+i*T1->strides[0]);
		//printf("%d %e %e %e %e %e %e %e %e %e %e %e %e\n", i, _G1[i], _T1[i], in[i].Voc,in[i].Isc,in[i].Vmpp,in[i].Impp,in[i].Vqmp1,in[i].Iqmp1,in[i].Vqmp2,in[i].Iqmp2,in[i].Gsc,in[i].Roc);
	}
	out=pv_craze_fit_IEC60891_v3(in, _G1, _T1, N, Tstd, Gstd, &alpha, &beta, &B1, &B2, &Rs, &kappa, fN,off);
	/* clean up*/
	free(in);
	free(_G1);	
	free(_T1);
	
	
	epars_out=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_DOUBLE);
	//build dictionary with std pars alpha, beta
	double_into_1darray(epars_out, 0, out.Voc);
	double_into_1darray(epars_out, 1, out.Isc);
	double_into_1darray(epars_out, 2, out.Vmpp);
	double_into_1darray(epars_out, 3, out.Impp);
	double_into_1darray(epars_out, 4, out.Vqmp1);
	double_into_1darray(epars_out, 5, out.Iqmp1);
	double_into_1darray(epars_out, 6, out.Vqmp2);
	double_into_1darray(epars_out, 7, out.Iqmp2);
	double_into_1darray(epars_out, 8, out.Gsc);
	double_into_1darray(epars_out, 9, out.Roc);
	
	
	ret=Py_BuildValue("{s:O,s:d,s:d,s:d,s:d,s:d,s:d}",
		"esp_std", epars_out, "alpha", alpha,
		"beta", beta, "B1", B1, "B2", B2, 
		"Rs", Rs, "kappa", kappa);
	Py_DECREF(epars_out);
	return ret;
} 

static PyObject *py_pv_craze_translate_diode_GT(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *mpars_in, *mpars_out;
	PyArrayObject *G, *T;
	double Tstd=-300, Gstd=-1, alpha=0, Ea=0;
	double *_G1, *_T1;
	int N, n, i, j;
	generic_pars in, *out;
	npy_intp ArrayDim[2]={0};
	char unit ='K', *sunit=NULL;
	
	static char *kwlist[] = {
		"modelpar",
		"Tstd",
		"Gstd",
		"T",
		"G",
		"alpha",
		"Ea",
		"Tunit",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!ddO!O!dd|s",kwlist,&PyArray_Type, &mpars_in,&Tstd, &Gstd, &PyArray_Type, &T,&PyArray_Type, &G,&alpha, &Ea, &sunit)) 
		return NULL;
	if(NULL==mpars_in) 
		return NULL;
	if(NULL==G) 
		return NULL;
	if(NULL==T) 
		return NULL;
	if (sunit)
	{
		int cc=0;
		cc+=(strncmp("C",		sunit, 10)==0);
		cc+=(strncmp("Celcius",	sunit, 10)==0);
		if (cc)
			unit='C';
		else 
		{
			cc+=(strncmp("K",		sunit, 10)==0);
			cc+=(strncmp("Kelvin",	sunit, 10)==0);
			if (cc)
				unit='K';
			else
			{
				PyErr_SetString(PyExc_ValueError, "Error: expect a temperature unit(\"C\"/\"Celcius\"/\"K\"/\"Kelvin\")");
				return NULL;	
			}
		}
	}
	switch(unit)
	{
		case 'C':					
			if (Tstd<=-273.15)
			{
				PyErr_SetString(PyExc_ValueError, "Error: temperature value out of range");
				return NULL;
			}
		default:		
			if (Tstd<=0)
			{
				PyErr_SetString(PyExc_ValueError, "Error: temperature value out of range");
				return NULL;
			}
	}
	if (Gstd<=0)
	{
		PyErr_SetString(PyExc_ValueError, "Error: irradiance value out of range");
		return NULL;
	}
	
	if (!IsDoubleArray(mpars_in,11,-1))
		return NULL;
	if (!IsDoubleArray(G,-1,-1))
		return NULL;
	if (!IsDoubleArray(T,-1,-1))
		return NULL;
	
	N=G->dimensions[0];
	n=T->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: T input array must be of length n (number of rows in G input array)");
		return NULL;	
	}
	
	{ // read input model parameters
		double a;
		a=*(double *)(mpars_in->data);
		if (fabs(a-round(a))>1e-15)
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, first number must be a round number");
			return NULL;	
		}
		if (((MODEL_TYPE)a<0)||((MODEL_TYPE)a>=M_NONE))
		{
			PyErr_SetString(PyExc_ValueError, "Error: wrong parameter array format, unknown model index");
			return NULL;	
		}
		in.TYPE=(MODEL_TYPE)a;
		for(i=0; i<10; i++)
			in.pars[i]=*(double *)(mpars_in->data+(i+1)*mpars_in->strides[0]);
	}
	
	if ((_G1=malloc(N*sizeof(double)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_T1=malloc(N*sizeof(double)))==NULL)
	{
		free(_G1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
	{
		_G1[i]=*(double *)(G->data+i*G->strides[0]);
		_T1[i]=*(double *)(T->data+i*T->strides[0]);
	}
	out=pv_craze_translate_diode_GT(in,Gstd,Tstd, _G1, _T1, N, alpha, Ea, unit);
	free(_G1);
	free(_T1);
	
	
	ArrayDim[0]=N;
	ArrayDim[1]=11;
	mpars_out=(PyArrayObject *) PyArray_SimpleNew(2, ArrayDim, NPY_DOUBLE);
	// output diode model parameters
	for(i=0; i<N; i++)
	{
		double_into_2darray(mpars_out, i, 0, (double)out[i].TYPE);
		for(j=0; j<10; j++)
			double_into_2darray(mpars_out, i, j+1, out[i].pars[j]);
	}
	pv_craze_free(out);
	return (PyObject *)mpars_out;
} 

static PyObject *py_pv_craze_fit_diode_GT(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *epars_in, *mpars_out;
	PyObject *ret;
	PyArrayObject *G1, *T1;
	double alpha=0, Ea=0, Tstd=298.15, Gstd=1000;
	double *_G1, *_T1;
	int N, n, i, fN=0, off=0;
	curvepars *in;
	generic_pars out;
	npy_intp ArrayDim[1];
	char unit ='K', *sunit=NULL;
				
			
	
	static char *kwlist[] = {
		"esp",
		"T",
		"G",
		"Tstd",
		"Gstd",
		"Tunit",
		"furthest_N",
		"offset",
		NULL
	};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!|ddsii",kwlist,&PyArray_Type, &epars_in,&PyArray_Type, &T1,&PyArray_Type, &G1,&Tstd, &Gstd,&sunit,&fN,&off)) 
		return NULL;
	if(NULL==epars_in) 
		return NULL;
	if(NULL==G1) 
		return NULL;
	if(NULL==T1) 
		return NULL;
	if (sunit)
	{
		int cc=0;
		cc+=(strncmp("C",		sunit, 10)==0);
		cc+=(strncmp("Celcius",	sunit, 10)==0);
		if (cc)
			unit='C';
		else 
		{
			cc+=(strncmp("K",		sunit, 10)==0);
			cc+=(strncmp("Kelvin",	sunit, 10)==0);
			if (cc)
				unit='K';
			else
			{
				PyErr_SetString(PyExc_ValueError, "Error: expect a temperature unit(\"C\"/\"Celcius\"/\"K\"/\"Kelvin\")");
				return NULL;	
			}
		}
	}
	if (fN<0)
		fN=0;	
	
	if (!IsDoubleArray(epars_in,-1,10))
		return NULL;
	if (!IsDoubleArray(G1,-1,-1))
		return NULL;
	if (!IsDoubleArray(T1,-1,-1))
		return NULL;
	n=epars_in->dimensions[0];
	N=n;
	n=G1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: G input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	n=T1->dimensions[0];
	if (n!=N)
	{
		PyErr_SetString(PyExc_ValueError, "Error: T input array must be of length n (number of rows in epar input array)");
		return NULL;	
	}
	if ((in=malloc(N*sizeof(curvepars)))==NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_G1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	if ((_T1=malloc(N*sizeof(double)))==NULL)
	{
		free(in);
		free(_G1);
		PyErr_SetString(PyExc_ValueError, "Error: failed to allocate memory");
		return NULL;	
	}
	for(i=0; i<N; i++)
	{
		in[i].Voc=*(double *)(epars_in->data+i*epars_in->strides[0]+0*epars_in->strides[1]);
		in[i].Isc=*(double *)(epars_in->data+i*epars_in->strides[0]+1*epars_in->strides[1]);
		in[i].Vmpp=*(double *)(epars_in->data+i*epars_in->strides[0]+2*epars_in->strides[1]);
		in[i].Impp=*(double *)(epars_in->data+i*epars_in->strides[0]+3*epars_in->strides[1]);
		in[i].Vqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+4*epars_in->strides[1]);
		in[i].Iqmp1=*(double *)(epars_in->data+i*epars_in->strides[0]+5*epars_in->strides[1]);
		in[i].Vqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+6*epars_in->strides[1]);
		in[i].Iqmp2=*(double *)(epars_in->data+i*epars_in->strides[0]+7*epars_in->strides[1]);
		in[i].Gsc=*(double *)(epars_in->data+i*epars_in->strides[0]+8*epars_in->strides[1]);
		in[i].Roc=*(double *)(epars_in->data+i*epars_in->strides[0]+9*epars_in->strides[1]);
		_G1[i]=*(double *)(G1->data+i*G1->strides[0]);
		_T1[i]=*(double *)(T1->data+i*T1->strides[0]);
	}
	out=pv_craze_fit_diode_GT(in, _G1, _T1, N, &alpha, &Ea,Tstd,Gstd, unit,fN,off);
	/* clean up*/
	free(in);
	free(_G1);	
	free(_T1);
	
	
	//build dictionary with model pars alpha, Ea	
	ArrayDim[0]=11;
	mpars_out=(PyArrayObject *) PyArray_SimpleNew(1, ArrayDim, NPY_DOUBLE);
	
	
	double_into_1darray(mpars_out, 0, (double)out.TYPE);
	for(i=0; i<10; i++)
		double_into_1darray(mpars_out, i+1, out.pars[i]);
	
	ret=Py_BuildValue("{s:O,s:d,s:d}","modelpar", mpars_out, "alpha", alpha,"Ea", Ea);
	Py_DECREF(mpars_out);
	return ret;
} 
