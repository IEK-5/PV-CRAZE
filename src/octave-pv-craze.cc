/*
    octave interface to PV-CRAZE
    Copyright (C) 2024  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich
*/
/*****************************************************************              
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
 +  IEK-5 PHOTOVOLTAIK                                           *              
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

/*****************************************************************           
 * FUNCTION:                                                     *                
 * octave interface for IV analysis                              *         
 *****************************************************************/  
#include <octave/oct.h>
#include <octave/variables.h>
#include <octave/Cell.h>
#include <math.h>
extern "C"
{
	#include <pvcraze.h>
}

// default model
MODEL_TYPE PVC_def_type=M_PVBEZIER;
typedef enum {IFILTER,SFILTER, MPPMETHOD, MPPDV, MPPDI, SCDV, OCDI, MODEL, RESAMPLE, NONE} OPT_PARS;


/* monolithic c++ file defining several octave commands
 * to allow multiple function definitions from one oct file we use 
 * PKG_ADD directives to call autoload. Here follows one directive per 
 * octave function defined in this file
 */
// PKG_ADD: autoload("pvcraze_quiet","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_verbose","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_verbose_setreport","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_verbose_setprogress","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_verbose_settiming","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_verbose_sethead","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_lookupmodel","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_reverselookup","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_flagtable","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_analyze","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_modelfit","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_curvemodel","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_modelpars","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_modelerror","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_errorcurve","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_iv","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_vi","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_fit_IEC60891_v2","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_fit_IEC60891_v3","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_translate_IEC60891_v2","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_translate_IEC60891_v3","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_fit_diode_GT","pv_craze.oct");
// PKG_ADD: autoload("pvcraze_translate_diode_GT","pv_craze.oct");


DEFUN_DLD (pvcraze_quiet, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} pvcraze_quiet ()\n\
set pvcraze verbose level to quiet\n\
@end deftypefn")
{
	pv_craze_quiet();
	return octave_value_list ();
}
DEFUN_DLD (pvcraze_verbose, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} pvcraze_quiet ()\n\
set pvcraze verbose level to verbose\n\
@end deftypefn")
{
	pv_craze_verbose();
	return octave_value_list ();
}

DEFUN_DLD (pvcraze_verbose_setreport, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} pvcraze_verbose_setreport (@var{enable})\n\
Set verbose level to print a report. The argument @var{enable} should be 0 (print nothing) or 1 (print stuff).\n\
@end deftypefn")
{
  	int nargin = args.length ();
  	if (nargin != 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	if (args(0).int_value())
		pv_craze_verbosereport(1);
	else
		pv_craze_verbosereport(0);
	return octave_value_list ();
}
DEFUN_DLD (pvcraze_verbose_setprogress, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} pvcraze_verbose_setprogress (@var{enable})\n\
Set verbose level to print a progress bar. The argument @var{enable} should be 0 (print nothing) or 1 (print stuff).\n\
@end deftypefn")
{
  	int nargin = args.length ();
  	if (nargin != 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	if (args(0).int_value())
		pv_craze_verboseprogress(1);
	else
		pv_craze_verboseprogress(0);
	return octave_value_list ();
}
DEFUN_DLD (pvcraze_verbose_settiming, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} pvcraze_verbose_settiming (@var{enable})\n\
Set verbose level to print timing information. The argument @var{enable} should be 0 (print nothing) or 1 (print stuff).\n\
@end deftypefn")
{
  	int nargin = args.length ();
  	if (nargin != 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	if (args(0).int_value())
		pv_craze_verbosetiming(1);
	else
		pv_craze_verbosetiming(0);
	return octave_value_list ();
}
DEFUN_DLD (pvcraze_verbose_sethead, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} pvcraze_verbose_sethead (@var{enable})\n\
Set verbose level to print a header. The argument @var{enable} should be 0 (print nothing) or 1 (print stuff).\n\
@end deftypefn")
{
  	int nargin = args.length ();
  	if (nargin != 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	if (args(0).int_value())
		pv_craze_verbosehead(1);
	else
		pv_craze_verbosehead(0);
	return octave_value_list ();
}

DEFUN_DLD (pvcraze_lookupmodel, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {@var{Mid} =} pvcraze_lookupmodel (@var{model-name})\n\
pvcraze_lookupmodel searches for the model id corresponding to model-name.\n\n\
Input:\n\
@itemize \n\
@item @var{model-name}: Name of the model, e.g. PV-Bezier or Diode\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{Mid}: Model ID\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	MODEL_TYPE T;
	octave_value_list mid;
	std::string mname;
  	if (nargin < 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	mname=args(0).string_value();
	T=pv_craze_lookupmodel((char *)mname.c_str());
	if (T==M_NONE)
	{
		fprintf(stderr, "Error: model %s does not exist\n", mname.c_str());
		return octave_value_list ();
	}
	else
		return octave_value((int) T);
	
}

DEFUN_DLD (pvcraze_reverselookup, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {@var{model-name} =} pvcraze_lookupmodel (@var{Mid})\n\
pvcraze_reverselookup searches for the model-name corresponding to a model id.\n\n\
Input:\n\
@itemize \n\
@item @var{Mid}: Model ID\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{model-name}: Name of the model, e.g. PV-Bezier or Diode\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	MODEL_TYPE T;
	octave_value_list mid;
	charMatrix oname;
  	char *name;
  	int i,n;
  	if (nargin < 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	T=(MODEL_TYPE)args(0).int_value();
	name=pv_craze_reverselookup(T);
	if (name==NULL)
	{
		fprintf(stderr, "Error: model with id %d does not exist\n", (int)T);
		return octave_value_list ();
	}
		
	n=strlen(name);
	oname=charMatrix(1,(octave_idx_type)n);
	for (i=0;i<n;i++)
		oname(0,i)=name[i];
	
	return octave_value(oname);
	
}

DEFUN_DLD (pvcraze_flagtable, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{flags},@var{messages}, @var{err}]=} pvcraze_flagtable ( )\n\
pvcraze_flagtable generates a table of flag codes and messages for the\n\
installed pv-craze library version.\n\n\
Input:\n\
None\n\
Output:\n\
@itemize \n\
@item @var{flags}: array with principle flag values (powers of two, actual flags may be combinations of principle flags)\n\
@item @var{messages}: Message (array with messages) corresponding to flags\n\
@item @var{err}: Boolean array whether the flag represents an error\n\
@end itemize\n\
@end deftypefn")
{
	AnalyzeFlag *T;
	int i,j,N, n, nm;
	charMatrix message;
	int32NDArray FF;
	boolMatrix err;
	octave_value_list retval;
	
	T=pv_craze_flagtable(&N);
	if (!T)
	{
		fprintf(stderr, "Error: could not generate a flag table, pv-crze appears to be broken!\n");
		return retval;	
	}
	dim_vector dim(N,1);
	FF=int32NDArray(dim);
	err=boolMatrix(N,1);
	nm=0;
	for (i=0;i<N;i++)
	{
		n=strlen(T[i].message);
		if (n>nm)
			nm=n;
	}
	message=charMatrix(N,nm);
	
	
	for (i=0;i<N;i++)
	{
		FF(i)=(int)T[i].FLAG;
		
		n=strlen(T[i].message);	
		for (j=0;j<n;j++)
			message(i,j)=T[i].message[j];
		for (;j<nm;j++)
			message(i,j)='\0';
		
		err(i,0)=(bool)T[i].ERR;
	}
	pv_craze_free(T); // allocated in pv_craze
	i=0;
	if (nargout>i)
		retval(i++)=octave_value(FF);
	if (nargout>i)
		retval(i++)=octave_value(message);
	if (nargout>i) 
		retval(i++)=octave_value(err);
	return retval;	
}

DEFUN_DLD (pvcraze_analyze, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{esp},@var{flags},@var{modelpar}, @var{err}] =} pvcraze_analyze (@var{IV}, [optional-parameters])\n\
pvcraze_analyze reads solar cell IV characteristics, determines an extended set of solar cell parameters \
and does sanity checks on the data.\n\n\
Input:\n\
@itemize \n\
@item @var{IV}: Cell array with arrays with two columns (voltage,current)\n\
@item \'models\':  (Optional) Takes a model string (this option may be passed multiple times)\n\
@item \'mpp-method\': (Optional) Select mpp method, can be \'f\' (fast, not robust), \'r\' (robust), 'p' (parametric), or 's' Savitsky-Golay \n\
@item \'S_filter\': (Optional) Takes boolean, Apply a 3rd order median filter on the data\n\
@item \'I_filter\':  (Optional) Takes positive float. Filter out subsequent identical voltage and current values. Specify minimum deviation. default 1e-12\n\
@item \'mpp-dv\':  (Optional) Takes float [0-1]. The mpp fitting range starts at (1-mpp-dv) Vmpp. default 0.12\n\
@item \'mpp-di\':  (Optional) Takes float [0-1]. The mpp fitting range ends at (1-mpp-di) Impp. default  0.12\n\
@item \'sc-dv\':  (Optional) Takes float [0-1]. Voltage range to determine sc slope between 0 and sc-dv Vqmp1, default 0.1 \n\
@item \'oc-di\':  (Optional) Takes float [0-1]. Current range to determine oc slope between 0 and oc-di Iqmp2, default 0.1 \n\
@item \'resample\': (Optional) Takes an integer >20. Resample the IV characteristics using piece-wise linear interpolation. \n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{esp}: Extended Solar Cell Parameter Matrix:\n\
\t\t[Voc;Isc;Vmpp;Impp;Vqmpp-,Iqmpp-,Vqmpp+,Iqmpp+,Gsc,Roc]\n\
@item @var{flags}: Flags for the analyzed curves\n\
@item @var{modelpar}: Cell array with model parameters (one column per model, one row per IV)\n\
@item @var{err}: Model fit errors (one column per model one row per IV)\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Cell IV;
	Matrix iv, EPar, FErr, PErr;
	int32NDArray Flags; 
	Cell MPar, MEPar;
	OPT_PARS parse=NONE;
	ProcessSettings Set;
	CurveProp * Prop;
	curve *C=NULL;
	int Nc=0, nm=0;
	int i, j, k;
  	if (nargin < 1)
	{
    		print_usage ();
		return octave_value_list ();
	}
	
	Set=pv_craze_initsettings();
	IV=args(0).cell_value();
	
	for (i=1;i<nargin;i++)
	{
		if (parse==NONE)
		{
			std::string parname;
			parname=args(i).string_value();
			if (parname.compare("models") == 0)
				parse=MODEL;
			else if (parname.compare("mpp-method") == 0)
				parse=MPPMETHOD;
			else if (parname.compare("I_filter") == 0)
				parse=IFILTER;  			
			else if (parname.compare("S_filter") == 0)
				parse=SFILTER;
			else if (parname.compare("mpp-dv") == 0)
				parse=MPPDV;
			else if (parname.compare("mpp-di") == 0)
				parse=MPPDI;
			else if (parname.compare("sc-dv") == 0)
				parse=SCDV;
			else if (parname.compare("oc-di") == 0)
				parse=OCDI;
  			else if (parname.compare("resample") == 0)
				parse=RESAMPLE;
			else
				printf("Unknown argument\n");
		}
		else
		{
			switch (parse)
			{
				case MPPMETHOD:
					Set.analyze.am=args(i).string_value().at(0);
					parse=NONE;
					break;	
				case IFILTER:
					if (args(i).double_value()>0)
					{
						Set.analyze.Identicals=1;
						Set.analyze.IdentEps=args(i).double_value();
					}
					else
						Set.analyze.Identicals=0;	
					parse=NONE;					
					break;	
				case SFILTER:
					Set.analyze.SpikeFilter=(args(i).int_value()>0);	
					parse=NONE;				
					break;	
				case MPPDV:
					Set.analyze.MPP_f_V=args(i).double_value();
					parse=NONE;
					break;	
				case MPPDI:
					Set.analyze.MPP_f_I=args(i).double_value();
					parse=NONE;
					break;	
				case SCDV:
					Set.analyze.SC_f_V=args(i).double_value();
					parse=NONE;
					break;	
				case OCDI:
					Set.analyze.OC_f_I=args(i).double_value();
					parse=NONE;
					break;	
				case RESAMPLE:
					Set.analyze.Resample=args(i).int_value();
					if (Set.analyze.Resample<20)
					{
						fprintf(stderr, "Warning: invalid resample numer (should be larger than 20)\n");
						Set.analyze.Resample=0;
					}
					parse=NONE;
					break;	
				case MODEL:
				{	
					charMatrix ch = args(i).char_matrix_value ();
					std::string tmp = ch.row_as_string (0);			
					pv_craze_addmodeltosettings(&Set, (char *)tmp.c_str());
					nm++;
					parse=NONE;
					break;	
				}
				default:
					parse=NONE;
					break;			
			}
		}
	
	}
	if (parse!=NONE)
		fprintf(stderr, "Warning: Incomplete optional argument list\n");
	Nc=IV.numel();
	if (Nc==0)
	{
		pv_craze_freesettings(&Set);
		return octave_value_list ();
	}
	
	EPar=Matrix(Nc,10);
	dim_vector dim(Nc,1);
	Flags=int32NDArray(dim);
	if (Set.L.N)
	{
		MPar=Cell(Nc, Set.L.N);
		FErr=Matrix(Nc, Set.L.N);
		MEPar=Cell(Nc, Set.L.N);
		PErr=Matrix(Nc, Set.L.N);
	}
	else
	{
		MPar=Cell(0, 0);
		FErr=Matrix(0, 0);
		MEPar=Cell(0, 0);
		PErr=Matrix(0, 0);
	}
		

	// create a list of curves
	if ((C=(curve *)malloc(Nc*sizeof(curve)))==NULL)
	{
		pv_craze_freesettings(&Set);
		return octave_value_list ();
	}
	
	for (i=0;i<Nc;i++)
	{
		iv=IV(i,0).matrix_value();
		if ((iv.cols()!=2)||(iv.rows()==0))
		{
			C[i].P=NULL;
			C[i].N=0;
		}
		else
		{
			if ((C[i].P=(point*) malloc((int)iv.rows()*sizeof(point)))==NULL)
			{
				for (i=i-1;i>=0;i--)
					free(C[i].P);
				free(C);
				pv_craze_freesettings(&Set);
				return octave_value_list ();
			}				
			C[i].N=iv.rows();
			for(j=0;j<C[i].N;j++)
			{
				C[i].P[j].x=iv(j,0);
				C[i].P[j].y=iv(j,1);
			}
		}		
	}
	// analyze and construct a list of curve properties
	Prop=pv_craze_analyze(C, Nc, &Set);
	
	// get data out to matrices
	for (i=0;i<Nc;i++)
	{
		free(C[i].P); // free curve data
		
		EPar(i,0)=Prop[i].par.Voc;
		EPar(i,1)=Prop[i].par.Isc;
		EPar(i,2)=Prop[i].par.Vmpp;
		EPar(i,3)=Prop[i].par.Impp;
		EPar(i,4)=Prop[i].par.Vqmp1;
		EPar(i,5)=Prop[i].par.Iqmp1;
		EPar(i,6)=Prop[i].par.Vqmp2;
		EPar(i,7)=Prop[i].par.Iqmp2;
		EPar(i,8)=Prop[i].par.Gsc;
		EPar(i,9)=Prop[i].par.Roc;
		
		Flags(i)=Prop[i].A_FLAGS;
		
		if (Set.L.N)
		{			
			for (k=0;k<Set.L.N;k++)
			{
				Matrix par, solpar;
					
				par=Matrix(1,11);				
				par(0,0)=(double)Set.L.List[k].M;
				for (j=0;j<pv_craze_numfields(Set.L.List[k].M);j++)
					par(0,j+1)=Prop[i].P[k].pars[j];
				for (;j<10;j++)
					par(0,j+1)=0;				
				FErr(i,k)=Prop[i].E[k];
			}
		}
		pv_craze_freeprop(&Prop[i]);
	}
	free(C);
	pv_craze_free(Prop);
	pv_craze_freesettings(&Set);
	
	
	i=0;
	if (nargout>i)
		retval(i++)=octave_value(EPar);
	if (nargout>i)
		retval(i++)=octave_value(Flags);
	if (nargout>i) 
		retval(i++)=octave_value(MPar);
	if (nargout>i) 
		retval(i++)=octave_value(FErr);
		
	return retval;	
	
}


DEFUN_DLD (pvcraze_modelfit, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{modelpar},  @var{err}] =} pvcraze_modelfit (@var{esp}, @var{models}, ...)\n\
pvcraze_modelfit takes extended solar cell parameters and fits models to them.\n\n\
Input:\n\
@itemize \n\
@item @var{esp}: Array of extended solar cell parameters (Nx10)\n\
@item @var{models}:  One or more model specification strings (this option may be passed multiple times)\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{modelpar}: Cell array with model parameters (one column per model, one row per esp row)\n\
@item @var{err}: Fit error\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Matrix EPar, FErr;
	Cell MPar;
	ProcessSettings Set;
	CurveProp *P;
	curvepars *pars;
	int Nc=0, nm=0;
	int i, j, k;
  	if (nargin < 2)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Set=pv_craze_initsettings();
	EPar=args(0).matrix_value();
	
	for (i=1;i<nargin;i++)
	{
		charMatrix ch = args(i).char_matrix_value ();
		std::string tmp = ch.row_as_string (0);			
		pv_craze_addmodeltosettings(&Set, (char *)tmp.c_str());
		nm++;
	}
			
	if ((EPar.cols()!=10)||(EPar.rows()<1))
	{
		fprintf(stderr, "Error: pvcraze_modelfit requires a Nx10 extended solar cell parameter array\n");
		return octave_value_list ();
	}
	if (Set.L.N==0)
	{
		fprintf(stderr, "Error: pvcraze_modelfit requires at least one valid model specification\n");
		MPar=Cell(0, 0);
		FErr=Matrix(0,0);
	}
	else
	{
		Nc=EPar.rows();
		MPar=Cell(Nc, Set.L.N);
		FErr=Matrix(Nc, Set.L.N);
		if ((pars=(curvepars *)malloc(Nc*sizeof(curvepars)))==NULL)
		{
			pv_craze_freesettings(&Set);
			fprintf(stderr, "Error: pvcraze_modelfit failed to allocate memory\n");
			pv_craze_freesettings(&Set);
			return octave_value_list ();
		}
		for(i=0; i<Nc; i++)
		{
			pars[i].Voc  = EPar(i,0);
			pars[i].Isc  = EPar(i,1);
			pars[i].Vmpp = EPar(i,2);
			pars[i].Impp = EPar(i,3);
			pars[i].Vqmp1= EPar(i,4);
			pars[i].Iqmp1= EPar(i,5);
			pars[i].Vqmp2= EPar(i,6);
			pars[i].Iqmp2= EPar(i,7);
			pars[i].Gsc  = EPar(i,8);
			pars[i].Roc  = EPar(i,9);
		}
		P=pv_craze_modelfit(pars, Nc, &Set);
		
		// get data out to matrices
		for (i=0;i<Nc;i++)
		{
			for (k=0;k<Set.L.N;k++)
			{
				Matrix par, solpar;
				
				//par=Matrix(1,1+NumFields(Set.L.List[k].M));	
				par=Matrix(1,11);				
				par(0,0)=(double)Set.L.List[k].M;
				for (j=0;j<pv_craze_numfields(Set.L.List[k].M);j++)
					par(0,j+1)=P[i].P[k].pars[j];
				for (;j<10;j++)
					par(0,j+1)=0;
				MPar(i,k)=par;
				FErr(i,k)=P[i].E[k];				
			}
			pv_craze_freeprop(&P[i]);
		}
		free(pars);
		pv_craze_free(P);
	}
	pv_craze_freesettings(&Set);
	
	
	i=0;
	if (nargout>i) 
		retval(i++)=octave_value(MPar);
	if (nargout>i) 
		retval(i++)=octave_value(FErr);
		
	return retval;	
	
}

DEFUN_DLD (pvcraze_curvemodel, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{IV}] =} pvcraze_curvemodel (@var{modelpar}, N) \n\
pvcraze_curvemodel is part of the PV-CRAZE package and computes IV \
curves from PV-CRAZE generic parameters\n\n\
Input:\n\
@itemize \n\
@item @var{modelpar}: Cell array with model parametersEach cell must contain a 1x10 or 1x11 sized array. If an array is 1x11, it is treated as generic model parameters. If an 1x10 array is passed the PV-Bezier model is used.\n\
@item @var{N}: Positive Integer: nuber of points in the IV characteristics\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{IV}: Cell Array with IV characteristics:\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Cell IV, Pars;
	generic_pars *P;
	generic_pars NULPAR = {0,0,0,0,0,0,0,0,0,0, M_NONE};
	Matrix PP;
	int N, Nr,Nc;
	curve *C=NULL;
	int i, j;
  	if (nargin != 2)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pars=args(0).cell_value();
	N=args(1).int_value();
	if (N<=0)
	{
		fprintf(stderr, "Number of points in an IV curve cannot be smaller or equal to 0\n");
		return octave_value_list ();
	}
	
	Nr=Pars.rows();
	Nc=Pars.cols();
	if ((Nc==0)||(Nr==0))
	{
		retval(0)=Cell(Nr,Nc);
		return retval;
	}
		
	if ((P=(generic_pars *)malloc(Nc*Nr*sizeof(generic_pars)))==NULL)
	{
		return octave_value_list ();
	}
		
	for (i=0;i<Nc*Nr;i++)
	{
		PP=Pars(i/Nc,i%Nc).matrix_value();
		if (PP.rows()!=1)
			P[i]=NULPAR;
		else if(PP.cols()==11)
		{
			P[i].TYPE=(MODEL_TYPE)PP(0);
			P[i].pars[0]=PP(1);
			P[i].pars[1]=PP(2);
			P[i].pars[2]=PP(3);
			P[i].pars[3]=PP(4);
			P[i].pars[4]=PP(5);
			P[i].pars[5]=PP(6);
			P[i].pars[6]=PP(7);
			P[i].pars[7]=PP(8);
			P[i].pars[8]=PP(9);
			P[i].pars[9]=PP(10);
		}		
		else if(PP.cols()==10)
		{
			P[i].TYPE=M_PVBEZIER;
			P[i].pars[0]=PP(0);
			P[i].pars[1]=PP(1);
			P[i].pars[2]=PP(2);
			P[i].pars[3]=PP(3);
			P[i].pars[4]=PP(4);
			P[i].pars[5]=PP(5);
			P[i].pars[6]=PP(6);
			P[i].pars[7]=PP(7);
			P[i].pars[8]=PP(8);
			P[i].pars[9]=PP(9);
		}
		else
			P[i]=NULPAR;
	}
	
	C=pv_craze_curvemodel(P, Nc*Nr, N);
	free(P);
	
	if (C==NULL)
	{
		fprintf(stderr, "PV-CRAZE returned a NULL list\n");
		return octave_value_list ();	
	}
	
	IV=Cell(Nr, Nc);
	for (i=0;i<Nr*Nc;i++)
	{
		Matrix iv;
		iv=Matrix(C[i].N,2);
		for (j=0;j<C[i].N;j++)
		{
			iv(j,0)=C[i].P[j].x;
			iv(j,1)=C[i].P[j].y;
		}
		IV(i/Nc,i%Nc)=iv;
		if (C[i].P)
			free(C[i].P);
	}
	pv_craze_free(C);

	retval(0)=octave_value(IV);
	return retval;	
	
}


DEFUN_DLD (pvcraze_modelpars, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{esp}] =} pvcraze_modelpars (@var{modelpar})\n\
pvcraze_modelpars is part of the PV-CRAZE package and computes extended solar cell parameters from \
PV-CRAZE generic parameters\n\n\
Input:\n\
@itemize \n\
@item @var{modelpar}: Cell array with model parameters. Each cell must contain a 1x10 or 1x11 sized array. If an array is 1x11, it is treated as generic model parameters. If an 1x10 array is passed the PV-Bezier model is used.\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{esp}: Cell Array with extended solar cell parameters:\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Cell Pars, MPars;
	generic_pars P;
	generic_pars NULPAR = {0,0,0,0,0,0,0,0,0,0, M_NONE};
	Matrix PP;
	int Nr,Nc;
	curvepars solpars;
	int i;
  	if (nargin != 1)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pars=args(0).cell_value();
	
	Nr=Pars.rows();
	Nc=Pars.cols();
	if ((Nc==0)||(Nr==0))
	{
		retval(0)=Cell(Nr,Nc);
		return retval;
	}
	MPars=Cell(Nr, Nc);
		
	for (i=0;i<Nc*Nr;i++)
	{
		Matrix solpar;
		PP=Pars(i/Nc,i%Nc).matrix_value();
		if (PP.rows()!=1)
			P=NULPAR;
		else if(PP.cols()==11)
		{
			P.TYPE=(MODEL_TYPE)PP(0);
			P.pars[0]=PP(1);
			P.pars[1]=PP(2);
			P.pars[2]=PP(3);
			P.pars[3]=PP(4);
			P.pars[4]=PP(5);
			P.pars[5]=PP(6);
			P.pars[6]=PP(7);
			P.pars[7]=PP(8);
			P.pars[8]=PP(9);
			P.pars[9]=PP(10);
		}		
		else if(PP.cols()==10)
		{
			P.TYPE=M_PVBEZIER;
			P.pars[0]=PP(0);
			P.pars[1]=PP(1);
			P.pars[2]=PP(2);
			P.pars[3]=PP(3);
			P.pars[4]=PP(4);
			P.pars[5]=PP(5);
			P.pars[6]=PP(6);
			P.pars[7]=PP(7);
			P.pars[8]=PP(8);
			P.pars[9]=PP(9);
		}
		else
			P=NULPAR;
		solpars=pv_craze_g_pars(P);
		solpar=Matrix(1,10);
		solpar(0,0)=solpars.Voc;
		solpar(0,1)=solpars.Isc;
		solpar(0,2)=solpars.Vmpp;
		solpar(0,3)=solpars.Impp;
		solpar(0,4)=solpars.Vqmp1;
		solpar(0,5)=solpars.Iqmp1;
		solpar(0,6)=solpars.Vqmp2;
		solpar(0,7)=solpars.Iqmp2;
		solpar(0,8)=solpars.Gsc;
		solpar(0,9)=solpars.Roc;
		MPars(i/Nc,i%Nc)=solpar;
	}
	retval(0)=octave_value(MPars);
	return retval;
}

DEFUN_DLD (pvcraze_modelerror, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{err}] =} pvcraze_modelerror (@var{modelpar}, @var{IV})\n\
pvcraze_modelerror is part of the PV-CRAZE package and computes a RMS distance IV curve and a model given and IV's\n\
PV-CRAZE generic parameters\n\n\
Input:\n\
@itemize \n\
@item @var{modelpar}: Cell array with model parameters. Each cell must contain a 1x10 or 1x11 sized array. If an array is 1x11, it is treated as generic model parameters. If an 1x10 array is passed the PV-Bezier model is used.\n\
@item @var{IV}: Cell array with IV curves corresponding to the modelpar variable\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{err}: Array with RMS distances for each IV an model parameter set\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Cell IV, Pars;
	generic_pars *P;
	generic_pars NULPAR = {0,0,0,0,0,0,0,0,0,0, M_NONE};
	Matrix PP, D;
	int Nr,Nc;
	double *DD;
	curve *C=NULL;
	int i, j;
  	if (nargin != 2)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pars=args(0).cell_value();
	IV=args(1).cell_value();
	
	Nr=Pars.rows();
	Nc=Pars.cols();
	
	if ((Nr!=IV.rows())||(Nc!=IV.cols()))	
	{
		fprintf(stderr, "Size of the modelpar cell array must be identical to the IV cell array\n");
		return octave_value_list ();
	}
	if ((Nc==0)||(Nr==0))
	{
		retval(0)=Cell(Nr,Nc);
		return retval;
	}
	

	// create a list of curves
	if ((C=(curve *)malloc(Nc*Nr*sizeof(curve)))==NULL)
	{
		fprintf(stderr, "Could not allocate curve list\n");
		return octave_value_list ();
	}
		
	if ((P=(generic_pars *)malloc(Nc*Nr*sizeof(generic_pars)))==NULL)
	{
		free(C);
		fprintf(stderr, "Could not allocate parameters list\n");
		return octave_value_list ();
	}
	
	for (i=0;i<Nc*Nr;i++)
	{
		Matrix iv;
		iv=IV(i%Nr,i/Nr).matrix_value();
		if ((iv.cols()!=2)||(iv.rows()==0))
		{
			C[i].P=NULL;
			C[i].N=0;
		}
		else
		{
			C[i].P=(point*) malloc((int)iv.rows()*sizeof(point));
			C[i].N=iv.rows();
			for(j=0;j<C[i].N;j++)
			{
				C[i].P[j].x=iv(j,0);
				C[i].P[j].y=iv(j,1);
			}
		}
		PP=Pars(i/Nc,i%Nc).matrix_value();
		if (PP.rows()!=1)
			P[i]=NULPAR;
		else if(PP.cols()==11)
		{
			P[i].TYPE=(MODEL_TYPE)PP(0);
			P[i].pars[0]=PP(1);
			P[i].pars[1]=PP(2);
			P[i].pars[2]=PP(3);
			P[i].pars[3]=PP(4);
			P[i].pars[4]=PP(5);
			P[i].pars[5]=PP(6);
			P[i].pars[6]=PP(7);
			P[i].pars[7]=PP(8);
			P[i].pars[8]=PP(9);
			P[i].pars[9]=PP(10);
		}		
		else if(PP.cols()==10)
		{
			P[i].TYPE=M_PVBEZIER;
			P[i].pars[0]=PP(0);
			P[i].pars[1]=PP(1);
			P[i].pars[2]=PP(2);
			P[i].pars[3]=PP(3);
			P[i].pars[4]=PP(4);
			P[i].pars[5]=PP(5);
			P[i].pars[6]=PP(6);
			P[i].pars[7]=PP(7);
			P[i].pars[8]=PP(8);
			P[i].pars[9]=PP(9);
		}
		else
			P[i]=NULPAR;
	}	
	DD=pv_craze_error(P, C, Nc*Nr);
	for (i=0;i<Nc*Nr;i++)
		free(C[i].P);
	free(C);
	free(P);
	
	if (!DD)
	{
		fprintf(stderr, "PV-CRAZE returned a NULL list\n");
		return octave_value_list ();			
	}
	D=Matrix(Nr, Nc);
	for (i=0;i<Nr*Nc;i++)
	{
		D(i%Nr,i/Nr)=DD[i];
	}
	pv_craze_free(DD);

	retval(0)=octave_value(D);
	return retval;	
	
}


DEFUN_DLD (pvcraze_errorcurve, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{errcurve}] =} pvcraze_errorcurve (@var{modelpar}, @var{IV})\n\
pvcraze_errorcurve is part of the PV-CRAZE package and computes a distance curve of an IV curve and a model given and IV's \
model parameters\n\n\
Input:\n\
@itemize \n\
@item @var{modelpar}: Cell array with model parameters. Each cell must contain a 1x10 or 1x11 sized array. If an array is 1x11, it is treated as generic model parameters. If an 1x10 array is passed the PV-Bezier model is used.\n\
@item @var{IV}: Cell array with IV curves corresponding to the modelpar variable\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{errcurve}: Cell Array with distance curves (cartesian distance in a normalized current-voltage plane between curve and model w.r.t. cell voltage\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Cell IV, Pars, IVo;
	generic_pars *P;
	generic_pars NULPAR = {0,0,0,0,0,0,0,0,0,0, M_NONE};
	Matrix PP;
	int Nr,Nc;
	curve *C=NULL, *Co=NULL;
	int i, j;
  	if (nargin != 2)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pars=args(0).cell_value();
	IV=args(1).cell_value();
	
	Nr=Pars.rows();
	Nc=Pars.cols();
	
	if ((Nr!=IV.rows())||(Nc!=IV.cols()))	
	{
		fprintf(stderr, "Size of the modelpar cell array must be identical to the IV cell array\n");
		return octave_value_list ();
	}
	if ((Nc==0)||(Nr==0))
	{
		retval(0)=Cell(Nr,Nc);
		return retval;
	}
	

	// create a list of curves
	if ((C=(curve *)malloc(Nc*Nr*sizeof(curve)))==NULL)
	{
		fprintf(stderr, "Could not allocate curve list\n");
		return octave_value_list ();
	}
		
	if ((P=(generic_pars *)malloc(Nc*Nr*sizeof(generic_pars)))==NULL)
	{
		free(C);
		fprintf(stderr, "Could not allocate parameters list\n");
		return octave_value_list ();
	}
	
	
	for (i=0;i<Nc*Nr;i++)
	{
		Matrix iv;
		iv=IV(i%Nr,i/Nr).matrix_value();
		if ((iv.cols()!=2)||(iv.rows()==0))
		{
			C[i].P=NULL;
			C[i].N=0;
		}
		else
		{
			C[i].P=(point*) malloc((int)iv.rows()*sizeof(point));
			C[i].N=iv.rows();
			for(j=0;j<C[i].N;j++)
			{
				C[i].P[j].x=iv(j,0);
				C[i].P[j].y=iv(j,1);
			}
		}
		PP=Pars(i%Nr,i/Nr).matrix_value();
		if (PP.rows()!=1)
			P[i]=NULPAR;
		else if(PP.cols()==11)
		{
			P[i].TYPE=(MODEL_TYPE)PP(0);
			P[i].pars[0]=PP(1);
			P[i].pars[1]=PP(2);
			P[i].pars[2]=PP(3);
			P[i].pars[3]=PP(4);
			P[i].pars[4]=PP(5);
			P[i].pars[5]=PP(6);
			P[i].pars[6]=PP(7);
			P[i].pars[7]=PP(8);
			P[i].pars[8]=PP(9);
			P[i].pars[9]=PP(10);
		}		
		else if(PP.cols()==10)
		{
			P[i].TYPE=M_PVBEZIER;
			P[i].pars[0]=PP(0);
			P[i].pars[1]=PP(1);
			P[i].pars[2]=PP(2);
			P[i].pars[3]=PP(3);
			P[i].pars[4]=PP(4);
			P[i].pars[5]=PP(5);
			P[i].pars[6]=PP(6);
			P[i].pars[7]=PP(7);
			P[i].pars[8]=PP(8);
			P[i].pars[9]=PP(9);
		}
		else
			P[i]=NULPAR;
	}
	
	
	
	Co=pv_craze_errorcurve(P, C, Nc*Nr);
	free(P);
	
	if (!Co)
	{
		free(C);
		fprintf(stderr, "PV-CRAZE returned a NULL list\n");
		return octave_value_list ();			
	}
		
	IVo=Cell(Nr, Nc);
	for (i=0;i<Nr*Nc;i++)
	{
		Matrix iv;
		iv=Matrix(Co[i].N,2);
		for (j=0;j<Co[i].N;j++)
		{
			iv(j,0)=Co[i].P[j].x;
			iv(j,1)=Co[i].P[j].y;
		}
		IVo(i%Nr,i/Nr)=iv;
		if (Co[i].P)
			free(Co[i].P);
		free(C[i].P);
	}
	free(C);
	pv_craze_free(Co);

	retval(0)=octave_value(IVo);
	return retval;	
	
}

DEFUN_DLD (pvcraze_iv, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{I}] =} pvcraze_iv (@var{V},@var{modelpar})\n\
pvcraze_iv is part of the PV-CRAZE package and computes c current given a voltage for\n\
PV-CRAZE generic parameters\n\n\
Input:\n\
@itemize \n\
@item @var{V}: Array of voltages\n\
@item @var{modelpar}: generic model parameters (sized 1x10 or 1x11). If an array is 1x11, it is treated as generic model parameters. If an 1x10 array is passed the PV-Bezier model is used.\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{I}: Array with currents corresponing to the given voltages\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Matrix I, V, Pars;
	generic_pars P;
	int Nr,Nc;
	int i, j;
  	if (nargin != 2)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pars=args(1).matrix_value();
	V=args(0).matrix_value();
	
	Nr=V.rows();
	Nc=V.cols();
	
	if ((Pars.rows()!=1)||((Pars.cols()!=10)&&(Pars.cols()!=11)))	
	{
		fprintf(stderr, "Size of the modelpar array must be sized 1x10 or 1x11\n");
		return octave_value_list ();
	}
	if ((Nc==0)||(Nr==0))
	{
		retval(0)=Matrix(Nr,Nc);
		return retval;
	}
	if (Pars.cols()==11)
	{
		P.TYPE=(MODEL_TYPE)Pars(0);
		P.pars[0]=Pars(1);
		P.pars[1]=Pars(2);
		P.pars[2]=Pars(3);
		P.pars[3]=Pars(4);
		P.pars[4]=Pars(5);
		P.pars[5]=Pars(6);
		P.pars[6]=Pars(7);
		P.pars[7]=Pars(8);
		P.pars[8]=Pars(9);
		P.pars[9]=Pars(10);
	}
	else
	{
		P.TYPE=M_PVBEZIER;
		P.pars[0]=Pars(0);
		P.pars[1]=Pars(1);
		P.pars[2]=Pars(2);
		P.pars[3]=Pars(3);
		P.pars[4]=Pars(4);
		P.pars[5]=Pars(5);
		P.pars[6]=Pars(6);
		P.pars[7]=Pars(7);
		P.pars[8]=Pars(8);
		P.pars[9]=Pars(9);
	}
	
	I=Matrix(Nr,Nc);
	
	for (i=0;i<Nr;i++)
		for (j=0;j<Nc;j++)
			I(i,j)=pv_craze_IV(V(i,j), P);
	
	retval(0)=octave_value(I);
	return retval;	
	
}

DEFUN_DLD (pvcraze_vi, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{V}] =} pvcraze_iv (@var{I},@var{modelpar})\n\
pvcraze_vi is part of the PV-CRAZE package and computes a voltage given a current for\n\
PV-CRAZE generic parameters\n\n\
Input:\n\
@itemize \n\
@item @var{I}: Array of currents\n\
@item @var{modelpar}: generic model parameters (sized 1x10 or 1x11). If an array is 1x11, it is treated as generic model parameters. If an 1x10 array is passed the PV-Bezier model is used.\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{V}: Array with voltages corresponing to the given currents\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	Matrix I, V, Pars;
	generic_pars P;
	int Nr,Nc;
	int i, j;
  	if (nargin != 2)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pars=args(1).matrix_value();
	I=args(0).matrix_value();
	
	Nr=I.rows();
	Nc=I.cols();
	
	if ((Pars.rows()!=1)||((Pars.cols()!=10)&&(Pars.cols()!=11)))
	{
		fprintf(stderr, "Size of the modelpar array must be sized 1x10 or 1x11\n");
		return octave_value_list ();
	}
	if ((Nc==0)||(Nr==0))
	{
		retval(0)=Matrix(Nr,Nc);
		return retval;
	}
	if (Pars.cols()==11)
	{
		P.TYPE=(MODEL_TYPE)Pars(0);
		P.pars[0]=Pars(1);
		P.pars[1]=Pars(2);
		P.pars[2]=Pars(3);
		P.pars[3]=Pars(4);
		P.pars[4]=Pars(5);
		P.pars[5]=Pars(6);
		P.pars[6]=Pars(7);
		P.pars[7]=Pars(8);
		P.pars[8]=Pars(9);
		P.pars[9]=Pars(10);
	}
	else
	{
		P.TYPE=M_PVBEZIER;
		P.pars[0]=Pars(0);
		P.pars[1]=Pars(1);
		P.pars[2]=Pars(2);
		P.pars[3]=Pars(3);
		P.pars[4]=Pars(4);
		P.pars[5]=Pars(5);
		P.pars[6]=Pars(6);
		P.pars[7]=Pars(7);
		P.pars[8]=Pars(8);
		P.pars[9]=Pars(9);
	}
	
	V=Matrix(Nr,Nc);
	
	for (i=0;i<Nr;i++)
		for (j=0;j<Nc;j++)
			V(i,j)=pv_craze_VI(I(i,j), P);
	
	retval(0)=octave_value(V);
	return retval;	
	
}

DEFUN_DLD (pvcraze_fit_IEC60891_v2, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{esp_std}, alpha, beta, a, Rs, kappa] =} pvcraze_fit_IEC60891_v2 (@var{esp}, T, G, Tstd, Gstd [,furthest_N[,offset]])\n\
pvcraze_fit_IEC60891_v2 is part of the PV-CRAZE package and fits IEC60891 translation\n\
coefficients and determines standard test condition curve parameters to best describe \n\
provided temperature and irradiance dependent curve parameters. Uses the IEC60891 version 2 (2009)\n\n\
Input:\n\
@itemize \n\
@item @var{esp}: array with extended solar cell parameters, one row per conditions\n\
@item @var{T}: Array with temperatures, one temperature per row in esp\n\
@item @var{G}: Array with irradiances, one irradiance per row in esp\n\
@item @var{Tstd}: Scalar value for the standard test conditions Temperature\n\
@item @var{Gstd}: Scalar value for the standard test conditions Irradiance\n\
@item @var{furthest_N}: Number of samples to select using furthest sampling\n\
@item @var{offset}: If furthest sampling is used, the offset to use (changes which samples are chosen)\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{esp_std}: Parameters at standard test conditions\n\
@item @var{alpha}:     IEC60891 alpha\n\
@item @var{beta}:      IEC60891 beta\n\
@item @var{a}:         IEC60891 a\n\
@item @var{Rs}:        IEC60891 Rs\n\
@item @var{kappa}:     IEC60891 kappa\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	curvepars *Cin;
	curvepars Cout;
	Matrix T,G,Pin;
	RowVector Pout;
	double *TT, *GG;
	double alpha, beta, a, Rs, kappa, Tstd, Gstd;
	int N, Nr,Nc;
	int FN=0,off=0;
	int i;
  	if ((nargin < 5)||(nargin > 7))
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pin=args(0).matrix_value();
	T=args(1).matrix_value();
	G=args(2).matrix_value();
	Tstd=args(3).double_value();
	Gstd=args(4).double_value();
	
	if((nargin==6)||(nargin==7))
		FN=args(5).int_value();
	if(nargin==7)
		off=args(6).int_value();
		
	Nr=Pin.rows();
	Nc=Pin.cols();
	if ((Nr<3)||(Nc!=10))
	{
		fprintf(stderr, "Number of rows in esp may not less than 3\n");
		fprintf(stderr, "and number of columns must be 10\n");
		return octave_value_list ();
	}
	N=T.rows();
	Nc=T.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in T must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	N=G.rows();
	Nc=G.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in G must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	
	
	// create a list of esp
	if ((Cin=(curvepars *)malloc(N*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	
	for (i=0;i<N;i++)
	{
		Cin[i].Voc=Pin(i,0);
		Cin[i].Isc=Pin(i,1);
		Cin[i].Vmpp=Pin(i,2);
		Cin[i].Impp=Pin(i,3);
		Cin[i].Vqmp1=Pin(i,4);
		Cin[i].Iqmp1=Pin(i,5);
		Cin[i].Vqmp2=Pin(i,6);
		Cin[i].Iqmp2=Pin(i,7);
		Cin[i].Gsc=Pin(i,8);
		Cin[i].Roc=Pin(i,9);	
		TT[i]=T(i,0);	
		GG[i]=G(i,0);			
	}
	Cout=pv_craze_fit_IEC60891_v2(Cin, GG, TT, N, Tstd,Gstd, &alpha, &beta, &a, &Rs, &kappa,FN,off);
	
	free(GG);
	free(TT);
	free(Cin);
	
	Pout=RowVector(10);
	Pout(0)=Cout.Voc;
	Pout(1)=Cout.Isc;
	Pout(2)=Cout.Vmpp;
	Pout(3)=Cout.Impp;
	Pout(4)=Cout.Vqmp1;
	Pout(5)=Cout.Iqmp1;
	Pout(6)=Cout.Vqmp2;
	Pout(7)=Cout.Iqmp2;
	Pout(8)=Cout.Gsc;
	Pout(9)=Cout.Roc;	

	i=0;
	if (nargout>i)
		retval(i++)=octave_value(Pout);
	if (nargout>i)
		retval(i++)=octave_value(alpha);
	if (nargout>i) 
		retval(i++)=octave_value(beta);
	if (nargout>i) 
		retval(i++)=octave_value(a);
	if (nargout>i) 
		retval(i++)=octave_value(Rs);
	if (nargout>i) 
		retval(i++)=octave_value(kappa);
	return retval;	
}


DEFUN_DLD (pvcraze_fit_IEC60891_v3, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{esp_std}, alpha, beta, B1, B2, Rs, kappa] =} pvcraze_fit_IEC60891_v3 (@var{esp}, T, G, Tstd, Gstd [,furthest_N[,offset]])\n\
pvcraze_fit_IEC60891_v3 is part of the PV-CRAZE package and fits IEC60891 translation\n\
coefficients and determines standard test condition curve parameters to best describe \n\
provided temperature and irradiance dependent curve parameters Uses the IEC60891 version 3 (2021)\n\n\
Input:\n\
@itemize \n\
@item @var{esp}: array with extended solar cell parameters, one row per conditions\n\
@item @var{T}: Array with temperatures, one temperature per row in esp\n\
@item @var{G}: Array with irradiances, one irradiance per row in esp\n\
@item @var{Tstd}: Scalar value for the standard test conditions Temperature\n\
@item @var{Gstd}: Scalar value for the standard test conditions Irradiance\n\
@item @var{furthest_N}: Number of samples to select using furthest sampling\n\
@item @var{offset}: If furthest sampling is used, the offset to use (changes which samples are chosen)\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{esp_std}: Parameters at standard test conditions\n\
@item @var{alpha}:     IEC60891 alpha\n\
@item @var{beta}:      IEC60891 beta\n\
@item @var{B1}:        IEC60891 B1\n\
@item @var{B2}:        IEC60891 B2\n\
@item @var{Rs}:        IEC60891 Rs\n\
@item @var{kappa}:     IEC60891 kappa\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	curvepars *Cin;
	curvepars Cout;
	Matrix T,G,Pin;
	RowVector Pout;
	double *TT, *GG;
	double alpha, beta, B1, B2, Rs, kappa, Tstd, Gstd;
	int N, Nr,Nc;
	int FN=0, off=0;
	int i;
  	if ((nargin < 5)||(nargin > 7))
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pin=args(0).matrix_value();
	T=args(1).matrix_value();
	G=args(2).matrix_value();
	Tstd=args(3).double_value();
	Gstd=args(4).double_value();
	
	if((nargin==6)||(nargin==7))
		FN=args(5).int_value();
	if(nargin==7)
		off=args(6).int_value();
		
	Nr=Pin.rows();
	Nc=Pin.cols();
	if ((Nr<3)||(Nc!=10))
	{
		fprintf(stderr, "Number of rows in esp may not less than 3\n");
		fprintf(stderr, "and number of columns must be 10\n");
		return octave_value_list ();
	}
	N=T.rows();
	Nc=T.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in T must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	N=G.rows();
	Nc=G.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in G must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	
	
	// create a list of esp	
	if ((Cin=(curvepars *)malloc(N*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	
	for (i=0;i<N;i++)
	{
		Cin[i].Voc=Pin(i,0);
		Cin[i].Isc=Pin(i,1);
		Cin[i].Vmpp=Pin(i,2);
		Cin[i].Impp=Pin(i,3);
		Cin[i].Vqmp1=Pin(i,4);
		Cin[i].Iqmp1=Pin(i,5);
		Cin[i].Vqmp2=Pin(i,6);
		Cin[i].Iqmp2=Pin(i,7);
		Cin[i].Gsc=Pin(i,8);
		Cin[i].Roc=Pin(i,9);	
		TT[i]=T(i,0);	
		GG[i]=G(i,0);
	
	}
	
	Cout=pv_craze_fit_IEC60891_v3(Cin, GG, TT, N, Tstd,Gstd, &alpha, &beta, &B1, &B2, &Rs, &kappa, FN, off);
	
	free(GG);
	free(TT);
	free(Cin);
	
	Pout=RowVector(10);
	Pout(0)=Cout.Voc;
	Pout(1)=Cout.Isc;
	Pout(2)=Cout.Vmpp;
	Pout(3)=Cout.Impp;
	Pout(4)=Cout.Vqmp1;
	Pout(5)=Cout.Iqmp1;
	Pout(6)=Cout.Vqmp2;
	Pout(7)=Cout.Iqmp2;
	Pout(8)=Cout.Gsc;
	Pout(9)=Cout.Roc;	

	i=0;
	if (nargout>i)
		retval(i++)=octave_value(Pout);
	if (nargout>i)
		retval(i++)=octave_value(alpha);
	if (nargout>i) 
		retval(i++)=octave_value(beta);
	if (nargout>i) 
		retval(i++)=octave_value(B1);
	if (nargout>i) 
		retval(i++)=octave_value(B2);
	if (nargout>i) 
		retval(i++)=octave_value(Rs);
	if (nargout>i) 
		retval(i++)=octave_value(kappa);
	return retval;	
}


DEFUN_DLD (pvcraze_translate_IEC60891_v2, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{esp}] =} pvcraze_translate_IEC60891_v2 (@var{esp}, T1, G1, T2, G2, alpha, beta, a, Rs, kappa)\n\
pvcraze_translate_IEC60891_v2 is part of the PV-CRAZE package and translates curve parameters from \n\
conditions G1, T1 to G2, T2, using procedure 2 from IEC60891 norm version 2 (2009)\n\n\
Input:\n\
@itemize \n\
@item @var{esp}: array with extended solar cell parameters, one row per conditions\n\
@item @var{T1}: Array with temperatures, one temperature per row in esp\n\
@item @var{G1}: Array with irradiances, one irradiance per row in esp\n\
@item @var{T2}: Array with output temperatures\n\
@item @var{G2}: Array with output irradiances\n\
@item @var{alpha}: IEC60891 alpha\n\
@item @var{beta}:  IEC60891 beta\n\
@item @var{a}:     IEC60891 a\n\
@item @var{Rs}:    IEC60891 Rs\n\
@item @var{kappa}: IEC60891 kappa\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{esp}: Array with translated extended parameters. The length of the array is the maximum of the length of the input (T1/G1 or G2/T2)\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	curvepars *Cin;
	curvepars *Cout;
	Matrix T1,G1,T2,G2,Pin,Pout;
	double *TT1, *GG1, *TT2, *GG2;
	double alpha, beta, a, Rs, kappa;
	int N, Nr,Nc, M, Mr, Mc;
	int i;
  	if (nargin != 10)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pin=args(0).matrix_value();
	T1=args(1).matrix_value();
	G1=args(2).matrix_value();
	T2=args(3).matrix_value();
	G2=args(4).matrix_value();
	alpha=args(5).double_value();
	beta=args(6).double_value();
	a=args(7).double_value();
	Rs=args(8).double_value();
	kappa=args(9).double_value();
	
	Nr=Pin.rows();
	Nc=Pin.cols();
	if ((Nr==0)||(Nc!=10))
	{
		fprintf(stderr, "Number of rows in esp may not be 0\n");
		fprintf(stderr, "and number of columns must be 10\n");
		return octave_value_list ();
	}
	N=T1.rows();
	Nc=T1.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in T1 must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	N=G1.rows();
	Nc=G1.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in G1 must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	Mr=T2.rows();
	Mc=T2.cols();
	if ((Mr==0)||(Mc!=1))
	{
		fprintf(stderr, "Number of rows in G2 must belarger than 0\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	M=G2.rows();
	Mc=G2.cols();
	if ((M!=Mr)||(Mc!=1))
	{
		fprintf(stderr, "Number of rows in G2 must be equal to rows in T2\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	
	
	// create a list of esp
	if ((Cin=(curvepars *)malloc(N*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT1=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG1=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT1);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT2=(double *)malloc(M*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT1);
		free(GG1);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG2=(double *)malloc(M*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT1);
		free(GG1);
		free(TT2);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	
	for (i=0;i<N;i++)
	{
		Cin[i].Voc=Pin(i,0);
		Cin[i].Isc=Pin(i,1);
		Cin[i].Vmpp=Pin(i,2);
		Cin[i].Impp=Pin(i,3);
		Cin[i].Vqmp1=Pin(i,4);
		Cin[i].Iqmp1=Pin(i,5);
		Cin[i].Vqmp2=Pin(i,6);
		Cin[i].Iqmp2=Pin(i,7);
		Cin[i].Gsc=Pin(i,8);
		Cin[i].Roc=Pin(i,9);	
		TT1[i]=T1(i,0);	
		GG1[i]=G1(i,0);			
	}
	for (i=0;i<M;i++)
	{
		TT2[i]=T2(i,0);	
		GG2[i]=G2(i,0);		
	}
	Cout=pv_craze_translate_IEC60891_v2(Cin, GG1,TT1, N, GG2, TT2, M, alpha, beta, a, Rs, kappa);
	free(GG1);
	free(GG2);
	free(TT1);
	free(TT2);
	free(Cin);
	
	if (!Cout)
	{
		fprintf(stderr, "PV-CRAZE returned a NULL list\n");
		return octave_value_list ();		
	}
	
	if (M>N)
		N=M;
	Pout=Matrix(N,10);
	for (i=0;i<N;i++)
	{
		Pout(i,0)=Cout[i].Voc;
		Pout(i,1)=Cout[i].Isc;
		Pout(i,2)=Cout[i].Vmpp;
		Pout(i,3)=Cout[i].Impp;
		Pout(i,4)=Cout[i].Vqmp1;
		Pout(i,5)=Cout[i].Iqmp1;
		Pout(i,6)=Cout[i].Vqmp2;
		Pout(i,7)=Cout[i].Iqmp2;
		Pout(i,8)=Cout[i].Gsc;
		Pout(i,9)=Cout[i].Roc;	
	}	
	pv_craze_free(Cout);	

	retval(0)=octave_value(Pout);
	return retval;	
	
}

DEFUN_DLD (pvcraze_translate_IEC60891_v3, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{esp}] =} pvcraze_translate_IEC60891_v3 (@var{esp}, T1, G1, T2, G2, Vocstd, Tstd, Gstd, alpha, beta, B1, B2, Rs, kappa)\n\
pvcraze_translate_IEC60891_v3 is part of the PV-CRAZE package and translates curve parameters from \n\
conditions G1, T1 to G2, T2, using procedure 2 from IEC60891 norm version 3 (2021)\n\n\
Input:\n\
@itemize \n\
@item @var{esp}: array with extended solar cell parameters, one row per conditions\n\
@item @var{T1}: Array with temperatures, one temperature per row in esp\n\
@item @var{G1}: Array with irradiances, one irradiance per row in esp\n\
@item @var{T2}: Array with output temperatures\n\
@item @var{G2}: Array with output irradiances\n\
@item @var{Vocstd}: Open-Circuit voltage at standard test conditions\n\
@item @var{Tstd}:   Standard test conditions temperature\n\
@item @var{Gstd}:   Standard test conditions irradiance\n\
@item @var{alpha}:  IEC60891 alpha\n\
@item @var{beta}:   IEC60891 beta\n\
@item @var{B1}:     IEC60891 B1\n\
@item @var{B2}:     IEC60891 B2\n\
@item @var{Rs}:     IEC60891 Rs\n\
@item @var{kappa}:  IEC60891 kappa\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{esp}: Array with translated extended parameters. The length of the array is the maximum of the length of the input (T1/G1 or G2/T2)\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	curvepars *Cin;
	curvepars *Cout;
	Matrix T1,G1,T2,G2,Pin,Pout;
	double *TT1, *GG1, *TT2, *GG2;
	double alpha, beta, B1, B2, Rs, kappa, Vocstd, Gstd, Tstd;
	int N, Nr,Nc, M, Mr, Mc;
	int i;
  	if (nargin != 14)
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pin=args(0).matrix_value();
	T1=args(1).matrix_value();
	G1=args(2).matrix_value();
	T2=args(3).matrix_value();
	G2=args(4).matrix_value();
	Vocstd=args(5).double_value();
	Tstd=args(6).double_value();
	Gstd=args(7).double_value();
	alpha=args(8).double_value();
	beta=args(9).double_value();
	B1=args(10).double_value();
	B2=args(11).double_value();
	Rs=args(12).double_value();
	kappa=args(13).double_value();
	
		
	Nr=Pin.rows();
	Nc=Pin.cols();
	if ((Nr==0)||(Nc!=10))
	{
		fprintf(stderr, "Number of rows in esp may not be 0\n");
		fprintf(stderr, "and number of columns must be 10\n");
		return octave_value_list ();
	}
	N=T1.rows();
	Nc=T1.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in T1 must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	N=G1.rows();
	Nc=G1.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in G1 must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	Mr=T2.rows();
	Mc=T2.cols();
	if ((Mr==0)||(Mc!=1))
	{
		fprintf(stderr, "Number of rows in G2 must belarger than 0\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	M=G2.rows();
	Mc=G2.cols();
	if ((M!=Mr)||(Mc!=1))
	{
		fprintf(stderr, "Number of rows in G2 must be equal to rows in T2\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	
	
	// create a list of esp
	if ((Cin=(curvepars *)malloc(N*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT1=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG1=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT1);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT2=(double *)malloc(M*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT1);
		free(GG1);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG2=(double *)malloc(M*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT1);
		free(GG1);
		free(TT2);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	
	for (i=0;i<N;i++)
	{
		Cin[i].Voc=Pin(i,0);
		Cin[i].Isc=Pin(i,1);
		Cin[i].Vmpp=Pin(i,2);
		Cin[i].Impp=Pin(i,3);
		Cin[i].Vqmp1=Pin(i,4);
		Cin[i].Iqmp1=Pin(i,5);
		Cin[i].Vqmp2=Pin(i,6);
		Cin[i].Iqmp2=Pin(i,7);
		Cin[i].Gsc=Pin(i,8);
		Cin[i].Roc=Pin(i,9);	
		TT1[i]=T1(i,0);	
		GG1[i]=G1(i,0);			
	}
	for (i=0;i<M;i++)
	{
		TT2[i]=T2(i,0);	
		GG2[i]=G2(i,0);		
	}
	Cout=pv_craze_translate_IEC60891_v3(Cin, GG1,TT1, N, GG2, TT2, M, Vocstd, Tstd, Gstd, alpha, beta, B1, B2, Rs, kappa);
	free(GG1);
	free(GG2);
	free(TT1);
	free(TT2);
	free(Cin);
	
	if (!Cout)
	{
		fprintf(stderr, "PV-CRAZE returned a NULL list\n");
		return octave_value_list ();		
	}
	
	if (M>N)
		N=M;
	Pout=Matrix(N,10);
	for (i=0;i<N;i++)
	{
		Pout(i,0)=Cout[i].Voc;
		Pout(i,1)=Cout[i].Isc;
		Pout(i,2)=Cout[i].Vmpp;
		Pout(i,3)=Cout[i].Impp;
		Pout(i,4)=Cout[i].Vqmp1;
		Pout(i,5)=Cout[i].Iqmp1;
		Pout(i,6)=Cout[i].Vqmp2;
		Pout(i,7)=Cout[i].Iqmp2;
		Pout(i,8)=Cout[i].Gsc;
		Pout(i,9)=Cout[i].Roc;	
	}	
	pv_craze_free(Cout);	

	retval(0)=octave_value(Pout);
	return retval;	
	
}


DEFUN_DLD (pvcraze_fit_diode_GT, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{modelpar_std}, alpha, Ea] =} pvcraze_fit_diode_GT (@var{esp}, T, G, Tstd, Gstd [,unit][,furthest_N [,offset]])\n\
pvcraze_fit_diode_GT is part of the PV-CRAZE package and fits a diode model \
to temperature and irradiance dependent extended solar cell parameters \
It returns diode model parameters at some specified standard temperatre \
and irradiance conditions. Furthermore, it returns two coefficients, alpha, \
the temperature coefficient of the short circuit current, and Ea, the \
activation energy of the saturation current.\n\n\
Input:\n\
@itemize \n\
@item @var{esp}: array with extended solar cell parameters, one row per conditions\n\
@item @var{T}: Array with temperatures, one temperature per row in esp\n\
@item @var{G}: Array with irradiances, one irradiance per row in esp\n\
@item @var{Tstd}: Scalar value for the standard test conditions Temperature\n\
@item @var{Gstd}: Scalar value for the standard test conditions Irradiance\n\
@item @var{unit}: Temperature unit (\"C\"/\"Celcius\" or \"K\"/\"Kelvin\")\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{modelpar_std}: Diode parameters at standard test conditions\n\
@item @var{alpha}:     Temperature coefficient of the photo-current.\n\
@item @var{Ea}:        Activation energy of the saturation current\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	curvepars *Cin;
	generic_pars Pstd;
	Matrix T,G,Pin;
	RowVector Pout;
	double *TT, *GG;
	double alpha, Ea, Tstd, Gstd;
	int N, Nr,Nc;
	int i;
	char unit ='K';
	int FN=0,off=0,ud=0;
  	if ((nargin < 5)&&(nargin >7))
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pin=args(0).matrix_value();
	T=args(1).matrix_value();
	G=args(2).matrix_value();
	Tstd=args(3).double_value();
	Gstd=args(4).double_value();
	for (i=5;i<nargin;i++)
  	{
		if (args(i).is_string ())
		{
			std::string Aunit;
			Aunit=args(i).string_value();
			int cc=0;
			if (ud)
			{
				fprintf(stderr,"Warning: cannot define temperature unit more than once!\n");
				print_usage ();
				return octave_value_list ();
			}
				
			cc+=(strncmp("C",		(char *) Aunit.c_str(), 10)==0);
			cc+=(strncmp("Celcius",	(char *) Aunit.c_str(), 10)==0);
			if (cc)
			{
				unit='C';
				ud=1;
			}
			else 
			{
				cc+=(strncmp("K",		(char *) Aunit.c_str(), 10)==0);
				cc+=(strncmp("Kelvin",	(char *) Aunit.c_str(), 10)==0);
				if (cc)
					ud=1;
			}
		}
		if (args(i).is_scalar_type ())
		{
			if (FN)
			{
				off=args(i).int_value ();
			}
			else
			{
				FN = args(i).int_value ();
				if (FN<0)
					FN=0;
			}
		}
	}
				
			
	
	Nr=Pin.rows();
	Nc=Pin.cols();
	if ((Nr<3)||(Nc!=10))
	{
		fprintf(stderr, "Number of rows in esp may not less than 3\n");
		fprintf(stderr, "and number of columns must be 10\n");
		return octave_value_list ();
	}
	N=T.rows();
	Nc=T.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in T must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	N=G.rows();
	Nc=G.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in G must be equal to rows in esp\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	
	
	// create a list of esp	
	if ((Cin=(curvepars *)malloc(N*sizeof(curvepars)))==NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((TT=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(Cin);
		free(TT);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	
	for (i=0;i<N;i++)
	{
		Cin[i].Voc=Pin(i,0);
		Cin[i].Isc=Pin(i,1);
		Cin[i].Vmpp=Pin(i,2);
		Cin[i].Impp=Pin(i,3);
		Cin[i].Vqmp1=Pin(i,4);
		Cin[i].Iqmp1=Pin(i,5);
		Cin[i].Vqmp2=Pin(i,6);
		Cin[i].Iqmp2=Pin(i,7);
		Cin[i].Gsc=Pin(i,8);
		Cin[i].Roc=Pin(i,9);	
		TT[i]=T(i,0);	
		GG[i]=G(i,0);	
	}
	
	Pstd=pv_craze_fit_diode_GT(Cin, GG,TT, N, &alpha, &Ea,Tstd,Gstd, unit,FN,off);
	
	free(GG);
	free(TT);
	free(Cin);
	
	Pout=RowVector(11);
	Pout(0) =(double)Pstd.TYPE;
	Pout(1) =Pstd.pars[0];
	Pout(2) =Pstd.pars[1];
	Pout(3) =Pstd.pars[2];
	Pout(4) =Pstd.pars[3];
	Pout(5) =Pstd.pars[4];
	Pout(6) =Pstd.pars[5];
	Pout(7) =Pstd.pars[6];
	Pout(8) =Pstd.pars[7];
	Pout(9) =Pstd.pars[8];	
	Pout(10)=Pstd.pars[9];	


	i=0;
	if (nargout>i)
		retval(i++)=octave_value(Pout);
	if (nargout>i)
		retval(i++)=octave_value(alpha);
	if (nargout>i) 
		retval(i++)=octave_value(Ea);
	return retval;	
}


DEFUN_DLD (pvcraze_translate_diode_GT, args, nargout, "-*- texinfo -*-\n\
@deftypefn {} {[@var{modelpar}] =} pvcraze_translate_diode_GT (@var{modelpar_std}, Tstd, Gstd, T, G, alpha, Ea [, unit])\n\
pvcraze_translate_diode_GT is part of the PV-CRAZE package and translates \
diode model parameters from Tstd,Gstd to T,G, using the temperature \
coeefficient alpha and the activation energy Ea.\n\n\
Input:\n\
@itemize \n\
@item @var{modelpar_std}: model parameters at Tstd,Gstd\n\
@item @var{Tstd}: Standard tempersture\n\
@item @var{Gstd}: Standard irradiance\n\
@item @var{T}: Array with output temperatures\n\
@item @var{G}: Array with output irradiances\n\
@item @var{alpha}:  temperature coefficient of the photocurrent\n\
@item @var{Ea}:     activation energy of the saturation current\n\
@item @var{unit}:   Temperature unit (\"C\"/\"Celcius\" or \"K\"/\"Kelvin\")\n\
@end itemize\n\
Output:\n\
@itemize \n\
@item @var{modelpar}: Array with translated model parameters, one per G,T pair\n\
@end itemize\n\
@end deftypefn")
{
  	int nargin = args.length ();
	octave_value_list retval;
	generic_pars P;
	generic_pars *Pt;
	Matrix T,G,Pin,Pout;
	double *TT, *GG;
	double alpha, Ea, Gstd, Tstd;
	int N, Nr,Nc;
	int i;
	char unit ='K';
  	if ((nargin != 7)&&(nargin != 8))
	{
    	print_usage ();
		return octave_value_list ();
	}
	
	Pin=args(0).matrix_value();
	Tstd=args(1).double_value();
	Gstd=args(2).double_value();	
	T=args(3).matrix_value();
	G=args(4).matrix_value();
	alpha=args(5).double_value();
	Ea=args(6).double_value();
	
  	if (nargin==8)
  	{
		std::string Aunit;
		Aunit=args(7).string_value();
		int cc=0;
		cc+=(strncmp("C",		(char *) Aunit.c_str(), 10)==0);
		cc+=(strncmp("Celcius",	(char *) Aunit.c_str(), 10)==0);
		if (cc)
			unit='C';
		else 
		{
			cc+=(strncmp("K",		(char *) Aunit.c_str(), 10)==0);
			cc+=(strncmp("Kelvin",	(char *) Aunit.c_str(), 10)==0);
			if (!cc)
			{
				fprintf(stderr, "Unknown temperature unit \"%s\"\n", (char *) Aunit.c_str());
				fprintf(stderr, "Using default unit Kelvin\n");
			}
		}
	}
	Nr=Pin.rows();
	Nc=Pin.cols();
	if ((Nr!=1)||(Nc!=11))
	{
		fprintf(stderr, "Number of rows in modelpar_std must be 1\n");
		fprintf(stderr, "and number of columns must be 11\n");
		return octave_value_list ();
	}
	N=T.rows();
	Nc=T.cols();
	if (Nc!=1)
	{
		fprintf(stderr, "Number of columns in T must be 1\n");
		return octave_value_list ();
	}
	Nr=G.rows();
	Nc=G.cols();
	if ((N!=Nr)||(Nc!=1))
	{
		fprintf(stderr, "Number of rows in G must be equal to rows in T\n");
		fprintf(stderr, "and number of columns must be 1\n");
		return octave_value_list ();
	}
	
	P.TYPE =(MODEL_TYPE)Pin(0);
	P.pars[0] = Pin(1);
	P.pars[1] = Pin(2);
	P.pars[2] = Pin(3);
	P.pars[3] = Pin(4);
	P.pars[4] = Pin(5);
	P.pars[5] = Pin(6);
	P.pars[6] = Pin(7);
	P.pars[7] = Pin(8);
	P.pars[8] = Pin(9);
	P.pars[9] = Pin(10);
	
	
	
	if ((TT=(double *)malloc(N*sizeof(double)))==NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	if ((GG=(double *)malloc(N*sizeof(double)))==NULL)
	{
		free(TT);
		fprintf(stderr, "Memory allocation failed\n");
		return octave_value_list ();
	}
	
	for (i=0;i<N;i++)
	{
		TT[i]=T(i,0);	
		GG[i]=G(i,0);			
	}
	Pt=pv_craze_translate_diode_GT(P,Gstd,Tstd, GG, TT, N, alpha, Ea, unit);
	free(GG);
	free(TT);
	
	if (!Pt)
	{
		fprintf(stderr, "PV-CRAZE returned a NULL list\n");
		return octave_value_list ();		
	}
	
	Pout=Matrix(N,11);
	for (i=0;i<N;i++)
	{
		Pout(i,0) =(double)Pt[i].TYPE;
		Pout(i,1) =Pt[i].pars[0];
		Pout(i,2) =Pt[i].pars[1];
		Pout(i,3) =Pt[i].pars[2];
		Pout(i,4) =Pt[i].pars[3];
		Pout(i,5) =Pt[i].pars[4];
		Pout(i,6) =Pt[i].pars[5];
		Pout(i,7) =Pt[i].pars[6];
		Pout(i,8) =Pt[i].pars[7];
		Pout(i,9) =Pt[i].pars[8];	
		Pout(i,10)=Pt[i].pars[9];	
	}	
	pv_craze_free(Pt);	

	retval(0)=octave_value(Pout);
	return retval;	
	
}

