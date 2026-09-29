# PV-CRAZE Documentation
The PV-CRAZE library is not a standalone tool but rather a library with bindings to Python and GNU Octave. In the General Section we discuss the library without the details of the Python and GNU Octave bindings. The specifics of the bindings are discussed in separate Sections.

## General
The PV-CRAZE library provides a collection of routines:
 
### pvcraze_analyze
This routine is the core of PV-CRAZE. It analyzes a list of IV curves and returns the corresponding ESP. In addition the routine returns flags. The flags are binary encoded and can signal various problem with the IV. The flags distinguish between warnings (results are a-typical or inaccurate) and errors (not possible to determine meaningful parameters). Optionally this routine may be used to directly fit a set of curve models to the ESP and compute the RMS error w.r.t. the provided IV curves. See the section on pvcraze_modelfit for a description of supported models. If models are provided, the routine also returns generic parameters (a unified parameter array for all supported curve models), and fit errors. Furthermore, one can specify various numerical settings to specify how the IV curves are analyzed.

The numerical settings for the curve analysis are:
- Spike Filter

> Boolean value. If true a 3rd order median filter is used to filter the IV prior to analyzing the data.

- Identical Filter

> Boolean value. Remove identical points in IV data. The threshold of what is considered identical is set using IdentEps

- Identical threshold

> Float value. Threshold value used to determine whether values are identical. 

- Resample

> Nicely distribute points in the IV along the curve. Often IV's are measured with equidistant voltages. With increasing voltage the curve tends to become ever steeper. This leads to a reduction in resolution at higher voltages. If you measure the IV with a small step you will be fine. However, for larger voltage steps you may find curves are flagged as having insufficient data to determine various parameters. The best is obviously to re-measure the IV with a finer voltage step. However, to salvage a measurement we can do resampling, where the poorly resolved sections get more points using piece-wise linear interpolation.

- Minimum Relative Range

> The minimum relative range is a constant between 0 and 1 used to check whether the voltage and current ranges in the IV are sufficient to do any analysis at all. This is formalized as

    xmin
    ---- < 1-R
    xmax

> where xmin and xmax are the minimum and maximum value of the current or voltage and R is the minimum relative range. 

- MPP_f_V

> the mpp fit uses data in the curve starting at the voltage (1-MPP_f_V) Vmpp'. Here Vmpp' is the first guess maximum power voltage.

- MPP_f_I;

> the mpp fit uses data in the curve up to the current (1-MPP_f_I) Impp', where Impp' is the first guess maximum power current.


- OC_f_I
> To extract the open circuit voltage we first use a third order polynomial to approximate the IV beyond the qmp+ point. Using this polynomial directly for extrapolating the Voc may become ill conditioned due to Runge's phenomenon. To overcome this we extract two points from the polynomial fit and use a linear extrapolation through these points. The two poinst should be chosen close enough to the Voc end of the curve to allow a linear approximation. The used poinst are V1, the last point in the IV with a current larger than 0. The second point, V2, is defined as V2 = V1 + OC_f_I Vqmp+. The default value for OC_f_I=0.12.

- SC_f_V
> To extract the short circuit current we again use a third order polynomial to approximate the IV up to the qmp- point. Using this polynomial directly for extrapolating the Isc may become ill conditioned due to Runge's phenomenon. To overcome this we extract two points from the polynomial fit and use a linear extrapolation through these points. The two poinst should be chosen close enough to the Isc end of the curve to allow a linear approximation. The used poinst are I1, the first point in the IV with a voltage larger than 0. The second point, I2, is defined as I2 = I1 + SC_f_V Iqmp-. The default value for SC_f_V=0.12.

- mpp-method

> The mpp is determined using 3rd order polynomial fits. The mpp method is one of `f` no noise assumption, fast but not robust. `r` robust OLS fit. `p`: parametric. The parametric option treats the voltage, current, and power all as polynomial functions of time. The optimum time is extracted from the third oder power curve. Subsequently the voltage and current are extracted from the time dependent voltage and current curves. `s` Savitsky-Golay filter. Provides a faster model which is still robust (Savitsky-Golay also uses the parametric approach).


### pvcraze_modelfit

Fit models to ESP data. Currently PV-CRAZE supports the following curve models

- B2: Quadratic Bézier interpolation model for the ESP
- KH: Karmalka-Haneefa model
- PV-spline: A cubic spline interpolation model for the ESP
- Spline: A more primitive interpolation model. It may result in non monotonous interpollations where a monotonous interpolation is possible. In general this model is not recommended. However, contrary to the other models it is possible to fit this model to the actual IV data (all other models are fit to the ESP alone). This way one can obtain a set of ESP which is optimized to represent the IV as a cubic spline.
- PV-Bezier: A cubic Bézier interpolation model for the ESP (the recommended interpolation model)
- Diode: The one diode model
- PC: Pindado-Cubas IV model
- PWL: Plain piece wise linear interpolation for the ESP. Note that the key ESP points are the knots, and hence, the derivatives at the key points are not defined.
- APWL: Derivative respecting piece wise linear interpolation for the ESP. This model computes knots such that the piece wise linear interpolation passes through each key ESP point with the correct derivative.

Model fitting will produce "generic parameters", an array of values, the first of which is a "model ID", followed by the actual model parameters.

Models are specified as a list of strings, where each string sets one model. The string for a model consist of keywords with a values specified between brackets. At the very least a model name must be specified using the ``model()`` keyword. For example "model(Diode)" specifies the diode model. The model names identical to the ones in the item list above (case sensitive). In addition to the model keyword various numerical settings may be specified. Note that the impact of these settings may vary from model to model, and some settings may have no effect with some models. The keywords are:

- "model(model-name)": specify model, where model-name is one of the listed models above
- "fitit(max-iterations)": specify the maximum number of iterations used in the parametrization routine
- "npop(population-size)": specify the population size (if applicable). Some parametrization routines may use a population based optimizer
- "fiteps(epsilon)": specify the stop criterion, where epsilon is some small error value.
- "mit(max-iterations)": Apart from the iterations to fit a model, the model may also emplot iterations for model evaluations (i.e. compute I(V) or V(I))
- "meps(epsilon)": The stop criterion for model evaluations.
- "oeps(epsilon)": Criterion to perform expensive optimizations (if applicable). 
- "refine": applicable only to the Spline model (the plain Spline model, not the PV-Spline model), fit the spline parameters to the IV data
- "rebase": obsolute option

### pvcraze_lookupmodel

Looks up a model ID from a model name.

### pvcraze_curvemodel

Reconstruct an IV from generic parameters. Takes an array with generic parameters and returns a list of IV

### pvcraze_modelerror

Compute RMS error between a modeled curve and and experimental curve. 

### pvcraze_errorcurve

Compute distance between modeled curve and and experimental curve as a function of voltage

### pvcraze_fit_IEC60891_v2

Fit the IEC60891 v2 norm model. This routine takes arrays with ESP, and the corresponding temperature and irradiance. The routine then fits the IEC60891_v2 and returns ESP at standard test conditions, and the various IEC60891_v2 coefficients.

### pvcraze_fit_IEC60891_v3

Fit the IEC60891 v3 norm model. This routine takes arrays with ESP, and the corresponding temperature and irradiance. The routine then fits the IEC60891_v3 and returns ESP at standard test conditions, and the various IEC60891_v3 coefficients.

### pvcraze_translate_IEC60891_v2

Translate ESP using the IEC60891 v2 norm model. Takes as input ESP, input temperature and irradiance, and output temperature and irradiance. The routine then computes for each input ESP the output ESP.

### pvcraze_translate_IEC60891_v3

Translate ESP using the IEC60891 v3 norm model. Takes as input ESP, input temperature and irradiance, and output temperature and irradiance. The routine then computes for each input ESP the output ESP.

## Python Interface
The Python package is a separate, object based interface; it does not mirror
the Octave routines below. See the section [Python](#python) at the end of
this document.

## GNU Octave Interface
In Octave all routine names start with ```pvcraze_<function>```. Mandatory arguments are positional. Optional arguments are passed as name value pairs, e.g. given a function with the mandatory (positional) arguments ```X``` and ```Y``` and the optional argument ```foo```, we can call as ```<function>(1,2,"foo", 3)```, which would make ```X=1```, ```Y=2```, and ```foo=3```.

To load pvcraze use in Octave:

    pkg load pv-craze

The routine decription is built in and accessable using the ```help``` command, e.g. ```help pvcraze_analyze```.

You may look at the scripts in the test directory for examples.

### pvcraze_analyze:
     [esp,flags,modelpar,err] = pvcraze_analyze (IV,[optional-parameters])
     pvcraze_analyze reads solar cell IV characteristics, determines an
     extended set of solar cell parameters and does sanity checks on the
     data.
     Input:
        • IV: Cell array with arrays with two columns (voltage,current)
        • ’models’: (Optional) Takes a model string (this option may 
          be passed multiple times)
        • ’mpp-method’: (Optional) Select mpp method, can be ’f’ (fast,
          not robust), ’r’ (robust), ’p’ (parametric), or ’s’
          Savitsky-Golay
        • ’S_filter’: (Optional) Takes boolean, apply a 3rd order
          median filter on the data to remove spikes
        • ’I_filter’: (Optional) Takes positive float.  Filter out
          subsequent identical voltage and current values.  Specify
          minimum deviation.  default 1e-12
        • ’mpp-dv’: (Optional) Takes float [0-1].  The mpp fitting range
          starts at (1-mpp-dv) Vmpp.  default 0.12
        • ’mpp-di’: (Optional) Takes float [0-1].  The mpp fitting range
          ends at (1-mpp-di) Impp.  default 0.12
        • ’sc-dv’: (Optional) Takes float [0-1].  Voltage range to
          determine sc slope between 0 and sc-dv Vqmp1, default 0.1
        • ’oc-di’: (Optional) Takes float [0-1].  Current range to
          determine oc slope between 0 and oc-di Iqmp2, default 0.1
        • ’resample’: (Optional) Takes an integer >20. Resample the IV
          to nicely distribute the IV points along the curve.
     Output:
        • esp: Extended Solar Cell Parameter Matrix:
          [Voc;Isc;Vmpp;Impp;Vqmpp-,Iqmpp-,Vqmpp+,Iqmpp+,Gsc,Roc]
        • flags: Flags for the analyzed curves
        • modelpar: Cell array with model parameters (one column per model,
          one row per IV)
        • err: Model fit errors w.r.t. the input IV's (one column 
          per model one row per IV)

### pvcraze_modelfit: 
     [modelpar,err] = pvcraze_modelfit (esp, models)
     pvcraze_modelfit takes extended solar cell parameters and fits
     models to them.
     Input:
        • esp: Array of extended solar cell parameters (Nx10)
        • ’models’: One or more model specification strings (multiple models 
          may be passed). 
          Note, that unlike the analysis function the model 
          specifications are not optional, and hence in Octave the 
          models are positional arguments (i.e. do not pass string 
          value pairs).
     Output:
        • modelpar: Cell array with model parameters (one column per model,
          one row per Epar row)
        • err: Fit error array (one column per model, one row per esp 
          row). Note that, as no IV characteristics are provided, the 
          fit error is computed with respect to the 5 key points of the 
          esp.
          
### pvcraze_curvemodel
     [IV] = pvcraze_curvemodel (modelpar, N)
     pvcraze_curvemodel computes IV curves from PV-CRAZE generic 
     parameters

     Input:
        • modelpar: Cell array with model parameters. Each cell must 
          contain a 1x10 or 1x11 sized array. If an array is 1x11, it is
          treated as generic model parameters. If an 1x10 array is 
          passed the PV-Bezier model is used.
        • N: Positive Integer: nuber of points in the IV characteristics
     Output:
        • IV: Cell Array with IV characteristics:
        
### pvcraze_modelpars
     [esp] = pvcraze_modelpars (modelpar)
     pvcraze_modelpars computes esp from generic model parameters

     Input:
        • modelpar: Cell array with model parameters. Each cell must 
          contain a 1x10 or 1x11 sized array. If an array is 1x11, it is
          treated as generic model parameters. If an 1x10 array is 
          passed the PV-Bezier model is used.
     Output:
        • esp: Extended Solar Cell Parameters

### pvcraze_modelerror
     [err] = pvcraze_modelerror (modelpar, IV)
     pvcraze_modelerror computes a RMS distance between IV curves and 
     corresponding generic model parameters

     Input:
        • modelpar: Cell array with model parameters. Each cell must 
          contain a 1x10 or 1x11 sized array. If an array is 1x11, it is
          treated as generic model parameters. If an 1x10 array is 
          passed the PV-Bezier model is used.
        • IV: Cell array with IV curves corresponding to the MPar
          variable
     Output:
        • D: Array with RMS distances for each IV an model parameter set



### pvcraze_errorcurve
     [errcurve] = pvcraze_errorcurve (modelpar, IV)
     pvcraze_errorcurve is part of the PV-CRAZE package and computes 
     the distance curve between IV curves and corresponding generic model 
     parameters. Currently only available in the Octave interface.

     Input:
        • modelpar: Cell array with model parameters. Each cell must 
          contain a 1x10 or 1x11 sized array. If an array is 1x11, it is
          treated as generic model parameters. If an 1x10 array is 
          passed the PV-Bezier model is used.
        • IV: Cell array with IV curves corresponding to the modelpar
          variable
     Output:
        • errcurve: Cell Array with distance curves (cartesian
          distance in a normalized current-voltage plane between curve
          and model w.r.t.  cell voltage
          
### pvcraze_fit_IEC60891_v2
     [esp_std, alpha, beta, a, Rs, kappa] = pvcraze_fit_IEC60891_v2 
     (esp, T, G, Tstd, Gstd [, unit,furthest_N])
     pvcraze_fit_IEC60891_v2 is part of the PV-CRAZE package and fits
     IEC60891 translation coefficients and determines standard test
     condition curve parameters to best describe provided temperature
     and irradiance dependent curve parameters.  Uses the IEC60891
     version 2 (2009). 

     Input:
        • esp: array with extended solar cell parameters, one row
          per conditions
        • T: Array with temperatures, one temperature per row in
          esp
        • G: Array with irradiances, one irradiance per row in esp
        • Tstd: Scalar value for the standard test conditions
          Temperature
        • Gstd: Scalar value for the standard test conditions Irradiance
        • furthest_N: If this value is 0 all esp data is used. If this 
          value is a positive integer, we will use furthest point 
          sampling in G-T space to reduce the number of datapoints used 
          for fitting. This ensures the full range of irradiance and 
          temperature is used but the number of used data points is 
          limited to a maximum of furthest_N.			
     Output:
        • esp_std: Parameters at standard test conditions
        • alpha: IEC60891 alpha
        • beta: IEC60891 beta
        • a: IEC60891 a
        • Rs: IEC60891 Rs
        • kappa: IEC60891 kappa
        
### pvcraze_fit_IEC60891_v3
	 [esp_std, alpha, beta, B1, B2, Rs, kappa] = pvcraze_fit_IEC60891_v3 
	 (esp, T, G, Tstd, Gstd [,furthest_N])
	 pvcraze_fit_IEC60891_v3 is part of the PV-CRAZE package and fits
	 IEC60891 translation coefficients and determines standard test
	 condition curve parameters to best describe provided temperature
	 and irradiance dependent curve parameters Uses the IEC60891 version
	 3 (2021). 

     Input:
        • esp: array with extended solar cell parameters, one row
          per conditions
        • T: Array with temperatures, one temperature per row in
          esp
        • G: Array with irradiances, one irradiance per row in esp
        • Tstd: Scalar value for the standard test conditions
          Temperature
        • Gstd: Scalar value for the standard test conditions Irradiance
        • furthest_N: If this value is 0 all esp data is used. If this 
          value is a positive integer, we will use furthest point 
          sampling in G-T space to reduce the number of datapoints used 
          for fitting. This ensures the full range of irradiance and 
          temperature is used but the number of used data points is 
          limited to a maximum of furthest_N.			
     Output:
        • esp_std: Parameters at standard test conditions
        • alpha: IEC60891 alpha
        • beta: IEC60891 beta
        • B1: IEC60891 B1
        • B2: IEC60891 B2
        • Rs: IEC60891 Rs
        • kappa: IEC60891 kappa
        
### pvcraze_translate_IEC60891_v2
     esp = pvcraze_translate_IEC60891_v2 (esp, T1, G1, T2, G2, alpha, 
     beta, a, Rs, kappa)
     pvcraze_translate_IEC60891_v2 is part of the PV-CRAZE package and
     translates curve parameters from conditions G1, T1 to G2, T2, using
     procedure 2 from IEC60891 norm version 2 (2009)

     Input:
        • esp: array with extended solar cell parameters, one row
          per conditions
        • T1: Array with temperatures, one temperature per row in
          esp
        • G1: Array with irradiances, one irradiance per row in esp
        • T2: Array with output temperatures
        • G2: Array with output irradiances
        • alpha: IEC60891 alpha
        • beta: IEC60891 beta
        • a: IEC60891 a
        • Rs: IEC60891 Rs
        • kappa: IEC60891 kappa
     Output:
        • esp: Array with translated extended parameters.  The
          length of the array is the maximum of the length of the input
          (T1/G1 or G2/T2)

### pvcraze_translate_IEC60891_v3
     esp = pvcraze_translate_IEC60891_v3 (esp, T1, G1, T2, G2, Vocstd, 
     Tstd, Gstd, alpha, beta, B1, B2, Rs, kappa)
     pvcraze_translate_IEC60891_v3 is part of the PV-CRAZE package and
     translates curve parameters from conditions G1, T1 to G2, T2, using
     procedure 2 from IEC60891 norm version 3 (2021)

     Input:
        • esp: array with extended solar cell parameters, one row
          per conditions
        • T1: Array with temperatures, one temperature per row in
          esp
        • G1: Array with irradiances, one irradiance per row in esp
        • T2: Array with output temperatures
        • G2: Array with output irradiances
        • Vocstd: Open-Circuit voltage at standard test conditions
        • Tstd: Standard test conditions temperature
        • Gstd: Standard test conditions irradiance
        • alpha: IEC60891 alpha
        • beta: IEC60891 beta
        • B1: IEC60891 B1
        • B2: IEC60891 B2
        • Rs: IEC60891 Rs
        • kappa: IEC60891 kappa
     Output:
        • esp: Array with translated extended parameters.  The
          length of the array is the maximum of the length of the input
          (T1/G1 or G2/T2)
 
### pvcraze_fit_diode_GT
	 [modelpar, alpha, Ea] = pvcraze_fit_diode_GT 
	 (esp, T, G, Tstd, Gstd [, unit,furthest_N])
	 pvcraze_fit_diode_GT is part of the PV-CRAZE package and fits
	 a temperature and irradiance dependent one diode model. The 
	 result is a set of model parameters at the given standard 
	 irradiance and temperature conditions. In addition the temperature 
	 coefficient of the short circuit current as well as the activation 
	 energy of the saturation current are returned. These values may be 
	 used to translate the standard modelpar parameters to a given 
	 temperature and irradiance condition. In Octave, the last two 
	 optional arguments are, as an exception to the general rule, not 
	 passed as key value pairs. As the unit variable has few valid 
	 string values and furthest_N is an integer, the meaning of the 
	 last two optional arguments is inferred from its contents.

     Input:
        • esp: array with extended solar cell parameters, one row
          per conditions
        • T: Array with temperatures, one temperature per row in
          esp
        • G: Array with irradiances, one irradiance per row in esp
        • Tstd: Scalar value for the standard test conditions
          Temperature
        • Gstd: Scalar value for the standard test conditions Irradiance        
        • unit: Temperature unit string. May have the values:
			- "C" or "Celcius"  : Temperature is in Celcius
			- "K" of "Kelvin"   : Temperature is in Kelvin
        • furthest_N: If this value is 0 all esp data is used. If this 
          value is a positive integer, we will use furthest point 
          sampling in G-T space to reduce the number of datapoints used 
          for fitting. This ensures the full range of irradiance and 
          temperature is used but the number of used data points is 
          limited to a maximum of furthest_N.			
     Output:
        • modelpar: Diode model parameters at standard test conditions
        • alpha: Temperature coefficient for the photo-current
        • Ea: Activation energy of the saturation current.
        
### pvcraze_translate_diode_GT
	 [modelpar] = pvcraze_translate_diode_GT 
	 (modelpar, Tstd, Gstd, T, G alpha, Ea [, unit])
	 pvcraze_translate_diode_GT is part of the PV-CRAZE package and 
	 translates a set of diode model parameters to a different 
	 temperature and irradiance, given a temperature coefficient for 
	 the photo-current and an activation energy for the saturation 
	 current. In Octave, the last optional argument is, as an exception 
	 to the general rule, not passed as a key value pair.

     Input:
        • modelpar: Diode model parameters at standard test conditions
        • Tstd: Scalar value for the standard test conditions
          Temperature
        • Gstd: Scalar value for the standard test conditions Irradiance    
        • T: Array with temperatures
        • G: Array with irradiances, one irradiance per temperature value  
        • unit: Temperature unit string. May have the values:
			- "C" or "Celcius"  : Temperature is in Celcius
			- "K" of "Kelvin"   : Temperature is in Kelvin
     Output:
        • modelpar: Diode model parameters at the given temperature and 
          irradiance conditions
             
### pvcraze_lookupmodel
     Mid = pvcraze_lookupmodel (modelname)
     pvcraze_lookupmodel searches for the model id corresponding to
     the string modelname.

     Input:
        • modelname: Name of the model (string), e.g. PV-Bezier or Diode
     Output:
        • Mid: Model ID
        
### pvcraze_reverselookup
     modelname = pvcraze_lookupmodel (Mid)
     pvcraze_reverselookup searches for the model name corresponding to
     the model id.

     Input:
        • Mid: id number of the model, e.g. 0
     Output:
        • modelname: Name of the model (string), e.g. PV-Bezier or Diode
   
### pvcraze_flagtable
     [flags,messages,err] = pvcraze_flagtable ()
     pvcraze_flagtable generates a table of principle flags,
     corresponding messages, and a boolean value whether the flag 
     constitutes an error. Note that this table may change between 
     pvcraze versions so you may want to store this data together with 
     analysis results.

     Input: (none)
     Output:
        • flags: Principle flags (powers of two, actual flags may be 
                 the sum of several principle flags)
        • messages: messages corresponding to principle flags 
        • err: some messages are merely warnings and no errors. Err is 
               a boolean array with True for error flags.
  
### pvcraze_quiet
     pvcraze_quiet()
     Set pv-craze verbose level to quiet.
     
### pvcraze_verbose
     pvcraze_verbose()
     Set pv-craze verbose level to verbose.
     
### pvcraze_verbose_setreport
     pvcraze_verbose_setreport(set-value)
     Enable/Disable analysis report printing. If set-value evaluates to true a report is printed.
     
### pvcraze_verbose_setprogress
     pvcraze_verbose_setprogress(set-value)
     Enable/Disable progress bar. If set-value evaluates to true a progress bar is displayed.
     
### pvcraze_verbose_settiming
     pvcraze_verbose_settiming(set-value)
     Enable/Disable timing report. If set-value evaluates to true a timing information is displayed.
     
### pvcraze_verbose_sethead
     pvcraze_verbose_sethead(set-value)
     Enable/Disable header output. If set-value evaluates to true a header is displayed.
	
## C
PV-CRAZE is written in C, so naturally it also provides a C interface (upon which the GNU Octave interface builds; the Python package compiles the same sources into its extension module). Installing the library also installs the library header, "pvcraze.h". This single header file contains all necessary data structures and declarations of the library symbols. 

The most important data structures of the library are:

typedef struct {
	double x,y;
} point;

typedef struct {
	point *P;
	int N;
} curve;

    typedef struct {
	    double Voc, Isc;     // sc oc
	    double Vmpp, Impp;   // mpp
	    double Vqmp1, Iqmp1; // quasi mpp-
	    double Vqmp2, Iqmp2; // quasi mpp+
	    double Gsc, Roc;
    } curvepars;

The ```curvepars``` data structure contains the extended solar cell parameters. Apart for esp's, pv-craze also deals with various other models. These models are handled with a "generic" model parameter struct. It is defined as:

     typedef struct {
	     double pars[10];
     	 MODEL_TYPE TYPE;
     } generic_pars;

The ```generic_pars``` data structure contains a fixed size array of 10 parameters. The meaning of the parameters, however, depends on the ```MODEL_TYPE```. The ```MODEL_TYPE``` is an enumerator type and is automatically generated during compilation and depends on the models which are registered during compilation. This automated registering of models it easier to add new models but may have as a side effect that a new version of the library may change the value of the ```MODEL_TYPE``` refering to a specific model (use ```pv_craze_lookupmodel``` for this, see below). 

The ```CurveProp``` data structure contains the esp (curvepars-struct), flags for the curve, and lists of ```generic_pars``` and fit errors for the various models which are fitted to the curve.

     typedef struct {
     	curvepars par;								// curve parameters
     	int A_FLAGS;								// analysis flags
     	double *E;									// fit errors
     	generic_pars *P;							// parameters for the various fit models
     	int N;										// number of fit models
     } CurveProp;



### Routine Index	 
#### pv_craze_initsettings
     ProcessSettings pv_craze_initsettings();
     Initialize a "ProcessSettings" struct needed to analyze IV Curves.

#### pv_craze_addmodeltosettings
     void pv_craze_addmodeltosettings(ProcessSettings *S, const char *str);
     Add a model specified by a string to a ProcessSettings-struct (see modelfit section)

#### pv_craze_freesettings
     void pv_craze_freesettings(ProcessSettings *S);
     Clean up a ProcessSettings-struct.

#### pv_craze_analyze
     CurveProp * pv_craze_analyze(curve *C, int N, ProcessSettings *S);
     Process a list of curves in C, with N elements and ProcessSettings S. It returns an allocated list of corresponding curve properties. Thus, pv_craze_analyze determined the esp of the IV data and fits the models which were added to the ProcessSettings (i.e. if no model were added only the esp are determined).
     
#### pv_craze_modelfit
     CurveProp * pv_craze_modelfit(curvepars *par, int N, ProcessSettings *S);
     Given a list of esp's and a ProcessSettings struct, fit the models which were added to the ProcessSettings to the esp's. It returns an allocated list of corresponding curve properties. 
     
#### pv_craze_freeprop
     void pv_craze_freeprop(CurveProp *P);
     Free the results produced by pv_craze_analyze or pv_craze_modelfit.

#### pv_craze_curvemodel
curve *pv_craze_curvemodel(generic_pars *pars, int N, int M);

#### pv_craze_error
double *pv_craze_error(generic_pars *pars, curve *C, int N);

#### pv_craze_errorcurve
curve *pv_craze_errorcurve(generic_pars *pars, curve *C, int N);

#### pv_craze_fit_IEC60891_v2

#### pv_craze_fit_IEC60891_v3

#### pv_craze_translate_IEC60891_v2

#### pv_craze_translate_IEC60891_v3


#### pv_craze_g_distcurve

#### pv_craze_g_fit

#### pv_craze_g_pars

#### pv_craze_g_rmse

#### pv_craze_g_tocurve

#### pv_craze_lookupmodel
    MODEL_TYPE pv_craze_lookupmodel(char *word);    
	This routine takes an identifyer string (e.g. "Diode") and returns the corresponding ```MODEL_TYPE```. 
	
#### pv_craze_reverselookup	
char *pv_craze_reverselookup(MODEL_TYPE T);


#### pv_craze_numfields

#### pv_craze_flagtable
AnalyzeFlag *pv_craze_flagtable(int *N);

#### pv_craze_quiet
     void pv_craze_quiet()
     Set pv-craze verbose level to quiet.

#### pv_craze_verbose
     void pv_craze_verbose()
     Set pv-craze verbose level to verbose.

#### pv_craze_verbosehead
     void pv_craze_verbosehead(int set);
     Enable/Disable header output. If set-value evaluates to true a header is displayed.

#### pv_craze_verboseprogress
     void pv_craze_verboseprogress(int set);
     Enable/Disable progress bar. If set-value evaluates to true a progress bar is displayed.  

#### pv_craze_verbosereport
     void pv_craze_verbosereport(int set);
     Enable/Disable analysis report printing. If set-value evaluates to true a report is printed.

#### pv_craze_verbosetiming
     void pv_craze_verbosetiming(int set);
     Enable/Disable timing report. If set-value evaluates to true a timing information is displayed.


## Python
### Installation
The Python package compiles the C sources into its extension module, so it
does not need an installed libpvcraze. It needs a C compiler and the
development files of GSL and LAPACK (found via pkg-config). From the
repository root:

    pip install .              # add [pandas] for DataFrame export

Build and install use [meson-python](https://mesonpy.readthedocs.io); the
build description is `meson.build` in the repository root. The C library and
the Octave package are still built with autotools.

### Overview
```python
import numpy as np
import pvcraze as pc

res = pc.analyze(V, I, models=["PV-Bezier", pc.FitSpec("Diode", population=200)])
res.ok                            # False if an error flag is set
res.flag().messages               # human readable warnings/errors
res.esp                           # ESP(voc=..., isc=..., ...)
res.esp.voc, res.esp.ff           # fields and derived quantities
bez = res.fits["PV-Bezier"]       # FitResult: .params, .error, .spec
bez.params.iv([0.0, 0.3, 0.5])    # evaluate I(V)
bez.params.curve(200)             # (n, 2) array (V, I)

fit = pc.fit_model(res.esp, "Diode")   # fit from ESP only
fit.params.rs, fit.params.rsh          # named parameters
```

Everything has docstrings; `help(pc)`, `help(pc.analyze)`, `help(pc.ESP)`
and so on are the reference. Conventions:

* **Units** are SI (V, A, S, Ohm); temperatures are in **degrees Celsius**,
  irradiances in W/m2. IV curves use the generator convention (positive
  current between short circuit and open circuit).
* **Single vs batch.** `ESP`, `ModelParams`, `AnalysisResult` and
  `FitResult` are either *single* (fields are floats) or a *batch* (fields are
  1-D arrays). The form follows the input: one curve in gives single
  results, a list of curves gives batch results. Batches support `len`,
  iteration, integer indexing (giving a single object), slicing and boolean
  masks. `np.asarray(esp)` gives the raw `(10,)` / `(N, 10)` array.
* **Models** are given by name (case and punctuation insensitive, e.g.
  `"pv-bezier"`), as a `pc.Model` member, or as a `pc.FitSpec` with fitting
  options. Unknown names raise `ValueError` listing the valid ones.
* **Objects are immutable**; `to_array()`, `to_dict()` and `to_dataframe()`
  return copies.
* The C library is **silent** by default. Use `with pc.verbose(): ...` or
  `pc.set_verbosity(report=True)` to see its console output.

### Curve input
`pc.analyze` accepts

| call | meaning |
|---|---|
| `analyze(iv)` | one curve, `iv` an `(n, 2)` array with columns (V, I) |
| `analyze([iv0, iv1, ...])` | a batch; curves may differ in length |
| `analyze(V, I)` | one curve given as 1-D arrays |
| `analyze([V0, V1, ...], [I0, I1, ...])` | a batch given as 1-D arrays |

Any real numeric dtype is accepted and converted to float64.

### ESP fields
| field | unit | meaning |
|---|---|---|
| `voc` | V | open-circuit voltage |
| `isc` | A | short-circuit current |
| `vmpp`, `impp` | V, A | maximum power point |
| `vqmp_minus`, `iqmp_minus` | V, A | qmp- point, between short circuit and the MPP |
| `vqmp_plus`, `iqmp_plus` | V, A | qmp+ point, between the MPP and open circuit |
| `gsc` | S | -dI/dV at short circuit |
| `roc` | Ohm | -dV/dI at open circuit |

plus the derived `pmpp` and `ff`.

### Models and parameters
`pc.Model` lists the models (`APWL`, `B2`, `DIODE`, `KH`, `PC`, `PV_BEZIER`,
`PV_SPLINE`, `PWL`, `SPLINE`); `Model.X.param_names` gives the parameter
names used as attributes of `ModelParams`:

| model | parameters |
|---|---|
| PV-Bezier, PV-spline, B2, Spline, PWL, APWL | the ESP fields (`voc`, `isc`, ...) |
| Diode | `i0` [A], `nkt` [V], `iph` [A], `rs` [Ohm], `rsh` [Ohm] |
| KH | `voc`, `isc`, `gamma`, `m` |
| PC | `voc`, `isc`, `vmpp`, `impp`, `eta` |

`ModelParams` can also be built by hand, e.g.
`pc.ModelParams("Diode", i0=1e-10, nkt=0.03, iph=0.035, rs=0.1, rsh=700)`,
and provides `iv(V)`, `vi(I)`, `curve(npoints)`, `esp()`, `rmse(curves)` and
`distance(curves)`.

`FitSpec(model, max_iter=, tol=, population=, eval_max_iter=, eval_tol=,
refine_tol=, refine=)` replaces the former model specification strings
(`fitit`, `fiteps`, `npop`, `mit`, `meps`, `oeps`, `refine`).

### Flags
`AnalysisResult.flags` holds the raw integer codes (an int64 array for a
batch); `res.ok` / `res.has_warning` are vectorised checks and
`res.flag(i)` returns a `pc.Flag` (an `enum.IntFlag`) with `.messages`,
`.has_error`, `.errors` and `.warnings`. Test individual flags by name, e.g.
`pc.Flag.NOMPP in res.flag(i)`.

### Irradiance and temperature translation
```python
from pvcraze import iec60891, diode_gt

coef = iec60891.fit(esp, T=T, G=G, version=3)   # IEC60891v3 (or version=2)
coef.esp_ref                                    # ESP at t_ref, g_ref (25 degC, 1000 W/m2)
coef.translate(T=60, G=800)                     # reference -> other conditions
coef.translate(T=25, G=1000, esp=esp, T_from=T, G_from=G)   # any -> any
iec60891.IEC60891v2(alpha=..., beta=..., a=..., rs=..., kappa=...)  # known coefficients

gt = diode_gt.fit(esp, T=T, G=G)                # DiodeGT
gt.at(T=60, G=800)                              # diode ModelParams at those conditions
gt.esp(T=T, G=G)
```
Condition arguments broadcast against each other and the ESP rows. The
`max_points` / `seed_index` options of the fitting routines enable furthest
point sampling in the G-T plane to limit the number of conditions used.

### Migrating from the pvcraze 1.x Python API
| old | new |
|---|---|
| `pvcraze.analyze(IV, models=["model(PV-Bezier)"])['esp']` | `pc.analyze(IV, models="PV-Bezier").esp` |
| `B['flags'][:, 0]` | `res.flags`, `res.ok`, `res.flag(i).messages` |
| `B['modelpar'][k]` (N x 11, model id in column 0) | `res.fits[name].params` |
| `B['err'][:, k]` | `res.fits[name].error` |
| `modelfit(esp, [spec])` | `fit_model(esp, "Diode", population=200)` |
| `modelpars([mp])[0]` | `params.esp()` |
| `curvemodel([mp], N)[0]` | `params.curve(N)` |
| `modelerror([mp], IV)` | `params.rmse(IV)` |
| `iv(V, mp)`, `vi(I, mp)` | `params.iv(V)`, `params.vi(I)` |
| `lookupmodel`, `reverselookup` | `pc.Model("name")`, `Model.X.id` |
| `flagtable()` | `pc.Flag`, `pc.describe_flags(code)` |
| `quiet()`, `verbose()`, `verbose_set_*` | silent by default, `pc.verbose()`, `pc.set_verbosity()` |
| `fit_IEC60891_v2/v3`, `translate_IEC60891_v2/v3` | `iec60891.fit(..., version=)`, `coef.translate(...)` |
| `fit_diode_GT`, `translate_diode_GT` | `diode_gt.fit(...)`, `gt.at(T, G)` |
| `'Tunit'` arguments | temperatures are always degrees Celsius |
