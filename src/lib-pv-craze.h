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
#include "curve.h"
#include "fitmodels.h"
#include "lib-pv-craze-structs.h"
//START-HEADER-EXPORT mark the start of a block that describes the C api
/* **************************************************
 * io and data/model management tools 
 * **************************************************/
ProcessSettings pv_craze_initsettings();
void pv_craze_addmodeltosettings(ProcessSettings *S, const char *str);
void pv_craze_freesettings(ProcessSettings *S);
void pv_craze_freeprop(CurveProp *P);

// verbosity controll
void pv_craze_verbose(); // all verbose
void pv_craze_quiet();	// all quiet
// set/unset single output types
void pv_craze_verbosetiming(int set);
void pv_craze_verbosereport(int set);
void pv_craze_verboseprogress(int set);
void pv_craze_verbosehead(int set);


MODEL_TYPE pv_craze_lookupmodel(char *word);
char *pv_craze_reverselookup(MODEL_TYPE T);
AnalyzeFlag *pv_craze_flagtable(int *N);

/* **************************************************
 * curve modeling and analysis routines 
 * **************************************************/
// analyze and optionally fit
CurveProp pv_craze_analyze(curve C, ProcessSettings *S, double *tana);
CurveProp * pv_craze_analyze_list(curve *C, int N, ProcessSettings *S);
// only fit
CurveProp pv_craze_modelfit(curvepars par, ProcessSettings *S);
CurveProp * pv_craze_modelfit_list(curvepars *par, int N, ProcessSettings *S);
// evaluate model and model errors
curve *pv_craze_curvemodel(generic_pars *pars, int N, int M);
double *pv_craze_error(generic_pars *pars, curve *C, int N);
curve *pv_craze_errorcurve(generic_pars *pars, curve *C, int N);

/* **************************************************
 * IEC60891 (version 2 and 3) procedure 2 routines
 * **************************************************/
// IEC60891 version 2 (2009)
curvepars *pv_craze_translate_IEC60891_v2(curvepars *pars, double *G1, double *T1, int N, double *G2, double *T2, int M, 
                   double alpha, double beta, double a, double Rs, double kappa);
curvepars pv_craze_fit_IEC60891_v2(curvepars *pars, double *G, double *T, int N, double Tstd, double Gstd,
                   double *alpha, double *beta, double *a, double *Rs, double *kappa, int FURTHEST_N);
                   
// IEC60891 version 3 (2021)
curvepars *pv_craze_translate_IEC60891_v3(curvepars *pars, double *G1, double *T1, int N, double *G2, double *T2, int M, 
                   double Vocstd, double Tstd, double Gstd, double alpha, double beta, double B1, double B2, double Rs, double kappa);
curvepars pv_craze_fit_IEC60891_v3(curvepars *pars, double *G, double *T, int N, double Tstd, double Gstd,
                   double *alpha, double *beta, double *B1, double *B2, double *Rs, double *kappa, int FURTHEST_N);
                   
/* **************************************************
 * Diode model G-T translation model
 * **************************************************/
// unit is either 'C' or 'K' and refers to the temperature unit.
generic_pars pv_craze_fit_diode_GT(curvepars *pars, double *G, 
                double *T, int N, double *alpha, double *Ea, 
                double Tstd, double Gstd, char unit, int FURTHEST_N);
generic_pars pv_craze_translate_single_diode_GT(generic_pars pars, 
                double Gstd, double Tstd, double G, double T, 
                double alpha, double Ea, char unit);
generic_pars *pv_craze_translate_diode_GT(generic_pars pars, double Gstd, 
                double Tstd, double *G, double *T, int N, double alpha, 
                double Ea, char unit);
//END-HEADER-EXPORT mark the end of a block that describes the C api
