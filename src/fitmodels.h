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
#include "models.h"
#include <stdbool.h>
//START-HEADER-EXPORT mark the start of a block that describes the C api
/* typedef for registering model parameter names */
typedef struct {
	MODEL_TYPE TYPE;
	int N;
	char *labels[10];
} model_pars;
/* typedef for a generic parameter struct. Parameters are passes as a flat array*/
typedef struct {
	double pars[10];
	MODEL_TYPE TYPE;
} generic_pars;
/* Lists all parameters for all models */
extern model_pars MODEL_PROP[];

typedef struct {
	MODEL_TYPE M;
	int fititer;		// max iteration for basic fitting algorithm (if applicale)
	double fiteps;		// eps for basic fitting optimizer (if applicale)
	int miter;			// iterations for model evaluation (if applicale)
	int npop;			// population size (if applicale)
	double meps;		// model evaluation eps (if applicale)
	double oeps;		// eps threshold for expensive optim,ization (if applicable)
	bool RefinePars;	// boolean to change the solar cell parameters such that the curve fit is better (not all models have to support this)
} MODELSET;
extern MODELSET DefaultMODSET;

typedef struct {
	MODELSET *List;
	generic_pars *pars;
	double *Efit;
	double *tm; // cumulative time counter for each model (never resets only adds)
	int N;
} MODEL_LIST;


/* make a lookup model function to return a MODEL struct
 * Then make generic routines which take
 * generic_pars and whatever else to create the 4 basic model routines
 * for the generic pars type
 */

curve 			pv_craze_g_tocurve(generic_pars P, int N);
double 			pv_craze_g_rmse(generic_pars P, curve C);
curve 			pv_craze_g_distcurve(generic_pars P, curve C);
generic_pars 	pv_craze_g_fit(MODEL_TYPE T, curvepars par, curve C, double *Efit);
curvepars 		pv_craze_g_pars(generic_pars);
int 			pv_craze_numfields(MODEL_TYPE T);

//END-HEADER-EXPORT mark the end of a block that describes the C api

MODEL_TYPE LookupModel(char *word);
void AddModel(MODEL_LIST *L, const char *spec);
void FitModels(MODEL_LIST *L, curvepars par, curve C, int R);
int PropIndex(MODEL_TYPE T);
void FreeList(MODEL_LIST *L);
