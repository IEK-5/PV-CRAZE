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
//START-HEADER-EXPORT mark the start of a block that describes the C api
typedef struct {
	curvepars par;								// curve parameters
	int A_FLAGS;								// analysis flags
	double *E;									// fit errors
	generic_pars *P;							// parameters for the various fit models
	int N;										// number of fit models
} CurveProp;

typedef struct {
	analyze_settings analyze;
	MODEL_LIST L;
} ProcessSettings;

typedef struct {
	int FLAG;
	char *message;
	int ERR;
} AnalyzeFlag;

//END-HEADER-EXPORT mark the end of a block that describes the C api
