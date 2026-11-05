/*
	Routines to aid building new curve models for PV-CRAZE.

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

curve 		ModelToCurve	(void * pars, point (*model)(double, point *, void *), int N, MODELSET *S);
curve 		DistCurve		(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), curve C, MODELSET *S);
/* number of points on which a model is sampled for the area distance */
#define AREA_NSAMPLE 1024
double 		AreaErr			(void * pars, point (*model)(double, point *, void *), curve C, MODELSET *S);
curve 		SampleModel		(void * pars, point (*model)(double, point *, void *), int N);
double 		CurveAreaDist	(curve R, curve C);
double 		Solve_IV		(double V, double Voc, point *dp, void * pars, point (*model)(double, point *, void *), MODELSET *S);
double 		Solve_VI		(double I, double Isc,  point *dp, void * pars, point (*model)(double, point *, void *), MODELSET *S);
int 		RootIllinois	(double a, double fa, double b, double fb, int (*f)(double, void *, double *), void *ctx, double *root);
point 		FindMp			(point off, void * pars, point (*model)(double, point *, void *), MODELSET *S);
curvepars 	SolPar			(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), MODELSET *S);
