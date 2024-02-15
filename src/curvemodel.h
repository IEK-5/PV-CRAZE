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
double 		RMSErr			(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), curve C, MODELSET *S);
double 		Solve_IV		(double V, double Voc, point *dp, void * pars, point (*model)(double, point *, void *), MODELSET *S);
double 		Solve_VI		(double I, double Isc,  point *dp, void * pars, point (*model)(double, point *, void *), MODELSET *S);
point 		FindMp			(point off, void * pars, point (*model)(double, point *, void *), MODELSET *S);
curvepars 	SolPar			(double Voc, double Isc, void * pars, point (*model)(double, point *, void *), MODELSET *S);
