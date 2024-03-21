/*
    Karmalkar-Haneefa IV model.
    Copyright (C) 2021  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2021                                   *
 *                                                               *             
 *****************************************************************/ 

/* Use some crafted comments with the register_models script to generate the headers
 * to integrate this model in the pv-craze routines
 * I use two flags, The first to register the model, the model name, and the model routines
 * The model is registered in models.h. The rest is registerd in model_table.h which is used
 * by pv-craze internally.
 * The sedond flag is to register the arguments of the model. This is passed in a table to pv-craze
 * but it is also exported beyond the pv_craze scope via fitmodels.h
 */
 
//MODELFLAG KH KH KHToCurve KHRMSErr KHDistCurve KHFit KHPars KHIV KHVI
//MODELPARS KH Voc Isc gamma m

curve 			KHToCurve(generic_pars P, int N, MODELSET *S);
double 			KHRMSErr(generic_pars P, curve C, MODELSET *S);
curve 			KHDistCurve(generic_pars P, curve C, MODELSET *S);
generic_pars 	KHFit(curvepars par, curve C, double *Efit, MODELSET *S);
curvepars 		KHPars(generic_pars P, MODELSET *S);
double 			KHIV(double V, generic_pars P, MODELSET *S);
double 			KHVI(double I, generic_pars P, MODELSET *S);
