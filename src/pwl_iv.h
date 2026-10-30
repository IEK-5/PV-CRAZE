/*
    Piece-Wise Linear IV model.
    Copyright (C) 2022  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2022                                   *
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
 
//MODELFLAG PWL PWL PWLToCurve PWLRMSErr PWLDistCurve PWLFit PWLPars PWLIV PWLVI
//MODELPARS PWL Voc Isc Vmpp Impp Vqmp- Iqmp- Vqmp+ Iqmp+ Gsc Roc

curve 			PWLToCurve(generic_pars P, int N, MODELSET *S);
double 			PWLRMSErr(generic_pars P, curve C, MODELSET *S);
curve 			PWLDistCurve(generic_pars P, curve C, MODELSET *S);
generic_pars 	PWLFit(curvepars par, curve C, double *Efit, MODELSET *S);
curvepars 		PWLPars(generic_pars P, MODELSET *S);
double 			PWLIV(double V, generic_pars P, MODELSET *S);
double 			PWLVI(double I, generic_pars P, MODELSET *S);
