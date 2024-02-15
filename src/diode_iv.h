
/*
    one diode IV model.
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
 
//MODELFLAG DIODE Diode DiodeToCurve DiodeRMSErr DiodeDistCurve DiodeFit DiodePars
//MODELPARS DIODE I0 nkT Iph Rs Rsh
 
// disabled model, remove X if you want to re-enable
//MXODELFLAG DIODE DiodeQ DiodeToCurve DiodeRMSErr DiodeDistCurve DiodeFitGqmp DiodePars
//MXODELPARS DIODE I0 nkT Iph Rs Rsh

// disabled model, remove X if you want to re-enable
//MXODELFLAG DIODEDE DioDE DiodeToCurve DiodeRMSErr DiodeDistCurve DiodeFitDE DiodePars
//MXODELPARS DIODEDE I0 nkT Iph Rs Rsh

curve 			DiodeToCurve(generic_pars P, int N, MODELSET *S);
double 			DiodeRMSErr(generic_pars P, curve C, MODELSET *S);
curve 			DiodeDistCurve(generic_pars P, curve C, MODELSET *S);
generic_pars 	DiodeFit(curvepars par, curve C, double *Efit, MODELSET *S);
generic_pars 	DiodeFitGqmp(curvepars par, curve C, double *Efit, MODELSET *S);
// generic_pars 	DiodeFit_red(curvepars par, curve C, double *Efit);
// generic_pars 	DiodeFitDE(curvepars par, curve C, double *Efit);
curvepars 		DiodePars(generic_pars P, MODELSET *S);


