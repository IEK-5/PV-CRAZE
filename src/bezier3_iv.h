/*
    Cubic Bézier IV model.
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
 
//MODELFLAG PVBEZIER PV-Bezier B3ToCurve B3RMSErr B3DistCurve B3Fit B3Pars B3IV B3VI
//MODELPARS PVBEZIER Voc Isc Vmpp Impp Vqmp- Iqmp- Vqmp+ Iqmp+ Gsc Roc


curve 			B3ToCurve(generic_pars P, int N, MODELSET *S);
double 			B3RMSErr(generic_pars P, curve C, MODELSET *S);
curve 			B3DistCurve(generic_pars P, curve C, MODELSET *S);
generic_pars 	B3Fit(curvepars par, curve C, double *Efit, MODELSET *S);
curvepars 		B3Pars(generic_pars P, MODELSET *S);
double			B3IV(double V, generic_pars P, MODELSET *S);
double			B3VI(double I, generic_pars P, MODELSET *S);

// IEC60891 Procedure 2 Translation of a PV-Bezier
generic_pars B3IEC60891_v2(generic_pars P, double alpha, double beta, double a, double Rs, double kappa, double T1, double T2, double G1, double G2, MODELSET *S);
generic_pars B3IEC60891_v3(generic_pars P, double alpha, double beta, double B1, double B2, double Rs, double kappa, double Vocstd, double T1, double T2, double Tstd, double G1, double G2, double Gstd, MODELSET *S);


