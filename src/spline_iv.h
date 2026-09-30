/*
    Cubic Spline IV model.
    Copyright (C) 2024  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2024                                   *
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

// provide two models, the plain Spline and the PV-Spline (monotone spline by inserting dummy knots)
//MODELFLAG SPLINE Spline SplineToCurve SplineRMSErr SplineDistCurve PlainSplineFit SplinePars SplineIV SplineVI
//MODELPARS SPLINE Voc Isc Vmpp Impp Vqmp- Iqmp- Vqmp+ Iqmp+ Gsc Roc
//MODELFLAG PVSPLINE PV-spline SplineToCurve SplineRMSErr SplineDistCurve PVSplineFit SplinePars SplineIV SplineVI
//MODELPARS PVSPLINE Voc Isc Vmpp Impp Vqmp- Iqmp- Vqmp+ Iqmp+ Gsc Roc

curve 			SplineToCurve(generic_pars P, int N, MODELSET *S);
double 			SplineRMSErr(generic_pars P, curve C, MODELSET *S);
curve 			SplineDistCurve(generic_pars P, curve C, MODELSET *S);
generic_pars 	PlainSplineFit(curvepars par, curve C, double *Efit, MODELSET *S);
generic_pars 	PVSplineFit(curvepars par, curve C, double *Efit, MODELSET *S);
curvepars 		SplinePars(generic_pars P, MODELSET *S);
double 			SplineIV(double V, generic_pars P, MODELSET *S);
double 			SplineVI(double I, generic_pars P, MODELSET *S);

// IEC60891 Procedure 2 Translation of an ESP
// The interpolation model is selected by P.TYPE: M_PVSPLINE (or M_SPLINE) or M_PVBEZIER
generic_pars IEC60891_v2(generic_pars P, double alpha, double beta, double a, double Rs, double kappa, double T1, double T2, double G1, double G2, MODELSET *S);
generic_pars IEC60891_v3(generic_pars P, double alpha, double beta, double B1, double B2, double Rs, double kappa, double Vocstd, double T1, double T2, double Tstd, double G1, double G2, double Gstd, MODELSET *S);
// fit the translation coefficients, M is the interpolation model (M_PVSPLINE or M_PVBEZIER)
generic_pars IEC60891Fit_v2(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double *alpha, double *beta, double *a, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S);
generic_pars IEC60891Fit_v3(curvepars *C, double *G, double *T, int N, double Tstd, double Gstd, double *alpha, double *beta, double *B1, double *B2, double *Rs, double *kappa, MODEL_TYPE M, MODELSET *S);
// translate a single point (V1, I1) of a curve (v2: Voc1 is the Voc of that curve)
void IEC_v2(double V1, double I1, double Voc1, double alpha, double beta, double a, double Rs, double kappa,
            double T1, double T2, double G1, double G2, double *V2o, double *I2o);
void IEC_v3(double V1, double Vs, double I1, double dIdV1,
            double G1, double G2, double Gs, double T1, double T2, double Ts,
            double alpha, double beta, double B1, double B2, double Rs, double kappa,
            double *V2o, double *I2o, double *dIdV2o);

// test curvepars
void Convexity(curvepars par, double *convex1, double *convex2, double *convex3, double *convex4);
