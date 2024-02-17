/*
    Curve analysis routines. Extract ESP and Flag faulty curves
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

#ifndef _CURVE_H
#define _CURVE_H
//START-HEADER-EXPORT mark the start of a block that describes the C api
typedef struct {
	double x,y;
} point;

typedef struct {
	point *P;
	int N;
} curve;

typedef struct {
	int SpikeFilter;
	int Identicals;
	int Resample;			/* if 0, do nothing, otherwise it is the number of resampling points */
	double IdentEps;
	double MinRelRange;
	double MPP_f_V, MPP_f_I;	/* MPP is fitted to a range left and right of MPP. To the left we go down to (1-MPP_f_V)*Vmpp. To ghe right we go to (1-MPP_f_I)*Impp */
	double OC_f_I, SC_f_V;	/* slopes at sc and oc are fitted to a range. Slope sc from V=0 to (SC_f_V)*Vmpp. Slope oc from I=0 to (OC_f_I)*Impp */
	char am;				/* mpp-method 
							 * am='f': no noise assumption, fast but not robust
							 * am='r': robust
							 * am='s': Savitsky-Golay
							 * am='p': parametric */
} analyze_settings;

typedef struct {
	double Voc, Isc;   // sc oc
	double Vmpp, Impp; // mpp
	double Vqmp1, Iqmp1; // meta mpp1
	double Vqmp2, Iqmp2; // meta mpp2
	double Gsc, Roc;
} curvepars;

//END-HEADER-EXPORT mark the end of a block that describes the C api
extern curvepars nullpars;

double dist2(point p1, point p2, double norm_x, double norm_y);
curve PreprocessCurve(curve C, int ext);
curve Resample(curve C, int N, double Voc, double Isc);
/* character am:
 * am='f': no noise assumption, fast but not robust
 * am='r': robust
 * am='s': Savitsky-Golay
 * am='p': parametric */
curvepars AnalyzeLSQ(curve C, curve *Cout, int *FLAGCODE, analyze_settings settings);
int FlagMono(curvepars par);
double CurveParDist(curvepars p, curvepars ref, double Voc0, double Isc0);
# endif
