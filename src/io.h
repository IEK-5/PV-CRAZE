/*
    io module for the commandline PV-CRAZE interface
    Copyright (C) 2020  B. E. Pieters, 
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
 *    Dr. Bart E. Pieters 2020                                   *
 *                                                               *             
 *****************************************************************/   
typedef struct {
	double **data;
	int n,m;
} MatrixData;

void FreeMat(MatrixData *M);
MatrixData ReadMatrix(char *fn, char sep);
MatrixData ReadSpecMatrix(char *fn, char sep, const char *cspec, const char *rspec);

curve MatrixRows2Curve(MatrixData *V, MatrixData *I, int rv, int ri);
curve MatrixCols2Curve(MatrixData *V, MatrixData *I, int cv, int ci);

void WriteExtrParMatrix(char *fn, MatrixData par);
void WriteParMatrix(char *fn, MatrixData par, MODEL_LIST L);
void WriteModelParsMatrix(char *fn, MatrixData par, MODEL_LIST L);
void WriteMatrix(char *fn, MatrixData M);

MatrixData ReadModelParMatrix(char *fn, char sep, int *M);
generic_pars GetNextPars(MatrixData par, int r /* row */, int *c /* column */);

curve FileToCurve(char *fn, int colx, int coly);
void WriteCurve(char *fn, curve C);
