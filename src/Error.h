/*
    Error/Flag handling routines
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
 
#include "ivanalyzeflags.h"

/* error macros, flag codes are unsigned int bit masks */
#define FLAG_CODE(N)  (((N)==FLAG_NOFLAG) ? 0u : (1u<<((N)-1)))
#define ISFLAG(FLAGCODE,N) (((FLAGCODE) & (FLAG_CODE(N))) ? (1) : (0))


char * Fmessage(int N);
char * Fname(int N);
int Fiserror(int N);
void F_Message(unsigned int FLAGCODE);
int IsError(unsigned int FLAGCODE);
