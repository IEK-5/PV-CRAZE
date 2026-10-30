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
   
#include <stdio.h>
#include "ivanalyzemessages.h"
#include "Error.h"

char * Fmessage(int N)
{
	if ((N>=0)&&(N<NIVF))
		return iv_flag_messages[N];
	else
		return NULL;
}

char * Fname(int N)
{
	if ((N>=0)&&(N<NIVF))
		return iv_flag_names[N];
	else
		return NULL;
}

/* 1 if flag number N is an error, 0 if it is a warning (or out of range) */
int Fiserror(int N)
{
	if ((N>0)&&(N<NIVF))
		return iv_flag_is_error[N];
	else
		return 0;
}

void F_Message(unsigned int FLAGCODE)
{
	if (FLAGCODE)
	{
		int i;
		for (i=1;i<NIVF;i++)
			if (ISFLAG(FLAGCODE,i))
				fprintf(stderr,"\t%s\n",  Fmessage(i));
	}	
}

int IsError(unsigned int FLAGCODE)
{
	return (FLAGCODE&IV_ERROR_MASK)!=0;
}
