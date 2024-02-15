/*
    malloc wrapper
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
 *   Copyright (C) Dr. Bart E. Pieters 2023                      *
 *                                                               *             
 *****************************************************************/ 
#include <stdlib.h>
#include <stdio.h>


#ifdef EMALLOC_DEBUG
// note that malloc fail tests are not thread safe!
int TEST_MALLOC_FAIL=0; // simulate a malloc fail for testing purposes
#endif // EMALLOC_DEBUG


/* the malloc fail flag allows one to signal malloc errors to higher up functions */
void *Emalloc(size_t size)
{
	void *p;
#ifdef EMALLOC_DEBUG
	if (TEST_MALLOC_FAIL)
	{
		TEST_MALLOC_FAIL=0; // reset
		return NULL;
	}
#endif // EMALLOC_DEBUG
    if ( (p = malloc(size)) == NULL )
		fprintf(stderr,"Error: failed to allocate memory for PV-CRAZE\n");
	return p;
	
}
void *Ecalloc(size_t nmemb, size_t size)
{
	void *p;
#ifdef EMALLOC_DEBUG
	if (TEST_MALLOC_FAIL)
	{
		TEST_MALLOC_FAIL=0; // reset
		return NULL;
	}
#endif // EMALLOC_DEBUG
    if ( (p = calloc(nmemb,size)) == NULL )
		fprintf(stderr,"Error: failed to allocate memory for PV-CRAZE\n");
	return p;
	
}
void *Erealloc(void *ptr, size_t size)
{
	void *p;
#ifdef EMALLOC_DEBUG
	if (TEST_MALLOC_FAIL)
	{
		TEST_MALLOC_FAIL=0; // reset
		return NULL;
	}
#endif // EMALLOC_DEBUG
    if ( (p = realloc(ptr,size)) == NULL )
		fprintf(stderr,"Error: failed to re-allocate memory for PV-CRAZE\n");
	return p;
	
}
