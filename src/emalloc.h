/*
    malloc wrapper
    Copyright (C) 2023  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich
*/
/*****************************************************************              
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
 +  IEK-5 PHOTOVOLTAIK                                           *              
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

#ifdef EMALLOC_DEBUG
/* simulate a malloc fail for testing purposes
 * Set TEST_MALLOC_FAIL before a memory allocation routine and the 
 * memory allocation will fail. The TEST_MALLOC_FAIL variable is
 * automatically reset, i.e. after the failed call the next call should 
 * succeed.
 * 
 * note that this construct is not thread safe
 */
extern int TEST_MALLOC_FAIL;
#endif // EMALLOC_DEBUG
void *Emalloc(size_t size);
void *Ecalloc(size_t nmemb, size_t size);
void *Erealloc(void *ptr, size_t size);
