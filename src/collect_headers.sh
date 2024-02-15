#!/bin/bash

# Collect the headers for the library
# this script combines curve.h models.h fitmodels.h lib-pv-craze-structs.h lib-pv-craze.h
# in that order

function getapi(){
	awk '/\/\/START-HEADER-EXPORT/{flag=1;next}/\/\/END-HEADER-EXPORT/{flag=0}flag' $1
}

if [ -z "$1" ]
then
	echo gen_ivanalyzeflags.sh needs arguments!
	exit 1
fi
srcdir="$1"

FILE="pvcraze.h"
echo "// PV-CRAZE library header" >$FILE

echo "/*****************************************************************">>$FILE          
echo " *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *">>$FILE           
echo " *  IEK-5 PHOTOVOLTAIK                                           *">>$FILE          
echo " *                                                               *">>$FILE           
echo " *        ########                _   _                          *">>$FILE        
echo " *     ##########                |_| |_|                         *">>$FILE       
echo " *    ##########     ##         _ _   _ _     ___ ____ _   _     *">>$FILE        
echo " *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *">>$FILE        
echo " *   #########     #####    _  | | | | | |    | | |   | |_| |    *">>$FILE        
echo " *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *">>$FILE        
echo " *    #          ######     \___/ \___/|_____|___\____|_| |_|    *">>$FILE        
echo " *     ##      #######      F o r s c h u n g s z e n t r u m    *">>$FILE        
echo " *       ##########                                              *">>$FILE       
echo " *                                                               *">>$FILE        
echo " *   http://www.fz-juelich.de/iek/iek-5/DE/Home/home_node.html   *">>$FILE         
echo " *****************************************************************">>$FILE
echo " *                                                               *">>$FILE
echo " *    Dr. Bart E. Pieters 2023                                   *">>$FILE
echo " *                                                               *">>$FILE      
echo " *****************************************************************/">>$FILE
echo "#include <stdbool.h>">>$FILE
getapi $srcdir/curve.h >>$FILE
cat $srcdir/models.h >>$FILE
getapi $srcdir/fitmodels.h >>$FILE
getapi $srcdir/lib-pv-craze-structs.h >>$FILE
getapi $srcdir/lib-pv-craze.h >>$FILE


           
