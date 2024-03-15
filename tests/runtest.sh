#!/bin/bash
F=0
T=0
# color escape codes
RED='\033[0;31m'		# red
GREEN='\033[0;32m'		# green
NC='\033[0m'			# no color (default color)

function runtest() {
	T=$((T+1))
	./$i &> /dev/null
	if [ $? -gt 0 ]
	then
		echo -e "${RED}FAIL${NC}	$i"
		F=$((F+1))
	else
		echo -e "${GREEN}PASS${NC}	$i"
	fi
}
		
for i in *.m
do
	runtest $i
done
		
for i in *.py
do
	runtest $i
done
if [ "$F" -gt 0 ]
then
	echo -e "${RED}$F/$T	failed${NC}"
else
	echo -e "${GREEN}$T/$T	success${NC}"
fi
exit $F 
