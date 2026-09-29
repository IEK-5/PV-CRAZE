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
		
# python tests (pytest, see tests/python)
T=$((T+1))
if python3 -m pytest -q python &> /dev/null
then
	echo -e "${GREEN}PASS${NC}	python (pytest)"
else
	echo -e "${RED}FAIL${NC}	python (pytest), run 'python3 -m pytest tests/python' for details"
	F=$((F+1))
fi
if [ "$F" -gt 0 ]
then
	echo -e "${RED}$F/$T	failed${NC}"
else
	echo -e "${GREEN}$T/$T	success${NC}"
fi
exit $F 
