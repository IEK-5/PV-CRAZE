#!/usr/bin/python
import numpy as np
import pvcraze
import sys
refmodels=["APWL", "B2", "Diode", "KH", "PC", "PV-Bezier", "PV-spline", "PWL", "Spline"]
models=[]

go=1
i=0
while (go==1):
	try:
		a=pvcraze.reverselookup(i)
	except:
		go=0
	if (go==1):
		models.append(a);
	i=i+1

if (len(models)!=len(refmodels)):
	print(f"Error: number of models has changed", file=sys.stderr)
	exit(1)
# we have (at least) the same number of models


for i in range(len(refmodels)):
	try:
		a=pvcraze.lookupmodel(refmodels[i]);
	except:
		print(f"Error: model name {refmodels[i]} does not exist anymore, name change?", file=sys.stderr)
		exit(1)
# each reference model still exists


for i in range(len(models)):
	if (refmodels[i] != models[i]):
		print(f"Error: model order has changed", file=sys.stderr)
		exit(1)
# order has not changed

# test forward backward search
for i in range(len(refmodels)):
	j=pvcraze.lookupmodel(pvcraze.reverselookup(i))
	if (i!=j):
		print(f"Error: reverse lookup and lookup mismatch", file=sys.stderr)
		exit(1)

print(f"All OK", file=sys.stderr)
