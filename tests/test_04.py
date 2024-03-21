#!/usr/bin/python
import numpy as np
import pvcraze
import sys
# test pvcraze_modelfit pvcraze_modelpars
ESP=np.loadtxt('data/refdiodeGT.dat')
ESP=ESP[:,2:]
# models
models=["model(PV-Bezier)","model(PV-spline)","model(KH)","model(Diode) oeps(0.6e-3) npop(200) fitit(4000)"]
# name of test
name  =["PV-Bezier","PV-spline","KH","Diode"]

TH=1e-2
file="data/ref{name}.dat"
for i in range(len(models)):
	# model fit
	mpar = pvcraze.modelfit (ESP, [models[i]])
	# convert to esp (some model parameters may change alot without affecting the curve so much)
	mesp=pvcraze.modelpars(mpar['modelpar'])[0];
	
	D=np.loadtxt(file.format(name=name[i]))
	# compute distance from reference
	dd=D-mesp
	x=np.diag(1/(np.mean(mesp,axis=0)))
	
	dd=np.matmul(dd,x)
	dd=np.abs(dd)
	ds=np.max(dd, axis=1)
	f=(ds>TH)
	r=sum(f)
	if (r>0):
		print(f"Error: Deviations detected in {r}/{np.prod(np.size(G))} results", file=sys.stderr)
		exit(1)
	else:
		print(f"All OK", file=sys.stderr)
	
	
