#!/usr/bin/python
import numpy as np
import pvcraze
import sys

# test pvcraze_curvemodel pvcraze_modelfit pvcraze_modelerror
ESP=np.loadtxt('data/refdiodeGT.dat')
ESP=ESP[:,2:]
mpar = pvcraze.modelfit (ESP, ["model(PV-Bezier)"]);
IV=pvcraze.curvemodel(mpar['modelpar'],100)[0];

# model specifiers
models=["model(PV-Bezier)","model(PV-spline)","model(KH)","model(Diode) oeps(0.6e-3) npop(200) fitit(4000)"]
# name of test
name  =["PV-Bezier","PV-spline","KH","Diode"]

TH=1e-2
for i in range(len(models)):
	# model fit
	mpar = pvcraze.modelfit (ESP, [models[i]])
	# reconstruct IV's
	dd = pvcraze.modelerror (mpar['modelpar'], IV);
	f=(dd>TH)
	r=sum(f)
	if (r>0):
		print(f"ESP Error: Deviations detected in {r}/{len(IV)} results", file=sys.stderr)
		exit(1)
	else:
		print(f"All OK", file=sys.stderr)
		exit(0)
