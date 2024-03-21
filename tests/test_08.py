#!/usr/bin/python
import numpy as np
import pvcraze
import sys
pvcraze.quiet()

# test pvcraze_iv pvcraze_vi pvcraze_modelfit pvcraze_curvemodel
ESP=np.loadtxt('data/refdiodeGT.dat')
ESP=ESP[:,2:]

# model specifiers
models=["model(PV-Bezier)","model(PV-spline)","model(KH)","model(Diode) oeps(0.6e-3) npop(200) fitit(4000)"]
# name of test
name  =["PV-Bezier","PV-spline","KH","Diode"]

TH=1e-6
for i in range(len(models)):
	# model fit
	mpar = pvcraze.modelfit (ESP, [models[i]])
	# reconstruct IV's
	IV=pvcraze.curvemodel(mpar['modelpar'],100);
	Ei=0;
	Ev=0;
	IV=IV[0]
	N=len(IV);
	mpar=mpar['modelpar'][0]
	for j in range(N):
		V=IV[j][:,0]
		I=IV[j][:,1]
		Io=pvcraze.iv(V,mpar[j])
		Vo=pvcraze.vi(I,mpar[j])
		Ei+=np.sqrt(np.mean((I-Io)**2))
		Ev+=np.sqrt(np.mean((V-Vo)**2))

	Ei=Ei/N
	Ev=Ev/N
	
	if (Ei>TH):
		print(f"I Error: pvcraze_iv deviates for {name[i]}")
		exit(1)
	if (Ev>TH):
		print(f"V Error: pvcraze_iv deviates for {name[i]}")
		exit(1)
		
print(f"All OK", file=sys.stderr)
exit(0)
