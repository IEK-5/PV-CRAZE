#!/usr/bin/python
import numpy as np
import pvcraze
import sys
# translating diode model parameters in the irradiance-temperature plane
TH=1e-2;

D=np.loadtxt('data/refdiodeGT.dat')

# first create some curves as a function of G and T
# diode parameters
I0=1e-8
nkT=1.5*8.6173899875662609e-05*298.15
Iph=0.035
Rs=0.1
Rsh=700
Ea=0.8
alpha=1e-3
stdmodelpar=np.array([pvcraze.lookupmodel('Diode'), I0, nkT, Iph, Rs, Rsh,0,0,0,0,0]).reshape(1,-1)
# irradiance/temperature grid
G=np.array([100,200,400,800,1000]);
T=np.array([0,10,20,40,60,80]);

# make flat list of G,T combinations
R=np.meshgrid(T,G);
G=R[1].reshape(-1,1);
T=R[0].reshape(-1,1);


MP=pvcraze.translate_diode_GT (stdmodelpar, 25, 1000, D[:,1], D[:,0], alpha, Ea, 'C');
mesp=pvcraze.modelpars([MP])[0]

E=mesp-D[:,2:]
En=np.abs(E/np.mean(D[:,2:],axis=0))
f=(En>TH)

r=np.sum(f,axis=1)
f=(r>0)
r=sum(f)
if (r>0):
	print(f"Error: Deviations detected in {r}/{np.prod(np.size(G))} results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)
	exit(0)



