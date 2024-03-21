#!/usr/bin/python
import numpy as np
import pvcraze
import sys
THd=2e-2
THiec=15e-2

# test pvcraze_curvemodel pvcraze_modelfit pvcraze_modelerror

def dist(par1,par2):
	Voc0=(par1[:,1]+par2[:,1])/2
	Isc0=(par1[:,2]+par2[:,2])/2
	dv=(par1[:,1]-par2[:,1])/Voc0
	E=dv*dv
	di=(par1[:,2]-par2[:,2])/Isc0
	E+=di*di
	dv=(par1[:,3]-par2[:,3])/Voc0
	di=(par1[:,4]-par2[:,4])/Isc0
	E+=dv*dv+di*di
	dv=(par1[:,5]-par2[:,5])/Voc0
	di=(par1[:,6]-par2[:,6])/Isc0
	E+=dv*dv+di*di
	dv=(par1[:,7]-par2[:,7])/Voc0
	di=(par1[:,8]-par2[:,8])/Isc0
	E+=dv*dv+di*di
	d=np.mean(np.sqrt(E))
	return d

ESP=np.loadtxt('data/refdiodeGT.dat')
# have to add type conversions and checks!
G=ESP[:,0]
T=ESP[:,1]
ESP=ESP[:,2:]


MF=pvcraze.fit_diode_GT (ESP, T, G, 25, 1000, 'C')
MP=pvcraze.translate_diode_GT (MF['modelpar'], 25.0, 1000.0, T, G, MF['alpha'], MF['Ea'], 'C')
ESPn=pvcraze.modelpars([MP])[0]


dd=dist(ESP, ESPn)
if (np.isfinite(dd)==False) or (dd>THd):
	print(f"Diode Translation Error: Deviation detected in results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)

# fit the ESP's
MF = pvcraze.fit_IEC60891_v2(ESP, T, G, 25.0, 1000.0)
# translate again
ESPn = pvcraze.translate_IEC60891_v2 (MF['esp_std'].reshape(-1,10), np.array([25.0]).reshape(-1,1), np.array([1000.0]).reshape(-1,1), T, G, MF['alpha'], MF['beta'], MF['a'], MF['Rs'], MF['kappa'])

dd=dist(ESP, ESPn)
if (np.isfinite(dd)==False) or (dd>THiec):
	print(f"IEC60891_v2 Translation Error: Deviation detected in results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)

# fit the ESP's
MF = pvcraze.fit_IEC60891_v3(ESP, T, G, 25, 1000)
# translate again
ESPn = pvcraze.translate_IEC60891_v3 (MF['esp_std'].reshape(-1,10), np.array([25.0]).reshape(-1,1), np.array([1000.0]).reshape(-1,1), T, G, MF['esp_std'][0], MF['alpha'], MF['beta'], MF['B1'], MF['B2'], MF['Rs'], MF['kappa'])

dd=dist(ESP, ESPn)
if  (np.isfinite(dd)==False) or(dd>THiec):
	print(f"IEC60891_v3 Translation Error: Deviation detected in results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)

