#!/usr/bin/octave
pkg load pv-craze
THd=2e-2;
THiec=15e-2;
# test pvcraze_fit_diode_GT
pvcraze_verbose_setreport(false)

makeref=false; % rebase the reference on the current state of the library
% load diode model ESP
ESP=load('data/refdiodeGT.dat');
G=ESP(:,1);
T=ESP(:,2);
ESP=ESP(:,3:end);
function d=dist(par1, par2)
	Voc0=(par1(:,1)+par2(:,1))./2;
	Isc0=(par1(:,2)+par2(:,2))./2;
	dv=(par1(:,1)-par2(:,1))./Voc0;
	E=dv.*dv;
	di=(par1(:,2)-par2(:,2))./Isc0;
	E+=di.*di;	
	dv=(par1(:,3)-par2(:,3))./Voc0;
	di=(par1(:,4)-par2(:,4))./Isc0;
	E+=dv.*dv+di.*di;
	dv=(par1(:,5)-par2(:,5))./Voc0;
	di=(par1(:,6)-par2(:,6))./Isc0;
	E+=dv.*dv+di.*di;
	dv=(par1(:,7)-par2(:,7))./Voc0;
	di=(par1(:,8)-par2(:,8))./Isc0;
	E+=dv.*dv+di.*di;
	d=mean(sqrt(E));
endfunction

% fit the ESP's
[modelstd, alpha, Ea] = pvcraze_fit_diode_GT(ESP, T, G, 25,1000, 'C');

% translate again
modelpar = pvcraze_translate_diode_GT (modelstd, 25, 1000, T, G, alpha, Ea, 'C');
% recompute the ESP
ESPn=cell2mat(pvcraze_modelpars(num2cell(modelpar,2)));
dd=dist(ESP, ESPn);

if (dd>THd)
	fprintf(stderr, "Diode Translation Error: Deviation detected in results\n");
	exit(1);
endif

% fit the ESP's
[ESPstd, alpha, beta, a, Rs, kappa] = pvcraze_fit_IEC60891_v2(ESP, T, G, 25, 1000);
% translate again
ESPn = pvcraze_translate_IEC60891_v2 (ESPstd, 25, 1000, T, G, alpha, beta, a, Rs, kappa);
dd=dist(ESP, ESPn);
if (dd>THiec)
	fprintf(stderr, "IEC60891_v2 Translation Error: Deviation detected in results\n");
	exit(1);
endif

% fit the ESP's
[ESPstd, alpha, beta, B1, B2, Rs, kappa] = pvcraze_fit_IEC60891_v3 (ESP, T, G, 25, 1000);
% translate again
ESPn = pvcraze_translate_IEC60891_v2 (ESPstd, 25, 1000, T, G, alpha, beta, a, Rs, kappa);
dd=dist(ESP, ESPn);

if (dd>THiec)
	fprintf(stderr, "IEC60891_v3 Translation Error: Deviation detected in results\n");
	exit(1);
endif
