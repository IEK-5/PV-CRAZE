#!/usr/bin/octave
pkg load pv-craze
% test pvcraze_lookupmodel pvcraze_modelpars pvcraze_translate_diode_GT
TH=1e-2;

makeref=false; % rebase the reference on the current state of the library

% first create some curves as a function of G and T
% diode parameters
I0=1e-8
nkT=1.5*8.6173899875662609e-05*298.15
Iph=0.035
Rs=0.1
Rsh=700
Ea=0.8
alpha=1e-3
modelstd=[pvcraze_lookupmodel('Diode'), I0, nkT, Iph, Rs, Rsh,0,0,0,0,0]

% irradiance/temperature grid
G=[100,200,400,800,1000];
T=[0,10,20,40,60,80];

% make flat list of G,T combinations
[GG,TT]=meshgrid(G,T);
G=reshape(GG,prod(size(GG)),1);
T=reshape(TT,prod(size(TT)),1);


% translate std parameters
modelpar = pvcraze_translate_diode_GT (modelstd, 25, 1000, T, G, alpha, Ea, 'C');
ESP=cell2mat(pvcraze_modelpars(num2cell(modelpar,2)));
if (makeref==true)
	D=[G,T,ESP];
	save "data/refdiodeGT.dat" D;
	system("sed -i '/^#.*$/d' data/refdiodeGT.dat");
endif


Dr=load('data/refdiodeGT.dat');

dd=Dr(:,3:end)-ESP;
x=diag(1./(mean(ESP)));
dd=dd*x;
dd=abs(dd);
ii=(dd>TH);
r=sum(ii,2);
ii=(r>0);
r=sum(ii);

if (r>0)
	fprintf(stderr, "ESP Error: Deviations detected in %d/%d results\n",r, length(IV));
	exit(1);
endif
