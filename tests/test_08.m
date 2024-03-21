#!/usr/bin/octave
pkg load pv-craze
TH=1e-6;
# test pvcraze_iv pvcraze_vi pvcraze_modelfit pvcraze_curvemodel
pvcraze_quiet()

makeref=true; % rebase the reference on the current state of the library

ESP=load('data/refdiodeGT.dat');
ESP=ESP(:,3:end);

% model specifiers
models={"model(PV-Bezier)","model(PV-spline)","model(KH)","model(Diode) oeps(0.6e-3) npop(200) fitit(4000)"};
% name of test
name  ={"PV-Bezier","PV-spline","KH","Diode"};

for i=1:length(models)
	% model fit
	[mpar] = pvcraze_modelfit (ESP, sprintf("%s",models{i}));
	% reconstruct IV's
	IV=pvcraze_curvemodel(mpar,100);
	Ei=0;
	Ev=0;
	N=length(IV);
	for j=1:N
		V=IV{j}(:,1);
		I=IV{j}(:,2);
		Io=pvcraze_iv(V,mpar{j});
		Vo=pvcraze_vi(I,mpar{j});
		Ei+=sqrt(mean((I-Io).^2));
		Ev+=sqrt(mean((V-Vo).^2));
	endfor
	Ei/=N
	Ev/=N
	
	if (Ei>TH)
		fprintf(stderr, "I Error: pvcraze_iv deviates for %s\n", name{i});
		exit(1);
	endif
	if (Ev>TH)
		fprintf(stderr, "I Error: pvcraze_vi deviates for %s\n", name{i});
		exit(1);
	endif
	
endfor
fprintf(stderr, "ALL OK\n");



