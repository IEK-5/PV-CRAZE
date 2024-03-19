#!/usr/bin/octave
pkg load pv-craze
% test pvcraze_modelfit pvcraze_modelpars

makeref=false; % rebase the reference on the current state of the library
% load diode model ESP
ESP=load('data/refdiodeGT.dat');
ESP=ESP(:,3:end);

% model specifiers
models={"model(PV-Bezier)","model(PV-spline)","model(KH)","model(Diode) oeps(0.6e-3) npop(200) fitit(4000)"}
% name of test
name  ={"PV-Bezier","PV-spline","KH","Diode"}
% model specific threshold (note that fitting routines may be randomly initialized!)
TH    =1e-2;

for i=1:length(models)
	% model fit
	mpar = pvcraze_modelfit (ESP, sprintf("%s",models{i}));
	% convert to esp (some model parameters may change alot without affecting the curve so much)
	% do not trust the modelfit error, pass it through modelpars and an independent reference
	mpar=cell2mat(mpar);
	mesp=cell2mat(pvcraze_modelpars(num2cell(mpar,2)));
	
	file=sprintf("data/ref%s.dat",name{i});
	if (makeref)
		save tmp.dat mesp;
		system(sprintf("sed '/^#.*$/d' tmp.dat > %s",file));
	endif
	% load reference data
	Dr=load(file);	
	
	% compute distance from reference
	dd=Dr(:,1:end)-mesp(:,1:end);
	x=diag(1./(mean(mesp(:,1:end))));
	
	dd=dd*x;
	dd=abs(dd);
	ds=max(dd, [], 2);
	
	ii=(ds>TH);
	r=sum(ii,2);
	ii=(r>0);
	r=sum(ii);
	
	
	if (r>0)
		fprintf(stderr, "%s Error: Deviations detected in %d/%d results\n",models{i}, r, length(mesp(:,1)));
		exit(1);
	endif
endfor
fprintf(stderr, "ALL OK\n");





