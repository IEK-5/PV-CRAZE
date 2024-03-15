#!/usr/bin/octave
fn="data/IVs.csv";
pkg load pv-craze
TH=1e-2;
# test pvcraze_curvemodel pvcraze_modelfit pvcraze_modelerror
pvcraze_verbose_setreport(false)

makeref=true; % rebase the reference on the current state of the library
function IV=ReadIV(fn)
	f=fopen(fn,'r');
	IV={};
	while(!feof(f))
		l=fgets(f);
		if (l(1)!='#')		
			[v,c]=sscanf(l,"%e,",Inf);
			n=v(1);
			v=v(2:end);
			iv=reshape(v,n,2);
			IV=[IV,iv];
		endif
	endwhile
endfunction

ESP=load('data/refdiodeGT.dat');
ESP=ESP(:,3:end);
[mpar] = pvcraze_modelfit (ESP, "model(PV-Bezier)");
IV=pvcraze_curvemodel(mpar,100);

% model specifiers
models={"model(PV-Bezier)","model(PV-spline)","model(KH)","model(Diode) oeps(0.6e-3) npop(200) fitit(4000)"}
% name of test
name  ={"PV-Bezier","PV-spline","KH","Diode"}

for i=1:length(models)
	% model fit
	[mpar] = pvcraze_modelfit (ESP, sprintf("%s",models{i}));
	% reconstruct IV's
	dd = pvcraze_modelerror (mpar, IV);
	ii=(dd>TH);
	r=sum(ii,2);
	ii=(r>0);
	r=sum(ii);

	if (r>0)
		fprintf(stderr, "ESP Error: Deviations detected in %d/%d results\n",r, length(IV));
		exit(1);
	endif
endfor
fprintf(stderr, "ALL OK\n");



