#!/usr/bin/octave
pkg load pv-craze
% test pvcraze_modelerror on a 2x3 cell array: parameters and IV curves must
% be paired element by element (column major), not mixed up
pvcraze_quiet();
ESP=load('data/refdiodeGT.dat');
ESP=ESP(1:6,3:end);
mpar=reshape(pvcraze_modelfit(ESP, "model(PV-Bezier)"), 2, 3);
IV=pvcraze_curvemodel(mpar, 1000);
err=pvcraze_modelerror(mpar, IV);
if (any(size(err)!=[2 3]) || !(max(err(:))<1e-4))
	disp(err);
	fprintf(stderr, "Error: model parameters and IV curves are not paired element by element\n");
	exit(1);
endif
fprintf(stderr, "ALL OK\n");
