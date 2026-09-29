#!/usr/bin/octave
pkg load pv-craze
% test analyze with some bad curves
IV=cell(1,1);
IV=[IV;ones(5,1)];
IV=[IV;ones(50,2)];
t=1:100;
t=t';
iv=[t, cos(t*pi/50)];
IV=[IV;iv];
iv=[t, sin(t*pi/50)];
IV=[IV;iv];
iv=[t, 1-exp(t*pi/50)];
IV=[IV;iv];
iv=[t, exp(t*pi/50)];
IV=[IV;iv];


[ESP,F]=pvcraze_analyze(IV);
% the last three currents increase with voltage: not an IV curve (NOISY_IV)
ref=[128;128;256;0;2147483648;2147483648;2147483648];
E=ref-double(F(:,1));
ii=(abs(E)>0.5);
r=sum(ii);
if (r>0)
	fprintf(stderr, "Error: Deviations detected in %d/%d results\n",r, length(IV));
	exit(1);
else
	fprintf(stderr, "ALL OK\n");
	exit(0);
endif




