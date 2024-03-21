#!/usr/bin/octave
pkg load pv-craze
refmodels={"APWL", "B2", "Diode", "KH", "PC", "PV-Bezier", "PV-spline", "PWL", "Spline"};

go=1;
i=0;
models={};
while (go==1)
	try a=pvcraze_reverselookup(i);
	catch
		go=0;
	end_try_catch
	if (go)
		models(i+1)=a;
	endif
	i++;
endwhile

if (length(models)!=length(refmodels))
	fprintf(stderr, "Error: number of models has changed\n");
	exit(1);
endif	
% we have the same number of models
for i=1:length(refmodels)
	try a=pvcraze_lookupmodel(refmodels{i});
	catch
		fprintf(stderr, "Error: model names %s does not exist anymore, name change?\n", refmodels{i});
		exit(1);
	end_try_catch
endfor
% each reference model still exists

for i=1:length(models)	
	if (strcmp(refmodels{i},models{i})!=1)
		fprintf(stderr, "Error: model order has changed\n");
		exit(1);
	endif		
endfor
% order has not changed

% test forward backward search
for i=0:length(refmodels)-1	
	j=pvcraze_lookupmodel(pvcraze_reverselookup(i));
	if (i!=j)
		fprintf(stderr, "Error: reverse lookup and lookup mismatch\n");
		exit(1);
	endif
endfor
fprintf(stderr, "All OK!\n");
