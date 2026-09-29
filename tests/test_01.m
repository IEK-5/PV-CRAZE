#!/usr/bin/octave
fn="data/IVs.csv";
pkg load pv-craze
TH=1e-2;
# test basic analyze functionality
pvcraze_verbose_setreport(false)

makeref=false; % rebase the reference on the current state of the library
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
IV=ReadIV(fn);
[ESP,F]=pvcraze_analyze(IV);
if (makeref)
	save data/refesp1.dat ESP;
	save data/refflags1.dat F;
	system("sed -i '/^#.*$/d' data/refesp1.dat");
	system("sed -i '/^#.*$/d' data/refflags1.dat");
endif
# test ESP	
Dr=load('data/refesp1.dat');
dd=Dr(:,1:10)-ESP;
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

# test FLAGS
Fr=load('data/refflags1.dat');
r=double(F(:,1))-Fr(:,1);
ii=(abs(r)>0.5);
r=sum(ii);

if (r>0)
	fprintf(stderr, "FLAGS Error: Deviations detected in %d/%d results\n",r, length(IV));
	exit(1);
else
	fprintf(stderr, "ALL OK\n");
endif



