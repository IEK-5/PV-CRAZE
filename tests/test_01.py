#!/usr/bin/python
import numpy as np
import pvcraze
import sys
# Simple test extracting ESP
TH=1e-2;

with open('data/IVs.csv', 'r') as f:
	Lines = f.readlines()
	N=len(Lines)-1
	IV=[]
	count=0
	s=1;
	# somehow this is ugly and unreadable compared to the good old scanf
	# one day I'll learn python...
	for line in Lines:
		if (s>0):
			s=s-1
		else:
			substrings = line.split(',')
			n=int(substrings[0])
			V=np.empty(n)
			I=np.empty(n)
			nn=0
			for v in substrings[1:1+n]:
				V[nn]=float(v)
				nn=nn+1
			nn=0
			for i in substrings[1+n:1+2*n]:
				I[nn]=float(i)
				nn=nn+1
			IV.append(np.concatenate((V.reshape(-1,1),I.reshape(-1,1)),axis=1))
			count=count+1
f.close()
B=pvcraze.analyze(IV)

Dr=np.loadtxt('data/refesp1.dat', comments='#')

E=B['esp']-Dr
En=np.abs(E/np.mean(B['esp'],axis=0))
f=(En>TH)

r=np.sum(f,axis=1)
f=(r>0)
r=sum(f)
if (r>0):
	print(f"ESP Error: Deviations detected in {r}/{len(IV)} results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)
	exit(0)


# test FLAGS
Fr=np.loadtxt('data/refflags1.dat');
r=B['flags'][:,1]-Fr[:,1];
f=(r>0.5)
r=sum(f)

if (r>0):
	print(f"FLAG Error: Deviations detected in {r}/{len(IV)} results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)
	exit(0)







