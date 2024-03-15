#!/usr/bin/python
import numpy as np
import math 
import pvcraze
import sys
# test analyze with some bad curves
IV=[];
IV.append(np.ones((50,2)))
t=np.linspace(1,100,num=100);
iv=np.concatenate((t.reshape(-1,1),np.cos(t*math.pi/50).reshape(-1,1)),axis=1)
IV.append(iv)
iv=np.concatenate((t.reshape(-1,1),np.sin(t*math.pi/50).reshape(-1,1)),axis=1)
IV.append(iv)
iv=np.concatenate((t.reshape(-1,1),(1-np.exp(t*math.pi/50)).reshape(-1,1)),axis=1)
IV.append(iv)
iv=np.concatenate((t.reshape(-1,1),np.exp(t*math.pi/50).reshape(-1,1)),axis=1)
IV.append(iv)
print(IV)
B=pvcraze.analyze(IV)
ref=np.array([256,0,2097264,136356976,136356976]).reshape(-1,1)
r=B['flags']-ref
f=(r>0.5)
r=sum(f)

if (r>0):
	print(f"Error: Deviations detected in {r}/{len(IV)} results", file=sys.stderr)
	exit(1)
else:
	print(f"All OK", file=sys.stderr)
	exit(0)
