import sys
from operator import itemgetter
import re


xlink_file = 'Ranish_Kornberg_thiih_xlinks_human.txt'
viols_file = 'violation_counts_alltopscoring_34A.txt'
rest_viol = 'restraint_violations_34A.txt'
DIC = {
	'P06839': ['P06839', 'RAD3_HUMAN', 'Rad3', 'YER171W', 760], #TODO
	'Q04673': ['Q04673', 'SSL1_HUMAN', 'Ssl1', 'YLR005W', 395], 
	'Q12004': ['Q12004', 'TFB4_HUMAN', 'Tfb4', 'YPR056W', 308], 
	'Q3E7C1': ['Q3E7C1', 'TFB5_HUMAN', 'Tfb5', 'YDR079C-A', 71], 
	'Q00578': ['Q00578', 'SSL2_HUMAN', 'Ssl2', 'YIL143C', 782], 
	'Q02939': ['Q02939', 'TFB2_HUMAN', 'Tfb2', 'YPL122C', 462], 
	'P32776': ['P32776', 'TFB1_HUMAN', 'Tfb1', 'YDR311W', 548],
	'P37366': ['P37366', 'CCL1_HUMAN', 'Ccl1', 'YPR025C', 323], 
	'P06242': ['P06242', 'KIN28_HUMAN', 'Kin28', 'YDL108W', 346], 
	'Q03290': ['Q03290', 'TFB3_HUMAN', 'Tfb3', 'YDR460W',309]
}


Proteins = ['P06242','P37366','Q03290','P06839','P32776','Q04673','Q12004','Q02939','Q3E7C1','Q00578']

Lenghts = [(DIC[i][2].lower(), DIC[i][-1]) for i in Proteins]
L = Lenghts #sorted(Lenghts, key=itemgetter(1), reverse=True)





data = open(xlink_file)
D = data.readlines()
data.close()

Xlinks = []
for d in D:
	p1,p2,i1,i2,sr = d.strip().split()
	Xlinks.append((p1,p2,i1,i2,sr))

data = open(rest_viol)
D = data.readlines()
data.close()

resviol = []
for vr in D:
	v = vr.strip().split()[0].split('_')
	if v[1]=="Score":
	    continue
	r1,p1 = v[-2].split(':')
	r2,p2 = v[-1].split(':')
	resviol.append((p1,p2,r1,r2))



X,Y=[],[]
for link in resviol:
	p1,p2,r1,r2=link
	r1=int(r1)
	r2=int(r2)
	x = 0
	for i in L:
		if i[0]==p1: break
		else: x+=i[1]
	x += r1
	
	y = 0
        for i in L:
                if i[0]==p2: break
                else: y+=i[1]
        y += r2

	X.append(x)
	Y.append(y)


	



import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib
import numpy as np



fig = plt.figure(figsize = (11,11))
ax = fig.add_subplot(111)

t = 0
for l in L:
	t += l[1]
	ax.plot([t,t], [0,sum([i[1] for i in L])], 'k-',lw=.5, alpha=.5)
	ax.plot([0,sum([i[1] for i in L])], [t,t], 'k-',lw=.5, alpha=.5)


Xt,Yt = [],[]
for i in xrange(len(X)):
	if X[i]>Y[i]:
		Xt.append(X[i])
		Yt.append(Y[i])
	else:
		Xt.append(Y[i])
		Yt.append(X[i])




data = open(viols_file)

#VR = data.readlines()
VR = eval(data.read())
data.close()
Xlinkst,DICT = [],{}
for xl in Xlinks:
	Xlinkst.append(tuple(list(xl)[:-1]))
	DICT[tuple(list(xl)[:-1])] = xl
Xlinks= Xlinkst

VRC = {}
for vr in VR:
	#if 'SimplifiedCrossLinkMS_Score_standar' not in vr: continue
	#if 'ssl2' in vr and 'tfb2' in vr: print vr, VR[vr]
	v = vr.split('_')
	r1,p1 = v[-2].split(':')
	r2,p2 = v[-1].split(':')
	if (p1,p2,r1,r2) in resviol: VRC[(p1,p2,r1,r2)] = VR[vr] #int(vr.strip().split()[1])
	else: print '!!!!!!',(p1,p2,r1,r2)
C = []
for xl in resviol:
	if xl in VRC:
		C.append(VRC[xl])
		# print all violated crosslinks
		#print DICT[xl][-1] #,VRC[xl]
		if VRC[xl]>68:  # violated in more than half of the structures
		  #print DICT[xl][-1] #,VRC[xl]
		  continue
		else: # violated in less than half of the structures
		  print DICT[xl][-1] #,VRC[xl]
	else:
		C.append(0)
		print DICT[xl][-1] #,0

#print len(C), len([i for i in C if i<max(C)/2.])

#ax.scatter(Xt,Yt,c=['green' if i<max(C)/2. else 'red' for i in C],s=[i*5+40 for i in C],linewidth=0,alpha=0.5)
ax.scatter(Xt,Yt,c=['green' if i<=max(C)/2. else 'red' for i in C],s=133,linewidth=0,alpha=0.5)
plt.ylim((0,sum([i[1] for i in L])))
plt.xlim((0,sum([i[1] for i in L])))
plt.savefig('finalFigures/XL_violations_human_withkinase.pdf', format='PDF')
plt.show()



