from __future__ import print_function, absolute_import
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import pickle
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import KMeans
import os,sys,string

run_number = 1
inpPrefix='cluster_run%d' % run_number

data = open(inpPrefix+'.pkl', 'rb')
K,M = pickle.load(data)
data.close()

print(M.max())

# result from hierarchical clustering
ax = pl.subplot(311)
R = dendrogram(linkage(M, method='complete'), color_threshold=100)

# heat map corresponding to hierarchical clustering
ax=pl.subplot(312)
l=ax.imshow(M[R['leaves']][:,R['leaves']], interpolation='nearest')
pl.colorbar(l)

km = KMeans(3) #TODO check based on hierarchical clustering and fit! 
clusters = km.fit_predict(M)

# heat map corresponding to kmeans clustering
R = list(np.argsort(clusters))
ax=pl.subplot(313)
l=ax.imshow(M[R][:,R], interpolation='nearest')

clusters=list(clusters)
Clusters = {}
for c in set(clusters):
	print(c,clusters.count(c))
	r = list(np.argwhere(clusters==c).T[0])
	m = M[r][:,r] 
	cls= [K[i] for i in r]
	if clusters.count(c)>1:
	  print(np.sum(m) / (len(m)**2 - len(m)))
	  Clusters[c] = [r,cls,np.sum(m) / (len(m)**2 - len(m))]	
	else:
	  Clusters[c] = [r,cls,0.0]
	  print("Single member cluster: how on earth do you expect me to calculate precision? Are you in your senses?")
pl.savefig('clusters.pdf')


# will this give results equivalent to just clustering 
out = open(inpPrefix+'.out','w')
out.write(str(Clusters))
out.close
