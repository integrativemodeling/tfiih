import pickle
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


CM_file = 'ContactMap_all_Matrix_CM.pkl_Matrix_CM.pkl'

data = open(CM_file)
CM = pickle.load(data)
data.close()


costum_ranges ={  # these domains are from the paper #TODO modify for human
                'kin28' :[['kin28_nlobe',1,96],['kin28_clobe',97,-1]],

                'ccl1' :[['ccl1_n',1,155],['ccl1_c',156,-1]],

                'tfb3' :[['tfb3_n',1,132],['tfb3_latch',133,-1]],
                
                'rad3' :[['rad3d1',1,245],['rad3d2',246,438],['rad3d3',439,701],['rad3d4',702,-1]],
                'ssl1' :[['ssl1d1',1,240],['ssl1d2',241,-1]],
                
                'ssl2' :[['ssl2d1',1,304],['ssl2d2',305,471],['ssl2d3',472,674],['ssl2d4',675,-1]],
           
		             
                'tfb1' :[['tfb1d1',1,107],['tfb1d2',108,254],['tfb1d3',255,457],['tfb1d4',458,-1]],
                
                'tfb2' :[['tfb2d1',1,274],['tfb2d2',275,381],['tfb2d3',382,-1]],
                
                'tfb4' :[['tfb4d1',1,230],['tfb4d2',231,-1]], # in keeping with the yeast case which is 8 residues before the VWA domain actually ends
                
                'tfb5' :[['tfb5',1,-1]]}

order = [[['kin28_nlobe','kin28_clobe'],\
	 ['ccl1_n','ccl1_c'],\
	 ['tfb3_n','tfb3_latch']],\

	 [['rad3d1','rad3d2','rad3d3','rad3d4']],\

	 [['tfb1d1','tfb1d2','tfb1d3','tfb1d4'],\
	 ['ssl1d1','ssl1d2'],\
	 ['tfb4d1','tfb4d2'],\
	 ['tfb2d1','tfb2d2','tfb2d3'],\
	 ['tfb5']],\
	 [['ssl2d1','ssl2d2','ssl2d3','ssl2d4']]]



def make_plot(groups, num_rmf, edges, out_fn='bla_temp'):
	'''
	plot the interaction matrix
	@param groups is the list of groups of domains, eg,
	[["protA_1-10","prot1A_11-100"],["protB"]....]
	it will plot a space between different groups
	'''

	ax=plt.gca()
	ax.set_aspect('equal', 'box')
	ax.xaxis.set_major_locator(plt.NullLocator())
	ax.yaxis.set_major_locator(plt.NullLocator())

	largespace=0.6
	smallspace=0.5
	squaredistance=1.0
	squaresize=0.99

	domain_xlocations={}
	domain_ylocations={}

	xoffset=squaredistance
	yoffset=squaredistance
	xlabels=[]
	ylabels=[]

	for group in groups:
		xoffset+=largespace
		yoffset+=largespace
		for subgroup in group:
			xoffset+=smallspace
			yoffset+=smallspace
			for domain in subgroup:
				domain_xlocations[domain]=xoffset
				domain_ylocations[domain]=yoffset
				#rect = plt.Rectangle([xoffset- squaresize / 2, yoffset - squaresize / 2], squaresize, squaresize,
				# facecolor=(1,1,1), edgecolor=(0.1,0.1,0.1))
				#ax.add_patch(rect)
				#ax.text(xoffset , yoffset ,domain,horizontalalignment='left',verticalalignment='center',rotation=-45.0)
				xoffset+=squaredistance
				yoffset+=squaredistance

	for edge,count in edges.iteritems():
		
		if edge[0]!=edge[1]: density=(1.0-float(count)/num_rmf)
		else: density=(0.0)
		color=(density,density,1.0)
		x=domain_xlocations[edge[0]]
		y=domain_ylocations[edge[1]]
		if x>y: xtmp=y; ytmp=x; x=xtmp; y=ytmp
		rect = plt.Rectangle([x - squaresize / 2, y - squaresize / 2], squaresize, squaresize,
		facecolor=color, edgecolor='Gray', linewidth=0.1)
		ax.add_patch(rect)
		rect = plt.Rectangle([y - squaresize / 2, x - squaresize / 2], squaresize, squaresize,
		facecolor=color, edgecolor='Gray', linewidth=0.1)

	ax.add_patch(rect)
	ax.autoscale_view()
	#plt.savefig(out_fn)
	plt.show()



lws = []
for j in order:
	lws += [0.]*(len(j)-1)
	lws += [1]
lws = lws[:-1] 

ot = sum(order,[])
subunits = costum_ranges.keys()

Domains = {}
num_max=0
for pair in CM.keys():
	for domain1 in costum_ranges[pair[0]]:
		for domain2 in costum_ranges[pair[1]]:
			d1,s1,e1 = domain1
			d2,s2,e2 = domain2
			
			# print pair, domain1, domain2 , np.average(CM[pair][s1:e1,s2:e2])
			
			Domains[ tuple(sorted([d1, d2])) ] = np.power(np.average(CM[pair][s1:e1,s2:e2]),.25)
			if np.power(CM[pair].max(),.25)>num_max: num_max = np.power(CM[pair].max(),.25)
			#Domains[ tuple(sorted([d1, d2])) ] = np.average(CM[pair][s1:e1,s2:e2])
			#if CM[pair].max()>num_max: num_max = CM[pair].max()


make_plot(order, num_max, Domains, out_fn='bla_temp')


