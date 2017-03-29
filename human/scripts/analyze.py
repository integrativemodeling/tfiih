#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
import os,sys
import IMP.rmf
import RMF

import IMP.pmi.restraints as restraints
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output
import IMP.pmi.analysis as analysis


nframes=5000 
xlinkAmplitude=17  

ncycl=100  #number of montecarlo steps cycles
rbmaxtrans=0.7
fbmaxtrans=0.7

outputobjects=[]
sampleobjects=[]

'''
#####################################################
#read EM maps and centroid coordinates
#####################################################
map_kinase = IMP.em.read_map('kinase.mrc', IMP.em.MRCReaderWriter())
map_rad3 = IMP.em.read_map('rad3.mrc', IMP.em.MRCReaderWriter())
map_ssl2 = IMP.em.read_map('ssl2.mrc', IMP.em.MRCReaderWriter())
map_tfb = IMP.em.read_map('tfb.mrc', IMP.em.MRCReaderWriter())
'''

#####################################################
#create hierarchies and rigid bodies and flexible parts
#####################################################

m=IMP.Model()
simo=representation.SimplifiedModel(m, upperharmonic=False, disorderedlength=True)

# --- Get atomic models available
# Sub component 1: Kinase 
simo.add_component_name("ccl1")
simo.add_pdb_and_intervening_beads("ccl1",'KIN28_CCL1.pdb',"B",resolutions=[1,30],resrange=(1,323),beadsize=30,color=60./360,attachbeads=True)
simo.setup_component_sequence_connectivity("ccl1", resolution=30)

simo.add_component_name("kin28")
simo.add_pdb_and_intervening_beads("kin28",'KIN28_CCL1.pdb',"A",resolutions=[1,30],resrange=(1,346), beadsize=30,color=185./360,attachbeads=True)
simo.setup_component_sequence_connectivity("kin28", resolution=30)

simo.add_component_name("tfb3")
simo.add_pdb_and_intervening_beads("tfb3",'TFB3.pdb',"A",resolutions=[1,30],resrange=(1,309), beadsize=30,color=285./360,attachbeads=True)
simo.setup_component_sequence_connectivity("tfb3", resolution=30)

# Sub component 2: Rad3
simo.add_component_name("rad3")
simo.add_pdb_and_intervening_beads("rad3",'RAD3.pdb'," ",resolutions=[1,30],resrange=(1,760), beadsize=30,color=0.4,attachbeads=True)
simo.setup_component_sequence_connectivity("rad3", resolution=30)

# Sub component 3: Ssl2 and Tfiicore
simo.add_component_name("ssl2")
simo.add_pdb_and_intervening_beads("ssl2",'SSL2.pdb'," ",resolutions=[1,30],resrange=(1,782), beadsize=30,color=0.0,attachbeads=True)
simo.setup_component_sequence_connectivity("ssl2", resolution=30)
  
simo.add_component_name("tfb1")
simo.add_pdb_and_intervening_beads("tfb1",'TFB1.pdb'," ",resolutions=[1,30],resrange=(1,548), beadsize=30,color=1.0,attachbeads=True)
simo.setup_component_sequence_connectivity("tfb1", resolution=30)

simo.add_component_name("tfb2")
simo.add_pdb_and_intervening_beads("tfb2",'TFB2_TFB5.pdb',"A",resolutions=[1,30],resrange=(1,462), beadsize=30,color=0.9,attachbeads=True)
simo.setup_component_sequence_connectivity("tfb2", resolution=30)

simo.add_component_name("tfb4")
simo.add_pdb_and_intervening_beads("tfb4",'TFB4.pdb'," ",resolutions=[1,30],resrange=(1,308), beadsize=30,color=0.8,attachbeads=True)
simo.setup_component_sequence_connectivity("tfb4", resolution=30)

simo.add_component_name("tfb5")
simo.add_pdb_and_intervening_beads("tfb5",'TFB2_TFB5.pdb',"B",resolutions=[1,30],resrange=(1,71), beadsize=30,color=0.,attachbeads=True)
simo.setup_component_sequence_connectivity("tfb5", resolution=30)

simo.add_component_name("ssl1")
simo.add_pdb_and_intervening_beads("ssl1",'SSL1.pdb'," ",resolutions=[1,30],resrange=(1,395), beadsize=30,color=0.6,attachbeads=True)
simo.setup_component_sequence_connectivity("ssl1", resolution=30)

# --- set rigid bodies
# Sub component 1: Kinase
emxk,emyk,emzk = 80.35, -23.66, 73.14   
simo.set_rigid_bodies(["ccl1","kin28"],(emxk,emyk,emzk))

# Sub component 2: Tfb3
emxk,emyk,emzk = 17.25, -41.81, -21.93   
simo.set_rigid_bodies(["tfb3"],(emxk,emyk,emzk))

# Sub component 2: Rad3
emxr,emyr,emzr = 37.25, 14.39, -1.38
simo.set_rigid_bodies(["rad3"],(emxr,emyr,emzr))

# Sub component 3: Ssl2 + Tfiihcore
emxt,emyt,emzt = -29.65, -1.27, -0.32
simo.set_rigid_bodies([("tfb2"),("tfb5")], (emxt,emyt,emzt))
simo.set_rigid_bodies(["tfb1"], (emxt,emyt,emzt))
simo.set_rigid_bodies(["tfb4"],(emxt,emyt,emzt))
simo.set_rigid_bodies(["ssl1"],(emxt,emyt,emzt))
simo.set_rigid_bodies(["ssl2"],(emxt,emyt,emzt))

# --- set simulation params
simo.set_floppy_bodies()
d=simo.get_particles_to_sample()
simo.set_rigid_bodies_max_trans(rbmaxtrans)
simo.set_floppy_bodies_max_trans(fbmaxtrans)

#--- re-orient initial orientation only
simo.shuffle_configuration(translate=False)

prot=simo.get_hierarchy()
outputobjects.append(simo)
sampleobjects.append(simo)

##simo.draw_hierarchy_composition()
#####################################################
#analyze RMFs
#####################################################

import pickle
import numpy as np
import glob
import pylab as pl

allStatFiles = glob.glob('/salilab/park1/shruthi/tfiih/tfb3changed_humanCAK/stat/stat*.dat') 

print allStatFiles
scoreCutoff=-4230 # for top 10% or about 100 models #TODO change this as needed

#allStatFiles =glob.glob('/salilab/park1/shruthi/tfiih/humanAnalysis/stats/stat_390.dat')
files=allStatFiles

"""
#####################################################
### --- Get list of xlinks
#####################################################
for z, fil in enumerate(files[:1]):

    print fil
    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    
    # get the names of crosslinks
    count=0
    for d in H:
      if H[d].startswith('SigmoidCrossLinkMS_Distance'):
	  print H[d]+'\t34'
	  count=count+1
    print count
exit()
"""

"""
#####################################################
### --- Get stat files and models belonging to half of best cluster
#####################################################

data = open('cluster_run'+seed+'.pkl')

data.close()

exit()
"""

"""
#####################################################
### --- Test RMF hierarchy
#####################################################

for z, fil in enumerate(files[:1]):

    print fil
    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    
    H = dict([(H[c],c) for c in H])
    
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
        print i,score

    # load the frame
    rx = int(fil.split('_')[-1].split('.')[0])
    rh= RMF.open_rmf_file('rmfs/models_%i.0.rmf' % (rx))
    IMP.rmf.link_hierarchies(rh, [prot])

    # save the frame
    frame_number = scores[min(scores.keys())]
    print
    print scores[min(scores.keys())], min(scores.keys())
    IMP.rmf.load_frame(rh, frame_number)
    m.update()
    output = output.Output()
    rmffile="models_tst_run"+seed+".rmf"
    output.init_rmf(rmffile, prot)
    output.write_rmf(rmffile,0)
exit()
"""
"""
#####################################################
### --- Score distribution
#####################################################

S = []
out = open('scores_run'+seed+'.out','w')
for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H]) # creates reverse dictionary

    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
    S.append(min(scores.keys()))
    print fil, min(scores.keys())
    out.write(fil+'\t'+str(min(scores.keys()))+'\n')
out.close()
print S, len(S)

pl.hist(S,bins=100,normed=0,linewidth=0)
pl.savefig('scores_run'+seed+'.png')
pl.show()

exit()
"""

#####################################################
### --- Analyze violations
#####################################################
import copy
Analysis = analysis.Violations('restraint_violations_34A.txt') # from where is this file coming?? 
X,Y=[],[]
Scores=[]
for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])
    otherscores={}
   
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
        otherscores[score] = [float(dct[H['SimplifiedModel_Linker_Score_None']]),\
                              float(dct[H['ExcludedVolumeSphere_None']]),\
                              float(dct[H['GaussianEMRestraint_kinase_em']]),\
                              float(dct[H['GaussianEMRestraint_rad3_em']]),\
                              float(dct[H['GaussianEMRestraint_tfb3_em']]),\
                              float(dct[H['GaussianEMRestraint_tfiicore_em']])]

    for score in sorted(scores.keys())[:1]:
        if score>scoreCutoff: continue  # consider only top 10% of models
        frame_number = scores[score]

        violrst = Analysis.get_number_violated_restraints( eval(D[frame_number].strip('\n')), eval(D[0].strip('\n')) )
        X.append(score); Y.append(violrst)
        print z,fil, score, violrst, otherscores[score]
        Scores.append(otherscores[score])
    

#pl.plot(X,Y,'k.')
#pl.savefig('score_vs_number_viol_restraints_run'+seed+'.png')
#pl.show()

V = Analysis.violation_counts  # only for cross-links, plot violated restraints over all structures
print
print Scores
#pl.bar(range(len(V)),sorted(V.values()))
#print V,len(V),'1111111'
#pl.savefig('bar_xlink_violation_counts_run'+seed+'.png')
#pl.show()

ofile=open('violation_counts_alltopscoring_34A.txt','w')
ofile.write(str(V))
ofile.close()

#exit()

"""
#####################################################
### --- Get clustered heat map  (Manual)
#####################################################
# Note that this only computes the all vs all RMSD matrix for top scoring models 
# To do the actual clustering we store the output pkl file and use clustering.py
# Clustering.py is using kmeans as it is more robust
Clusters = analysis.Clustering(output_filename='cluster_run'+seed) 
for cnt,fil in enumerate(files): 

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])
    scores = {}
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
    if sorted(scores.keys())[0]>scoreCutoff: continue #TODO

    # load the frame
    rx = int(fil.split('_')[-1].split('.')[0])
    rh= RMF.open_rmf_file('rmf/models_%i.0.rmf' % (rx ))  
   
    IMP.rmf.link_hierarchies(rh, [prot])

    for frm,score in enumerate(sorted(scores.keys())[:1]):
        frame_number = scores[score]
        print frame_number, scores[score]
        IMP.rmf.load_frame(rh, frame_number)
        m.update()
        #output.write_rmf(rmffile,1); output.close_rmf(rmffile); exit()

        Clusters.set_prot(prot)

        # set particles to calculate RMSDs on 
        # (if global, list all proteins, or just a few for local)    
        Coords = {}
        for pr in [p.get_name() for p in prot.get_children()]:  # for each subunit of each protein
            parts = IMP.atom.Selection(prot,molecule=pr).get_selected_particles()
            coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in parts])
            Coords[pr] = coords 
            
            #if pr=='tfb2' or pr=='ssl2': #TODO
	    #for i in parts: print i, IMP.core.XYZ(i).get_coordinates()

        Clusters.fill(fil+'.'+str(frame_number), Coords, alignment=0)

        print cnt,fil,score,frame_number
  
print Clusters.all_coords.keys()
print len(Clusters.all_coords.keys())

print "Global clustered heat map, not aligned"
Clusters.dist_matrix()
# these two were removed because we made analysis take in the filename for cluster file or pickle file
#os.rename("tmp_clustering.pdf", "Seh1_Sec13_global.pdf")
#os.remove("tmp_cluster_493.pkl") #info for plotting, later on. Also contains list of clusters. 
exit()
"""

"""
# This is needed for the analysis steps that depend on clustering
infile = open('cluster_run'+seed+'.out')
dominantClusterNumber=0  # this is true for both independent runs
Clusters = eval(infile.read())
infile.close()
import sys

files = [i.rsplit('.',1)[0] for i in Clusters[dominantClusterNumber][1]]

# just taking half of this cluster hoping to get better resolved localizations
files= files[:len(files)/2]
print len(files)
"""

"""
#####################################################
### --- Save as PDBs
#####################################################
for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i

    # load the frame
    rx = int(fil.split('_')[-1].split('.')[0])
    rh= RMF.open_rmf_file('%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
    IMP.rmf.link_hierarchies(rh, [prot])

    for frm,score in enumerate(sorted(scores.keys())[:1]):
        frame_number = scores[score]
        IMP.rmf.load_frame(rh, frame_number)
        m.update()

        outp = output.Output()
        outp.write_pdb_from_model(prot,name='cluster_0_%i.pdb' % z)
    #exit()

    print z,fil
exit()
"""
"""
#####################################################
### --- Get contact map
#####################################################
ContactMap = analysis.GetContactMap(distance=5.)

numModels=0

for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
    if sorted(scores.keys())[0]>scoreCutoff: continue 

    # load the frame
    rx = int(fil.split('_')[-1].split('.')[0])
    rh= RMF.open_rmf_file('rmf/models_%i.0.rmf' % (rx ))
    IMP.rmf.link_hierarchies(rh, [prot])

    for frm,score in enumerate(sorted(scores.keys())[:1]):
        frame_number = scores[score]
        IMP.rmf.load_frame(rh, frame_number)
        m.update()
        
        ContactMap.set_prot(prot)
    
        ContactMap.get_subunit_coords(fil+'.'+str(frm))
        
               
        print fil,frame_number
        
        numModels=numModels+1
        
        
print "Total number of models",numModels

ContactMap.add_xlinks('Ranish_Kornberg_thiih_xlinks_human.txt')  
ContactMap.dist_matrix(skip_cmap=0, skip_xl=0, outname='ContactMap_all_Matrix_CM.pkl')
exit()
"""
"""
#####################################################
### ---- Get global or local density map
#####################################################
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
                
                'tfb5' :[['tfb5',2,-1]]}

DensModule = analysis.GetModelDensity(margin=50., voxel=5.)  # margin is for making box bigger
#DensModule = analysis.GetModelDensity(margin=50., voxel=5., costum_ranges=costum_ranges)


for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])
    scores = {}
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
    if sorted(scores.keys())[0]>scoreCutoff: continue #TODO


    # load the frame
    rx = int(fil.split('_')[-1].split('.')[0])
    rh= RMF.open_rmf_file('rmf/models_%i.0.rmf' % (rx ))
    IMP.rmf.link_hierarchies(rh, [prot])

    for frm,score in enumerate(sorted(scores.keys())[:1]):
        frame_number = scores[score]
        print fil, frame_number
        IMP.rmf.load_frame(rh, frame_number)
        m.update()

        DensModule.fill(prot, alignment=0)
    print z,fil,score #,max(DensModule.densities[DensModule.densities.keys()[0]])

#print len(DensModule.densities[DensModule.densities.keys()[0]])

# for all top scoring runs (NOT the best model, but all top scoring)
#DensModule.write_mrc('tfiih_human_kinase_cluster_run'+seed)    # required

# for all top scoring models of both runs with domain splits
#DensModule.write_mrc('tfiih_human_kinase_domain_split_all_top')

# for all top scoring models of both runs without domain splits
#DensModule.write_mrc('tfiih_human_kinase_all_top')

# for dominant half of single cluster runs
#DensModule.write_mrc('tfiih_human_kinase_dominant_cluster_half_run'+seed)

# for dominant half of single cluster, with domain splits
#DensModule.write_mrc('tfiih_human_kinase_domain_split_dominant_cluster_run'+seed)
exit()
"""

