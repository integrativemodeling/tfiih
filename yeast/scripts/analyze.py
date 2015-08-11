"""
import RMF
import sys

def _do_add(node, key, rbi):
    print 'node: ',node.get_name(), len(node.get_children())
    if node.get_has_value(key):
        rbi = node.get_value(key)
    elif len(node.get_children()) == 0:
        assert(rbi >= 0)
        print "adding", node.get_name(), "to", rbi
        node.set_value(key, rbi)
    else:
        for c in node.get_children():
            rbi = _do_add(c, key, rbi)
    return rbi

def _add_file(name):
    print "processing", name
    f = RMF.open_rmf_file(name)
    imp_cat = f.get_category("IMP")
    rbk = f.get_index_key(imp_cat, "rigid body")
    f.set_current_frame(RMF.ALL_FRAMES)
    _do_add(f.get_root_node(), rbk, None)

for f in sys.argv[1:]:
    _add_file(f)
    print f


exit()
"""


#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
#import IMP.isd2
import random
import os
import IMP.rmf
import RMF
import inspect

from numpy import random
#from numpy import *

from math import cos
from math import sqrt

import IMP.pmi.restraints as restraints
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output
import IMP.pmi.analysis as analysis


'''
#####################################################
#read EM maps and centroid coordinates
#####################################################
map_kinase = IMP.em.read_map('kinase.mrc', IMP.em.MRCReaderWriter())
map_rad3 = IMP.em.read_map('rad3.mrc', IMP.em.MRCReaderWriter())
map_ssl2 = IMP.em.read_map('ssl2.mrc', IMP.em.MRCReaderWriter())
map_tfb = IMP.em.read_map('tfb.mrc', IMP.em.MRCReaderWriter())
'''

ncycl=100  #number of montecarlo steps cycles
rbmaxtrans=0.5
fbmaxtrans=0.5
outputobjects=[]
sampleobjects=[]


#####################################################
#create hierarchies and rigid bodies and flexible parts
#####################################################

m=IMP.Model()
simo=representation.SimplifiedModel(m, upperharmonic=False, disorderedlength=True)



activated_components=[1,2,3,4]


if 1 in activated_components:
  
  simo.add_component_name("ccl1")
  simo.add_pdb_and_intervening_beads("ccl1",'../inputs/CCL1_48-393.pdb',"A",resolutions=[1,30],resrange=(1,393),beadsize=30,color=1.0,attachbeads=True)
  simo.setup_component_sequence_connectivity("ccl1", resolution=30)

  simo.add_component_name("kin28")
  simo.add_pdb_and_intervening_beads("kin28",'../inputs/KIN28_5-299.pdb',"A",resolutions=[1,30],resrange=(1,306), beadsize=30,color=0.9,attachbeads=True)
  simo.setup_component_sequence_connectivity("kin28", resolution=30)

  simo.add_component_name("tfb3")
  simo.add_pdb_and_intervening_beads("tfb3",'../inputs/TFB3_8-142.pdb',"A",resolutions=[1,30],resrange=(1,321), beadsize=30,color=0.8,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb3", resolution=30)


if 2 in activated_components:

  simo.add_component_name("rad3")
  simo.add_pdb_and_intervening_beads("rad3",'../inputs/RAD3_14-725.pdb',"A",resolutions=[1,30],resrange=(1,778), beadsize=30,color=150./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("rad3", resolution=30)


if 4 in activated_components:

  simo.add_component_name("tfb1")
  simo.add_pdb_and_intervening_beads("tfb1",'../inputs/TFB1_2-115.pdb',"A",resolutions=[1,30],resrange=(1,642), beadsize=30,color=60./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb1", resolution=30)

  simo.add_component_name("tfb2")
  #simo.add_pdb_and_intervening_beads("tfb2",'../inputs/TFB2_1-168.pdb',"A",resolutions=[1,30],resrange=(1,170), beadsize=30,color=185./360,attachbeads=True)
  #simo.add_pdb_and_intervening_beads("tfb2",'../inputs/TFB2_186-417.pdb',"A",resolutions=[1,30],resrange=(171,417), beadsize=30,color=185./360,attachbeads=True)
  simo.add_pdb_and_intervening_beads("tfb2",'../inputs/TFB2_392-513.pdb',"A",resolutions=[1,30],resrange=(418,513), beadsize=30,color=185./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb2", resolution=30)

  simo.add_component_name("tfb4")
  simo.add_pdb_and_intervening_beads("tfb4",'../inputs/TFB4_24-250.pdb',"A",resolutions=[1,30],resrange=(1,338), beadsize=30,color=210./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb4", resolution=30)

  simo.add_component_name("tfb5")
  simo.add_pdb_and_intervening_beads("tfb5",'../inputs/TFB5.pdb',"B",resolutions=[1,30],resrange=(1,72), beadsize=30,color=0.,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb5", resolution=30)

  simo.add_component_name("ssl1")
  simo.add_pdb_and_intervening_beads("ssl1",'../inputs/SSL1_123-302.pdb',"A",resolutions=[1,30],resrange=(1,302), beadsize=30,color=285./360,attachbeads=True)
  simo.add_pdb_and_intervening_beads("ssl1",'../inputs/SSL1_386-455.pdb',"A",resolutions=[1,30],resrange=(303,461), beadsize=30,color=285./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("ssl1", resolution=30)


if 1 in activated_components:
  simo.add_component_name("ssl2")
  simo.add_pdb_and_intervening_beads("ssl2",'../inputs/fit_gtfs_3.pdb',"C",resolutions=[1,30],resrange=(1,538), beadsize=30,color=0.0,attachbeads=True)
  simo.add_pdb_and_intervening_beads("ssl2",'../inputs/fit_gtfs_3.pdb',"D",resolutions=[1,30],resrange=(539,843), beadsize=30,color=0.0,attachbeads=True)
  simo.setup_component_sequence_connectivity("ssl2", resolution=30)






# --- set rigid bodies
if 1 in activated_components:
  emxk,emyk,emzk = 215.249,220.731,167.764 #21
  simo.set_rigid_bodies(["tfb3"],(emxk,emyk,emzk))
  simo.set_rigid_bodies(["ccl1","kin28"],(emxk,emyk,emzk))

if 2 in activated_components:
  emxr,emyr,emzr = 147.133,219.683,188.062 #0
  simo.set_rigid_bodies(["rad3"],(emxr,emyr,emzr))

if 4 in activated_components:
  emxt,emyt,emzt = 126.646,130.588,192.099 #90
  simo.set_rigid_bodies([("tfb2",(437,513)),("tfb5",(1,72))], (emxt,emyt,emzt))
  simo.set_rigid_bodies([("tfb2",(1,436))], (emxt,emyt,emzt))
  #simo.set_rigid_bodies([("tfb2",(1,168))], (emxt,emyt,emzt))
  #simo.set_rigid_bodies([("tfb2",(169,436))], (emxt,emyt,emzt))
  simo.set_rigid_bodies(["tfb1"], (emxt,emyt,emzt))
  simo.set_rigid_bodies(["tfb4"],(emxt,emyt,emzt))
  simo.set_rigid_bodies([("ssl1",(1,385))],(emxt,emyt,emzt))
  simo.set_rigid_bodies([("ssl1",(386,461))],(emxt,emyt,emzt))


if 3 in activated_components:
  emxk,emyk,emzk = 212.272,129.032,178.206 #80
  simo.set_rigid_bodies(["ssl2"],(emxk,emyk,emzk))


simo.set_floppy_bodies()
d=simo.get_particles_to_sample()
simo.set_rigid_bodies_max_trans(rbmaxtrans)
simo.set_floppy_bodies_max_trans(fbmaxtrans)


#re-orient initial positions
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
import glob, random


files = glob.glob('../outputs/stat*.dat')


"""
#####################################################
### --- Test
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
    rh= RMF.open_rmf_file('%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
    IMP.rmf.link_hierarchies(rh, [prot])

    # save the frame
    frame_number = scores[min(scores.keys())]
    print
    print scores[min(scores.keys())], min(scores.keys())
    IMP.rmf.load_frame(rh, frame_number)
    m.update()
    output = output.Output()
    rmffile="models_tst.rmf"
    output.init_rmf(rmffile, prot)
    output.write_rmf(rmffile,0)
exit()
"""


#####################################################
### --- Score distribution
#####################################################

S = []
out = open('scores.out','w')
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
    S.append(min(scores.keys()))
    print fil, min(scores.keys())
    out.write(fil+'\t'+str(min(scores.keys()))+'\n')
out.close()
print S, len(S)
import pylab as pl
pl.hist(S,bins=100,normed=0,linewidth=0)
pl.show()




#####################################################
### --- Analyze violations
#####################################################
import copy
Analysis = analysis.Violations('../inputs/restraint_violations.txt')
X,Y=[],[]
Scores = []
for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    scores = {}
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])

    otherscores = {}
   
    for i,d in enumerate(D[1:]):
        dct = eval(d.strip('\n'))
        score = float(dct[H['SimplifiedModel_Total_Score_None']])
        scores[score] = i
        otherscores[score] = [float(dct[H['SimplifiedModel_Linker_Score_None']]),\
                              float(dct[H['ExcludedVolumeSphere_None']]),\
                              float(dct[H['GaussianEMRestraint_kinase_em']]),\
                              float(dct[H['GaussianEMRestraint_rad3_em']]),\
                              float(dct[H['GaussianEMRestraint_ssl2_em']]),\
                              float(dct[H['GaussianEMRestraint_tfiihcore_em']])]

    for score in sorted(scores.keys())[:1]:
        if score>1000.: continue
        frame_number = scores[score]

        violrst = Analysis.get_number_violated_restraints( eval(D[frame_number].strip('\n')), eval(D[0].strip('\n')) )
        X.append(score); Y.append(violrst)
        print z,fil, score, violrst, otherscores[score]
        Scores.append( otherscores[score] )
    
import pylab as pl
pl.plot(X,Y,'k.')
pl.show()

V = Analysis.violation_counts
pl.bar(range(len(V)),sorted(V.values()))
print V,len(V),'1111111'
print
print Scores
pl.show()





#####################################################
### --- Get clustered heat map  (Manual)
#####################################################
Clusters = analysis.Clustering()
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
    if sorted(scores.keys())[0]>1000.: continue

    # load the frame
    rx = int(fil.split('_')[-1].split('.')[0])
    rh= RMF.open_rmf_file('%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
    IMP.rmf.link_hierarchies(rh, [prot])

    for frm,score in enumerate(sorted(scores.keys())[:1]):
        frame_number = scores[score]
        print frame_number, scores[score]
        IMP.rmf.load_frame(rh, frame_number)
        m.update()
        #output.write_rmf(rmffile,1); output.close_rmf(rmffile); exit()

        Clusters.set_prot(prot)

        # set alignment template (if global, list all proteins, or just a few for local)   
        ''' 
        template = {}
        for pr in ['Seh1..1', 'Sec13']:
            parts = IMP.atom.Selection(prot,molecule=pr).get_selected_particles()
            coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in parts])
            template[pr] = coords

        if len(Clusters.all_coords)==0:
            Clusters.set_template(template)
        '''
        # set particles to calculate RMSDs on 
        # (if global, list all proteins, or just a few for local)    
        Coords = {}
        for pr in [p.get_name() for p in prot.get_children()]:
            parts = IMP.atom.Selection(prot,molecule=pr).get_selected_particles()
            coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in parts])
            Coords[pr] = coords  

        Clusters.fill(fil+'.'+str(frame_number), Coords, alignment=0)

        print cnt,fil,score,frame_number
    #if cnt==1: break
print Clusters.all_coords.keys()
print len(Clusters.all_coords.keys())

print "Global clustered heat map, not aligned"
Clusters.dist_matrix()
#os.rename("tmp_clustering.pdf", "Seh1_Sec13_global.pdf")
#os.remove("tmp_cluster_493.pkl")




# use a subset of output files corresponding to a given cluster of solutions
files = files # will use all this time


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


#####################################################
### --- Get contact map
#####################################################
ContactMap = analysis.GetContactMap(distance=5.)
print files
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
        
        ContactMap.set_prot(prot)
    
        ContactMap.get_subunit_coords(fil+'.'+str(frm))
    print z,fil

ContactMap.add_xlinks('Ranish_Kornberg_thiih_xlinks.txt')  
ContactMap.dist_matrix(skip_cmap=0, skip_xl=0, outname='ContactMap_cluster_%i' % int(sys.argv[-1]))




#####################################################
### ---- Get global or local density map
#####################################################
costum_ranges ={'kin28':[['kin28',1,-1]],
                'ccl1' :[['ccl1',1,-1]],
                'tfb3' :[['tfb3d1',1,142],['tfb3d2',143,-1]],
                'rad3' :[['rad3d1',1,235],['rad3d2',236,452],['rad3d3',453,-1]],
                'ssl2' :[['ssl2d1',1,297],['ssl2d2',298,538],['ssl2d3',538,712],['ssl2d4',713,-1]],
                'ssl1' :[['ssl1d1',1,122],['ssl1d2',123,302],['ssl1d3',386,-1]],
                'tfb1' :[['tfb1d1',2,115],['tfb1d2',116,400],['tfb1d3',401,-1]],
                'tfb2' :[['tfb2d1',1,168],['tfb2d2',186,417],['tfb2d3',418,-1]],
                'tfb4' :[['tfb4d1',1,250],['tfb4d2',251,-1]],
                'tfb5' :[['tfb5',2,-1]]}


costum_ranges ={  # these domains are from the paper
                'kin28' :[['kin28_nlobe',1,88],['kin28_clobe',89,-1]],

                'ccl1' :[['ccl1_n',1,206],['ccl1_c',207,-1]],

                'tfb3' :[['tfb3_n',1,142],['tfb3_latch',142,-1]],

                'rad3' :[['rad3d1',1,246],['rad3d2',247,440],['rad3d3',441,701],['rad3d4',702,-1]],
                
                'ssl2' :[['ssl2d1',1,350],['ssl2d2',351,518],['ssl2d3',519,714],['ssl2d4',715,-1]],
           
		'ssl1' :[['ssl1d1',1,309],['ssl1d2',310,-1]],
                
                'tfb1' :[['tfb1d1',2,115],['tfb1d2',116,332],['tfb1d3',333,546],['tfb1d4',547,-1]],
                
                'tfb2' :[['tfb2d1',1,262],['tfb2d2',263,432],['tfb2d3',433,-1]],
                
                'tfb4' :[['tfb4d1',1,250],['tfb4d2',251,-1]],
                
                'tfb5' :[['tfb5',2,-1]]}


DensModule = analysis.GetModelDensity(margin=50., voxel=5., costum_ranges=costum_ranges)
for z, fil in enumerate(files):

    data = open(fil)
    D = data.readlines()
    data.close()

    # find frame with the lowest score
    #scores = {}
    H = eval(D[0].strip())
    H = dict([(H[c],c) for c in H])
    
    scores = {}
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

        '''
        Coords = {}
        for pr in prot.get_children():
            #parts = IMP.atom.Selection(prot,molecule=pr.get_name()).get_selected_particles()
            parts = IMP.pmi.tools.get_particles_by_resolution(pr,1.)
            coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in parts])
            Coords[pr.get_name()] = coords  
        '''
        DensModule.fill(prot, alignment=0)
    print z,fil,score#,max(DensModule.densities[DensModule.densities.keys()[0]])

#print len(DensModule.densities[DensModule.densities.keys()[0]])
DensModule.write_mrc('tfiih_cluster_1')






