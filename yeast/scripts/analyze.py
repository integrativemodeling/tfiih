#!/usr/bin/env python

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import os
import IMP.rmf
import RMF

import topology
import IMP.pmi.restraints as restraints
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output
import IMP.pmi.analysis as analysis

#####################################################
#analyze RMFs
#####################################################


import pickle
import numpy as np
import glob

# Get stats files from all independent runs
files = glob.glob('../outputs/stat*.dat')

m, simo = topology.make_topology()
prot = simo.prot

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
    rh= RMF.open_rmf_file_read_only(
                    '%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
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
import matplotlib
matplotlib.use('Agg')
import pylab as pl
pl.hist(S,bins=100,normed=0,linewidth=0)
pl.savefig('score_distribution.pdf')

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
        dct = eval(D[frame_number].strip('\n'))
        score_dct = {}
        for k in H.keys():
            if isinstance(H[k], int):
                score_dct[k] = dct[H[k]]

        violrst = Analysis.get_number_violated_restraints(score_dct)
        X.append(score); Y.append(violrst)
        print z,fil, score, violrst, otherscores[score]
        Scores.append( otherscores[score] )
    
pl.plot(X,Y,'k.')
pl.savefig('violations.pdf')

V = Analysis.violation_counts
pl.bar(range(len(V)),sorted(V.values()))
print V,len(V),'1111111'
print
print Scores
pl.savefig('violation_counts.pdf')


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
    rh= RMF.open_rmf_file_read_only(
                        '%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
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
    rh= RMF.open_rmf_file_read_only(
                      '%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
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
    rh= RMF.open_rmf_file_read_only(
                        '%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
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
    rh= RMF.open_rmf_file_read_only(
                         '%s/models_%i.0.rmf' % (fil.rsplit('/',1)[0], rx ))
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
