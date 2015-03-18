#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
#import IMP.isd2
import random
import os,sys
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

seed=sys.argv[1]

nframes=5000
nframes=10
xlinkAmplitude=17  

ncycl=100  #number of montecarlo steps cycles
rbmaxtrans=1.0
fbmaxtrans=1.0

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
#simo.add_component_name("ccl1")
#simo.add_pdb_and_intervening_beads("ccl1",'CCL1.pdb',"A",resolutions=[1,30],resrange=(?,?),beadsize=30,color=60./360,attachbeads=True)
#simo.setup_component_sequence_connectivity("ccl1", resolution=30)

#simo.add_component_name("kin28")
#simo.add_pdb_and_intervening_beads("kin28",'KIN28.pdb',"A",resolutions=[1,30],resrange=(?,?), beadsize=30,color=185./360,attachbeads=True)
#simo.setup_component_sequence_connectivity("kin28", resolution=30)

#simo.add_component_name("tfb3")
#simo.add_pdb_and_intervening_beads("tfb3",'TFB3.pdb',"A",resolutions=[1,30],resrange=(?,?), beadsize=30,color=285./360,attachbeads=True)
#simo.setup_component_sequence_connectivity("tfb3", resolution=30)

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
#emxk,emyk,emzk = ?, ?, ? #99
#simo.set_rigid_bodies(["tfb3"],(emxk,emyk,emzk))
#simo.set_rigid_bodies(["ccl1","kin28"],(emxk,emyk,emzk))

# Sub component 2: Rad3
emxr,emyr,emzr =  39.118, 7.312, -11.998 # 75
simo.set_rigid_bodies(["rad3"],(emxr,emyr,emzr))

# Sub component 3: Ssl2 + Tfiihcore
emxt,emyt,emzt = -42.520,  -2.832, -10.375 # 123
simo.set_rigid_bodies([("tfb2"),("tfb5")], (emxt,emyt,emzt))

# No need to set extra rigid bodies: do so only when you have the structure.
#simo.set_rigid_bodies([("tfb2",(395,453)),("tfb5")], (emxt,emyt,emzt))

#simo.set_rigid_bodies([("tfb2",(1,394))], (emxt,emyt,emzt))
#simo.set_rigid_bodies([("tfb2",(454,462))], (emxt,emyt,emzt))

simo.set_rigid_bodies(["tfb1"], (emxt,emyt,emzt))
simo.set_rigid_bodies(["tfb4"],(emxt,emyt,emzt))
simo.set_rigid_bodies(["ssl1"],(emxt,emyt,emzt))
#simo.set_rigid_bodies([("ssl1",(386,461))],(emxt,emyt,emzt))
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
#restraint setup
#####################################################

#EV restraint
ev=restraints.ExcludedVolumeSphere(prot, resolution=30)
ev.add_to_model()
outputobjects.append(ev)

#cross-link restraint
xl=restraints.SigmoidCrossLinkMS(prot,'Ranish_Kornberg_thiih_xlinks_human.txt',inflection=25.0,slope=5.0,amplitude=xlinkAmplitude,linear_slope=0.05,resolution=1)
#xl=restraints.SimplifiedCrossLinkMS(prot,'XLINKS_Ranish_Kornberg.txt',expdistance=17.,strength=0.2,resolution=1)
xl.add_to_model()
outputobjects.append(xl)
#xl.plot_restraint(5.,5.,maxdist=100,npoints=20)

#GaussianEMRestraint
# Kinase = []
Rad3 = [2,3,4,13,14,15,25,27,28,29,30,31,32,33,34,35,36,38,39,40,41,42,43,44,45,51,52, \
   53,64,65,66,67,68,69,74,75,76,77,78,79,80,81,82,83,91,92,93,94,95,98,99,100,106,107,108,111, \
   112,113,114,115,116,117,118,119,120,121,122,123,124,125]

Ssl2Tfiihcore = [i for i in range(133) if i not in Rad3]
 
# Sub component 1: Kinase
#em1=restraints.GaussianEMRestraint(prot,'TFIIH_core.txt',segment_anchors=Kinase,
			      #segment_parts=['tfb3','ccl1','kin28'],resolution=30)
#em1.set_label('kinase_em')
#em1.add_to_model()
#sampleobjects.append(em1)
#outputobjects.append(em1)

#Sub component 2: Rad3
em2=restraints.GaussianEMRestraint(prot,'TFIIH_core_match_to_mrc.txt',segment_anchors=Rad3,segment_parts=['rad3'],resolution=30)
em2.set_label('rad3_em')
em2.add_to_model()
sampleobjects.append(em2)
outputobjects.append(em2)

## Sub component 3:Ssl2+Tfiihcore
em3=restraints.GaussianEMRestraint(prot,'TFIIH_core_match_to_mrc.txt',segment_anchors=Ssl2Tfiihcore,
		segment_parts=['ssl2','tfb1','tfb2','tfb4','tfb5','ssl1'],resolution=30)
em3.set_label('ssl2_tfiicore_em')
em3.add_to_model()
sampleobjects.append(em3)
outputobjects.append(em3)


####################################################
# Monte Carlo
####################################################

mc=samplers.MonteCarlo(m,sampleobjects,5.)
mc.set_simulated_annealing(5.,15.,50,10) # 50 steps at 15, 10 steps at 45. 
#mc.set_simulated_annealing(25.,75.,50,10)
#cg=samplers.ConjugateGradients(m,sampleobjects)
outputobjects.append(mc)
#outputobjects.append(cg)

sw = tools.Stopwatch()
outputobjects.append(sw)

output = output.Output()
output.init_stat2("stat_%s.dat" % seed, outputobjects,extralabels=["rmf_file","rmf_frame_index"])
# earlier stat file which had restraint scores etc was different from the RMF file. 


#####################################################
#running simulation
#####################################################
nrmffiles=1
# nframes=500
bestscore, step = 1000000000000000, 0
for k in range(nrmffiles):
  rmffile="models_%s.%d.rmf" % (seed,k)
  output.init_rmf(rmffile, prot)
  output.add_restraints_to_rmf(rmffile,[xl])

  for i in range(nframes):
    mc.run(ncycl)
    print mc.get_frame_number()
    score = m.evaluate(False)

    output.set_output_entry("rmf_file",rmffile)
    output.set_output_entry("rmf_frame_index",step)
    
    if score < bestscore:
      output.write_stats2()
      print '\tBest score: ',score
      output.write_rmf(rmffile,step)
      bestscore = score
      step += 1
      
  output.close_rmf(rmffile)

