#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
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
xlinkAmplitude=17  

ncycl=100  #number of montecarlo steps cycles
rbmaxtrans=0.7
fbmaxtrans=0.7

outputobjects=[]
sampleobjects=[]


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

#GaussianEMRestraint
Kinase = [0,8,14,16,24,39,45,49,51,53,62,63,70,72,73,75,86,89,93,98,101,103,110,115,117,118,132,133,137,139,145,146,154,159,164,168,170]
Tfb3 = [3,9,20,21,36,59,67,91,94,97,121,126,156,167]
#Tfb3 = [2,3,9,13,20,21,29,36,59,61,67,79,82,91,94,97,107,114,121,124,126,149,156,166,167]
Rad3 = [4,7,18,23,26,30,31,33,37,41,43,47,52,57,60,65,84,85,99,102,108,112,119,120,134,147,148,149,152,155,157,161,162,163,165,166,169,171]
Tfiihcore = [1,2,5,6,10,11,12,13,15,17,19,22,25,26,28,29,31,32,34,35,36,37,38,40,42,44,46,48,50,52,54,55,56,59,60,61,66,67,68,69,71,74,76,77,
    78,79,80,81,82,83,88,90,91,95,96,100,104,105,106,107,109,111,112,113,114,116,122,123,124,125,127,128,129,130,131,135,136,138,140,141,143,
    144,149,150,151,153,157,158,160,162,165,166,171]
 
# Sub component 1: Kinase
em1=restraints.GaussianEMRestraint(prot,'tfiih_apo_multifit_moved.txt',segment_anchors=Kinase,
			      segment_parts=['ccl1','kin28'],resolution=30)
em1.set_label('kinase_em')
em1.add_to_model()
sampleobjects.append(em1)
outputobjects.append(em1)

# Sub component 2: Tfb3
em2=restraints.GaussianEMRestraint(prot,'tfiih_apo_multifit_moved.txt',segment_anchors=Tfb3,
			      segment_parts=['tfb3'],resolution=30)
em2.set_label('tfb3_em')
em2.add_to_model()
sampleobjects.append(em2)
outputobjects.append(em2)

#Sub component 2: Rad3
em3=restraints.GaussianEMRestraint(prot,'tfiih_apo_multifit_moved.txt',segment_anchors=Rad3,segment_parts=['rad3'],resolution=30)
em3.set_label('rad3_em')
em3.add_to_model()
sampleobjects.append(em3)
outputobjects.append(em3)

## Sub component 3:Ssl2+Tfiihcore
em4=restraints.GaussianEMRestraint(prot,'tfiih_apo_multifit_moved.txt',segment_anchors=Tfiihcore,
		segment_parts=['ssl2','tfb1','tfb2','tfb4','tfb5','ssl1'],resolution=30)
em4.set_label('tfiicore_em')
em4.add_to_model()
sampleobjects.append(em4)
outputobjects.append(em4)


####################################################
# Monte Carlo
####################################################

mc=samplers.MonteCarlo(m,sampleobjects,5.)
mc.set_simulated_annealing(5.,15.,50,10) # 50 steps at 5, 10 steps at 15. 
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

