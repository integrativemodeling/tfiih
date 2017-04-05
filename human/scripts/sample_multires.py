#!/usr/bin/env python
from __future__ import print_function, absolute_import
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import os,sys
import IMP.rmf
import RMF
import restraints
import topology
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output

run_number = 0
xlinkAmplitude=17  

ncycl=100  #number of montecarlo steps cycles
outputobjects=[]
sampleobjects=[]

m, simo = topology.make_topology()

prot=simo.prot
outputobjects.append(simo)
sampleobjects.append(simo)

##simo.draw_hierarchy_composition()

#####################################################
#restraint setup
#####################################################

#EV restraint
ev=IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution=30)
ev.add_to_model()
outputobjects.append(ev)

#cross-link restraint
xl=IMP.pmi.restraints.crosslinking.SigmoidalCrossLinkMS(simo,
                     '../inputs/Ranish_Kornberg_thiih_xlinks_human.txt',
                     inflection=25.0,slope=5.0,amplitude=xlinkAmplitude,
                     linear_slope=0.05,resolution=1)
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
em1=restraints.GaussianEMRestraint(prot,
           '../inputs/tfiih_apo_multifit_moved.txt',segment_anchors=Kinase,
	   segment_parts=['ccl1','kin28'],resolution=30)
em1.set_label('kinase_em')
em1.add_to_model()
sampleobjects.append(em1)
outputobjects.append(em1)

# Sub component 2: Tfb3
em2=restraints.GaussianEMRestraint(prot,
        '../inputs/tfiih_apo_multifit_moved.txt',segment_anchors=Tfb3,
        segment_parts=['tfb3'],resolution=30)
em2.set_label('tfb3_em')
em2.add_to_model()
sampleobjects.append(em2)
outputobjects.append(em2)

#Sub component 2: Rad3
em3=restraints.GaussianEMRestraint(prot,
        '../inputs/tfiih_apo_multifit_moved.txt',segment_anchors=Rad3,
        segment_parts=['rad3'],resolution=30)
em3.set_label('rad3_em')
em3.add_to_model()
sampleobjects.append(em3)
outputobjects.append(em3)

## Sub component 3:Ssl2+Tfiihcore
em4=restraints.GaussianEMRestraint(prot,
        '../inputs/tfiih_apo_multifit_moved.txt',segment_anchors=Tfiihcore,
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
output.init_stat2("../outputs/stat_%d.dat" % run_number,
                  outputobjects,extralabels=["rmf_file","rmf_frame_index"])
# earlier stat file which had restraint scores etc was different from the RMF file. 

#####################################################
#running simulation
#####################################################
nrmffiles=1
bestscore, step = 1000000000000000, 0
nframes=20 if '--test' in sys.argv else 5000
for k in range(nrmffiles):
  rmffile="../outputs/models_%d.%d.rmf" % (run_number, k)
  output.init_rmf(rmffile, [prot])
  output.add_restraints_to_rmf(rmffile,[xl])
  rset = IMP.pmi.tools.get_restraint_set(m)

  for i in range(nframes):
    mc.optimize(ncycl)
    print(mc.get_frame_number())
    m.update()
    score = rset.evaluate(False)

    output.set_output_entry("rmf_file",rmffile)
    output.set_output_entry("rmf_frame_index",step)
    
    if score < bestscore:
      output.write_stats2()
      print('\tBest score: ',score)
      output.write_rmf(rmffile)
      bestscore = score
      step += 1
      
  output.close_rmf(rmffile)
